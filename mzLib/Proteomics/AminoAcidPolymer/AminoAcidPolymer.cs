// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (AminoAcidPolymer.cs) is part of Proteomics.
//
// Proteomics is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Proteomics is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with Proteomics. If not, see <http://www.gnu.org/licenses/>.

using Chemistry;
using Easy.Common.Extensions;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Text;

namespace Proteomics.AminoAcidPolymer
{
    /// <summary>
    /// A linear polymer of amino acids
    /// </summary>
    public abstract class AminoAcidPolymer : IEquatable<AminoAcidPolymer>, IHasMass
    {
        /// <summary>
        /// The C-terminus chemical formula cap. This is different from the C-Terminus modification.
        /// </summary>
        private IHasChemicalFormula _cTerminus;

        /// <summary>
        /// The N-terminus chemical formula cap. This is different from the N-Terminus modification.
        /// </summary>
        private IHasChemicalFormula _nTerminus;

        /// <summary>
        /// All of the modifications indexed by position from N to C. This array is 2 bigger than the amino acid array
        /// as index 0 and Count - 1 represent the N and C terminus, respectively
        /// </summary>
        private IHasMass[] _modifications;

        /// <summary>
        /// All of the amino acid residues indexed by position from N to C.
        /// </summary>
        private Residue[] residues;

        protected AminoAcidPolymer()
            : this(string.Empty, new ChemicalFormulaTerminus(ChemicalFormula.ParseFormula("H")), new ChemicalFormulaTerminus(ChemicalFormula.ParseFormula("OH")))
        {
        }

        protected AminoAcidPolymer(string sequence)
            : this(sequence, new ChemicalFormulaTerminus(ChemicalFormula.ParseFormula("H")), new ChemicalFormulaTerminus(ChemicalFormula.ParseFormula("OH")))
        {
        }

        protected AminoAcidPolymer(string sequence, IHasChemicalFormula nTerm, IHasChemicalFormula cTerm)
        {
            MonoisotopicMass = 0;
            Length = sequence.Length;
            residues = new Residue[Length];
            NTerminus = nTerm;
            CTerminus = cTerm;
            ParseSequence(sequence);
        }

        protected AminoAcidPolymer(AminoAcidPolymer aminoAcidPolymer, bool includeModifications)
            : this(aminoAcidPolymer, 0, aminoAcidPolymer.Length, includeModifications)
        {
        }

        protected AminoAcidPolymer(AminoAcidPolymer aminoAcidPolymer, int firstResidue, int length, bool includeModifications)
        {
            Length = length;
            residues = new Residue[length];

            bool isNterm = firstResidue == 0;
            bool isCterm = length + firstResidue == aminoAcidPolymer.Length;

            _nTerminus = isNterm ? aminoAcidPolymer.NTerminus : new ChemicalFormulaTerminus(ChemicalFormula.ParseFormula("H"));
            _cTerminus = isCterm ? aminoAcidPolymer.CTerminus : new ChemicalFormulaTerminus(ChemicalFormula.ParseFormula("OH"));

            double monoMass = _nTerminus.MonoisotopicMass + _cTerminus.MonoisotopicMass;

            Residue[] otherAminoAcids = aminoAcidPolymer.residues;

            if (includeModifications && aminoAcidPolymer.ContainsModifications())
            {
                _modifications = new IHasMass[length + 2];
                for (int i = 0; i < length; i++)
                {
                    var aa = otherAminoAcids[i + firstResidue];
                    residues[i] = aa;
                    monoMass += aa.MonoisotopicMass;

                    IHasMass mod = aminoAcidPolymer._modifications[i + firstResidue + 1];
                    if (mod == null)
                        continue;

                    _modifications[i + 1] = mod;
                    monoMass += mod.MonoisotopicMass;
                }
            }
            else
            {
                for (int i = 0, j = firstResidue; i < length; i++, j++)
                {
                    var aa = otherAminoAcids[j];
                    residues[i] = aa;
                    monoMass += aa.MonoisotopicMass;
                }
            }

            MonoisotopicMass = monoMass;

            if (includeModifications)
            {
                if (isNterm)
                    NTerminusModification = aminoAcidPolymer.NTerminusModification;

                if (isCterm)
                    CTerminusModification = aminoAcidPolymer.CTerminusModification;
            }
        }

        public ReadOnlyCollection<IHasMass> Modifications
        {
            get
            {
                return new ReadOnlyCollection<IHasMass>(_modifications);
            }
        }

        /// <summary>
        /// Gets or sets the C terminus of this amino acid polymer
        /// </summary>
        public IHasChemicalFormula CTerminus
        {
            get { return _cTerminus; }
            set { ReplaceTerminus(ref _cTerminus, value); }
        }

        /// <summary>
        /// Gets or sets the N terminus of this amino acid polymer
        /// </summary>
        public IHasChemicalFormula NTerminus
        {
            get { return _nTerminus; }
            set { ReplaceTerminus(ref _nTerminus, value); }
        }

        /// <summary>
        /// Gets the number of amino acids in this amino acid polymer
        /// </summary>
        public int Length { get; private set; }

        /// <summary>
        /// The total monoisotopic mass of this peptide and all of its modifications
        /// </summary>
        public double MonoisotopicMass { get; private set; }

        /// <summary>
        /// Returns the amino acid sequence with all isoleucines (I) replaced with leucines (L);
        /// </summary>
        /// <returns>The amino acid sequence with all I's into L's</returns>
        public virtual string BaseLeucineSequence
        {
            get
            {
                return BaseSequence.Replace('I', 'L');
            }
        }

        /// <summary>
        /// Gets the base amino acid sequence
        /// </summary>
        public string BaseSequence
        {
            get
            {
                return new string(residues.Select(aa => aa.Letter).ToArray());
            }
        }

        /// <summary>
        /// Gets or sets the modification of the C terminus on this amino acid polymer
        /// </summary>
        public IHasMass CTerminusModification
        {
            get { return GetModification(Length + 1); }
            set { ReplaceMod(Length + 1, value); }
        }

        /// <summary>
        /// Gets or sets the modification of the C terminus on this amino acid polymer
        /// </summary>
        public IHasMass NTerminusModification
        {
            get { return GetModification(0); }
            set { ReplaceMod(0, value); }
        }

        /// <summary>
        /// Returns all fragments that are present in either fragmentation of A or B, but not in both
        /// </summary>
        public static IEnumerable<Fragment> GetSiteDeterminingFragments(AminoAcidPolymer peptideA, AminoAcidPolymer peptideB, FragmentTypes types)
        {
            HashSet<Fragment> aFrags = new HashSet<Fragment>(peptideA.Fragment(types));
            aFrags.SymmetricExceptWith(peptideB.Fragment(types));
            return aFrags;
        }

        /// <summary>
        /// Gets the digestion points (starting index and length) of a amino acid sequence
        /// </summary>
        /// <param name="sequence">The sequence to cleave</param>
        /// <param name="proteases">The proteases to cleave with</param>
        /// <param name="maxMissedCleavages">The maximum number of missed clevages to allow</param>
        /// <param name="minLength">The minimum amino acid length of the peptides</param>
        /// <param name="maxLength">The maximum amino acid length of the peptides</param>
        /// <param name="methionineInitiator"></param>
        /// <param name="semiDigestion"></param>
        /// <returns>A collection of clevage points and the length of the cut (Item1 = index, Item2 = length)</returns>
        public static IEnumerable<DigestionPointAndLength> GetDigestionPointsAndLengths(string sequence, IEnumerable<IProtease> proteases, int maxMissedCleavages, int minLength, int maxLength, bool methionineInitiator, bool semiDigestion)
        {
            int[] indices = GetCleavageIndexes(sequence, proteases).ToArray();

            bool includeMethionineCut = methionineInitiator && sequence[0] == 'M';

            int indiciesCount = indices.Length - 1;

            for (int missedCleavages = 0; missedCleavages <= maxMissedCleavages; missedCleavages++)
            {
                int max = indiciesCount - missedCleavages;
                int offset = missedCleavages + 1;
                for (int i = 0; i < max; i++)
                {
                    int len = indices[i + offset] - indices[i];

                    // Case for initiator methionine
                    if (indices[i] == -1 && includeMethionineCut)
                    {
                        int newLength = len - 1;
                        if (newLength >= minLength && newLength <= maxLength)
                        {
                            yield return new DigestionPointAndLength(1, newLength); if (semiDigestion)
                            {
                                for (int j = 1; j < newLength; j++)
                                {
                                    if (j >= minLength && j <= maxLength)
                                    {
                                        yield return new DigestionPointAndLength(1, j);
                                    }
                                }
                            }
                        }
                    }

                    if (len < minLength || len > maxLength)
                        continue;

                    yield return new DigestionPointAndLength(indices[i] + 1, len);
                    if (semiDigestion)
                    {
                        for (int j = 1; j < len; j++)
                        {
                            if (len - j >= minLength && len - j <= maxLength)
                            {
                                yield return new DigestionPointAndLength(indices[i] + 1 + j, len - j);
                            }
                            if (j >= minLength && j <= maxLength)
                            {
                                yield return new DigestionPointAndLength(indices[i] + 1, j);
                            }
                        }
                    }
                }
            }
        }

        public static IEnumerable<int> GetCleavageIndexes(string sequence, IEnumerable<IProtease> proteases)
        {
            return GetCleavageIndexes(sequence, proteases, true);
        }

        /// <summary>
        /// Gets the location of all the possible cleavage points for a given sequence and set of proteases
        /// </summary>
        /// <param name="sequence">The sequence to determine the cleavage points for</param>
        /// <param name="proteases">The proteases to cleave with</param>
        /// <param name="includeTermini">Include the N and C terminus (-1 and Length + 1)</param>
        /// <returns>A collection of all the sites where the proteases would cleave</returns>
        public static IEnumerable<int> GetCleavageIndexes(string sequence, IEnumerable<IProtease> proteases, bool includeTermini)
        {
            // Combine all the proteases digestion sites
            SortedSet<int> locations = new SortedSet<int>();
            foreach (IProtease protease in proteases.Where(protease => protease != null))
            {
                locations.UnionWith(protease.GetDigestionSites(sequence));
            }

            if (!includeTermini)
                return locations;

            locations.Add(-1);
            locations.Add(sequence.Length - 1);

            return locations;
        }

        public static IEnumerable<string> Digest(string sequence, IEnumerable<IProtease> proteases, int maxMissedCleavages, int minLength, int maxLength, bool methionineInitiator, bool semiDigestion)
        {
            return GetDigestionPointsAndLengths(sequence, proteases, maxMissedCleavages, minLength, maxLength, methionineInitiator, semiDigestion).Select(points => sequence.Substring(points.Index, points.Length));
        }

        public static IEnumerable<string> Digest(AminoAcidPolymer sequence, IProtease protease)
        {
            return Digest(sequence, protease, 3, 1, int.MaxValue, true, false);
        }

        public static IEnumerable<string> Digest(AminoAcidPolymer polymer, IProtease protease, int maxMissedCleavages, int minLength, int maxLength, bool methionineInitiator, bool semiDigestion)
        {
            return Digest(polymer.BaseSequence, new[] { protease }, maxMissedCleavages, minLength, maxLength, methionineInitiator, semiDigestion);
        }

        public Residue GetResidue(int position)
        {
            if (position < 0 || position >= Length)
                return null;
            return residues[position];
        }

        /// <summary>
        /// Checks if an amino acid residue with the value of 'residue' is contained in this polymer
        /// </summary>
        /// <param name="residue">The character code for the amino acid residue</param>
        /// <returns>True if any amino acid residue is the same as the specified character</returns>
        public bool Contains(char residue)
        {
            return residues.Any(aa => aa.Letter.Equals(residue));
        }

        /// <summary>
        /// Checks if the amino acid residue is contained in this polymer
        /// </summary>
        /// <param name="residue">The residue to check for</param>
        /// <returns>True if the polymer contains the specified residue, False otherwise</returns>
        public bool Contains(Residue residue)
        {
            return residues.Contains(residue);
        }

        public string GetSequenceWithModifications()
        {
            return GetSequenceWithModifications(false);
        }

        public string GetSequenceWithModifications(bool leucineSequence)
        {
            if (_modifications == null)
                return (leucineSequence) ? BaseLeucineSequence : BaseSequence;

            StringBuilder modSeqSb = new StringBuilder(Length);

            IHasMass mod;

            // Handle N-Terminus Modification
            if ((mod = _modifications[0]) != null && mod.MonoisotopicMass > 0)
            {
                modSeqSb.Append('[');
                modSeqSb.Append(mod);
                modSeqSb.Append("]-");
            }

            // Handle Amino Acid Residues
            for (int i = 0; i < Length; i++)
            {
                if (leucineSequence && residues[i].Letter == 'I')
                    modSeqSb.Append('L');
                else
                    modSeqSb.Append(residues[i].Letter);

                // Handle Amino Acid Modification (1-based)
                if ((mod = _modifications[i + 1]) != null && mod.MonoisotopicMass > 0)
                {
                    modSeqSb.Append('[');
                    modSeqSb.Append(mod);
                    modSeqSb.Append(']');
                }
            }

            // Handle C-Terminus Modification
            if ((mod = _modifications[Length + 1]) != null && mod.MonoisotopicMass > 0)
            {
                modSeqSb.Append("-[");
                modSeqSb.Append(mod);
                modSeqSb.Append(']');
            }

            return modSeqSb.ToString();
        }

        /// <summary>
        /// Gets the total number of amino acid residues in this amino acid polymer
        /// </summary>
        /// <returns>The number of amino acid residues</returns>
        public int ResidueCount()
        {
            return Length;
        }

        public int ResidueCount(Residue aminoAcid)
        {
            return aminoAcid == null ? 0 : residues.Count(aar => aar.Equals(aminoAcid));
        }

        /// <summary>
        /// Gets the number of amino acids residues in this amino acid polymer that
        /// has the specified residue letter
        /// </summary>
        /// <param name="residueLetter">The residue letter to search for</param>
        /// <returns>The number of amino acid residues that have the same letter in this polymer</returns>
        public int ResidueCount(char residueLetter)
        {
            return residues.Count(aar => aar.Letter.Equals(residueLetter));
        }

        public int ResidueCount(char residueLetter, int index, int length)
        {
            return residues.SubArray(index, length).Count(aar => aar.Letter.Equals(residueLetter));
        }

        public int ResidueCount(Residue aminoAcid, int index, int length)
        {
            return residues.SubArray(index, length).Count(aar => aar.Equals(aminoAcid));
        }

        public int ElementCountWithIsotopes(string element)
        {
            // Residues count
            int count = residues.Sum(aar => aar.ThisChemicalFormula.CountWithIsotopes(element));
            // Modifications count (if the mod is a IHasChemicalFormula)
            if (_modifications != null)
                count += _modifications.Where(mod => mod is IHasChemicalFormula).Cast<IHasChemicalFormula>().Sum(mod => mod.ThisChemicalFormula.CountWithIsotopes(element));

            count += ChemicalFormula.ParseFormula("H2O").CountWithIsotopes(element);
            return count;
        }

        public int SpecificIsotopeCount(Isotope isotope)
        {
            // Residues count
            int count = residues.Sum(aar => aar.ThisChemicalFormula.CountSpecificIsotopes(isotope));
            // Modifications count (if the mod is a IHasChemicalFormula)
            if (_modifications != null)
                count += _modifications.Where(mod => mod is IHasChemicalFormula).Cast<IHasChemicalFormula>().Sum(mod => mod.ThisChemicalFormula.CountSpecificIsotopes(isotope));
            return count;
        }

        /// <summary>
        /// Calculates the fragments that are different between this and another aminoacidpolymer
        /// </summary>
        /// <param name="other"></param>
        /// <param name="type"></param>
        /// <returns></returns>
        public IEnumerable<Fragment> GetSiteDeterminingFragments(AminoAcidPolymer other, FragmentTypes type)
        {
            return GetSiteDeterminingFragments(this, other, type);
        }

        public IEnumerable<Fragment> Fragment(FragmentTypes types)
        {
            return Fragment(types, false);
        }

        /// <summary>
        /// Calculates all the fragments of the types you specify
        /// </summary>
        /// <param name="types"></param>
        /// <param name="calculateChemicalFormula"></param>
        /// <returns></returns>
        public IEnumerable<Fragment> Fragment(FragmentTypes types, bool calculateChemicalFormula)
        {
            return Fragment(types, 1, Length - 1, calculateChemicalFormula);
        }

        public IEnumerable<Fragment> Fragment(FragmentTypes types, int number)
        {
            return Fragment(types, number, false);
        }

        public IEnumerable<Fragment> Fragment(FragmentTypes types, int number, bool calculateChemicalFormula)
        {
            return Fragment(types, number, number, calculateChemicalFormula);
        }

        public IEnumerable<Fragment> Fragment(FragmentTypes types, int minIndex, int maxIndex)
        {
            return Fragment(types, minIndex, maxIndex, false);
        }

        public IEnumerable<Fragment> Fragment(FragmentTypes types, int minIndex, int maxIndex, bool calculateChemicalFormula)
        {
            foreach (FragmentTypes type in types.GetIndividualFragmentTypes())
            {
                bool isChemicalFormula = calculateChemicalFormula;
                ChemicalFormula capFormula = type.GetIonCap();
                bool isCTerminal = type.GetTerminus() == Terminus.C;

                double monoMass = capFormula.MonoisotopicMass;
                ChemicalFormula formula = new ChemicalFormula(capFormula);

                IHasChemicalFormula terminus = isCTerminal ? CTerminus : NTerminus;
                monoMass += terminus.MonoisotopicMass;
                if (isChemicalFormula)
                    formula.Add(terminus);

                bool first = true;
                bool hasMod = _modifications != null;

                for (int i = 0; i <= maxIndex; i++)
                {
                    int aaIndex = isCTerminal ? Length - i : i - 1;

                    // Handle the terminus mods first in a special case
                    IHasMass mod;
                    if (first)
                    {
                        first = false;
                        if (hasMod)
                        {
                            mod = _modifications[aaIndex + 1];
                            if (mod != null)
                            {
                                monoMass += mod.MonoisotopicMass;
                                if (isChemicalFormula)
                                {
                                    if (mod is IHasChemicalFormula modFormula)
                                    {
                                        formula.Add(modFormula);
                                    }
                                    else
                                    {
                                        isChemicalFormula = false;
                                    }
                                }
                            }
                        }
                        continue;
                    }

                    monoMass += residues[aaIndex].MonoisotopicMass;
                    formula.Add(residues[aaIndex]);

                    if (hasMod)
                    {
                        mod = _modifications[aaIndex + 1];

                        if (mod != null)
                        {
                            monoMass += mod.MonoisotopicMass;
                            if (isChemicalFormula)
                            {
                                if (mod is IHasChemicalFormula modFormula)
                                {
                                    formula.Add(modFormula);
                                }
                                else
                                {
                                    isChemicalFormula = false;
                                }
                            }
                        }
                    }

                    if (i < minIndex)
                        continue;

                    if (isChemicalFormula)
                    {
                        yield return new ChemicalFormulaFragment(type, i, formula, this);
                    }
                    else
                    {
                        yield return new Fragment(type, i, monoMass, this);
                    }
                }
            }
        }

        public bool ContainsModifications()
        {
            return _modifications != null && _modifications.Any(m => m != null);
        }

        public ISet<T> GetUniqueModifications<T>() where T : IHasMass
        {
            HashSet<T> uniqueMods = new HashSet<T>();

            if (_modifications == null)
                return uniqueMods;

            foreach (IHasMass mod in _modifications)
            {
                if (mod is T)
                    uniqueMods.Add((T)mod);
            }
            return uniqueMods;
        }

        /// <summary>
        /// Counts the total number of modifications on this polymer that are not null
        /// </summary>
        /// <returns>The number of modifications</returns>
        public int ModificationCount()
        {
            return _modifications == null ? 0 : _modifications.Count(mod => mod != null);
        }

        /// <summary>
        /// Get the modification at the given residue number
        /// </summary>
        /// <param name="residueNumber">The amino acid residue number</param>
        /// <returns>The modification at the site, null if there isn't any modification present</returns>
        public IHasMass GetModification(int residueNumber)
        {
            return _modifications?[residueNumber];
        }

        /// <summary>
        /// Sets the modification at the terminus of this amino acid polymer
        /// </summary>
        /// <param name="modification">The modification to set</param>
        /// <param name="terminus">The termini to set the mod at</param>
        public virtual void SetModification(IHasMass modification, Terminus terminus)
        {
            if ((terminus & Terminus.N) == Terminus.N)
                NTerminusModification = modification;

            if ((terminus & Terminus.C) == Terminus.C)
                CTerminusModification = modification;
        }

        /// <summary>
        /// Sets the modification at specific sites on this amino acid polymer
        /// </summary>
        /// <param name="modification">The modification to set</param>
        /// <param name="sites">The sites to set the modification at</param>
        /// <returns>The number of modifications added to this amino acid polymer</returns>
        public virtual int SetModification(IHasMass modification, ModificationSites sites)
        {
            int count = 0;

            if ((sites & ModificationSites.NPep) == ModificationSites.NPep)
            {
                NTerminusModification = modification;
                count++;
            }

            for (int i = 0; i < Length; i++)
            {
                ModificationSites site = residues[i].Site;
                if ((sites & site) == site)
                {
                    ReplaceMod(i + 1, modification);
                    count++;
                }
            }

            if ((sites & ModificationSites.PepC) == ModificationSites.PepC)
            {
                CTerminusModification = modification;
                count++;
            }

            return count;
        }

        /// <summary>
        /// Sets the modification at specific sites on this amino acid polymer
        /// </summary>
        /// <param name="modification">The modification to set</param>
        /// <param name="letter">The residue character to set the modification at</param>
        /// <returns>The number of modifications added to this amino acid polymer</returns>
        public virtual int SetModification(IHasMass modification, char letter)
        {
            int count = 0;
            for (int i = 0; i < Length; i++)
            {
                if (!letter.Equals(residues[i].Letter))
                    continue;

                ReplaceMod(i + 1, modification);
                count++;
            }

            return count;
        }

        /// <summary>
        /// Sets the modification at specific sites on this amino acid polymer
        /// </summary>
        /// <param name="modification">The modification to set</param>
        /// <param name="residue">The residue to set the modification at</param>
        /// <returns>The number of modifications added to this amino acid polymer</returns>
        public virtual int SetModification(IHasMass modification, Residue residue)
        {
            int count = 0;
            for (int i = 0; i < Length; i++)
            {
                if (!residue.Letter.Equals(residues[i].Letter))
                    continue;
                ReplaceMod(i + 1, modification);
                count++;
            }
            return count;
        }

        /// <summary>
        /// Sets the modification at specific sites on this amino acid polymer
        /// </summary>
        /// <param name="modification">The modification to set</param>
        /// <param name="residueNumber">The residue number to set the modification at</param>
        public virtual void SetModification(IHasMass modification, int residueNumber)
        {
            if (residueNumber > Length || residueNumber < 1)
                throw new MzLibException(string.Format(CultureInfo.InvariantCulture, "Residue number not in the correct range: [{0}-{1}] you specified: {2}", 1, Length, residueNumber));

            ReplaceMod(residueNumber, modification);
        }

        public void SetModifications(IEnumerable<OldSchoolModification> modifications)
        {
            foreach (OldSchoolModification mod in modifications)
            {
                SetModification(mod, mod.Sites);
            }
        }

        public void SetModification(OldSchoolModification mod)
        {
            SetModification(mod, mod.Sites);
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="mod"></param>
        /// <param name="residueNumbers">(1-based) residue number</param>
        public void SetModification(IHasMass mod, params int[] residueNumbers)
        {
            foreach (int residueNumber in residueNumbers)
            {
                SetModification(mod, residueNumber);
            }
        }

        /// <summary>
        /// Replaces all instances of the old modification with the new modification in this polymer
        /// </summary>
        /// <param name="oldMod">The modification to remove</param>
        /// <param name="newMod">The modification to replace it with</param>
        /// <returns>The number of modifications added to this amino acid polymer</returns>
        public virtual int ReplaceModification(IHasMass oldMod, IHasMass newMod)
        {
            if (oldMod == null)
                throw new MzLibException("Cannot replace a null modification");

            int count = 0;
            for (int i = 0; i < Length + 2; i++)
            {
                IHasMass mod = GetModification(i);
                if (mod == null || !oldMod.Equals(mod))
                    continue;

                ReplaceMod(i, newMod);
                count++;
            }
            return count;
        }

        /// <summary>
        /// Adds the modification at the terminus of this amino acid polymer, combining modifications if a modification is already present
        /// </summary>
        /// <param name="modification">The modification to set</param>
        /// <param name="terminus">The termini to set the mod at</param>
        public virtual int AddModification(IHasMass modification, Terminus terminus)
        {
            IHasMass currentMod;
            int count = 0;

            if ((terminus & Terminus.N) == Terminus.N)
            {
                currentMod = NTerminusModification;
                NTerminusModification = currentMod == null ? modification : new ModificationCollection(currentMod, modification);
                count++;
            }

            if ((terminus & Terminus.C) == Terminus.C)
            {
                currentMod = CTerminusModification;
                CTerminusModification = currentMod == null ? modification : new ModificationCollection(currentMod, modification);
                count++;
            }
            return count;
        }

        public virtual int AddModification(OldSchoolModification modification)
        {
            return AddModification(modification, modification.Sites);
        }

        public virtual int AddModification(IHasMass modification, ModificationSites sites)
        {
            if (_modifications == null)
                _modifications = new IHasMass[Length + 2];

            int count = 0;
            IHasMass currentMod;
            if ((sites & ModificationSites.NPep) == ModificationSites.NPep)
            {
                currentMod = NTerminusModification;
                NTerminusModification = currentMod == null ? modification : new ModificationCollection(currentMod, modification);
                count++;
            }

            for (int i = 0; i < Length; i++)
            {
                ModificationSites site = residues[i].Site;
                if ((sites & site) == site)
                {
                    currentMod = _modifications[i + 1];
                    ReplaceMod(i + 1, currentMod == null ? modification : new ModificationCollection(currentMod, modification));
                    count++;
                }
            }

            if ((sites & ModificationSites.PepC) == ModificationSites.PepC)
            {
                currentMod = CTerminusModification;
                CTerminusModification = currentMod == null ? modification : new ModificationCollection(currentMod, modification);
                count++;
            }

            return count;
        }

        /// <summary>
        /// Adds the modification at specific sites on this amino acid polymer, combining modifications if a modification is already present
        /// </summary>
        /// <param name="modification">The modification to set</param>
        /// <param name="location">The location to set the modification at</param>
        public virtual void AddModification(IHasMass modification, int location)
        {
            IHasMass currentMod = GetModification(location);
            ReplaceMod(location, currentMod == null ? modification : new ModificationCollection(currentMod, modification));
        }

        /// <summary>
        /// Clears the modification set at the terminus of this amino acid polymer back
        /// to the default C or N modifications.
        /// </summary>
        /// <param name="terminus">The termini to clear the mod at</param>
        public void ClearModifications(Terminus terminus)
        {
            if ((terminus & Terminus.N) == Terminus.N)
                NTerminusModification = null;

            if ((terminus & Terminus.C) == Terminus.C)
                CTerminusModification = null;
        }

        /// <summary>
        /// Clear the modifications from the specified sites(s)
        /// </summary>
        /// <param name="sites">The sites to remove modifications from</param>
        public void ClearModifications(ModificationSites sites)
        {
            if ((sites & ModificationSites.NPep) == ModificationSites.NPep || (sites & ModificationSites.NProt) == ModificationSites.NProt)
            {
                ReplaceMod(0, null);
            }

            for (int i = 0; i < Length; i++)
            {
                int modIndex = i + 1;

                if (_modifications[modIndex] == null)
                    continue;

                ModificationSites curSite = residues[i].Site;

                if ((curSite & sites) == curSite)
                {
                    ReplaceMod(modIndex, null);
                }
            }

            if ((sites & ModificationSites.PepC) == ModificationSites.PepC || (sites & ModificationSites.ProtC) == ModificationSites.ProtC)
            {
                ReplaceMod(Length + 1, null);
            }
        }

        /// <summary>
        /// Clear all modifications from this amino acid polymer.
        /// Includes N and C terminus modifications.
        /// </summary>
        public void ClearModifications()
        {
            if (!ContainsModifications())
                return;

            for (int i = 0; i <= Length + 1; i++)
            {
                if (_modifications[i] == null)
                    continue;

                MonoisotopicMass -= _modifications[i].MonoisotopicMass;
                _modifications[i] = null;
            }
        }

        /// <summary>
        /// Removes the specified mod from all locations on this polymer
        /// </summary>
        /// <param name="mod">The modification to remove from this polymer</param>
        public void ClearModifications(IHasMass mod)
        {
            for (int i = 0; i <= Length + 1; i++)
            {
                if (!mod.Equals(_modifications[i]))
                    continue;

                MonoisotopicMass -= mod.MonoisotopicMass;
                _modifications[i] = null;
            }
        }

        /// <summary>
        /// Gets the chemical formula of this amino acid polymer.
        /// </summary>
        /// <returns></returns>
        public ChemicalFormula GetChemicalFormula()
        {
            var formula = new ChemicalFormula();

            // Handle Modifications
            if (ContainsModifications())
            {
                for (int i = 0; i < Length + 2; i++)
                {
                    if (_modifications[i] == null)
                        continue;

                    if (!(_modifications[i] is IHasChemicalFormula chemMod))
                        throw new MzLibException("Modification " + _modifications[i] + " does not have a chemical formula!");

                    formula.Add(chemMod.ThisChemicalFormula);
                }
            }

            // Handle N-Terminus
            formula.Add(NTerminus.ThisChemicalFormula);

            // Handle C-Terminus
            formula.Add(CTerminus.ThisChemicalFormula);

            // Handle Amino Acid Residues
            for (int i = 0; i < Length; i++)
            {
                formula.Add(residues[i].ThisChemicalFormula);
            }

            return formula;
        }

        public override string ToString()
        {
            return GetSequenceWithModifications();
        }

        public override int GetHashCode()
        {
            return BaseSequence.GetHashCode();
        }

        public override bool Equals(object obj)
        {
            AminoAcidPolymer aap = obj as AminoAcidPolymer;
            return aap != null && Equals(aap);
        }

        public bool Equals(AminoAcidPolymer other)
        {
            if (other == null ||
                Length != other.Length ||
                !NTerminus.ThisChemicalFormula.Equals(other.NTerminus.ThisChemicalFormula) ||
                !CTerminus.ThisChemicalFormula.Equals(other.CTerminus.ThisChemicalFormula))
                return false;

            bool containsMod = ContainsModifications();

            if (containsMod != other.ContainsModifications())
                return false;

            for (int i = 0; i <= Length + 1; i++)
            {
                if (containsMod && !Equals(_modifications[i], other._modifications[i]))
                {
                    return false;
                }

                if (i == 0 || i == Length + 1)
                {
                    continue; // uneven arrays, so skip these two conditions
                }

                if (!residues[i - 1].Equals(other.residues[i - 1]))
                {
                    return false;
                }
            }
            return true;
        }

        private void ReplaceTerminus(ref IHasChemicalFormula terminus, IHasChemicalFormula value)
        {
            if (terminus != null)
                MonoisotopicMass -= terminus.MonoisotopicMass;

            terminus = value;

            if (value != null)
                MonoisotopicMass += value.MonoisotopicMass;
        }

        /// <summary>
        /// Replaces a modification (if present) at the specific index in the residue (0-based for N and C termini)
        /// </summary>
        /// <param name="index">The residue index to replace at</param>
        /// <param name="mod">The modification to replace with</param>
        private void ReplaceMod(int index, IHasMass mod)
        {
            // No error checking here as all validation will occur before this method is call. This is to prevent
            // unneeded bounds checking

            if (_modifications == null)
            {
                _modifications = new IHasMass[Length + 2];
            }

            IHasMass oldMod = _modifications[index]; // Get the mod at the index, if present

            if (Equals(mod, oldMod))
                return; // Same modifications, no change is required

            if (oldMod != null)
                MonoisotopicMass -= oldMod.MonoisotopicMass; // remove the old mod mass

            _modifications[index] = mod;

            if (mod != null)
                MonoisotopicMass += mod.MonoisotopicMass; // add the new mod mass
        }

        /// <summary>
        /// Parses a string sequence of amino acids characters into a peptide object
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns></returns>
        private void ParseSequence(string sequence)
        {
            bool inMod = false;
            bool cterminalMod = false; // n or c terminal modification
            int index = 0;

            double monoMass = 0;

            StringBuilder modSb = new StringBuilder(10);
            foreach (char letter in sequence)
            {
                if (inMod)
                {
                    if (letter == ']')
                    {
                        inMod = false; // end the modification phase

                        string modString = modSb.ToString();
                        modSb.Clear();
                        IHasMass modification;
                        try
                        {
                            modification = new OldSchoolChemicalFormulaModification(ChemicalFormula.ParseFormula(modString));
                        }
                        catch (MzLibException)
                        {
                            if (double.TryParse(modString, out double mass))
                            {
                                modification = new ModWithOnlyMass(mass);
                            }
                            else
                            {
                                throw new MzLibException("Unable to correctly parse the following modification: " + modString);
                            }
                        }

                        monoMass += modification.MonoisotopicMass;

                        if (_modifications == null)
                            _modifications = new IHasMass[Length + 2];

                        if (cterminalMod)
                        {
                            _modifications[index + 1] = modification;
                        }
                        else
                        {
                            _modifications[index] = modification;
                        }

                        cterminalMod = false;
                    }
                    else
                    {
                        modSb.Append(letter);
                    }
                }
                else
                {
                    //char upperletter = char.ToUpper(letter); // moved to amino acid dictionary
                    if (Residue.TryGetResidue(letter, out Residue residue))
                    {
                        residues[index++] = residue;
                        monoMass += residue.MonoisotopicMass;
                    }
                    else
                    {
                        switch (letter)
                        {
                            case '[': // start of a modification
                                inMod = true;
                                break;

                            case '-': // End of an n-terminus mod or start of a c-terminus mod
                                cterminalMod = (index > 0);
                                break;

                            default:
                                throw new MzLibException(string.Format(CultureInfo.InvariantCulture, "Amino Acid Letter {0} does not exist in the Amino Acid Dictionary. {0} is also not a valid character", letter));
                        }
                    }
                }
            }

            if (inMod)
                throw new MzLibException("Couldn't find the closing ] for a modification in this sequence: " + sequence);

            Length = index;
            MonoisotopicMass += monoMass;
            Array.Resize(ref residues, Length);
            if (_modifications != null)
            {
                Array.Resize(ref _modifications, Length + 2);
            }
        }

        private class ModWithOnlyMass : IHasMass
        {
            private readonly double mass;

            public ModWithOnlyMass(double mass)
            {
                this.mass = mass;
            }

            public double MonoisotopicMass
            {
                get
                {
                    return mass;
                }
            }

            public override string ToString()
            {
                return mass.ToString(CultureInfo.InvariantCulture);
            }
        }
    }
}