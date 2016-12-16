// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
// 
// This file (AminoAcid.cs) is part of Proteomics.
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
using System;
using System.Collections.Generic;
using System.Globalization;

namespace Proteomics
{
    public class AminoAcid : IHasChemicalFormula
    {
        private static readonly Dictionary<string, AminoAcid> Residues = new Dictionary<string, AminoAcid>
        {
            {"Alanine",        new AminoAcid("Alanine",       'A', "Ala","C3H5NO",   ModificationSites.A)},
            {"Arginine",       new AminoAcid("Arginine",      'R', "Arg","C6H12N4O", ModificationSites.R)},
            {"Asparagine",     new AminoAcid("Asparagine",    'N', "Asn","C4H6N2O2", ModificationSites.N)},
            {"Aspartic Acid",  new AminoAcid("Aspartic Acid", 'D', "Asp","C4H5NO3",  ModificationSites.D)},
            {"Cysteine",       new AminoAcid("Cysteine",      'C', "Cys","C3H5NOS",  ModificationSites.C)},
            {"Glutamic Acid",  new AminoAcid("Glutamic Acid", 'E', "Glu","C5H7NO3",  ModificationSites.E)},
            {"Glutamine",      new AminoAcid("Glutamine", 'Q', "Gln","C5H8N2O2",  ModificationSites.Q)},
            {"Glycine",        new AminoAcid("Glycine", 'G', "Gly","C2H3NO",  ModificationSites.G)},
            {"Histidine",      new AminoAcid("Histidine", 'H', "His","C6H7N3O",  ModificationSites.H)},
            {"Isoleucine",     new AminoAcid("Isoleucine", 'I', "Ile","C6H11NO",  ModificationSites.I)},
            {"Leucine",        new AminoAcid("Leucine", 'L', "Leu","C6H11NO",  ModificationSites.L)},
            {"Lysine",         new AminoAcid("Lysine", 'K', "Lys","C6H12N2O",  ModificationSites.K)},
            {"Methionine",     new AminoAcid("Methionine", 'M', "Met","C5H9NOS",  ModificationSites.M)},
            {"Phenylalanine",  new AminoAcid("Phenylalanine", 'F', "Phe","C9H9NO",  ModificationSites.F)},
            {"Proline",        new AminoAcid("Proline", 'P', "Pro","C5H7NO",  ModificationSites.P)},
            {"Selenocysteine", new AminoAcid("Selenocysteine", 'U', "Sec","C3H5NOSe",  ModificationSites.U)},
            {"Serine",         new AminoAcid("Serine", 'S', "Ser","C3H5NO2",  ModificationSites.S)},
            {"Threonine",      new AminoAcid("Threonine", 'T', "Thr","C4H7NO2",  ModificationSites.T)},
            {"Tryptophan",     new AminoAcid("Tryptophan", 'W', "Trp","C11H10N2O",  ModificationSites.W)},
            {"Tyrosine",       new AminoAcid("Tyrosine", 'Y', "Try","C9H9NO2",  ModificationSites.Y)},
            {"Valine",         new AminoAcid("Valine", 'V', "Val","C5H9NO",  ModificationSites.V)},
        };

        private static readonly AminoAcid[] ResiduesByLetter = new AminoAcid['z' + 1]
        {
            null,null,null,null,null,null,null,null,null,null,null,null,null,
            null,null,null,null,null,null,null,null,null,null,null,null,null,
            null,null,null,null,null,null,null,null,null,null,null,null,null,
            null,null,null,null,null,null,null,null,null,null,null,null,null,
            null,null,null,null,null,null,null,null,null,null,null,null,null,
            Residues["Alanine"],
            null, // B
            Residues["Cysteine"],
            Residues["Aspartic Acid"],
            Residues["Glutamic Acid"],
            Residues["Phenylalanine"],
            Residues["Glycine"],
            Residues["Histidine"],
            Residues["Isoleucine"],
            null, // J
            Residues["Lysine"],
            Residues["Leucine"],
            Residues["Methionine"],
            Residues["Asparagine"],
            null, // O
            Residues["Proline"],
            Residues["Glutamine"],
            Residues["Arginine"],
            Residues["Serine"],
            Residues["Threonine"],
            Residues["Selenocysteine"],
            Residues["Valine"],
            Residues["Tryptophan"],
            null, // X
            Residues["Tyrosine"],
            null, // Z
            null,null,null,null,null,null,null,null,null,null,null,null,null,
            null,null,null,null,null,null,null,null,null,null,null,null,null,
            null,null,null,null,null,null,
        };

        /// <summary>
        /// Get the residue based on the residues's symbol
        /// </summary>
        /// <param name="symbol"></param>
        /// <returns></returns>
        public static AminoAcid GetResidue(string symbol)
        {
            if (symbol == null)
                throw new ArgumentNullException("symbol", "Cannot get residue of null parameter");
            return symbol.Length == 1 ? ResiduesByLetter[symbol[0]] : Residues[symbol];
        }

        /// <summary>
        /// Gets the resdiue based on the residue's one-character symbol
        /// </summary>
        /// <param name="letter"></param>
        /// <returns></returns>
        public static AminoAcid GetResidue(char letter)
        {
            return ResiduesByLetter[letter];
        }

        public static bool TryGetResidue(char letter, out AminoAcid residue)
        {
            residue = ResiduesByLetter[letter];
            return residue != null;
        }

        public static bool TryGetResidue(string name, out AminoAcid residue)
        {
            return Residues.TryGetValue(name, out residue);
        }

        internal AminoAcid(string name, char oneLetterAbbreviation, string threeLetterAbbreviation, ChemicalFormula chemicalFormula, ModificationSites site)
        {
            Name = name;
            Letter = oneLetterAbbreviation;
            Symbol = threeLetterAbbreviation;
            ThisChemicalFormula = chemicalFormula;
            MonoisotopicMass = ThisChemicalFormula.MonoisotopicMass;
            Site = site;
        }

        public ChemicalFormula ThisChemicalFormula { get; private set; }

        public char Letter { get; private set; }

        public ModificationSites Site { get; private set; }

        public double MonoisotopicMass { get; private set; }

        public string Name { get; private set; }

        public string Symbol { get; private set; }

        public override string ToString()
        {
            return string.Format(CultureInfo.InvariantCulture, "{0} {1} ({2})", Letter, Symbol, Name);
        }
    }
}