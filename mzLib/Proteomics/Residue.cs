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
    public class Residue : IHasChemicalFormula
    {
        private static readonly Dictionary<string, Residue> ResiduesDictionary = new Dictionary<string, Residue>
        {
            {"Alanine",        new Residue("Alanine",       'A', "Ala","C3H5NO",   ModificationSites.A)},
            {"Arginine",       new Residue("Arginine",      'R', "Arg","C6H12N4O", ModificationSites.R)},
            {"Asparagine",     new Residue("Asparagine",    'N', "Asn","C4H6N2O2", ModificationSites.N)},
            {"Aspartic Acid",  new Residue("Aspartic Acid", 'D', "Asp","C4H5NO3",  ModificationSites.D)},
            {"Cysteine",       new Residue("Cysteine",      'C', "Cys","C3H5NOS",  ModificationSites.C)},
            {"Glutamic Acid",  new Residue("Glutamic Acid", 'E', "Glu","C5H7NO3",  ModificationSites.E)},
            {"Glutamine",      new Residue("Glutamine", 'Q', "Gln","C5H8N2O2",  ModificationSites.Q)},
            {"Glycine",        new Residue("Glycine", 'G', "Gly","C2H3NO",  ModificationSites.G)},
            {"Histidine",      new Residue("Histidine", 'H', "His","C6H7N3O",  ModificationSites.H)},
            {"Isoleucine",     new Residue("Isoleucine", 'I', "Ile","C6H11NO",  ModificationSites.I)},
            {"Leucine",        new Residue("Leucine", 'L', "Leu","C6H11NO",  ModificationSites.L)},
            {"Lysine",         new Residue("Lysine", 'K', "Lys","C6H12N2O",  ModificationSites.K)},
            {"Methionine",     new Residue("Methionine", 'M', "Met","C5H9NOS",  ModificationSites.M)},
            {"Phenylalanine",  new Residue("Phenylalanine", 'F', "Phe","C9H9NO",  ModificationSites.F)},
            {"Proline",        new Residue("Proline", 'P', "Pro","C5H7NO",  ModificationSites.P)},
            {"Selenocysteine", new Residue("Selenocysteine", 'U', "Sec","C3H5NOSe",  ModificationSites.U)},
            {"Serine",         new Residue("Serine", 'S', "Ser","C3H5NO2",  ModificationSites.S)},
            {"Threonine",      new Residue("Threonine", 'T', "Thr","C4H7NO2",  ModificationSites.T)},
            {"Tryptophan",     new Residue("Tryptophan", 'W', "Trp","C11H10N2O",  ModificationSites.W)},
            {"Tyrosine",       new Residue("Tyrosine", 'Y', "Try","C9H9NO2",  ModificationSites.Y)},
            {"Valine",         new Residue("Valine", 'V', "Val","C5H9NO",  ModificationSites.V)},
        };

        private static readonly Residue[] ResiduesByLetter = new Residue['z' + 1]
        {
            null,null,null,null,null,null,null,null,null,null,null,null,null,
            null,null,null,null,null,null,null,null,null,null,null,null,null,
            null,null,null,null,null,null,null,null,null,null,null,null,null,
            null,null,null,null,null,null,null,null,null,null,null,null,null,
            null,null,null,null,null,null,null,null,null,null,null,null,null,
            ResiduesDictionary["Alanine"],
            null, // B
            ResiduesDictionary["Cysteine"],
            ResiduesDictionary["Aspartic Acid"],
            ResiduesDictionary["Glutamic Acid"],
            ResiduesDictionary["Phenylalanine"],
            ResiduesDictionary["Glycine"],
            ResiduesDictionary["Histidine"],
            ResiduesDictionary["Isoleucine"],
            null, // J
            ResiduesDictionary["Lysine"],
            ResiduesDictionary["Leucine"],
            ResiduesDictionary["Methionine"],
            ResiduesDictionary["Asparagine"],
            null, // O
            ResiduesDictionary["Proline"],
            ResiduesDictionary["Glutamine"],
            ResiduesDictionary["Arginine"],
            ResiduesDictionary["Serine"],
            ResiduesDictionary["Threonine"],
            ResiduesDictionary["Selenocysteine"],
            ResiduesDictionary["Valine"],
            ResiduesDictionary["Tryptophan"],
            null, // X
            ResiduesDictionary["Tyrosine"],
            null, // Z
            null,null,null,null,null,null,null,null,null,null,null,null,null,
            null,null,null,null,null,null,null,null,null,null,null,null,null,
            null,null,null,null,null,null,
        };

        public static readonly double[] ResidueMonoisotopicMass = new double['z' + 1]
        {
            0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,
            ResiduesDictionary["Alanine"].MonoisotopicMass,
            0, // B
            ResiduesDictionary["Cysteine"].MonoisotopicMass,
            ResiduesDictionary["Aspartic Acid"].MonoisotopicMass,
            ResiduesDictionary["Glutamic Acid"].MonoisotopicMass,
            ResiduesDictionary["Phenylalanine"].MonoisotopicMass,
            ResiduesDictionary["Glycine"].MonoisotopicMass,
            ResiduesDictionary["Histidine"].MonoisotopicMass,
            ResiduesDictionary["Isoleucine"].MonoisotopicMass,
            0, // J
            ResiduesDictionary["Lysine"].MonoisotopicMass,
            ResiduesDictionary["Leucine"].MonoisotopicMass,
            ResiduesDictionary["Methionine"].MonoisotopicMass,
            ResiduesDictionary["Asparagine"].MonoisotopicMass,
            0, // O
            ResiduesDictionary["Proline"].MonoisotopicMass,
            ResiduesDictionary["Glutamine"].MonoisotopicMass,
            ResiduesDictionary["Arginine"].MonoisotopicMass,
            ResiduesDictionary["Serine"].MonoisotopicMass,
            ResiduesDictionary["Threonine"].MonoisotopicMass,
            ResiduesDictionary["Selenocysteine"].MonoisotopicMass,
            ResiduesDictionary["Valine"].MonoisotopicMass,
            ResiduesDictionary["Tryptophan"].MonoisotopicMass,
            0, // X
            ResiduesDictionary["Tyrosine"].MonoisotopicMass,
            0, // Z
            0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,
        };

        /// <summary>
        /// Get the residue based on the residues's symbol
        /// </summary>
        /// <param name="symbol"></param>
        /// <returns></returns>
        public static Residue GetResidue(string symbol)
        {
            if (symbol == null)
                throw new ArgumentNullException("symbol", "Cannot get residue of null parameter");
            return symbol.Length == 1 ? ResiduesByLetter[symbol[0]] : ResiduesDictionary[symbol];
        }

        /// <summary>
        /// Gets the resdiue based on the residue's one-character symbol
        /// </summary>
        /// <param name="letter"></param>
        /// <returns></returns>
        public static Residue GetResidue(char letter)
        {
            return ResiduesByLetter[letter];
        }

        public static bool TryGetResidue(char letter, out Residue residue)
        {
            residue = ResiduesByLetter[letter];
            return residue != null;
        }

        public static bool TryGetResidue(string name, out Residue residue)
        {
            return ResiduesDictionary.TryGetValue(name, out residue);
        }

        internal Residue(string name, char oneLetterAbbreviation, string threeLetterAbbreviation, ChemicalFormula chemicalFormula, ModificationSites site)
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