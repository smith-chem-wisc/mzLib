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
using System.Collections.Generic;

namespace Proteomics.AminoAcidPolymer
{
    public class Residue : IHasChemicalFormula
    {
        public static readonly double[] ResidueMonoisotopicMass;

        public static readonly Dictionary<string, Residue> ResiduesDictionary;
        private static readonly Residue[] ResiduesByLetter;

        static Residue()
        {
            ResiduesDictionary = new Dictionary<string, Residue>
            {
            {"Alanine",        new Residue("Alanine",       'A', "Ala",ChemicalFormula.ParseFormula("C3H5NO"),   ModificationSites.A)},
            {"Arginine",       new Residue("Arginine",      'R', "Arg",ChemicalFormula.ParseFormula("C6H12N4O"), ModificationSites.R)},
            {"Asparagine",     new Residue("Asparagine",    'N', "Asn",ChemicalFormula.ParseFormula("C4H6N2O2"), ModificationSites.N)},
            {"Aspartic Acid",  new Residue("Aspartic Acid", 'D', "Asp",ChemicalFormula.ParseFormula("C4H5NO3"),  ModificationSites.D)},
            {"Cysteine",       new Residue("Cysteine",      'C', "Cys",ChemicalFormula.ParseFormula("C3H5NOS"),  ModificationSites.C)},
            {"Glutamic Acid",  new Residue("Glutamic Acid", 'E', "Glu",ChemicalFormula.ParseFormula("C5H7NO3"),  ModificationSites.E)},
            {"Glutamine",      new Residue("Glutamine", 'Q', "Gln",ChemicalFormula.ParseFormula("C5H8N2O2"),  ModificationSites.Q)},
            {"Glycine",        new Residue("Glycine", 'G', "Gly",ChemicalFormula.ParseFormula("C2H3NO"),  ModificationSites.G)},
            {"Histidine",      new Residue("Histidine", 'H', "His",ChemicalFormula.ParseFormula("C6H7N3O"),  ModificationSites.H)},
            {"Isoleucine",     new Residue("Isoleucine", 'I', "Ile",ChemicalFormula.ParseFormula("C6H11NO"),  ModificationSites.I)},
            {"Leucine",        new Residue("Leucine", 'L', "Leu",ChemicalFormula.ParseFormula("C6H11NO"),  ModificationSites.L)},
            {"Lysine",         new Residue("Lysine", 'K', "Lys",ChemicalFormula.ParseFormula("C6H12N2O"),  ModificationSites.K)},
            {"Methionine",     new Residue("Methionine", 'M', "Met",ChemicalFormula.ParseFormula("C5H9NOS"),  ModificationSites.M)},
            {"Phenylalanine",  new Residue("Phenylalanine", 'F', "Phe",ChemicalFormula.ParseFormula("C9H9NO"),  ModificationSites.F)},
            {"Proline",        new Residue("Proline", 'P', "Pro",ChemicalFormula.ParseFormula("C5H7NO"),  ModificationSites.P)},
            {"Pyrrolysine",     new Residue("Pyrrolysine", 'O', "Pyl",ChemicalFormula.ParseFormula("C12H19N3O2"),  ModificationSites.P)},
            {"Selenocysteine", new Residue("Selenocysteine", 'U', "Sec",ChemicalFormula.ParseFormula("C3H5NOSe"),  ModificationSites.U)},
            {"Serine",         new Residue("Serine", 'S', "Ser",ChemicalFormula.ParseFormula("C3H5NO2"),  ModificationSites.S)},
            {"Threonine",      new Residue("Threonine", 'T', "Thr",ChemicalFormula.ParseFormula("C4H7NO2"),  ModificationSites.T)},
            {"Tryptophan",     new Residue("Tryptophan", 'W', "Trp",ChemicalFormula.ParseFormula("C11H10N2O"),  ModificationSites.W)},
            {"Tyrosine",       new Residue("Tyrosine", 'Y', "Try",ChemicalFormula.ParseFormula("C9H9NO2"),  ModificationSites.Y)},
            {"Valine",         new Residue("Valine", 'V', "Val",ChemicalFormula.ParseFormula("C5H9NO"),  ModificationSites.V)}
        };

            ResiduesByLetter = new Residue[]
        {
            null,null,null,null,null,null,null,null,null,null,null,null,null, //12
            null,null,null,null,null,null,null,null,null,null,null,null,null, //25
            null,null,null,null,null,null,null,null,null,null,null,null,null, //38
            null,null,null,null,null,null,null,null,null,null,null,null,null, //51
            null,null,null,null,null,null,null,null,null,null,null,null,null, //64
            ResiduesDictionary["Alanine"], //65
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
            ResiduesDictionary["Pyrrolysine"], // O
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
            null, // Z  //90
            null,null,null,null,null,null,null,null,null,null,null,null,null, //103
            null,null,null,null,null,null,null,null,null,null,null,null,null, //116
            null,null,null,null,null,null //122
        };
            ResidueMonoisotopicMass = new double[]
        {
            double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,
            double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,
            double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,
            double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,
            double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,
            ResiduesDictionary["Alanine"].MonoisotopicMass,
            double.NaN, // B
            ResiduesDictionary["Cysteine"].MonoisotopicMass,
            ResiduesDictionary["Aspartic Acid"].MonoisotopicMass,
            ResiduesDictionary["Glutamic Acid"].MonoisotopicMass,
            ResiduesDictionary["Phenylalanine"].MonoisotopicMass,
            ResiduesDictionary["Glycine"].MonoisotopicMass,
            ResiduesDictionary["Histidine"].MonoisotopicMass,
            ResiduesDictionary["Isoleucine"].MonoisotopicMass,
            ResiduesDictionary["Isoleucine"].MonoisotopicMass, // J - SPECIAL CASE!!!
            ResiduesDictionary["Lysine"].MonoisotopicMass,
            ResiduesDictionary["Leucine"].MonoisotopicMass,
            ResiduesDictionary["Methionine"].MonoisotopicMass,
            ResiduesDictionary["Asparagine"].MonoisotopicMass,
            ResiduesDictionary["Pyrrolysine"].MonoisotopicMass, // O
            ResiduesDictionary["Proline"].MonoisotopicMass,
            ResiduesDictionary["Glutamine"].MonoisotopicMass,
            ResiduesDictionary["Arginine"].MonoisotopicMass,
            ResiduesDictionary["Serine"].MonoisotopicMass,
            ResiduesDictionary["Threonine"].MonoisotopicMass,
            ResiduesDictionary["Selenocysteine"].MonoisotopicMass,
            ResiduesDictionary["Valine"].MonoisotopicMass,
            ResiduesDictionary["Tryptophan"].MonoisotopicMass,
            double.NaN, // X
            ResiduesDictionary["Tyrosine"].MonoisotopicMass,
            double.NaN, // Z
            double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,
            double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,
            double.NaN,double.NaN,double.NaN,double.NaN,double.NaN,double.NaN
        };
        }

        /// <summary>
        /// Adds a list of new residues to the dictionary at their specified index.
        /// </summary>
        /// <param name="residuesToAdd"></param>
        /// <returns></returns>
        public static void AddNewResiduesToDictionary(List<Residue> residuesToAdd)
        {
            foreach (Residue residue in residuesToAdd)
            {
                ResiduesDictionary[residue.Name] = residue;
                ResiduesByLetter[residue.Letter] = residue;
                ResidueMonoisotopicMass[residue.Letter] = residue.MonoisotopicMass;
            }
        }


        public Residue(string name, char oneLetterAbbreviation, string threeLetterAbbreviation, ChemicalFormula chemicalFormula, ModificationSites site)
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

        /// <summary>
        /// Get the residue based on the residues's symbol
        /// </summary>
        /// <param name="symbol"></param>
        /// <returns></returns>
        public static Residue GetResidue(string symbol)
        {
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
            if (letter < ResiduesByLetter.Length && letter >= 0)
            {
                residue = ResiduesByLetter[letter];
            }
            else
            {
                residue = null;
            }

            return residue != null;
        }

        public static bool TryGetResidue(string name, out Residue residue)
        {
            return ResiduesDictionary.TryGetValue(name, out residue);
        }
    }
}