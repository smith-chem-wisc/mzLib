// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (ModificationFragmentationTests.cs) is part of Proteomics.
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
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.Modifications.IO;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UsefulProteomicsDatabases;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class ModificationFragmentationTests
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public static void TestFragmentationNoMod()
        {
            // First we're checking to see if the fragment masses of an unmodified peptide a calculated correctly
            var prot = new Protein("PEPTIDE", null);
            DigestionParams digestionParams = new DigestionParams(

                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            List<Modification> variableModifications = new List<Modification>();
            var ye = prot.Digest(digestionParams, new List<Modification>(), variableModifications).ToList();

            // check unmodified
            var unmodPeptide = ye.Where(p => p.AllModsOneIsNterminus.Count == 0).First();
            var fragments = new List<Product>();
            unmodPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragments);
            var myUnmodFragmentMasses = fragments.Select(v => (int)Math.Round(v.NeutralMass.ToMz(1), 1)).ToList();
            HashSet<int> expectedMzs = new HashSet<int> { 98, 227, 324, 425, 538, 653, 703, 574, 477, 376, 263, 148 };

            Assert.That(expectedMzs.SetEquals(myUnmodFragmentMasses));
        }

        [Test]
        public static void TestFragmentationModNoNeutralLoss()
        {
            // Now we'll check the mass of modified peptide with no neutral losses
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            Modification mod = new Modification(_originalId: "oxidation", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("O1"), _locationRestriction: "Anywhere.");
            List<Modification> modlist = new List<Modification> { mod };
            DigestionParams digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var prot = new Protein("PEPTIDE", null, oneBasedModifications: new Dictionary<int, List<Modification>> { { 4, modlist } });
            var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

            // check unmodified
            var unmodPeptide = ye.Where(p => p.AllModsOneIsNterminus.Count == 0).First();
            var myUnmodFragments = new List<Product>();
            unmodPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, myUnmodFragments);
            var neutralMasses = new List<double>();
            neutralMasses.AddRange(myUnmodFragments.Select(m => m.NeutralMass).ToList());
            var expectedMasses = new List<double> { 97, 226, 323, 424, 537, 652, 147, 262, 375, 476, 573, 702 };
            for (int i = 0; i < neutralMasses.Count; i++)
            {
                neutralMasses[i] = Chemistry.ClassExtensions.RoundedDouble(neutralMasses[i], 0).Value;
            }

            var firstNotSecond = neutralMasses.Except(expectedMasses).ToList();
            var secondNotFirst = expectedMasses.Except(neutralMasses).ToList();

            //this is the set without oxidation
            Assert.AreEqual(12, myUnmodFragments.Count);
            Assert.AreEqual(0, firstNotSecond.Count);
            Assert.AreEqual(0, secondNotFirst.Count);

            // with oxidation, no neutral loss
            var modPeptide = ye.Where(p => p.AllModsOneIsNterminus.Count == 1).First();

            var myModFragments = new List<Product>();
            modPeptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, myModFragments);
            neutralMasses = new List<double>();
            neutralMasses.AddRange(myModFragments.Select(m => m.NeutralMass).ToList());
            expectedMasses = new List<double> { 97, 226, 323, 440, 553, 668, 147, 262, 375, 492, 589, 718 };
            for (int i = 0; i < neutralMasses.Count; i++)
            {
                neutralMasses[i] = Chemistry.ClassExtensions.RoundedDouble(neutralMasses[i], 0).Value;
            }

            firstNotSecond = neutralMasses.Except(expectedMasses).ToList();
            secondNotFirst = expectedMasses.Except(neutralMasses).ToList();

            //this is the set with oxidation
            Assert.AreEqual(12, myUnmodFragments.Count);
            Assert.AreEqual(0, firstNotSecond.Count);
            Assert.AreEqual(0, secondNotFirst.Count);
        }

        [Test]
        public static void Test_FragmentationModNeutralLoss()
        {
            // Now we'll check the mass of modified peptide with no neutral losses
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            Modification mod = new Modification(_originalId: "phospho", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: new Dictionary<DissociationType, List<double>> { { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } } }, _locationRestriction: "Anywhere.");
            List<Modification> modlist = new List<Modification> { mod };
            DigestionParams digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var prot = new Protein("PEPTIDE", null, oneBasedModifications: new Dictionary<int, List<Modification>> { { 4, modlist } });
            var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

            var peptideWithNeutralMassMod = ye.Where(v => v.AllModsOneIsNterminus.Count > 0).First();

            var myModFragments = new List<Product>();
            peptideWithNeutralMassMod.Fragment(DissociationType.HCD, FragmentationTerminus.Both, myModFragments);
            HashSet<int> neutralMasses = new HashSet<int>(myModFragments.Select(m => (int)m.NeutralMass.ToMz(1)).ToList());
            HashSet<int> expectedMasses = new HashSet<int> { 98,227, 324, 407, 520, 635, 505, 618, 733, //b-ions with and without neutral loss
            148, 263, 376, 459, 556, 685, 557, 654, 783, //y-ions with and without neutral loss
            782}; //molecular ion with neutral loss

            CollectionAssert.AreEquivalent(neutralMasses, expectedMasses);
        }

        [Test]
        public static void Test_FragmentationTwoModNeutralLoss()
        {
            // Now we'll check the mass of modified peptide with 2 neutral loss mods
            ModificationMotif.TryGetMotif("Q", out ModificationMotif motifone);
            Modification modone = new Modification(_originalId: "ammonia", _modificationType: "testModType", _target: motifone, _monoisotopicMass: 0, _neutralLosses: new Dictionary<DissociationType,
                List<double>> { { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 N1").MonoisotopicMass } } }, _locationRestriction: "Anywhere.");

            ModificationMotif.TryGetMotif("T", out ModificationMotif motiftwo);
            Modification modtwo = new Modification(_originalId: "phospho", _modificationType: "testModType", _target: motiftwo, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: new Dictionary<DissociationType,
                List<double>> { { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } } }, _locationRestriction: "Anywhere.");

            List<Modification> modlistone = new List<Modification> { modone };
            List<Modification> modlisttwo = new List<Modification> { modtwo };

            DigestionParams digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var prot = new Protein("PEQTIDE", null, oneBasedModifications: new Dictionary<int, List<Modification>> { { 3, modlistone }, { 4, modlisttwo } });
            var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

            var peptideWithNeutralMassMod = ye.Where(v => v.AllModsOneIsNterminus.Count == 2).First();

            var myModFragments = new List<Product>();
            peptideWithNeutralMassMod.Fragment(DissociationType.HCD, FragmentationTerminus.Both, myModFragments);
            HashSet<int> neutralMasses = new HashSet<int>(myModFragments.Select(m => (int)m.NeutralMass.ToMz(1)).ToList());
            HashSet<int> expectedMasses = new HashSet<int> { 98, 227, 355, 536, 649, 764, 438, 551, 666, 338, 519, 632, 747, // b-ions with and without neutral losses
                                                             148, 263, 376, 557, 685, 814, 668, 797, 459, 587, 716, //y ions with and without neutral losses
                                                               813, 894, }; //molecular ion with neutral losses (phospho and ammonia respectively)

            CollectionAssert.AreEquivalent(neutralMasses, expectedMasses);
        }

        [Test]
        public static void Test_FragmentationTwoModNeutralLossTwoFragTypes()
        {
            // Now we'll check the mass of modified peptide with no neutral losses
            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);

            Dictionary<DissociationType, List<double>> myNeutralLosses = new Dictionary<DissociationType, List<double>>()
            {
                { DissociationType.HCD, new List<double> { ChemicalFormula.ParseFormula("H3 O4 P1").MonoisotopicMass } },
                { DissociationType.ETD, new List<double>() { ChemicalFormula.ParseFormula("H3 N1").MonoisotopicMass } } // this makes no sense in real life, it's just for a unit test
            };

            Modification mod = new Modification(_originalId: "phospho", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("H1 O3 P1"), _neutralLosses: myNeutralLosses, _locationRestriction: "Anywhere.");
            List<Modification> modlist = new List<Modification> { mod };
            DigestionParams digestionParams = new DigestionParams(
                protease: "trypsin",
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            var prot = new Protein("PEPTIDE", null, oneBasedModifications: new Dictionary<int, List<Modification>> { { 4, modlist } });
            var ye = prot.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

            var peptideWithNeutralMassMod = ye.Where(v => v.AllModsOneIsNterminus.Count == 1).First();

            var myModFragmentsHCD = new List<Product>();
            peptideWithNeutralMassMod.Fragment(DissociationType.HCD, FragmentationTerminus.Both, myModFragmentsHCD);

            var neutralMassesHCD = myModFragmentsHCD.Select(m => (int)m.NeutralMass.ToMz(1));
            var expectedMassesHCD = new HashSet<int> { 98, 227, 324, 407, 520, 635, 505, 618, 733,// b-ions with and without neutral loss
                                                         148, 263, 376, 459, 556, 685, 557, 654, 783,//y-ions with and without neutral loss
                                                         782};// molecular ion with neutral loss

            CollectionAssert.AreEquivalent(expectedMassesHCD, neutralMassesHCD);

            //Now try the other half
            var myModFragmentsETD = new List<Product>();
            peptideWithNeutralMassMod.Fragment(DissociationType.ETD, FragmentationTerminus.Both, myModFragmentsETD);

            var neutralMassesETD = myModFragmentsETD.Select(m => (int)m.NeutralMass.ToMz(1));
            var expectedMassesETD = new HashSet<int> { 115, 244, 341, 505, 618, 733, 522, 635, 750,  // c-ions and c-17 ions
            148, 263, 376, 540, 637, 766, 557, 654, 783,          // y and y-17 ions
            133, 248, 361, 525, 622, 751, 542, 639, 768,         // z+1 and z+1-17 ions
            863 };//Molecular ions minus ammonia

            CollectionAssert.AreEquivalent(expectedMassesHCD, neutralMassesHCD);
        }

        [Test]
        public static void TestFragmentNterminalModifiedPeptide()
        {
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            Modification nTermMod = new Modification(_originalId: "acetylation", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "N-terminal.");

            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>> { { 1, new List<Modification> { nTermMod } } };

            Protein protein = new Protein("PEPTIDE", "", oneBasedModifications: mods);
            PeptideWithSetModifications peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).Where(p => p.AllModsOneIsNterminus.Count == 1).First();
            Assert.That(peptide.FullSequence == "[testModType:acetylation on P]PEPTIDE");

            var fragments = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragments);

            var roundedFragments = fragments.Select(f => (int)f.NeutralMass).ToList();
            CollectionAssert.AreEquivalent(roundedFragments, new int[] { 139, 268, 365, 466, 579, 694, 147, 262, 375, 476, 573, 702 });
        }

        [Test]
        public static void TestFragmentCTerminalModifiedPeptide()
        {
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            Modification cTermMod = new Modification(_originalId: "acetylation", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "C-terminal.");

            Dictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>> { { 7, new List<Modification> { cTermMod } } };

            Protein protein = new Protein("PEPTIDE", "", oneBasedModifications: mods);
            PeptideWithSetModifications peptide = protein.Digest(new DigestionParams(), new List<Modification>(), new List<Modification>()).Where(p => p.AllModsOneIsNterminus.Count == 1).First();
            Assert.That(peptide.FullSequence == "PEPTIDE-[testModType:acetylation on E]");

            var fragments = new List<Product>();
            peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragments);

            var roundedFragments = fragments.Select(f => (int)f.NeutralMass).ToList();
            CollectionAssert.AreEquivalent(roundedFragments, new int[] { 97, 226, 323, 424, 537, 652, 189, 304, 417, 518, 615, 744 });
        }
    }
}
