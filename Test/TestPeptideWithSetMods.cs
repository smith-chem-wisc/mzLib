using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test
{
    [TestFixture]
    public static class TestPeptideWithSetMods
    {
        /// <summary>
        /// The purpose of this test is to ensure that two peptides digested from two different proteases are not equal even if their sequences are equal
        /// This is important for multiprotease parsimony in MetaMorpheus
        /// </summary>
        [Test]
        public static void TestDifferentProteaseEquals()
        {
            Protein myProtein = new Protein("SEQUENCEK", "accession");

            DigestionParams digest1 = new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            DigestionParams digest2 = new DigestionParams(protease: "Lys-C (cleave before proline)", maxMissedCleavages: 0, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            PeptideWithSetModifications pep1 = myProtein.Digest(digest1, new List<Modification>(), new List<Modification>()).First();
            PeptideWithSetModifications pep2 = myProtein.Digest(digest2, new List<Modification>(), new List<Modification>()).First();

            Assert.That(pep1.FullSequence.Equals(pep2.FullSequence));
            Assert.That(pep1.Protein.Equals(pep2.Protein));
            Assert.That(!pep1.DigestionParams.Protease.Equals(pep2.DigestionParams.Protease));
            Assert.That(!pep1.Equals(pep2));
            Assert.That(!pep1.GetHashCode().Equals(pep2.GetHashCode()));
        }

        [Test]
        public static void TestNonAndSemiSpecificDigests()
        {
            Protein fiveCleavages = new Protein("MAAKCCKDDKEEKFFKGG", "fiveCleavages");
            List<Tuple<string, FragmentationTerminus>> trypticSequencesInducingClevage = new List<Tuple<string, FragmentationTerminus>>
            {
                new Tuple<string, FragmentationTerminus>("K",FragmentationTerminus.C )
            };
            List<Tuple<string, FragmentationTerminus>> trypticSequencesPreventingClevage = new List<Tuple<string, FragmentationTerminus>>();

            Protease trypsinForTestNonAndSemiSpecificDigests = new Protease("trypsinForTestNonAndSemiSpecificDigests", trypticSequencesInducingClevage, trypticSequencesPreventingClevage, CleavageSpecificity.Full, "asdf", "asdf", "asdf");
            Protease semiTrypsinForTestNonAndSemiSpecificDigests = new Protease("semitrypsinForTestNonAndSemiSpecificDigests", trypticSequencesInducingClevage, trypticSequencesPreventingClevage, CleavageSpecificity.Semi, "asdf", "asdf", "asdf");

            ProteaseDictionary.Dictionary.Add(trypsinForTestNonAndSemiSpecificDigests.Name, trypsinForTestNonAndSemiSpecificDigests);
            ProteaseDictionary.Dictionary.Add(semiTrypsinForTestNonAndSemiSpecificDigests.Name, semiTrypsinForTestNonAndSemiSpecificDigests);

            DigestionParams fullyDigestParams = new DigestionParams(trypsinForTestNonAndSemiSpecificDigests.Name, 3, 2);
            var fiveCleavageProductsTrypsin = fiveCleavages.Digest(fullyDigestParams, null, null).ToList();

            Assert.AreEqual(22, fiveCleavageProductsTrypsin.Count);

            DigestionParams semiDigestionParams = new DigestionParams(semiTrypsinForTestNonAndSemiSpecificDigests.Name, 3, 2);
            var fiveCleavageProductsSemiTrypsin = fiveCleavages.Digest(semiDigestionParams, null, null).ToList();

            List<string> expectedProductsSemiFiveCleavages = new List<string>
            {
                "MAAKCCKDDKEEK",
                "AAKCCKDDKEEK",
                "CCKDDKEEKFFK",
                "DDKEEKFFKGG",
                "EEKFFKGG",
                "FFKGG",
                "AAK",
                "AAKCCK",
                "AAKCCKDDK",
            };
            //do cleavage cleave and retain
            foreach (string s in expectedProductsSemiFiveCleavages)
            {
                for (int i = 0; i < s.Length - semiDigestionParams.MinPeptideLength; i++)
                {
                    string sToFind = s.Substring(i);
                    var peps = fiveCleavageProductsSemiTrypsin.Where(x => x.BaseSequence.Equals(sToFind)).ToArray();
                    Assert.IsTrue(peps.Length == 1);
                    var pep = peps[0];
                    if ((pep.BaseSequence[0] == pep.BaseSequence[1] || pep.BaseSequence[0] == 'M') && (pep.BaseSequence[pep.BaseSequence.Length - 1] == 'K'
                        || (pep.BaseSequence[pep.BaseSequence.Length - 1] == 'G' && pep.BaseSequence[pep.BaseSequence.Length - 2] == 'G')))
                    {
                        Assert.IsTrue(pep.CleavageSpecificityForFdrCategory == CleavageSpecificity.Full);
                    }
                    else
                    {
                        Assert.IsTrue(pep.CleavageSpecificityForFdrCategory == CleavageSpecificity.Semi);
                    }
                    PeptideWithSetModifications pwsmRemake = new PeptideWithSetModifications(fiveCleavages, semiDigestionParams, pep.OneBasedStartResidueInProtein, pep.OneBasedEndResidueInProtein, CleavageSpecificity.Unknown, "", 3, pep.AllModsOneIsNterminus, 0);
                    Assert.IsTrue(pwsmRemake.CleavageSpecificityForFdrCategory == pep.CleavageSpecificityForFdrCategory);

                    sToFind = s.Substring(0, semiDigestionParams.MinPeptideLength + i);
                    peps = fiveCleavageProductsSemiTrypsin.Where(x => x.BaseSequence.Equals(sToFind)).ToArray();
                    Assert.IsTrue(peps.Length == 1);
                    pep = peps[0];
                    if ((pep.BaseSequence[0] == pep.BaseSequence[1] || pep.BaseSequence[0] == 'M') && pep.BaseSequence.Last() == 'K')
                    {
                        Assert.IsTrue(pep.CleavageSpecificityForFdrCategory == CleavageSpecificity.Full);
                    }
                    else
                    {
                        Assert.IsTrue(pep.CleavageSpecificityForFdrCategory == CleavageSpecificity.Semi);
                    }
                    pwsmRemake = new PeptideWithSetModifications(fiveCleavages, semiDigestionParams, pep.OneBasedStartResidueInProtein, pep.OneBasedEndResidueInProtein, CleavageSpecificity.Unknown, "", 3, pep.AllModsOneIsNterminus, 0);
                    Assert.IsTrue(pwsmRemake.CleavageSpecificityForFdrCategory == pep.CleavageSpecificityForFdrCategory);
                }
            }
            Assert.AreEqual(85, fiveCleavageProductsSemiTrypsin.Count);

            DigestionParams semiCleaveDigestionParams = new DigestionParams(semiTrypsinForTestNonAndSemiSpecificDigests.Name, 3, 2, initiatorMethionineBehavior: InitiatorMethionineBehavior.Cleave);
            var fiveCleavageProductsSemiTrypsinCleave = fiveCleavages.Digest(semiCleaveDigestionParams, null, null).ToList();
            int numVariableWithMet = fiveCleavageProductsSemiTrypsin.Where(x => x.BaseSequence[0] == 'M').Count();
            Assert.AreEqual(fiveCleavageProductsSemiTrypsin.Count, fiveCleavageProductsSemiTrypsinCleave.Count + numVariableWithMet);

            DigestionParams semiRetainDigestionParams = new DigestionParams(semiTrypsinForTestNonAndSemiSpecificDigests.Name, 3, 2, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);
            var fiveCleavageProductsSemiTrypsinRetain = fiveCleavages.Digest(semiRetainDigestionParams, null, null).ToList();
            int numNotRetained = fiveCleavageProductsSemiTrypsin.Where(x => x.BaseSequence[0] == 'A' && x.BaseSequence[1] == 'A'
            && (x.BaseSequence[x.BaseSequence.Length - 1] != 'K' && !(x.BaseSequence[x.BaseSequence.Length - 1] == 'G' && x.BaseSequence[x.BaseSequence.Length - 2] == 'G'))).Count();
            Assert.AreEqual(fiveCleavageProductsSemiTrypsinRetain.Count + numNotRetained, fiveCleavageProductsSemiTrypsin.Count);

            DigestionParams modernSemiDigestionParamsN = new DigestionParams(trypsinForTestNonAndSemiSpecificDigests.Name, 3, 2, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.N);
            var fiveCleavageProductsModernSemiTrypsinN = fiveCleavages.Digest(modernSemiDigestionParamsN, null, null).ToList();
            Assert.AreEqual(7, fiveCleavageProductsModernSemiTrypsinN.Count);

            DigestionParams modernSemiDigestionParamsC = new DigestionParams(trypsinForTestNonAndSemiSpecificDigests.Name, 3, 2, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.C);
            var fiveCleavageProductsModernSemiTrypsinC = fiveCleavages.Digest(modernSemiDigestionParamsC, null, null).ToList();
            Assert.AreEqual(6, fiveCleavageProductsModernSemiTrypsinC.Count);

            //Variable
            DigestionParams modernNonSpecificN = new DigestionParams("singleN", 50, 2, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.N);
            var fiveCleavageProductsModernNonSpecificN = fiveCleavages.Digest(modernNonSpecificN, null, null).ToList();
            Assert.AreEqual(17, fiveCleavageProductsModernNonSpecificN.Count);

            DigestionParams modernNonSpecificC = new DigestionParams("singleC", 50, 2, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.C);
            var fiveCleavageProductsModernNonSpecificC = fiveCleavages.Digest(modernNonSpecificC, null, null).ToList();
            Assert.AreEqual(17, fiveCleavageProductsModernNonSpecificC.Count);

            //maxRangeLimit
            modernNonSpecificN = new DigestionParams("singleN", 4, 2, 4, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.N);
            fiveCleavageProductsModernNonSpecificN = fiveCleavages.Digest(modernNonSpecificN, null, null).ToList();
            Assert.AreEqual(17, fiveCleavageProductsModernNonSpecificN.Count);
            foreach (var pep in fiveCleavageProductsModernNonSpecificN)
            {
                Assert.IsTrue(pep.BaseSequence.Length <= 4 && pep.BaseSequence.Length >= 2);
            }

            modernNonSpecificC = new DigestionParams("singleC", 4, 2, 4, searchModeType: CleavageSpecificity.None, fragmentationTerminus: FragmentationTerminus.C);
            fiveCleavageProductsModernNonSpecificC = fiveCleavages.Digest(modernNonSpecificC, null, null).ToList();
            Assert.AreEqual(17, fiveCleavageProductsModernNonSpecificC.Count);
            foreach (var pep in fiveCleavageProductsModernNonSpecificC)
            {
                Assert.IsTrue(pep.BaseSequence.Length <= 4 && pep.BaseSequence.Length >= 2);
            }

            //Cleave
            DigestionParams modernNonSpecificNCleave = new DigestionParams("singleN", 50, 2, searchModeType: CleavageSpecificity.None, initiatorMethionineBehavior: InitiatorMethionineBehavior.Cleave, fragmentationTerminus: FragmentationTerminus.N);
            fiveCleavageProductsModernNonSpecificN = fiveCleavages.Digest(modernNonSpecificNCleave, null, null).ToList();
            Assert.AreEqual(16, fiveCleavageProductsModernNonSpecificN.Count);

            DigestionParams modernNonSpecificCCleave = new DigestionParams("singleC", 50, 2, searchModeType: CleavageSpecificity.None, initiatorMethionineBehavior: InitiatorMethionineBehavior.Cleave, fragmentationTerminus: FragmentationTerminus.C);
            fiveCleavageProductsModernNonSpecificC = fiveCleavages.Digest(modernNonSpecificCCleave, null, null).ToList();
            Assert.AreEqual(16, fiveCleavageProductsModernNonSpecificC.Count);

            //Retain
            DigestionParams modernNonSpecificNRetain = new DigestionParams("singleN", 50, 2, searchModeType: CleavageSpecificity.None, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain, fragmentationTerminus: FragmentationTerminus.N);
            fiveCleavageProductsModernNonSpecificN = fiveCleavages.Digest(modernNonSpecificNRetain, null, null).ToList();
            Assert.AreEqual(17, fiveCleavageProductsModernNonSpecificN.Count);

            DigestionParams modernNonSpecificCRetain = new DigestionParams("singleC", 50, 2, searchModeType: CleavageSpecificity.None, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain, fragmentationTerminus: FragmentationTerminus.C);
            fiveCleavageProductsModernNonSpecificC = fiveCleavages.Digest(modernNonSpecificCRetain, null, null).ToList();
            Assert.AreEqual(17, fiveCleavageProductsModernNonSpecificC.Count);
        }

        [Test]
        public static void TestHardToParseModifiedSequence()
        {
            string fullSequence = "PE[Metal:Cation:Fe[III] on X]PTIDE";

            ModificationMotif.TryGetMotif("X", out var motif);

            Modification mod = new Modification(_originalId: "Cation:Fe[III]", _modificationType: "Metal",
                _monoisotopicMass: 1, _locationRestriction: "Anywhere.", _target: motif);

            Dictionary<string, Modification> mods = new Dictionary<string, Modification> { { "Cation:Fe[III] on X", mod } };

            PeptideWithSetModifications pep = new PeptideWithSetModifications(fullSequence, mods);

            Assert.That(pep.AllModsOneIsNterminus.Count == 1);
            var annotatedMod = pep.AllModsOneIsNterminus.First();
            Assert.That(annotatedMod.Key == 3);
            Assert.That(annotatedMod.Value.IdWithMotif == "Cation:Fe[III] on X");
            Assert.That(annotatedMod.Value.OriginalId == "Cation:Fe[III]");
            Assert.That(annotatedMod.Value.ModificationType == "Metal");

            fullSequence = "[Metal:Cation:Fe[III] on X]PE[Metal:Cation:Fe[III] on X]PTIDE[Metal:Cation:Fe[III] on X]";
            pep = new PeptideWithSetModifications(fullSequence, mods);
            Assert.That(pep.AllModsOneIsNterminus.Count == 3);
            Assert.That(pep.AllModsOneIsNterminus.Keys.ToList().SequenceEqual(new int[] { 1, 3, 8 }));
        }

        [Test]
        public static void TestCTermAndLastSideChainModParsing()
        {
            string fullSequence = "PEPTIDE[Mod:MyMod on E][PeptideCTermMod:MyCTermMod on E]";

            ModificationMotif.TryGetMotif("E", out var motif);

            Modification mod = new Modification(_originalId: "MyMod", _modificationType: "Mod",
                _monoisotopicMass: 1, _locationRestriction: "Anywhere.", _target: motif);

            Modification cTermMod = new Modification(_originalId: "MyCTermMod", _modificationType: "PeptideCTermMod",
                _monoisotopicMass: 1, _locationRestriction: "Peptide C-terminal.", _target: motif);

            Dictionary<string, Modification> mods = new Dictionary<string, Modification>
            {
                { "MyMod on E", mod },
                { "MyCTermMod on E", cTermMod }
            };

            PeptideWithSetModifications pep = new PeptideWithSetModifications(fullSequence, mods);

            Assert.That(pep.AllModsOneIsNterminus.Count == 2);
            Assert.That(pep.AllModsOneIsNterminus.Keys.SequenceEqual(new int[] { 8, 9 }));
        }
    }
}
