using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using MzLibUtil;
using Readers;
using System.Collections.Generic;
using FlashLFQ;
using System.Linq;
using Proteomics.AminoAcidPolymer;
using System;
using NUnit.Framework.Legacy;
using MzLibUtil.PositionFrequencyAnalysis;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestMzLibUtil
    {
        [Test]
        [TestCase(@"C:\Users\bubba\Documents\Projects\K562\K562_2\20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase("C:\\Users\bubba\\Documents\\Projects\\K562\\K562_2\\20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase("20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase("20100730_Velos1_TaGe_SA_K565_4", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase(@"C:\Users\bubba\Documents\Projects\K562\K562_2\20100730_Velos1_TaGe_SA_K565_4", "20100730_Velos1_TaGe_SA_K565_4")]
        //test extra period in folder name of path
        [TestCase(@"C:\Users\bubba\Documents.docs\Projects\K562\K562_2\20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        //test extra period in filename
        [TestCase(@"C:\Users\bubba\Documents\Projects\K562\K562_2\20100730_Velos1_.TaGe_SA_K565_4.raw", "20100730_Velos1_.TaGe_SA_K565_4")]
        [TestCase("/home/seth/Pictures/penguin.jpg","penguin")]
        [TestCase("/home/seth/Pictures/penguin", "penguin")]
        [TestCase("penguin.jpg", "penguin")]
        [TestCase("penguin", "penguin")]
        [TestCase("penguin.jpg.gz", "penguin")]
        [TestCase("penguin.jpg.zip", "penguin")]
        [TestCase("penguin.jpg.mzXML", "penguin.jpg")]
        public static void TestPeriodTolerantFilenameWithoutExtension(string filenameAndOrPath, string expectedResult)
        {
            string result = PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(filenameAndOrPath);
            string extensionResult = filenameAndOrPath.GetPeriodTolerantFilenameWithoutExtension();
            Assert.AreEqual(expectedResult, result);
            Assert.AreEqual(expectedResult, extensionResult);
        }

        [Test]
        public static void TestToEnum()
        {
            Assert.IsTrue(0.ToEnum<TimsTofMsMsType>(out var result));
            Assert.AreEqual(TimsTofMsMsType.MS, result);

            Assert.IsTrue(2.ToEnum<TimsTofMsMsType>(out result));
            Assert.AreEqual(TimsTofMsMsType.MSMSFragment, result);

            Assert.IsTrue(8.ToEnum<TimsTofMsMsType>(out result));
            Assert.AreEqual(TimsTofMsMsType.PASEF, result);

            Assert.IsTrue(9.ToEnum<TimsTofMsMsType>(out result));
            Assert.AreEqual(TimsTofMsMsType.DIA, result);

            Assert.IsTrue(10.ToEnum<TimsTofMsMsType>(out result));
            Assert.AreEqual(TimsTofMsMsType.PRM, result);

            Assert.IsTrue(0.ToEnum<TimsTofAcquisitionMode>(out var result2));
            Assert.AreEqual(TimsTofAcquisitionMode.MS, result2);

            Assert.IsFalse(1.ToEnum<TimsTofMsMsType>(out result));
            Assert.IsFalse(11.ToEnum<TimsTofMsMsType>(out result));
            Assert.IsFalse(7.ToEnum<TimsTofMsMsType>(out result));
            
        }

        [Test]
        public void IsNullOrEmpty_IEnumerable_Null_ReturnsTrue()
        {
            IEnumerable<int> nullEnumerable = null;
            Assert.IsTrue(nullEnumerable.IsNullOrEmpty());
        }

        [Test]
        public void IsNullOrEmpty_IEnumerable_Empty_ReturnsTrue()
        {
            IEnumerable<int> emptyEnumerable = new List<int>();
            Assert.IsTrue(emptyEnumerable.IsNullOrEmpty());
        }

        [Test]
        public void IsNullOrEmpty_IEnumerable_NotEmpty_ReturnsFalse()
        {
            IEnumerable<int> notEmptyEnumerable = new List<int> { 1 };
            Assert.IsFalse(notEmptyEnumerable.IsNullOrEmpty());
        }

        [Test]
        public void IsNullOrEmpty_IDictionary_Null_ReturnsTrue()
        {
            IDictionary<int, int> nullDictionary = null;
            Assert.IsTrue(nullDictionary.IsNullOrEmpty());
        }

        [Test]
        public void IsNullOrEmpty_IDictionary_Empty_ReturnsTrue()
        {
            IDictionary<int, int> emptyDictionary = new Dictionary<int, int>();
            Assert.IsTrue(emptyDictionary.IsNullOrEmpty());
        }

        [Test]
        public void IsNullOrEmpty_IDictionary_NotEmpty_ReturnsFalse()
        {
            IDictionary<int, int> notEmptyDictionary = new Dictionary<int, int> { { 1, 1 } };
            Assert.IsFalse(notEmptyDictionary.IsNullOrEmpty());
        }

        [Test]
        public void IsDefaultOrNull_WithNullReferenceType_ReturnsTrue()
        {
            string value = null;
            bool result = value.IsDefaultOrNull();
            Assert.IsTrue(result);
        }

        [Test]
        public void IsDefaultOrNull_WithNonNullReferenceType_ReturnsFalse()
        {
            string value = "test";
            bool result = value.IsDefaultOrNull();
            Assert.IsFalse(result);
        }

        [Test]
        public void IsDefaultOrNull_WithDefaultStruct_ReturnsTrue()
        {
            TestStruct value = default(TestStruct);
            bool result = value.IsDefaultOrNull();
            Assert.IsTrue(result);
        }

        [Test]
        public void IsDefaultOrNull_WithNonDefaultStruct_ReturnsFalse()
        {
            TestStruct value = new TestStruct { X = 1, Y = 2 };
            bool result = value.IsDefaultOrNull();
            Assert.IsFalse(result);
        }

        [Test]
        public void IsNotDefaultOrNull_WithNullReferenceType_ReturnsFalse()
        {
            string value = null;
            bool result = value.IsNotDefaultOrNull();
            Assert.IsFalse(result);
        }

        [Test]
        public void TestRemoveSpecialCharacters()
        {
            // Test default pipe removal
            string seqWithPipes = "PE|PTI|DE";
            string seqNoPipes = seqWithPipes.ToString();
            ClassExtensions.RemoveSpecialCharacters(ref seqNoPipes);
            Assert.AreEqual("PEPTIDE", seqNoPipes);


            // Test specified character replacement
            string seqWithHash = seqWithPipes.ToString();
            ClassExtensions.RemoveSpecialCharacters(ref seqWithHash, replacement: "#", specialCharacter: @"\|");
            Assert.AreEqual("PE#PTI#DE", seqWithHash);

            // Test specified character removal
            string cleanSeq = seqWithHash.ToString();
            ClassExtensions.RemoveSpecialCharacters(ref cleanSeq, specialCharacter: "#");
            Assert.AreEqual("PEPTIDE", cleanSeq);
        }

        [Test]
        public void TestQuantifiedModification()
        {
            var quantmod = new QuantifiedModification(idWithMotif: "TestMod: ModX on AAY", positionInPeptide: 1, positionInProtein: 2, intensity: 10);
            Assert.AreEqual(quantmod.IdWithMotif, "TestMod: ModX on AAY");
            Assert.AreEqual(quantmod.PeptidePositionZeroIsNTerminus, 1);
            Assert.AreEqual(quantmod.ProteinPositionZeroIsNTerminus, 2);
            Assert.AreEqual(quantmod.Intensity, 10);
            Assert.AreEqual(quantmod.ModificationLocalization, "Unknown");
        }

        [Test]
        public void TestQuantifiedPeptide()
        {
            var fullSeq1 = "[UniProt: N - palmitoyl glycine on G]G[UniProt: N - methylglycine on G]K[UniProt: O - linked(Hex) hydroxylysine on K]";
            var peptide1 = new QuantifiedPeptide(fullSeq1, intensity: 1);
            Assert.That(peptide1.FullSequences.Contains(fullSeq1));
            Assert.AreEqual(peptide1.BaseSequence, "GK");
            Assert.AreEqual(peptide1.Intensity, 1);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions.Count, 3);
            Assert.That(peptide1.ModifiedAminoAcidPositions.ContainsKey(0));
            Assert.That(peptide1.ModifiedAminoAcidPositions.ContainsKey(1));
            Assert.That(peptide1.ModifiedAminoAcidPositions.ContainsKey(2));
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[0].First().Value.IdWithMotif, "UniProt: N - palmitoyl glycine on G");
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[1].First().Value.IdWithMotif, "UniProt: N - methylglycine on G");
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[2].First().Value.IdWithMotif, "UniProt: O - linked(Hex) hydroxylysine on K");
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[0].First().Value.Intensity, 1);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[1].First().Value.Intensity, 1);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[2].First().Value.Intensity, 1);

            // Test MergePeptide method
            var fullSeq2 = "[UniProt: N - acetylglycine on G]G[UniProt: N - methylglycine on G]K[UniProt: O - linked(Hex) hydroxylysine on K]";
            var peptide2 = new QuantifiedPeptide(fullSeq2, intensity: 10);
            peptide1.MergePeptide(peptide2);

            Assert.That(peptide1.FullSequences.Contains(fullSeq2));
            Assert.AreEqual(peptide1.Intensity, 11);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions.Count, 3);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[0].Count, 2);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[1].Count, 1);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[2].Count, 1);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[0]["UniProt: N - palmitoyl glycine on G"].Intensity, 1);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[0]["UniProt: N - acetylglycine on G"].Intensity, 10);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[1].First().Value.Intensity, 11);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[2].First().Value.Intensity, 11);

            // Test AddFullSequence method
            var fullSeq3 = "GK[UniProt: O - linked(Hex) hydroxylysine on K]";
            peptide1.AddFullSequence(fullSeq3, intensity:100);

            Assert.That(peptide1.FullSequences.Contains(fullSeq3));
            Assert.AreEqual(peptide1.Intensity, 111);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions.Count, 3);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[0].Count, 2);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[1].Count, 1);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[2].Count, 1);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[0]["UniProt: N - palmitoyl glycine on G"].Intensity, 1);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[0]["UniProt: N - acetylglycine on G"].Intensity, 10);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[1].First().Value.Intensity, 11);
            Assert.AreEqual(peptide1.ModifiedAminoAcidPositions[2].First().Value.Intensity, 111);

            // Test failed merge due to base sequence mismatch
            var errorMessage = "The base sequence of the peptide being added does not match the base sequence of this peptide.";
            var exception1 = Assert.Throws<System.Exception>(() => peptide1.AddFullSequence("AK", intensity: 1));
            Assert.AreEqual(exception1.Message, errorMessage);

            var peptide3 = new QuantifiedPeptide("AK", intensity: 1);
            var exception2 = Assert.Throws<System.Exception>(() => peptide1.MergePeptide(peptide3));
            Assert.AreEqual(exception2.Message, errorMessage);
        }

        [Test]
        public void TestQuantifiedProtein()
        {

            var fullSeq1 = "[UniProt: N - palmitoyl glycine on G]G[UniProt: N - methylglycine on G]K[UniProt: O - linked(Hex) hydroxylysine on K]";
            var fullSeq2 = "[UniProt: N - acetylglycine on G]G[UniProt: N - methylglycine on G]K-[C-Terminal UniProt: Lysine Amide on K]";
            var fullSeq3 = "A[UniProt:N-methylalanine on A]K[UniProt: O - linked(Hex) hydroxylysine on K]-[C-Terminal UniProt: Lysine Amide on K]";

            var basePeptide1 = new QuantifiedPeptide(fullSeq1, intensity: 1);
            var basePeptide2 = new QuantifiedPeptide(fullSeq3, intensity: 100);

            basePeptide1.AddFullSequence(fullSeq2, intensity: 10);
            var peptides = new Dictionary<string, QuantifiedPeptide> {{ basePeptide1.BaseSequence, basePeptide1},
                                                                      { basePeptide2.BaseSequence, basePeptide2 }};

            var proteinSeq = "GKAAAAAAK";
            var protein = new QuantifiedProtein(accession: "TESTPROT", sequence: proteinSeq, peptides: peptides);
            var stoich = protein.GetModStoichiometryFromProteinMods();

            // Check object fields modified by SetProteinModsFromPeptides, which gets called first in the GetModStoichiometryFromProteinMods method. 
            Assert.AreEqual(protein.Accession, "TESTPROT");
            Assert.AreEqual(protein.Sequence, proteinSeq);
            Assert.AreEqual(protein.Peptides.Count, 2);
            Assert.AreEqual(protein.ModifiedAminoAcidPositionsInProtein.Count, 6);
            Assert.AreEqual(protein.ModifiedAminoAcidPositionsInProtein[0].Count, 2);
            Assert.AreEqual(protein.ModifiedAminoAcidPositionsInProtein[1].Count, 1);
            Assert.AreEqual(protein.ModifiedAminoAcidPositionsInProtein[2].Count, 1);
            Assert.AreEqual(protein.ModifiedAminoAcidPositionsInProtein[8].Count, 1);
            Assert.AreEqual(protein.ModifiedAminoAcidPositionsInProtein[9].Count, 1);
            Assert.AreEqual(protein.ModifiedAminoAcidPositionsInProtein[10].Count, 1);

            Assert.AreEqual(protein.ModifiedAminoAcidPositionsInProtein[0]["UniProt: N - palmitoyl glycine on G"].Intensity, 1);
            Assert.AreEqual(protein.ModifiedAminoAcidPositionsInProtein[0]["UniProt: N - acetylglycine on G"].Intensity, 10);
            Assert.AreEqual(protein.ModifiedAminoAcidPositionsInProtein[1]["UniProt: N - methylglycine on G"].Intensity, 11);
            Assert.AreEqual(protein.ModifiedAminoAcidPositionsInProtein[2]["UniProt: O - linked(Hex) hydroxylysine on K"].Intensity, 1);
            Assert.AreEqual(protein.ModifiedAminoAcidPositionsInProtein[8]["UniProt:N-methylalanine on A"].Intensity, 100);
            Assert.AreEqual(protein.ModifiedAminoAcidPositionsInProtein[9]["UniProt: O - linked(Hex) hydroxylysine on K"].Intensity, 100);
            Assert.AreEqual(protein.ModifiedAminoAcidPositionsInProtein[10]["C-Terminal UniProt: Lysine Amide on K"].Intensity, 100);

            // Check stoichiometry results
            Assert.AreEqual(stoich.Count, 6);
            Assert.AreEqual(stoich[0]["UniProt: N - palmitoyl glycine on G"], 1 / 11.0);
            Assert.AreEqual(stoich[0]["UniProt: N - acetylglycine on G"], 10 / 11.0);
            Assert.AreEqual(stoich[1]["UniProt: N - methylglycine on G"], 11 / 11.0);
            Assert.AreEqual(stoich[2]["UniProt: O - linked(Hex) hydroxylysine on K"], 1 / 11.0);
            Assert.AreEqual(stoich[8]["UniProt:N-methylalanine on A"], 1);
            Assert.AreEqual(stoich[9]["UniProt: O - linked(Hex) hydroxylysine on K"], 1);
            Assert.AreEqual(stoich[10]["C-Terminal UniProt: Lysine Amide on K"], 1);
        }

        [Test]
        public void TestQuantifiedProteinGroup()
        {   
            // Test correct arguments where protein group name contains the names of the proteins
            var protein1 = new QuantifiedProtein(accession: "PROT1", sequence: "AAAYYY", peptides: new Dictionary<string, QuantifiedPeptide>());
            var protein2 = new QuantifiedProtein(accession: "PROT2", sequence: "AAARRR", peptides: new Dictionary<string, QuantifiedPeptide>());
            var proteins = new Dictionary<string, QuantifiedProtein> { { protein1.Accession, protein1 },
                                                                       { protein2.Accession, protein2 } };
            var proteinGroup = new QuantifiedProteinGroup("PROT1|PROT2", proteins);
            Assert.AreEqual(proteinGroup.Proteins.Count, 2);
            Assert.AreEqual(proteinGroup.Proteins["PROT1"].Accession, "PROT1");
            Assert.AreEqual(proteinGroup.Proteins["PROT2"].Accession, "PROT2");

            // Test incorrect argument where protein group name does not contain the names of the proteins
            var errorMessage = "The number of proteins provided does not match the number of proteins in the protein group name.";
            var exception1 = Assert.Throws<System.Exception>(() => new QuantifiedProteinGroup("PROT1|PROT2", new Dictionary<string, QuantifiedProtein> { { protein1.Accession, protein1 } }));
            Assert.AreEqual(exception1.Message, errorMessage);

            var exception2 = Assert.Throws<System.Exception>(() => new QuantifiedProteinGroup("PROT1", proteins));
            Assert.AreEqual(exception2.Message, errorMessage);

            var exception3 = Assert.Throws<System.Exception>(() => new QuantifiedProteinGroup("PROT1|PROT2|PROT3", proteins));
            Assert.AreEqual(exception3.Message, errorMessage);
        }

        [Test]
        public void TestSetUpQuantificationObjects()
        {
            var fullSeq1 = "[UniProt: N - palmitoyl glycine on G]G[UniProt: N - methylglycine on G]K[UniProt: O - linked(Hex) hydroxylysine on K]";
            var fullSeq2 = "[UniProt: N - acetylglycine on G]G[UniProt: N - methylglycine on G]K-[C-Terminal UniProt: Lysine Amide on K]";
            var fullSequences = new List<string> { fullSeq1, fullSeq2 };
            var proteinGroups = new List<string> { "TESTPROT1|TESTPROT2", "TESTPROT3" };
            var proteinSequences = new Dictionary<string, string> { { "TESTPROT1", "GKAAAAAAK" },
                                                                    { "TESTPROT2", "AKAAAAAGK" },
                                                                    { "TESTPROT3", "AKGK"} };
            var intensities = new List<double> { 1, 5 };
            var sequenceInputs = new List<(string, List<string>, double)> { };
            for (int i = 0; i < 2; i++)
            {
                sequenceInputs.Add((fullSequences[i], proteinGroups, intensities[i]));
            }
            sequenceInputs.Add(("AAAA", new List<string> { "TESTPROT1|TESTPROT2" }, 10));

            var quant = new PositionFrequencyAnalysis();
            quant.SetUpQuantificationObjectsFromFullSequences(sequenceInputs, proteinSequences);
            Assert.AreEqual(quant.ProteinGroups.Count, 2);
            Assert.That(quant.ProteinGroups["TESTPROT1|TESTPROT2"].Proteins.Keys.Contains("TESTPROT1"));
            Assert.That(quant.ProteinGroups["TESTPROT1|TESTPROT2"].Proteins.Keys.Contains("TESTPROT2"));
            Assert.AreEqual(quant.ProteinGroups["TESTPROT1|TESTPROT2"].Proteins["TESTPROT1"].Accession, "TESTPROT1");
            Assert.AreEqual(quant.ProteinGroups["TESTPROT1|TESTPROT2"].Proteins["TESTPROT1"].Sequence, "GKAAAAAAK");
            Assert.AreEqual(quant.ProteinGroups["TESTPROT1|TESTPROT2"].Proteins["TESTPROT1"].Peptides.Count, 2);
            Assert.AreEqual(quant.ProteinGroups["TESTPROT1|TESTPROT2"].Proteins["TESTPROT1"].Peptides["GK"].FullSequences.Count, 2);
            Assert.AreEqual(quant.ProteinGroups["TESTPROT1|TESTPROT2"].Proteins["TESTPROT1"].Peptides["AAAA"].FullSequences.Count, 1);
            Assert.AreEqual(quant.ProteinGroups["TESTPROT1|TESTPROT2"].Proteins["TESTPROT2"].Accession, "TESTPROT2");
            Assert.AreEqual(quant.ProteinGroups["TESTPROT1|TESTPROT2"].Proteins["TESTPROT2"].Sequence, "AKAAAAAGK");
            Assert.AreEqual(quant.ProteinGroups["TESTPROT1|TESTPROT2"].Proteins["TESTPROT2"].Peptides.Count, 2);
            Assert.AreEqual(quant.ProteinGroups["TESTPROT1|TESTPROT2"].Proteins["TESTPROT2"].Peptides["GK"].FullSequences.Count, 2);
            Assert.AreEqual(quant.ProteinGroups["TESTPROT1|TESTPROT2"].Proteins["TESTPROT2"].Peptides["AAAA"].FullSequences.Count, 1);

            Assert.That(quant.ProteinGroups["TESTPROT3"].Proteins.Keys.Contains("TESTPROT3"));
            Assert.AreEqual(quant.ProteinGroups["TESTPROT3"].Proteins["TESTPROT3"].Accession, "TESTPROT3");
            Assert.AreEqual(quant.ProteinGroups["TESTPROT3"].Proteins["TESTPROT3"].Sequence, "AKGK");
            Assert.AreEqual(quant.ProteinGroups["TESTPROT3"].Proteins["TESTPROT3"].Peptides.Count, 1);
        }

        public struct TestStruct
        {
            public int X { get; set; }
            public int Y { get; set; }
        }
    }
}
