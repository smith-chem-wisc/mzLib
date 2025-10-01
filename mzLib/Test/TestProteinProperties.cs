using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics;
using Omics.BioPolymer;
using Omics.Modifications;
using Stopwatch = System.Diagnostics.Stopwatch;
using Transcriptomics;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestProteinProperties
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setup()
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
        public void TestHashAndEqualsProtein()
        {
            Protein p1 = new Protein("MSEQ", "accession");
            Protein p11 = new Protein("MSEQ", "accession");
            Assert.AreEqual(p1, p11); // default object hash and equals are used
        }
        [Test]
        public void TestHashAndEqualsSequenceVariation()
        {
            // Base modifications
            var modM1 = new Modification("m1");
            var modM1Clone = new Modification("m1"); // logically identical (same id)
            var modM2 = new Modification("m2");

            // Variant-specific modification dictionaries (post-variation coordinates)
            var modsPos11_M1 = new Dictionary<int, List<Modification>> { { 11, new() { modM1 } } };
            var modsPos11_M1Clone = new Dictionary<int, List<Modification>> { { 11, new() { modM1Clone } } }; // value-equal
            var modsPos11_M2 = new Dictionary<int, List<Modification>> { { 11, new() { modM2 } } };
            var modsPos12_M1 = new Dictionary<int, List<Modification>> { { 12, new() { modM1 } } };

            // Multiple mods at same site (order-insensitive)
            var modsMultiAB = new Dictionary<int, List<Modification>>
            {
                { 11, new() { new Modification("mA"), new Modification("mB") } }
            };
            var modsMultiBA = new Dictionary<int, List<Modification>>
            {
                { 11, new() { new Modification("mB"), new Modification("mA") } }
            };

            // Baseline valid synonymous (no-op) but WITH a variant-specific mod (required for validity)
            var svBase1 = new SequenceVariation(
                oneBasedBeginPosition: 10,
                oneBasedEndPosition: 12,
                originalSequence: "AAA",
                variantSequence: "AAA",
                description: "desc",
                variantCallFormatDataString: "VCF1",
                oneBasedModifications: modsPos11_M1);

            // Same logical content, different description (ignored in equality)
            var svBase2 = new SequenceVariation(
                10, 12, "AAA", "AAA",
                "different description",
                "VCF1",
                modsPos11_M1Clone);

            var svDiffDescription = new SequenceVariation(
                10, 12, "AAA", "AAA",
                "another annotation",
                "VCF1",
                modsPos11_M1); // still equal to svBase1

            // Different modification position
            var svDiffModSite = new SequenceVariation(10, 12, "AAA", "AAA", "desc", "VCF1", modsPos12_M1);
            // Different modification identity
            var svDiffModIdentity = new SequenceVariation(10, 12, "AAA", "AAA", "desc", "VCF1", modsPos11_M2);
            // Different VCF metadata
            var svDiffVcf = new SequenceVariation(10, 12, "AAA", "AAA", "desc", "VCF2", modsPos11_M1);
            // Different span
            var svDiffSpan = new SequenceVariation(11, 13, "AAA", "AAA", "desc", "VCF1", modsPos11_M1);
            // Different original sequence
            var svDiffOriginal = new SequenceVariation(10, 12, "AAB", "AAA", "desc", "VCF1", modsPos11_M1);
            // Different variant sequence
            var svDiffVariant = new SequenceVariation(10, 12, "AAA", "AAT", "desc", "VCF1", modsPos11_M1);

            // Multi-mod order-insensitivity
            var svMultiA = new SequenceVariation(10, 12, "AAA", "AAA", "multiA", "VCF1", modsMultiAB);
            var svMultiB = new SequenceVariation(10, 12, "AAA", "AAA", "multiB", "VCF1", modsMultiBA);

            // Insertion (expansion)
            var svInsertion1 = new SequenceVariation(
                5, 5, "A", "ATG",
                "insertion", "VCF_INS",
                new Dictionary<int, List<Modification>> { { 5, new() { new Modification("mI") } } });

            var svInsertion2 = new SequenceVariation(
                5, 5, "A", "ATG",
                "insertion alt desc", "VCF_INS",
                new Dictionary<int, List<Modification>> { { 5, new() { new Modification("mI") } } });

            // Deletion (contraction)
            var svDeletion1 = new SequenceVariation(
                7, 9, "ATG", "A",
                "deletion", "VCF_DEL",
                null);

            var svDeletion2 = new SequenceVariation(
                7, 9, "ATG", "A",
                "deletion alt", "VCF_DEL",
                null);

            // INVALID CASES (no-op without variant-specific modifications) should throw
            // 1. Synonymous without mods
            Assert.Throws<ArgumentException>(() => _ = new SequenceVariation(15, 15, "G", "G", "no_op", "VCF_SYN", null),
                "No-op variant without variant-specific modifications must be invalid.");
            // 2. Whole-span no-op without mods
            Assert.Throws<ArgumentException>(() => _ = new SequenceVariation(10, 12, "AAA", "AAA", "no_op2", "VCF1", null),
                "Whole-span no-op without mods must be invalid.");

            // Positive equality
            Assert.AreEqual(svBase1, svBase2, "Baseline synonymous with equivalent mods should be equal.");
            Assert.AreEqual(svBase1, svDiffDescription, "Description difference should be ignored.");
            Assert.AreEqual(svMultiA, svMultiB, "Modification order should not affect equality.");
            Assert.AreEqual(svInsertion1, svInsertion2, "Equivalent insertions should be equal.");
            Assert.AreEqual(svDeletion1, svDeletion2, "Equivalent deletions should be equal.");

            // Hash code parity for equal objects
            Assert.AreEqual(svBase1.GetHashCode(), svBase2.GetHashCode(), "Equal variations must share hash code.");
            Assert.AreEqual(svInsertion1.GetHashCode(), svInsertion2.GetHashCode(), "Equal insertions must share hash code.");
            Assert.AreEqual(svMultiA.GetHashCode(), svMultiB.GetHashCode(), "Equal multi-mod variants must share hash code.");
            Assert.AreEqual(svDeletion1.GetHashCode(), svDeletion2.GetHashCode(), "Equal deletions must share hash code.");

            // Negative equality
            Assert.AreNotEqual(svBase1, svDiffModSite, "Different modification site should differ.");
            Assert.AreNotEqual(svBase1, svDiffModIdentity, "Different modification identity should differ.");
            Assert.AreNotEqual(svBase1, svDiffVcf, "Different VCF metadata should differ.");
            Assert.AreNotEqual(svBase1, svDiffSpan, "Different span should differ.");
            Assert.AreNotEqual(svBase1, svDiffOriginal, "Different original sequence should differ.");
            Assert.AreNotEqual(svBase1, svDiffVariant, "Different variant sequence should differ.");
            Assert.AreNotEqual(svBase1, svMultiA, "Different modification sets (different content) should differ.");

            // Collapsed set (description ignored). Unique logical keys:
            // 1. (10-12 AAA->AAA, mod at 11 m1)
            // 2. (10-12 AAA->AAA, mods at 11 mA+mB)
            // 3. (5-5  A->ATG)
            // 4. (7-9  ATG->A)
            var collapsed = new HashSet<SequenceVariation>
            {
                svBase1, svBase2, svDiffDescription,
                svMultiA, svMultiB,
                svInsertion1, svInsertion2,
                svDeletion1, svDeletion2
            };
            Assert.AreEqual(4, collapsed.Count, "HashSet should collapse logically equivalent variants.");
            Assert.IsTrue(collapsed.Contains(svBase1));
            Assert.IsTrue(collapsed.Contains(svInsertion1));
            Assert.IsTrue(collapsed.Contains(svDeletion1));
            Assert.IsTrue(collapsed.Contains(svMultiA));
        }
        [Test]
        public void TestProteinVariantModMethods()
        {
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            Modification mp = new Modification("mod", null, "type", null, motifP, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);
            int mpModLocationInCanonicalProtein = 4;

            ModificationMotif.TryGetMotif("T", out ModificationMotif motifT);
            Modification mt = new Modification("mod", null, "type", null, motifT, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);
            int mtModLocationInCanonicalProtein = 5;

            Protein protein1 = new Protein("MPEPTIDE", "protein1",
                sequenceVariations: new List<SequenceVariation>
                {
                    new SequenceVariation(4, 4, "P", "PV",
                        @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30")
                },
                oneBasedModifications: new Dictionary<int, List<Modification>> 
                { 
                    { mpModLocationInCanonicalProtein, new[] { mp }.ToList() },
                    { mtModLocationInCanonicalProtein, new[] { mt }.ToList() }
                });

            // If a protein has a modification at a specific residue and that residue is a potential variant, we need to include a variant protein bearing the same
            // modification in the variant database. IsSequenceVariantModifiation is used in MetaMorpheus to facilitate that process.
            Assert.IsTrue(VariantApplication.IsSequenceVariantModification(protein1.SequenceVariations.First(), mpModLocationInCanonicalProtein));
            Assert.IsFalse(VariantApplication.IsSequenceVariantModification(protein1.SequenceVariations.First(), mtModLocationInCanonicalProtein));
            
            int mtModLocationInVariant = 7;
            Protein protein2 = new Protein("MPEPTIDE", "protein2",
                appliedSequenceVariations: new List<SequenceVariation>
                {
                    new SequenceVariation(4, 4, "P", "PPP",
                    "",
                        @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30",
                        new Dictionary<int, List<Modification>> { { mtModLocationInVariant, new[] { mt }.ToList() } })
                });

            // RestoreModificationIndex translates the position of a modification in a variant protein to the equivalent location on the canonical protein
            // Here, the variant is MPEPPPT[Modification on T: mt]IDE, with the mt mod at position 7
            // and the canonical is MPEPT[[Modification on T: mt]IDE, with the mt mod at position 5
            Assert.AreEqual(VariantApplication.RestoreModificationIndex(protein2, mtModLocationInVariant), mtModLocationInCanonicalProtein);
        }

        [Test]
        public void TestHashAndEqualsDbRef()
        {
            DatabaseReference db1 = new DatabaseReference("type", "id", new List<Tuple<string, string>> { new Tuple<string, string>("1", "2") });
            DatabaseReference db11 = new DatabaseReference("type", "id", new List<Tuple<string, string>> { new Tuple<string, string>("1", "2") });
            DatabaseReference db2 = new DatabaseReference("type", "id", new List<Tuple<string, string>> { new Tuple<string, string>("1", "3") });
            DatabaseReference db3 = new DatabaseReference("type", "id", null);
            DatabaseReference db4 = new DatabaseReference("type", null, new List<Tuple<string, string>> { new Tuple<string, string>("1", "2") });
            DatabaseReference db5 = new DatabaseReference(null, "id", new List<Tuple<string, string>> { new Tuple<string, string>("1", "2") });
            Assert.AreEqual(db1, db11);
            Assert.AreNotEqual(db1, db2);
            Assert.AreNotEqual(db1, db3);
            Assert.AreNotEqual(db1, db4);
            Assert.AreNotEqual(db1, db5);
        }

        [Test]
        public void TestHashAndEqualsSpliceSite()
        {
            SpliceSite ss1 = new SpliceSite(1, 2, "description");
            SpliceSite ss11 = new SpliceSite(1, 2, "description");
            SpliceSite ss2 = new SpliceSite(1, 2, null);
            SpliceSite ss3 = new SpliceSite(1, "description");
            SpliceSite ss4 = new SpliceSite(2, "description");
            Assert.AreEqual(ss1, ss11);
            Assert.AreNotEqual(ss1, ss2);
            Assert.AreNotEqual(ss1, ss3);
            Assert.AreNotEqual(ss1, ss4);
        }

        [Test]
        public void TestHashAndEqualsDisulfide()
        {
            DisulfideBond bond7 = new DisulfideBond(1, 2, "description");
            DisulfideBond bond007 = new DisulfideBond(1, 2, "description");
            DisulfideBond bond8 = new DisulfideBond(1, 2, null);
            DisulfideBond bond9 = new DisulfideBond(1, "description");
            DisulfideBond bond17 = new DisulfideBond(2, "description");
            Assert.AreEqual(bond7, bond007);
            Assert.AreNotEqual(bond007, bond8);
            Assert.AreNotEqual(bond007, bond9);
            Assert.AreNotEqual(bond007, bond17);
        }

        [Test]
        public void TestHashAndEqualsProteolysis()
        {
            TruncationProduct pp1 = new TruncationProduct(1, 2, "type");
            TruncationProduct pp11 = new TruncationProduct(1, 2, "type");
            TruncationProduct pp2 = new TruncationProduct(1, 2, null);
            TruncationProduct pp3 = new TruncationProduct(1, null, "type");
            TruncationProduct pp4 = new TruncationProduct(null, 2, "type");
            TruncationProduct pp5 = new TruncationProduct(1, 1, "type");
            TruncationProduct pp6 = new TruncationProduct(2, 2, "type");
            Assert.AreEqual(pp1, pp11);
            Assert.AreNotEqual(pp1, pp2);
            Assert.AreNotEqual(pp1, pp3);
            Assert.AreNotEqual(pp1, pp4);
            Assert.AreNotEqual(pp1, pp5);
            Assert.AreNotEqual(pp1, pp6);
        }
        [Test]
        public static void CompareProteinProperties()
        {
            DatabaseReference d = new DatabaseReference("asdf", "asdfg", new List<Tuple<string, string>> { new Tuple<string, string>("bbb", "ccc") });
            DatabaseReference dd = new DatabaseReference("asdf", "asdfg", new List<Tuple<string, string>> { new Tuple<string, string>("bbb", "ccc") });
            DatabaseReference de = new DatabaseReference("asdf", "asdefg", new List<Tuple<string, string>> { new Tuple<string, string>("bbb", "ccc") });
            DatabaseReference df = new DatabaseReference("asddf", "asdfg", new List<Tuple<string, string>> { new Tuple<string, string>("bbb", "ccc") });
            DatabaseReference dg = new DatabaseReference("asdf", "asdfg", new List<Tuple<string, string>> { new Tuple<string, string>("babb", "ccc") });
            DatabaseReference dh = new DatabaseReference("asdf", "asdfg", new List<Tuple<string, string>> { new Tuple<string, string>("bbb", "cccf") });
            Assert.True(dd.Equals(d));
            Assert.False(de.Equals(d));
            Assert.False(df.Equals(d));
            Assert.False(dg.Equals(d));
            Assert.False(dh.Equals(d));
            Assert.AreEqual(5, new HashSet<DatabaseReference> { d, dd, de, df, dg, dh }.Count);

            // SequenceVariation equality DOES NOT include Description (see SequenceVariation.Equals)
            // Only coordinates, original/variant sequences, VCF data, and modification dictionaries are compared.
            SequenceVariation s = new SequenceVariation(1, "hello", "hey", "hi");
            SequenceVariation sv = new SequenceVariation(1, "hello", "hey", "hi");   // identical
            SequenceVariation sss = new SequenceVariation(2, "hallo", "hey", "hi");  // different begin/original
            SequenceVariation ssss = new SequenceVariation(1, "hello", "heyy", "hi"); // different variant seq
            SequenceVariation sssss = new SequenceVariation(1, "hello", "hey", "hii"); // ONLY description differs -> equal to s

            Assert.True(s.Equals(sv));
            Assert.False(s.Equals(sss));
            Assert.False(s.Equals(ssss));
            Assert.True(s.Equals(sssss)); // updated: description difference alone does NOT affect equality

            // Unique set should collapse s, sv, sssss into one entry
            Assert.AreEqual(3, new HashSet<SequenceVariation> { s, sv, sss, ssss, sssss }.Count);

            DisulfideBond b = new DisulfideBond(1, "hello");
            DisulfideBond bb = new DisulfideBond(1, "hello");
            DisulfideBond bbb = new DisulfideBond(1, 2, "hello");
            DisulfideBond bbbb = new DisulfideBond(1, 2, "hello");
            DisulfideBond ba = new DisulfideBond(1, 3, "hello");
            DisulfideBond baa = new DisulfideBond(2, 2, "hello");
            DisulfideBond baaa = new DisulfideBond(1, 2, "hallo");
            Assert.AreEqual(b, bb);
            Assert.AreEqual(bbb, bbbb);
            Assert.AreNotEqual(b, bbb);
            Assert.AreNotEqual(ba, bbb);
            Assert.AreNotEqual(baa, bbb);
            Assert.AreNotEqual(baaa, bbb);
            Assert.AreEqual(5, new HashSet<DisulfideBond> { b, bb, bbb, bbbb, ba, baa, baaa }.Count);

            TruncationProduct pp = new TruncationProduct(1, 1, "hello");
            TruncationProduct paaa = new TruncationProduct(1, 1, "hello");
            TruncationProduct p = new TruncationProduct(null, null, "hello");
            TruncationProduct ppp = new TruncationProduct(1, 2, "hello");
            TruncationProduct pa = new TruncationProduct(2, 1, "hello");
            TruncationProduct paa = new TruncationProduct(1, 1, "hallo");
            Assert.AreEqual(pp, paaa);
            Assert.AreNotEqual(p, pp);
            Assert.AreNotEqual(pp, ppp);
            Assert.AreNotEqual(pp, pa);
            Assert.AreNotEqual(pp, paa);
            Assert.AreEqual(5, new HashSet<TruncationProduct> { p, pp, ppp, pa, paa, paaa }.Count);
        }
        [Test]
        public static void TestProteoformClassification()//string inputPath)
        {
            //Test classifier
            string inputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "ProteoformClassificationUnitTest.csv");
            string[] lines = File.ReadAllLines(inputPath);

            List<string> expectedLevels = new List<string> { "1", "2A", "2B", "2C", "2D", "3", "4", "5", "2C", "2B" }; //each of these should be identified in the vignette

            //iterate through each result, check if there's a header or not
            for (int i = 1; i < lines.Length; i++)
            {
                string[] line = lines[i].Split(',');
                //should be scan number, sequence(s), gene(s)
                string level = ProteoformLevelClassifier.ClassifyPrSM(line[1], line[2]);

                Assert.IsTrue(level.Equals(expectedLevels[i - 1]));
            }

            //Test weird case where residues and PTMs are the same, but flipped. Should be 2A.
            string fullSeq = "ARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVKKPHRYRPGTVALREIRRYQK[UniProt:N6-glutaryllysine on K]STELLIRKLPFQRLVREIAQDFK[UniProt:N6-succinyllysine on K]TDLRFQSAAIGALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVKKPHRYRPGTVALREIRRYQK[UniProt:N6-glutaryllysine on K]STELLIRKLPFQRLVREIAQDFK[UniProt:N6-succinyllysine on K]TDLRFQSAAIGALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVKKPHRYRPGTVALREIRRYQK[UniProt:N6-succinyllysine on K]STELLIRKLPFQRLVREIAQDFK[UniProt:N6-glutaryllysine on K]TDLRFQSAAIGALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVKKPHRYRPGTVALREIRRYQK[UniProt:N6-succinyllysine on K]STELLIRKLPFQRLVREIAQDFK[UniProt:N6-glutaryllysine on K]TDLRFQSAAIGALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA";
            string gene = "primary:H3-3A, synonym:H3.3A, synonym:H3F3, synonym:H3F3A, ORF:PP781, primary:H3-3B, synonym:H3.3B, synonym:H3F3B";
            string weirdLevel = ProteoformLevelClassifier.ClassifyPrSM(fullSeq, gene);
            Assert.IsTrue(weirdLevel.Equals("2A"));

            ///Explanation of weird case, which is test case 1
            ///Test case 1, which should be level 2A:
            ///…K[UniProt: N6 - glutaryllysine on K]…K[UniProt: N6 - succinyllysine on K]…
            ///…K[UniProt: N6 - succinyllysine on K]…K[UniProt: N6 - glutaryllysine on K]…

            ///Test case 2, which is level 2B:
            ///…K[Acetylation]…
            ///…K[Trimethylation]…

            ///Test case 3, which should be level 2B:
            ///…K[Acetylation]…K[Acetylation]…
            ///…K[Acetylation]…K[Trimethylation]…
            ///…K[Trimethylation]…K[Acetylation]…
            ///…K[Trimethylation]…K[Trimethylation]…

            ///Test case 3 is 2B (not level 3) because you've localized the mod, you just aren't sure what mod it is.
            ///In test case 1, you know what the mods are, but you're not sure where they belong.
        }

        [Test]
        public void TestProteinEquals()
        {
            string sequence = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQR";
            string accession = "P02768";
            Protein protein1 = new Protein(sequence, accession);
            Protein protein2 = new Protein(sequence, accession);
            
            NUnit.Framework.Assert.That(protein1.Equals(protein2), Is.True);
            NUnit.Framework.Assert.That(protein1.Equals((object)protein2), Is.True);
            NUnit.Framework.Assert.That(protein1.Equals(null), Is.False);
        }

        [Test]
        public void TestProteinGetHashCode()
        {
            string sequence = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQR";
            string accession = "P02768";
            Protein protein = new Protein(sequence, accession);

            NUnit.Framework.Assert.That(protein.GetHashCode(), Is.EqualTo(sequence.GetHashCode()));
        }

        [Test]
        public void TestProteinRnaEquality()
        {
            string sequence = "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLFFAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQR";
            string accession = "P02768";
            Protein protein1 = new Protein(sequence, accession);
            RNA rna = new RNA("GUACUG");


            Assert.That(!rna.Equals(protein1));
            Assert.That(!protein1.Equals(rna));
            Assert.That(!((IBioPolymer)rna).Equals(protein1));
            Assert.That(!((IBioPolymer)protein1).Equals(rna));
        }
    }
}