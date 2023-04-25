using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Stopwatch = System.Diagnostics.Stopwatch;

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
            SequenceVariation sv1 = new SequenceVariation(1, "MAA", "MAA", "description", new Dictionary<int, List<Modification>> { { 2, new[] { new Modification("mod") }.ToList() } });
            SequenceVariation sv2 = new SequenceVariation(1, "MAA", "MAA", "description", new Dictionary<int, List<Modification>> { { 2, new[] { new Modification("mod") }.ToList() } });
            SequenceVariation sv22 = new SequenceVariation(1, "MAA", "MAA", "description", new Dictionary<int, List<Modification>> { { 3, new[] { new Modification("mod") }.ToList() } });
            SequenceVariation sv222 = new SequenceVariation(1, "MAA", "MAA", "description", new Dictionary<int, List<Modification>> { { 2, new[] { new Modification("another") }.ToList() } });
            SequenceVariation sv3 = new SequenceVariation(1, "MAA", "MAA", "description", null);
            SequenceVariation sv4 = new SequenceVariation(1, "MAA", "MAA", null, new Dictionary<int, List<Modification>> { { 2, new[] { new Modification("mod") }.ToList() } });
            SequenceVariation sv5 = new SequenceVariation(1, null, null, "description", new Dictionary<int, List<Modification>> { { 2, new[] { new Modification("mod") }.ToList() } });
            SequenceVariation sv6 = new SequenceVariation(2, "MAA", "MAA", "description", new Dictionary<int, List<Modification>> { { 2, new[] { new Modification("mod") }.ToList() } });
            Assert.AreEqual(sv1, sv2);
            Assert.AreNotEqual(sv1, sv22);
            Assert.AreNotEqual(sv1, sv222);
            Assert.AreNotEqual(sv1, sv3);
            Assert.AreNotEqual(sv1, sv4);
            Assert.AreNotEqual(sv1, sv5);
            Assert.AreNotEqual(sv1, sv6);
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
            ProteolysisProduct pp1 = new ProteolysisProduct(1, 2, "type");
            ProteolysisProduct pp11 = new ProteolysisProduct(1, 2, "type");
            ProteolysisProduct pp2 = new ProteolysisProduct(1, 2, null);
            ProteolysisProduct pp3 = new ProteolysisProduct(1, null, "type");
            ProteolysisProduct pp4 = new ProteolysisProduct(null, 2, "type");
            ProteolysisProduct pp5 = new ProteolysisProduct(1, 1, "type");
            ProteolysisProduct pp6 = new ProteolysisProduct(2, 2, "type");
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

            SequenceVariation s = new SequenceVariation(1, "hello", "hey", "hi");
            SequenceVariation sv = new SequenceVariation(1, "hello", "hey", "hi");
            SequenceVariation sss = new SequenceVariation(2, "hallo", "hey", "hi");
            SequenceVariation ssss = new SequenceVariation(1, "hello", "heyy", "hi");
            SequenceVariation sssss = new SequenceVariation(1, "hello", "hey", "hii");
            Assert.True(s.Equals(sv));
            Assert.False(s.Equals(sss));
            Assert.False(s.Equals(ssss));
            Assert.False(s.Equals(sssss));
            Assert.AreEqual(4, new HashSet<SequenceVariation> { s, sv, sss, ssss, sssss }.Count);

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

            ProteolysisProduct pp = new ProteolysisProduct(1, 1, "hello");
            ProteolysisProduct paaa = new ProteolysisProduct(1, 1, "hello");
            ProteolysisProduct p = new ProteolysisProduct(null, null, "hello");
            ProteolysisProduct ppp = new ProteolysisProduct(1, 2, "hello");
            ProteolysisProduct pa = new ProteolysisProduct(2, 1, "hello");
            ProteolysisProduct paa = new ProteolysisProduct(1, 1, "hallo");
            Assert.AreEqual(pp, paaa);
            Assert.AreNotEqual(p, pp);
            Assert.AreNotEqual(pp, ppp);
            Assert.AreNotEqual(pp, pa);
            Assert.AreNotEqual(pp, paa);
            Assert.AreEqual(5, new HashSet<ProteolysisProduct> { p, pp, ppp, pa, paa, paaa }.Count);
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
    }
}