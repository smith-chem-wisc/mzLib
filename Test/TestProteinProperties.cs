using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    public class TestProteinProperties
    {
        private static Stopwatch Stopwatch { get; set; }

        [OneTimeSetUp]
        public static void Setup()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [OneTimeTearDown]
        public static void TearDown()
        {
            lock (FixtureSetUp.ConsoleLock)
                Console.WriteLine($"TestProteinProperties Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public void TestHashAndEqualsProtein()
        {
            Protein p1 = new Protein("MSEQ", "accession");
            Protein p11 = new Protein("MSEQ", "accession");
            Assert.AreNotEqual(p1, p11); // default object hash and equals are used
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
    }
}