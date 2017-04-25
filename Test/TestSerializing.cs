using Chemistry;
using NetSerializer;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    internal class TestSerializing
    {
        #region Public Methods

        [Test]
        public void basic_test_serialize_mod_and_protein()
        {
            something yup = new something();
            something yup2 = new something();
            Loaders.LoadElements(Path.Combine(TestContext.CurrentContext.TestDirectory, @"elements.txt"));
            Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, @"ptmlist.txt"));
            yup.nice = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, @"ptmlist.txt")).OfType<Modification>().ToList();
            yup.ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml2.xml"), false, yup.nice, false, null, out Dictionary<string, Modification> un);

            Serializer ser = new Serializer(
                new Type[] {
                    typeof(something),
                    typeof(Protein),
                    typeof(ModificationWithLocation),
                    typeof(ModificationWithMass),
                    typeof(ModificationWithMassAndCf),
                    typeof(List<DatabaseReference>),
                    typeof(List<Tuple<string,string>>),
                    typeof(Dictionary<int, List<Modification>>),
                    typeof(Dictionary<string, IList<string>>),
                    typeof(List<ProteolysisProduct>),
                    typeof(ChemicalFormulaTerminus),
                    typeof(List<double>),
                    typeof(List<string>),
                });
            using (var file = File.Create(Path.Combine(TestContext.CurrentContext.TestDirectory, @"serialtest")))
                ser.Serialize(file, yup);
            using (var file = File.OpenRead(Path.Combine(TestContext.CurrentContext.TestDirectory, @"serialtest")))
                yup2 = (something)ser.Deserialize(file);

            Assert.AreEqual(yup.nice[0], yup2.nice[0]);
            Assert.AreEqual(yup.ok.Count, yup2.ok.Count);
        }

        [Test]
        public void test_serialize_mod_and_protein()
        {
            something yup = new something();
            something yup2 = new something();
            Loaders.LoadElements(Path.Combine(TestContext.CurrentContext.TestDirectory, @"elements.txt"));
            Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, @"ptmlist.txt"));
            yup.nice = PtmListLoader.ReadModsFromFile(Path.Combine(TestContext.CurrentContext.TestDirectory, @"ptmlist.txt")).OfType<Modification>().ToList();
            yup.ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"xml2.xml"), false, yup.nice, false, null, out Dictionary<string, Modification> un);

            Serializer ser = new Serializer(
                new Type[] {
                    typeof(something),
                    typeof(Protein),
                    typeof(ModificationWithLocation),
                    typeof(ModificationWithMass),
                    typeof(ModificationWithMassAndCf),
                    typeof(List<DatabaseReference>),
                    typeof(List<Tuple<string,string>>),
                    typeof(Dictionary<int, List<Modification>>),
                    typeof(Dictionary<string, IList<string>>),
                    typeof(List<ProteolysisProduct>),
                    typeof(ChemicalFormulaTerminus),
                    typeof(List<double>),
                    typeof(List<string>),
                });
            using (var file = File.Create(Path.Combine(TestContext.CurrentContext.TestDirectory, @"serialtest")))
                ser.Serialize(file, yup);
            using (var file = File.OpenRead(Path.Combine(TestContext.CurrentContext.TestDirectory, @"serialtest")))
                yup2 = (something) ser.Deserialize(file);

            Assert.AreEqual(yup.nice[0], yup2.nice[0]);
            Assert.AreEqual(yup.ok.Count, yup2.ok.Count);
        }

        #endregion Public Methods

        [Serializable]
        class something
        {
            public List<Modification> nice;
            public List<Protein> ok;
        }
    }
}
