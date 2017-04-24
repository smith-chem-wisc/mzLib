using Chemistry;
using NetSerializer;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    internal class TestSerializing
    {
        #region Public Methods

        [Test]
        public void test_serialize_mod_and_protein()
        {
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("K", out motif);
            ModificationWithMass m = new ModificationWithMass("asdf", new Tuple<string, string>("", ""), motif, ModificationSites.K, 10.01, new Dictionary<string, IList<string>>(), new List<double>(), new List<double>(), "");
            ModificationWithMassAndCf mcf = new ModificationWithMassAndCf("asdf", new Tuple<string, string>("", ""), motif, ModificationSites.K, new ChemicalFormula(), 10.01, new Dictionary<string, IList<string>>(), new List<double>(), new List<double>(), "");

            something yup = new something();
            something yup2 = new something();
            yup.nice = new List<Modification>
            {
                new Modification("", ""),
                new ModificationWithLocation("fayk",null, null,ModificationSites.A,null,  null),
                m,
                mcf
            };
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
                    typeof(List<double>)
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
