using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    public sealed class TestProductMassesMightHaveDuplicates
    {
        private static Stopwatch Stopwatch { get; set; }

        [OneTimeSetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [OneTimeTearDown]
        public static void TearDown()
        {
            lock (FixtureSetUp.ConsoleLock)
                Console.WriteLine($"TestProductMassesMightHaveDuplicates Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public static void Test_UnmodifiedPeptide_AllProductType_fragmentMasses()
        {
            Protein p = new Protein("PET", "accession");

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);

            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();

            var allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.Unknown, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.CID, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { 227, 120, 249 }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.IRMPD, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { 98, 227, 120, 249 }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.ECD, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { 115, 244, 120, 249, 104, 233 }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.PQD, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.ETD, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { 115, 244, 120, 249, 104, 233 }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { 98, 227, 120, 249 }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.AnyActivationType, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { 98, 227, 120, 249 }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.EThcD, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { 98, 227, 115, 244, 120, 249, 104, 233 }));

            DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom] = new List<ProductType> { };
            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.Custom, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.ISCID, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { }));
        }
    }
}