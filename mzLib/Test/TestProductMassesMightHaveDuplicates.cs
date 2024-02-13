using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Omics.Fragmentation;
using Stopwatch = System.Diagnostics.Stopwatch;
using Omics.Fragmentation.Peptide;
using Omics.Modifications;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestProductMassesMightHaveDuplicates
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
        public static void Test_UnmodifiedPeptide_AllProductType_fragmentMasses()
        {
            Protein p = new Protein("PET", "accession");

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);

            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();
            var fragments = new List<Product>();
            aPeptideWithSetModifications.Fragment(DissociationType.Unknown, FragmentationTerminus.Both, fragments);

            var allFragmentIonMzs = new HashSet<int>(fragments.Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { }));

            aPeptideWithSetModifications.Fragment(DissociationType.CID, FragmentationTerminus.Both, fragments);
            allFragmentIonMzs = new HashSet<int>(fragments.Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { 227, 120, 249 }));

            aPeptideWithSetModifications.Fragment(DissociationType.IRMPD, FragmentationTerminus.Both, fragments);
            allFragmentIonMzs = new HashSet<int>(fragments.Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { 98, 227, 120, 249 }));

            aPeptideWithSetModifications.Fragment(DissociationType.ECD, FragmentationTerminus.Both, fragments);
            allFragmentIonMzs = new HashSet<int>(fragments.Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { 115, 244, 120, 249, 104, 233 }));

            aPeptideWithSetModifications.Fragment(DissociationType.PQD, FragmentationTerminus.Both, fragments);
            allFragmentIonMzs = new HashSet<int>(fragments.Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { }));

            aPeptideWithSetModifications.Fragment(DissociationType.ETD, FragmentationTerminus.Both, fragments);
            allFragmentIonMzs = new HashSet<int>(fragments.Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { 115, 244, 120, 249, 104, 233 }));

            aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both, fragments);
            allFragmentIonMzs = new HashSet<int>(fragments.Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { 98, 227, 120, 249 }));

            aPeptideWithSetModifications.Fragment(DissociationType.AnyActivationType, FragmentationTerminus.Both, fragments);
            allFragmentIonMzs = new HashSet<int>(fragments.Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { 98, 227, 120, 249 }));

            aPeptideWithSetModifications.Fragment(DissociationType.EThcD, FragmentationTerminus.Both, fragments);
            allFragmentIonMzs = new HashSet<int>(fragments.Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { 98, 227, 115, 244, 120, 249, 104, 233 }));

            DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom] = new List<ProductType> { };
            aPeptideWithSetModifications.Fragment(DissociationType.Custom, FragmentationTerminus.Both, fragments);
            allFragmentIonMzs = new HashSet<int>(fragments.Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { }));

            aPeptideWithSetModifications.Fragment(DissociationType.ISCID, FragmentationTerminus.Both, fragments);
            allFragmentIonMzs = new HashSet<int>(fragments.Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new HashSet<int> { }));
        }
    }
}