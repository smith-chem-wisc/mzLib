using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public sealed class TestProductMassesMightHaveDuplicates
    {
        [Test]
        public static void Test_UnmodifiedPeptide_AllProductType_fragmentMasses()
        {
            Protein p = new Protein("PET", "accession");

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 2);

            var aPeptideWithSetModifications = p.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();
            var aCompactPeptide = aPeptideWithSetModifications.CompactPeptide(FragmentationTerminus.Both, DissociationType.AnyActivationType);

            var allFragmentIonMzs = aCompactPeptide.ProductMassesMightHaveDuplicatesAndNaNsNew(DissociationType.Unknown).Select(i => (int)Math.Round(i.ToMz(1))).ToArray();
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { }));
            allFragmentIonMzs = aCompactPeptide.ProductMassesMightHaveDuplicatesAndNaNsNew(DissociationType.CID).Select(i => (int)Math.Round(i.ToMz(1))).ToArray();
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { 98, 227, 120, 249 }));
            allFragmentIonMzs = aCompactPeptide.ProductMassesMightHaveDuplicatesAndNaNsNew(DissociationType.MPD).Select(i => (int)Math.Round(i.ToMz(1))).ToArray();
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { }));
            allFragmentIonMzs = aCompactPeptide.ProductMassesMightHaveDuplicatesAndNaNsNew(DissociationType.ECD).Select(i => (int)Math.Round(i.ToMz(1))).ToArray();
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { 115, 244, 120, 103, 249, 232 }));
            allFragmentIonMzs = aCompactPeptide.ProductMassesMightHaveDuplicatesAndNaNsNew(DissociationType.PQD).Select(i => (int)Math.Round(i.ToMz(1))).ToArray();
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { }));
            allFragmentIonMzs = aCompactPeptide.ProductMassesMightHaveDuplicatesAndNaNsNew(DissociationType.ETD).Select(i => (int)Math.Round(i.ToMz(1))).ToArray();
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { 115, 244, 120, 103, 249, 232 }));
            allFragmentIonMzs = aCompactPeptide.ProductMassesMightHaveDuplicatesAndNaNsNew(DissociationType.HCD).Select(i => (int)Math.Round(i.ToMz(1))).ToArray();
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { 98, 227, 120, 249 }));
            allFragmentIonMzs = aCompactPeptide.ProductMassesMightHaveDuplicatesAndNaNsNew(DissociationType.AnyActivationType).Select(i => (int)Math.Round(i.ToMz(1))).ToArray();
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { 98, 227, 120, 249 }));
            allFragmentIonMzs = aCompactPeptide.ProductMassesMightHaveDuplicatesAndNaNsNew(DissociationType.EThCD).Select(i => (int)Math.Round(i.ToMz(1))).ToArray();
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { 98, 227, 120, 249 }));
            allFragmentIonMzs = aCompactPeptide.ProductMassesMightHaveDuplicatesAndNaNsNew(DissociationType.Custom).Select(i => (int)Math.Round(i.ToMz(1))).ToArray();
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { }));
            allFragmentIonMzs = aCompactPeptide.ProductMassesMightHaveDuplicatesAndNaNsNew(DissociationType.ISCID).Select(i => (int)Math.Round(i.ToMz(1))).ToArray();
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { }));
        }
    }
}