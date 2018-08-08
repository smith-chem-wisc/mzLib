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

            var allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.Unknown, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SetEquals(new int[] { }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.CID, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { 98, 227, 120, 249 }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.MPD, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.ECD, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { 115, 244, 120, 249, 105, 234 }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.PQD, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.ETD, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { 115, 244, 120, 249, 105, 234 }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.HCD, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { 98, 227, 120, 249 }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.AnyActivationType, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { 98, 227, 120, 249 }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.EThCD, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { 98, 227, 115, 244, 120, 249, 105, 234 }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.Custom, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { }));

            allFragmentIonMzs = new HashSet<int>(aPeptideWithSetModifications.Fragment(DissociationType.ISCID, FragmentationTerminus.Both).Select(i => (int)Math.Round(i.NeutralMass.ToMz(1))));
            Assert.IsTrue(allFragmentIonMzs.SequenceEqual(new int[] { }));
        }
    }
}