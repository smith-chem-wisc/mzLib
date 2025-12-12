using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using MassSpectrometry;
using NUnit.Framework;

namespace Test.Omics
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class CrossClassComparisonTests
    {
        #region SpectraFileInfo Equality Tests

        [Test]
        public static void SpectraFileInfo_Equality_AllCastingVariants()
        {
            var spectra1 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var spectra2 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);

            // Same objects - typed
            Assert.That(spectra1.Equals(spectra2));
            Assert.That(spectra1.Equals((SpectraFileInfo)spectra2));
            Assert.That(spectra1.Equals((ISampleInfo)spectra2));
            Assert.That(spectra1.Equals((object)spectra2));

            // Same reference
            Assert.That(spectra1.Equals(spectra1));
            Assert.That(spectra1.Equals((SpectraFileInfo)spectra1));
            Assert.That(spectra1.Equals((ISampleInfo)spectra1));
            Assert.That(spectra1.Equals((object)spectra1));

            // Hash codes match
            Assert.That(spectra1.GetHashCode(), Is.EqualTo(spectra2.GetHashCode()));

            // Null checks - all fail on null
            Assert.That(!spectra1.Equals(null));
            Assert.That(!spectra1.Equals((SpectraFileInfo?)null));
            Assert.That(!spectra1.Equals((ISampleInfo?)null));
            Assert.That(!spectra1.Equals((object?)null));
        }

        [Test]
        public static void SpectraFileInfo_Inequality_DifferentValues_AllCastingVariants()
        {
            var spectra1 = new SpectraFileInfo(@"C:\Data\sample1.raw", "Control", 1, 1, 0);
            var spectra2 = new SpectraFileInfo(@"C:\Data\sample2.raw", "Treatment", 2, 2, 1);

            // Different objects - typed
            Assert.That(!spectra1.Equals(spectra2));
            Assert.That(!spectra1.Equals((SpectraFileInfo)spectra2));
            Assert.That(!spectra1.Equals((ISampleInfo)spectra2));
            Assert.That(!spectra1.Equals((object)spectra2));

            // Hash codes differ
            Assert.That(spectra1.GetHashCode(), Is.Not.EqualTo(spectra2.GetHashCode()));
        }

        [Test]
        public static void SpectraFileInfo_Inequality_SingleFieldDiffers_AllCastingVariants()
        {
            var spectra1 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);

            // Different file path only
            var differentPath = new SpectraFileInfo(@"C:\Data\other.raw", "Control", 1, 1, 0);
            Assert.That(!spectra1.Equals(differentPath));
            Assert.That(!spectra1.Equals((SpectraFileInfo)differentPath));
            Assert.That(!spectra1.Equals((ISampleInfo)differentPath));
            Assert.That(!spectra1.Equals((object)differentPath));

            // Different condition only
            var differentCondition = new SpectraFileInfo(@"C:\Data\sample.raw", "Treatment", 1, 1, 0);
            Assert.That(!spectra1.Equals(differentCondition));
            Assert.That(!spectra1.Equals((SpectraFileInfo)differentCondition));
            Assert.That(!spectra1.Equals((ISampleInfo)differentCondition));
            Assert.That(!spectra1.Equals((object)differentCondition));

            // Different biological replicate only
            var differentBioRep = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 2, 1, 0);
            Assert.That(!spectra1.Equals(differentBioRep));
            Assert.That(!spectra1.Equals((SpectraFileInfo)differentBioRep));
            Assert.That(!spectra1.Equals((ISampleInfo)differentBioRep));
            Assert.That(!spectra1.Equals((object)differentBioRep));

            // Different technical replicate only
            var differentTechRep = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 2, 0);
            Assert.That(!spectra1.Equals(differentTechRep));
            Assert.That(!spectra1.Equals((SpectraFileInfo)differentTechRep));
            Assert.That(!spectra1.Equals((ISampleInfo)differentTechRep));
            Assert.That(!spectra1.Equals((object)differentTechRep));

            // Different fraction only
            var differentFraction = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 1);
            Assert.That(!spectra1.Equals(differentFraction));
            Assert.That(!spectra1.Equals((SpectraFileInfo)differentFraction));
            Assert.That(!spectra1.Equals((ISampleInfo)differentFraction));
            Assert.That(!spectra1.Equals((object)differentFraction));
        }

        #endregion

        #region IsobaricQuantSampleInfo Equality Tests

        [Test]
        public static void IsobaricQuantSampleInfo_Equality_AllCastingVariants()
        {
            var isobaric1 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);
            var isobaric2 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);

            // Same objects - typed
            Assert.That(isobaric1.Equals(isobaric2));
            Assert.That(isobaric1.Equals((IsobaricQuantSampleInfo)isobaric2));
            Assert.That(isobaric1.Equals((ISampleInfo)isobaric2));
            Assert.That(isobaric1.Equals((object)isobaric2));

            // Same reference
            Assert.That(isobaric1.Equals(isobaric1));
            Assert.That(isobaric1.Equals((IsobaricQuantSampleInfo)isobaric1));
            Assert.That(isobaric1.Equals((ISampleInfo)isobaric1));
            Assert.That(isobaric1.Equals((object)isobaric1));

            // Hash codes match
            Assert.That(isobaric1.GetHashCode(), Is.EqualTo(isobaric2.GetHashCode()));

            // Null checks - all fail on null
            Assert.That(!isobaric1.Equals(null));
            Assert.That(!isobaric1.Equals((IsobaricQuantSampleInfo?)null));
            Assert.That(!isobaric1.Equals((ISampleInfo?)null));
            Assert.That(!isobaric1.Equals((object?)null));
        }

        [Test]
        public static void IsobaricQuantSampleInfo_Equality_SameFilePathAndChannelLabel_DifferentOtherProperties()
        {
            // IsobaricQuantSampleInfo equality is based ONLY on FullFilePathWithExtension and ChannelLabel
            var isobaric1 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);
            var isobaric2 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Treatment", 99, 99, 99, 99, "126", 999.999, true);

            // Same FullFilePathWithExtension and ChannelLabel = equal, regardless of other properties
            Assert.That(isobaric1.Equals(isobaric2));
            Assert.That(isobaric1.Equals((IsobaricQuantSampleInfo)isobaric2));
            Assert.That(isobaric1.Equals((ISampleInfo)isobaric2));
            Assert.That(isobaric1.Equals((object)isobaric2));

            // Hash codes match
            Assert.That(isobaric1.GetHashCode(), Is.EqualTo(isobaric2.GetHashCode()));
        }

        [Test]
        public static void IsobaricQuantSampleInfo_Equality_DifferentPlexId_SameFilePathAndChannel_ReturnsTrue()
        {
            // PlexId is NOT part of equality - only FullFilePathWithExtension and ChannelLabel
            var isobaric1 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);
            var isobaric2 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 99, "126", 126.127726, false);

            // Same FullFilePathWithExtension and ChannelLabel = equal
            Assert.That(isobaric1.Equals(isobaric2));
            Assert.That(isobaric1.Equals((IsobaricQuantSampleInfo)isobaric2));
            Assert.That(isobaric1.Equals((ISampleInfo)isobaric2));
            Assert.That(isobaric1.Equals((object)isobaric2));

            // Hash codes match
            Assert.That(isobaric1.GetHashCode(), Is.EqualTo(isobaric2.GetHashCode()));
        }

        [Test]
        public static void IsobaricQuantSampleInfo_Inequality_DifferentFilePath_AllCastingVariants()
        {
            var isobaric1 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample1.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);
            var isobaric2 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample2.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);

            // Different FullFilePathWithExtension = not equal
            Assert.That(!isobaric1.Equals(isobaric2));
            Assert.That(!isobaric1.Equals((IsobaricQuantSampleInfo)isobaric2));
            Assert.That(!isobaric1.Equals((ISampleInfo)isobaric2));
            Assert.That(!isobaric1.Equals((object)isobaric2));

            // Hash codes differ
            Assert.That(isobaric1.GetHashCode(), Is.Not.EqualTo(isobaric2.GetHashCode()));
        }

        [Test]
        public static void IsobaricQuantSampleInfo_Inequality_DifferentChannelLabel_AllCastingVariants()
        {
            var isobaric1 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);
            var isobaric2 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "127N", 127.124761, false);

            // Different ChannelLabel = not equal
            Assert.That(!isobaric1.Equals(isobaric2));
            Assert.That(!isobaric1.Equals((IsobaricQuantSampleInfo)isobaric2));
            Assert.That(!isobaric1.Equals((ISampleInfo)isobaric2));
            Assert.That(!isobaric1.Equals((object)isobaric2));

            // Hash codes differ
            Assert.That(isobaric1.GetHashCode(), Is.Not.EqualTo(isobaric2.GetHashCode()));
        }

        [Test]
        public static void IsobaricQuantSampleInfo_Inequality_ChannelLabel_CaseSensitive()
        {
            var isobaric1 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false);
            var isobaric2 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "127n", 127.0, false);

            // Case-sensitive ChannelLabel comparison
            Assert.That(!isobaric1.Equals(isobaric2));
            Assert.That(!isobaric1.Equals((IsobaricQuantSampleInfo)isobaric2));
            Assert.That(!isobaric1.Equals((ISampleInfo)isobaric2));
            Assert.That(!isobaric1.Equals((object)isobaric2));
        }

        [Test]
        public static void IsobaricQuantSampleInfo_Inequality_FilePath_CaseSensitive()
        {
            var isobaric1 = new IsobaricQuantSampleInfo(
                @"C:\Data\Sample.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var isobaric2 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            // Case-sensitive FilePath comparison
            Assert.That(!isobaric1.Equals(isobaric2));
            Assert.That(!isobaric1.Equals((IsobaricQuantSampleInfo)isobaric2));
            Assert.That(!isobaric1.Equals((ISampleInfo)isobaric2));
            Assert.That(!isobaric1.Equals((object)isobaric2));
        }

        #endregion

        #region Cross-Type Equality Tests (SpectraFileInfo vs IsobaricQuantSampleInfo)

        [Test]
        public static void CrossType_SpectraFileInfo_NotEqualTo_IsobaricQuantSampleInfo_AllCastingVariants()
        {
            var spectra = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var isobaric = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);

            // SpectraFileInfo comparing to IsobaricQuantSampleInfo - never equal
            Assert.That(!spectra.Equals(isobaric));
            Assert.That(!spectra.Equals((ISampleInfo)isobaric));
            Assert.That(!spectra.Equals((object)isobaric));

            // IsobaricQuantSampleInfo comparing to SpectraFileInfo - never equal
            Assert.That(!isobaric.Equals(spectra));
            Assert.That(!isobaric.Equals((ISampleInfo)spectra));
            Assert.That(!isobaric.Equals((object)spectra));

            // Hash codes should differ (different types)
            Assert.That(spectra.GetHashCode(), Is.Not.EqualTo(isobaric.GetHashCode()));
        }

        [Test]
        public static void CrossType_InHashSet_DifferentTypes_TreatedAsDifferent()
        {
            var spectra = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var isobaric = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);

            // Both can be added to a HashSet<ISampleInfo>
            var set = new HashSet<ISampleInfo> { spectra, isobaric };
            Assert.That(set.Count, Is.EqualTo(2));
            Assert.That(set.Contains(spectra));
            Assert.That(set.Contains(isobaric));
        }

        [Test]
        public static void CrossType_InDictionary_DifferentTypes_TreatedAsDifferent()
        {
            var spectra = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var isobaric = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);

            // Both can be used as keys in a Dictionary<ISampleInfo, string>
            var dict = new Dictionary<ISampleInfo, string>
            {
                { spectra, "LabelFree" },
                { isobaric, "Isobaric" }
            };

            Assert.That(dict.Count, Is.EqualTo(2));
            Assert.That(dict[spectra], Is.EqualTo("LabelFree"));
            Assert.That(dict[isobaric], Is.EqualTo("Isobaric"));
        }

        [Test]
        public static void CrossType_Equals_WrongType_ReturnsFalse()
        {
            var spectra = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var isobaric = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);

            // Comparing to completely unrelated types
            Assert.That(!spectra.Equals("a string"));
            Assert.That(!spectra.Equals(123));
            Assert.That(!spectra.Equals(new object()));

            Assert.That(!isobaric.Equals("a string"));
            Assert.That(!isobaric.Equals(123));
            Assert.That(!isobaric.Equals(new object()));
        }

        #endregion

        #region Equality Operator Tests

        [Test]
        public static void IsobaricQuantSampleInfo_EqualityOperators_AllCastingVariants()
        {
            var isobaric1 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);
            var isobaric2 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);
            var isobaric3 = new IsobaricQuantSampleInfo(
                @"C:\Data\other.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false);

            // Equality operator
            Assert.That(isobaric1 == isobaric2, Is.True);
            Assert.That(isobaric1 == isobaric3, Is.False);

            // Inequality operator
            Assert.That(isobaric1 != isobaric2, Is.False);
            Assert.That(isobaric1 != isobaric3, Is.True);

            // Null handling
            IsobaricQuantSampleInfo? nullSample = null;
            Assert.That(isobaric1 == nullSample, Is.False);
            Assert.That(nullSample == isobaric1, Is.False);
            Assert.That(nullSample == nullSample, Is.True);

            Assert.That(isobaric1 != nullSample, Is.True);
            Assert.That(nullSample != isobaric1, Is.True);
            Assert.That(nullSample != nullSample, Is.False);
        }

        [Test]
        public static void IsobaricQuantSampleInfo_EqualityOperators_DifferentPlexId_SameFilePathAndChannel()
        {
            var isobaric1 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);
            var isobaric2 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 99, "126", 126.127726, false);

            // Same FilePath and ChannelLabel, different PlexId - should be equal
            Assert.That(isobaric1 == isobaric2, Is.True);
            Assert.That(isobaric1 != isobaric2, Is.False);
        }

        #endregion

        #region CompareTo Cross-Type Tests

        [Test]
        public static void CompareTo_SpectraFileInfo_ToIsobaricQuantSampleInfo_CrossTypeOrdering()
        {
            var spectra = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var isobaric = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);

            // SpectraFileInfo.CompareTo returns positive when comparing to non-SpectraFileInfo types
            // This means SpectraFileInfo sorts AFTER IsobaricQuantSampleInfo
            int spectraToIsobaric = spectra.CompareTo(isobaric);
            Assert.That(spectraToIsobaric, Is.GreaterThan(0));

            // IsobaricQuantSampleInfo.CompareTo returns positive when comparing to non-IsobaricQuantSampleInfo types
            // This means IsobaricQuantSampleInfo also sorts AFTER SpectraFileInfo
            int isobaricToSpectra = isobaric.CompareTo(spectra);
            Assert.That(isobaricToSpectra, Is.GreaterThan(0));

            // Note: This creates an inconsistent ordering between the two types.
            // When sorting a mixed collection, the behavior depends on which type's CompareTo is called.
            // This is acceptable since mixed-type collections are not a typical use case.
        }

        [Test]
        public static void CompareTo_SpectraFileInfo_ToIsobaricQuantSampleInfo_DifferentProperties()
        {
            var spectra = new SpectraFileInfo(@"A:\Data\sample.raw", "Alpha", 1, 1, 0);
            var isobaric = new IsobaricQuantSampleInfo(
                @"B:\Data\sample.raw", "Beta", 2, 2, 1, 1, "126", 126.127726, false);

            // Cross-type comparison always returns positive (sorts after the other type)
            int spectraToIsobaric = spectra.CompareTo(isobaric);
            Assert.That(spectraToIsobaric, Is.GreaterThan(0));

            int isobaricToSpectra = isobaric.CompareTo(spectra);
            Assert.That(isobaricToSpectra, Is.GreaterThan(0));
        }

        [Test]
        public static void CompareTo_Null_AllTypes()
        {
            var spectra = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var isobaric = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);

            // Both types return negative when comparing to null (non-null comes before null)
            Assert.That(spectra.CompareTo(null), Is.LessThan(0));
            Assert.That(isobaric.CompareTo(null), Is.LessThan(0));
        }

        [Test]
        public static void CompareTo_IsobaricQuantSampleInfo_SameFilePathAndChannel_DifferentPlexId_ReturnsZero()
        {
            var isobaric1 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);
            var isobaric2 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 99, "126", 126.127726, false);

            // Same FilePath and ChannelLabel - CompareTo returns 0
            Assert.That(isobaric1.CompareTo(isobaric2), Is.EqualTo(0));
        }

        [Test]
        public static void CompareTo_IsobaricQuantSampleInfo_DifferentFilePath_OrdersByFilePath()
        {
            var isobaricA = new IsobaricQuantSampleInfo(
                @"C:\a.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);
            var isobaricB = new IsobaricQuantSampleInfo(
                @"C:\b.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);

            // Different FilePath - orders by file path
            Assert.That(isobaricA.CompareTo(isobaricB), Is.LessThan(0));
            Assert.That(isobaricB.CompareTo(isobaricA), Is.GreaterThan(0));
        }

        [Test]
        public static void CompareTo_IsobaricQuantSampleInfo_SameFilePath_DifferentChannel_OrdersByChannel()
        {
            var isobaric126 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);
            var isobaric127 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "127N", 127.124761, false);

            // Same FilePath, different ChannelLabel - orders by channel
            Assert.That(isobaric126.CompareTo(isobaric127), Is.LessThan(0));
            Assert.That(isobaric127.CompareTo(isobaric126), Is.GreaterThan(0));
        }

        #endregion

        #region ISampleInfo Interface Equality Tests

        [Test]
        public static void ISampleInfo_InterfaceEquality_SpectraFileInfo()
        {
            ISampleInfo sample1 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            ISampleInfo sample2 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            ISampleInfo sample3 = new SpectraFileInfo(@"C:\Data\different.raw", "Treatment", 2, 2, 1);

            // Equality through interface
            Assert.That(sample1.Equals(sample2));
            Assert.That(sample1.Equals((object)sample2));
            Assert.That(!sample1.Equals(sample3));
            Assert.That(!sample1.Equals((object)sample3));

            // Hash codes
            Assert.That(sample1.GetHashCode(), Is.EqualTo(sample2.GetHashCode()));
            Assert.That(sample1.GetHashCode(), Is.Not.EqualTo(sample3.GetHashCode()));
        }

        [Test]
        public static void ISampleInfo_InterfaceEquality_IsobaricQuantSampleInfo()
        {
            ISampleInfo sample1 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            ISampleInfo sample2 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Treatment", 99, 99, 99, 99, "126", 999.0, true);
            ISampleInfo sample3 = new IsobaricQuantSampleInfo(
                @"C:\Data\different.raw", "Control", 1, 1, 0, 1, "127N", 127.0, false);

            // Equality through interface (only FullFilePathWithExtension and ChannelLabel matter)
            Assert.That(sample1.Equals(sample2)); // Same FilePath and ChannelLabel
            Assert.That(sample1.Equals((object)sample2));
            Assert.That(!sample1.Equals(sample3)); // Different FilePath or ChannelLabel
            Assert.That(!sample1.Equals((object)sample3));

            // Hash codes
            Assert.That(sample1.GetHashCode(), Is.EqualTo(sample2.GetHashCode()));
            Assert.That(sample1.GetHashCode(), Is.Not.EqualTo(sample3.GetHashCode()));
        }

        [Test]
        public static void ISampleInfo_InterfaceEquality_MixedTypes()
        {
            ISampleInfo spectra = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            ISampleInfo isobaric = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            // Different types are never equal, even through interface
            Assert.That(!spectra.Equals(isobaric));
            Assert.That(!spectra.Equals((object)isobaric));
            Assert.That(!isobaric.Equals(spectra));
            Assert.That(!isobaric.Equals((object)spectra));
        }

        #endregion

        #region Collection Behavior Tests

        [Test]
        public static void HashSet_SpectraFileInfo_DeduplicatesCorrectly()
        {
            var sample1 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var sample2 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var sample3 = new SpectraFileInfo(@"C:\Data\other.raw", "Control", 1, 1, 0);

            var set = new HashSet<SpectraFileInfo> { sample1, sample2, sample3 };

            Assert.That(set.Count, Is.EqualTo(2)); // sample1 and sample2 are duplicates
            Assert.That(set.Contains(sample1));
            Assert.That(set.Contains(sample2));
            Assert.That(set.Contains(sample3));
        }

        [Test]
        public static void HashSet_IsobaricQuantSampleInfo_DeduplicatesCorrectly()
        {
            var sample1 = new IsobaricQuantSampleInfo(@"C:\a.raw", "A", 1, 1, 0, 1, "126", 126.0, false);
            var sample2 = new IsobaricQuantSampleInfo(@"C:\a.raw", "B", 2, 2, 1, 99, "126", 127.0, true);
            var sample3 = new IsobaricQuantSampleInfo(@"C:\a.raw", "A", 1, 1, 0, 1, "127N", 127.0, false);

            var set = new HashSet<IsobaricQuantSampleInfo> { sample1, sample2, sample3 };

            // sample1 and sample2 have same FullFilePathWithExtension and ChannelLabel, so they're duplicates
            Assert.That(set.Count, Is.EqualTo(2));
        }

        [Test]
        public static void Dictionary_SpectraFileInfo_WorksAsKey()
        {
            var sample1 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var sample2 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);

            var dict = new Dictionary<SpectraFileInfo, string>
            {
                { sample1, "Value1" }
            };

            // sample2 should find the same entry as sample1
            Assert.That(dict.TryGetValue(sample2, out var value));
            Assert.That(value, Is.EqualTo("Value1"));
        }

        [Test]
        public static void Dictionary_IsobaricQuantSampleInfo_WorksAsKey()
        {
            var sample1 = new IsobaricQuantSampleInfo(@"C:\a.raw", "A", 1, 1, 0, 1, "126", 126.0, false);
            var sample2 = new IsobaricQuantSampleInfo(@"C:\a.raw", "B", 2, 2, 1, 99, "126", 127.0, true);

            var dict = new Dictionary<IsobaricQuantSampleInfo, string>
            {
                { sample1, "Value1" }
            };

            // sample2 has same FullFilePathWithExtension and ChannelLabel, should find same entry
            Assert.That(dict.TryGetValue(sample2, out var value));
            Assert.That(value, Is.EqualTo("Value1"));
        }

        [Test]
        public static void HashSet_ISampleInfo_MixedTypes()
        {
            var spectra1 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var spectra2 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var isobaric1 = new IsobaricQuantSampleInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var isobaric2 = new IsobaricQuantSampleInfo(@"C:\Data\sample.raw", "Treatment", 99, 99, 99, 99, "126", 999.0, true);

            var set = new HashSet<ISampleInfo> { spectra1, spectra2, isobaric1, isobaric2 };

            // spectra1 and spectra2 are duplicates
            // isobaric1 and isobaric2 are duplicates (same FullFilePathWithExtension and ChannelLabel)
            // spectra and isobaric are different types
            Assert.That(set.Count, Is.EqualTo(2));
        }

        #endregion
    }
}