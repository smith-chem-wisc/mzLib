using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using MassSpectrometry;
using MassSpectrometry.ExperimentalDesign;
using NUnit.Framework;

namespace Test.Omics
{
    /// <summary>
    /// Tests for ISampleInfo implementations (SpectraFileInfo and IsobaricQuantSampleInfo).
    /// Validates equality, hashing, and collection behavior for quantification sample tracking.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class CrossClassComparisonTests
    {
        #region SpectraFileInfo Tests

        /// <summary>
        /// Verifies SpectraFileInfo equality is based on all properties.
        /// Critical: Used as Dictionary keys for intensity mapping in quantification.
        /// </summary>
        [Test]
        public void SpectraFileInfo_Equality_BasedOnAllProperties()
        {
            var sample1 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var sample2 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var different = new SpectraFileInfo(@"C:\Data\other.raw", "Treatment", 2, 2, 1);

            Assert.Multiple(() =>
            {
                Assert.That(sample1.Equals(sample2), Is.True);
                Assert.That(sample1.Equals(different), Is.False);
                Assert.That(sample1.Equals(null), Is.False);
                Assert.That(sample1.GetHashCode(), Is.EqualTo(sample2.GetHashCode()));
                Assert.That(sample1.GetHashCode(), Is.Not.EqualTo(different.GetHashCode()));
            });
        }

        /// <summary>
        /// Verifies SpectraFileInfo works correctly as Dictionary key.
        /// Critical: IntensitiesBySample uses ISampleInfo as key.
        /// </summary>
        [Test]
        public void SpectraFileInfo_WorksAsDictionaryKey()
        {
            var sample1 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var sample2 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);

            var dict = new Dictionary<ISampleInfo, double> { { sample1, 1000.0 } };

            Assert.That(dict.TryGetValue(sample2, out var value), Is.True);
            Assert.That(value, Is.EqualTo(1000.0));
        }

        #endregion

        #region IsobaricQuantSampleInfo Tests

        /// <summary>
        /// Verifies IsobaricQuantSampleInfo equality is based ONLY on FilePath and ChannelLabel.
        /// Critical: Other properties (Condition, BioRep, etc.) are ignored for deduplication.
        /// This enables the same channel to be referenced with different metadata.
        /// </summary>
        [Test]
        public void IsobaricQuantSampleInfo_Equality_BasedOnFilePathAndChannelOnly()
        {
            var sample1 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);
            var sameChannel = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Treatment", 99, 99, 99, 99, "126", 999.999, true);
            var differentChannel = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "127N", 127.124761, false);
            var differentFile = new IsobaricQuantSampleInfo(
                @"C:\Data\other.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);

            Assert.Multiple(() =>
            {
                // Same FilePath + ChannelLabel = equal, regardless of other properties
                Assert.That(sample1.Equals(sameChannel), Is.True);
                Assert.That(sample1.GetHashCode(), Is.EqualTo(sameChannel.GetHashCode()));

                // Different ChannelLabel = not equal
                Assert.That(sample1.Equals(differentChannel), Is.False);

                // Different FilePath = not equal
                Assert.That(sample1.Equals(differentFile), Is.False);

                // Null handling
                Assert.That(sample1.Equals(null), Is.False);
            });
        }

        /// <summary>
        /// Verifies IsobaricQuantSampleInfo deduplicates correctly in HashSet.
        /// Critical: Prevents duplicate channel entries in quantification results.
        /// </summary>
        [Test]
        public void IsobaricQuantSampleInfo_HashSet_DeduplicatesByFilePathAndChannel()
        {
            var sample1 = new IsobaricQuantSampleInfo(@"C:\a.raw", "A", 1, 1, 0, 1, "126", 126.0, false);
            var duplicate = new IsobaricQuantSampleInfo(@"C:\a.raw", "B", 2, 2, 1, 99, "126", 127.0, true);
            var differentChannel = new IsobaricQuantSampleInfo(@"C:\a.raw", "A", 1, 1, 0, 1, "127N", 127.0, false);

            var set = new HashSet<IsobaricQuantSampleInfo> { sample1, duplicate, differentChannel };

            // sample1 and duplicate have same FilePath+Channel, so only 2 unique entries
            Assert.That(set.Count, Is.EqualTo(2));
        }

        #endregion

        #region Cross-Type Tests

        /// <summary>
        /// Verifies different ISampleInfo types are never equal, even with same properties.
        /// Critical: Prevents accidental mixing of label-free and isobaric samples.
        /// </summary>
        [Test]
        public void CrossType_DifferentTypesNeverEqual()
        {
            var spectra = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var isobaric = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);

            Assert.Multiple(() =>
            {
                Assert.That(spectra.Equals(isobaric), Is.False);
                Assert.That(isobaric.Equals(spectra), Is.False);
                Assert.That(spectra.Equals((ISampleInfo)isobaric), Is.False);
                Assert.That(isobaric.Equals((ISampleInfo)spectra), Is.False);
            });
        }

        /// <summary>
        /// Verifies mixed ISampleInfo types work correctly in same collection.
        /// Critical: BioPolymerGroup.IntensitiesBySample can contain both types.
        /// </summary>
        [Test]
        public void CrossType_MixedTypesInCollection()
        {
            var spectra1 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var spectra2 = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0); // Duplicate
            var isobaric1 = new IsobaricQuantSampleInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);
            var isobaric2 = new IsobaricQuantSampleInfo(@"C:\Data\sample.raw", "Treatment", 99, 99, 99, 99, "126", 999.0, true); // Duplicate

            var set = new HashSet<ISampleInfo> { spectra1, spectra2, isobaric1, isobaric2 };

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
        public static void CompareTo_IsobaricQuantSampleInfo_SameFilePath_DifferentChannel_ReturnsPositive()
        {
            // CompareTo does NOT compare by ChannelLabel - it only uses Equals check
            // When FilePath and all other comparison properties are the same but ChannelLabel differs,
            // Equals returns false and comparison falls through to return 1
            var isobaric126 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.127726, false);
            var isobaric127 = new IsobaricQuantSampleInfo(
                @"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "127N", 127.124761, false);

            // Both return positive (1) because they're not equal and all comparison properties match
            Assert.That(isobaric126.CompareTo(isobaric127), Is.GreaterThan(0));
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
            Assert.That(set.Contains(spectra1), Is.True);
            Assert.That(set.Contains(isobaric1), Is.True);
        }

        /// <summary>
        /// Verifies mixed types work as Dictionary keys without collision.
        /// Critical: Ensures correct intensity retrieval for both quantification methods.
        /// </summary>
        [Test]
        public void CrossType_MixedTypesAsDictionaryKeys()
        {
            var spectra = new SpectraFileInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0);
            var isobaric = new IsobaricQuantSampleInfo(@"C:\Data\sample.raw", "Control", 1, 1, 0, 1, "126", 126.0, false);

            var dict = new Dictionary<ISampleInfo, double>
            {
                { spectra, 1000.0 },
                { isobaric, 2000.0 }
            };

            Assert.Multiple(() =>
            {
                Assert.That(dict.Count, Is.EqualTo(2));
                Assert.That(dict[spectra], Is.EqualTo(1000.0));
                Assert.That(dict[isobaric], Is.EqualTo(2000.0));
            });
        }

        #endregion
    }
}