using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestResultFile
    {
        [Test]
        public static void TestEquals()
        {
            string path1 =
                @"FileReadingTests\ExternalFileTypes\Ms1Feature_TopFDv1.6.2_ms1.feature";
            string path2 =
                @"FileReadingTests\ExternalFileTypes\Ms1Feature_FlashDeconvOpenMs3.0.0_ms1.feature";
            string ms2FeaturePath =
                @"FileReadingTests\ExternalFileTypes\Ms2Feature_TopFDv1.6.2_ms2.feature";

            Ms1FeatureFile ms1Features = new Ms1FeatureFile(path1);
            Assert.That(ms1Features.Software == Software.TopFD);

            Ms1FeatureFile ms1FeatureCopy = FileReader.ReadFile<Ms1FeatureFile>(path1);
            Assert.That(ms1FeatureCopy.Software, Is.EqualTo(Software.TopFD));

            Ms2FeatureFile ms2Features = new Ms2FeatureFile(ms2FeaturePath);

            Assert.That(ms1Features.Equals(ms1FeatureCopy));
            Assert.That(ms1Features.Equals(ms1Features));
            Assert.That(!ms1Features.Equals(null));
            Assert.That(!ms1Features.Equals(ms2Features));
            Assert.That(ms1Features.Equals((object)ms1FeatureCopy));
            Assert.That(ms1Features.Equals((object)ms1Features));
            Assert.That(!ms1Features.Equals((object)null));
            Assert.That(!ms1Features.Equals((object)ms2Features));

            Ms1FeatureFile ms1Features2 = new Ms1FeatureFile(path2);
            Assert.That(ms1Features2.Software, Is.EqualTo(Software.FLASHDeconv));

            Ms1FeatureFile ms1FeaturesCopy2 = FileReader.ReadFile<Ms1FeatureFile>(path2);
            Assert.That(ms1FeaturesCopy2.Software, Is.EqualTo(Software.FLASHDeconv));

            
            Assert.That(ms1Features2.Equals(ms1FeaturesCopy2));
            Assert.That(ms1Features2.Equals((object)ms1FeaturesCopy2));
        }

        [Test]
        public static void TestGetEnumerator()
        {
            string path1 =
                @"FileReadingTests\ExternalFileTypes\Ms1Feature_TopFDv1.6.2_ms1.feature";

            Ms1FeatureFile ms1Features = new Ms1FeatureFile(path1);

            foreach (var feature in ms1Features)
            {
                Assert.That(feature.GetType() == typeof(Ms1Feature));
            }
        }

        [Test]
        public static void TestOperatorOverloads()
        {
            string path1 = @"FileReadingTests\ExternalFileTypes\Ms1Feature_TopFDv1.6.2_ms1.feature";
            string path2 = @"FileReadingTests\ExternalFileTypes\Ms1Feature_FlashDeconvOpenMs3.0.0_ms1.feature";

            Ms1FeatureFile fileToTest = new Ms1FeatureFile(path1);
            Ms1FeatureFile ms1Features2 = new Ms1FeatureFile(path2);

            Ms1Feature feature = new Ms1Feature();
            List<Ms1Feature> features = new List<Ms1Feature> { new Ms1Feature() };


            var initialCount = fileToTest.Count();
            var file2Count = ms1Features2.Count();

            // Test operator +
            var result1 = fileToTest + feature;
            Assert.That(result1.Results.Contains(feature));
            Assert.That(fileToTest.Count, Is.EqualTo(initialCount + 1));

            var result2 = fileToTest + features;
            Assert.That(result2.Results.Contains(feature));
            Assert.That(fileToTest.Count, Is.EqualTo(initialCount + 2));

            var result3 = fileToTest + ms1Features2;
            Assert.That(result3.Results.Contains(feature));
            Assert.That(fileToTest.Count, Is.EqualTo(initialCount + 2 + file2Count));

            // Test operator -
            var result4 = fileToTest - feature;
            Assert.That(!result4.Results.Contains(feature));
            Assert.That(fileToTest.Count, Is.EqualTo(initialCount + 1 + file2Count));

            var result5 = fileToTest - features;
            Assert.That(!result5.Results.Contains(feature));
            Assert.That(fileToTest.Count, Is.EqualTo(initialCount + file2Count));

            var result6 = fileToTest - ms1Features2;
            Assert.That(!result6.Results.Contains(feature));
            Assert.That(fileToTest.Count, Is.EqualTo(initialCount));
        }
    }
}
