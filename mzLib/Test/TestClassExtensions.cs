using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.Diagnostics.CodeAnalysis;
using MzLibUtil;
using System;
using MassSpectrometry;
using System.Collections.Generic;
using Chemistry;
using Transcriptomics;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestClassExtensions
    {

        [Test]
        public static void TestGetClosestIndex()
        {
            double[] sortedArray = {1, 2, 3, 4.5, 5 };
            // First, make sure default binary search works how I think it does
            Assert.AreEqual(Array.BinarySearch(sortedArray, 3.1), ~3);

            // Then, the new extension
            Assert.AreEqual(sortedArray.GetClosestIndex(3.1), 2);
            Assert.AreEqual(sortedArray.GetClosestValue(3.1), 3);

            // What happens when the value is greater than the max
            Assert.AreEqual(sortedArray.GetClosestIndex(5.06), 4);
            // Less than the min?
            Assert.AreEqual(sortedArray.GetClosestIndex(-1), 0);

            // Test different search options
            Assert.AreEqual(sortedArray.GetClosestIndex(4.75, ArraySearchOption.Next), 4);
            Assert.AreEqual(sortedArray.GetClosestIndex(4.75, ArraySearchOption.Previous), 3);

            double[] smallerSortedArray = { 1, 2 };

            Assert.AreEqual(smallerSortedArray.GetClosestIndex( -1), 0);
            Assert.AreEqual(smallerSortedArray.GetClosestIndex(0), 0);
            Assert.AreEqual(smallerSortedArray.GetClosestIndex(1), 0);
            Assert.AreEqual(smallerSortedArray.GetClosestIndex(1.9), 1);
            Assert.AreEqual(smallerSortedArray.GetClosestIndex(Double.MaxValue), 1);

            double[] smallestSortedArray = { 1 };
            Assert.AreEqual(smallestSortedArray.GetClosestIndex(-1), 0);
            Assert.AreEqual(smallestSortedArray.GetClosestIndex(0), 0);
            Assert.AreEqual(smallestSortedArray.GetClosestIndex(1), 0);
            Assert.AreEqual(smallestSortedArray.GetClosestIndex(Double.MaxValue), 0);
        }

        [Test]
        public static void TestBoxCarSmooth()
        {
            double[] inputData = new double[] { 0.19, 0.69, 0.03, 0.85, 0.84, 0.46, 0.09, 0.05, 0.11, 0.5, 0.6, 0.78, 0.48, 0.66, 0.61, 0.78, 0.82, 0.18, 0.77, 0.14, 0.97, 0.48, 0.54, 0.98, 0.01, 0.38, 0.26, 0.4, 0.31, 0.41, 0.03, 0.2, 0.98, 0.36, 0.24, 0.51, 0.14, 0.96, 0.32, 0.9, 0.36, 0.57, 0.97, 0.07, 0.12, 0.73, 0.92, 0.51, 0.04, 0.2, 0.39, 0.32, 0.33, 0.62, 0.32, 0.68, 0.91, 0.3, 0.68, 0.22, 0.89, 0.27, 0.68, 0.08, 0.61, 0.25, 0.82, 0.73, 0.49, 0.76, 0.01, 0.15, 0.13, 0.96, 0.57, 0.58, 0.96, 0.93, 0.5, 0.45, 0.89, 0.44, 0.59, 0.68, 0.71, 0.85, 0.16, 0.18, 0.68, 0.37, 0.22, 0.81, 0.53, 0.26, 0.94, 0.52, 0.66, 0.55, 0.51, 0.14 };
            double[] mySmoothedArray = inputData.BoxCarSmooth(3);
            double[] expectedOutput = new[] { 0.3, 0.52, 0.57, 0.72, 0.46, 0.2, 0.08, 0.22, 0.4, 0.63, 0.62, 0.64, 0.58, 0.68, 
                0.74, 0.59, 0.59, 0.36, 0.63, 0.53, 0.66, 0.67, 0.51, 0.46, 0.22, 0.35, 0.32, 0.37, 0.25, 
                0.21, 0.4, 0.51, 0.53, 0.37, 0.3, 0.54, 0.47, 0.73, 0.53, 0.61, 0.63, 0.54, 0.39, 0.31, 0.59, 
                0.72, 0.49, 0.25, 0.21, 0.3, 0.35, 0.42, 0.42, 0.54, 0.64, 0.63, 0.63, 0.4, 0.6, 0.46, 0.61, 0.34, 0.46, 
                0.31, 0.56, 0.6, 0.68, 0.66, 0.42, 0.31, 0.1, 0.41, 0.55, 0.7, 0.7, 0.82, 0.8, 0.63, 0.61, 0.59, 0.64, 0.57, 0.66, 0.75, 0.57, 0.4, 0.34, 0.41, 0.42, 0.47, 0.52, 0.53, 0.58, 0.57, 0.71, 0.58, 0.57, 0.4 };
            Assert.That(expectedOutput, Is.EqualTo(mySmoothedArray).Within(0.1));
        }

        [Test]
        public static void TestScrambledEquals()
        {
            List<int> list1 = new() { 1, 2, 3, 4, 5, 6 };
            List<int> list2 = new() { 1, 3, 5, 6, 4, 2 };
            List<int> list3 = new() { 1, 2, 3, 4, 6, 7 };
            bool expectedTrue = list1.ScrambledEquals(list2);
            bool expectedFalse = list1.ScrambledEquals(list3); 
            
            Assert.True(expectedTrue);
            Assert.False(expectedFalse);
        }

        [Test]
        public static void TestAllSame()
        {
            MzSpectrum spectrum1 = new MzSpectrum(new[] { 2.0, 2.0 }, new[] { 1.0, 1.0 }, false);
            MzSpectrum spectrum2 = new MzSpectrum(new[] { 5.0, 5.0 }, new[] { 2.0, 2.0 }, false);

            var sameInt = new[] { 2, 2 };
            var sameDouble = new[] { 2.0, 2.0 };
            var sameSpectrum = new[] { spectrum1, spectrum1 };
            Assert.That(sameInt.AllSame());
            Assert.That(sameDouble.AllSame());
            Assert.That(sameSpectrum.AllSame());

            var differentInt = new[] { 2, 5 };
            var differentDouble = new[] { 2.0, 5.0 };
            var differentSpectrum = new[] { spectrum1, spectrum2 };
            Assert.That(!differentInt.AllSame());
            Assert.That(!differentDouble.AllSame());
            Assert.That(!differentSpectrum.AllSame());
        }

        [Test]
        [TestCase(1874.28, 373.8487, -5)]
        [TestCase(1874.28, 467.5627, -4)]
        [TestCase(1874.28, 623.7527, -3)]
        [TestCase(1874.28, 936.1327, -2)]
        [TestCase(1874.28, 1873.273, -1)]
        [TestCase(1874.28, 375.8633, 5)]
        [TestCase(1874.28, 469.5773, 4)]
        [TestCase(1874.28, 625.7673, 3)]
        [TestCase(1874.28, 938.1473, 2)]
        [TestCase(1874.28, 1875.287, 1)]

        public static void TestToMzAndMass(double mass, double mz, int charge)
        {
            Assert.That(mass, Is.EqualTo(mz.ToMass(charge)).Within(0.01));
            Assert.That(mz, Is.EqualTo(mass.ToMz(charge)).Within(0.01));
        }

        [Test]
        [TestCase("ATCG", "AUCG", true)]
        [TestCase("ATCG", "UAGC", false)]
        [TestCase("ATCGZ", "AUCGZ", true)]
        [TestCase("ATCGZ", "UAGCZ", false)]
        [TestCase("ATCGACGAATCACGATCAGTCATGCATTGCTAACT", "AUCGACGAAUCACGAUCAGUCAUGCAUUGCUAACU", true)]
        [TestCase("ATCGACGAATCACGATCAGTCATGCATTGCTAACT", "UAGCUGCUUAGUGCUAGUCAGUACGUAACGAUUGA", false)]
        [TestCase("ATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACT", "AUCGACGAAUCACGAUCAGUCAUGCAUUGCUAACUAUCGACGAAUCACGAUCAGUCAUGCAUUGCUAACUAUCGACGAAUCACGAUCAGUCAUGCAUUGCUAACUAUCGACGAAUCACGAUCAGUCAUGCAUUGCUAACUAUCGACGAAUCACGAUCAGUCAUGCAUUGCUAACUAUCGACGAAUCACGAUCAGUCAUGCAUUGCUAACU", true)]
        [TestCase("ATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACTATCGACGAATCACGATCAGTCATGCATTGCTAACT", "UAGCUGCUUAGUGCUAGUCAGUACGUAACGAUUGAUAGCUGCUUAGUGCUAGUCAGUACGUAACGAUUGAUAGCUGCUUAGUGCUAGUCAGUACGUAACGAUUGAUAGCUGCUUAGUGCUAGUCAGUACGUAACGAUUGAUAGCUGCUUAGUGCUAGUCAGUACGUAACGAUUGAUAGCUGCUUAGUGCUAGUCAGUACGUAACGAUUGA", false)]
        public static void TestTranscribe(string input, string expected, bool isCodingStrand)
        {
            Assert.That(input.Transcribe(isCodingStrand), Is.EqualTo(expected));
        }
    }
}