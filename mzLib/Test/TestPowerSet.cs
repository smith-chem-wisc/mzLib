using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Omics.Fragmentation;

namespace Test
{
    public class TestPowerSet
    {
        [Test]
        public void TestPowerSetGeneration_EmptyList()
        {
            List<double> numberList = new List<double>();
            var uniqueSumOfSubsets = PowerSet.UniqueSubsetSums(numberList, 3);
            List<double> expectedSums = new List<double> { 0.0 };
            Assert.AreEqual(expectedSums.Count, uniqueSumOfSubsets.Count);
            Assert.IsTrue(expectedSums.SequenceEqual(uniqueSumOfSubsets));
        }

        [Test]
        public void TestPowerSetGenration_StackOverFlow()
        {
            List<double> numberList = new List<double> {
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0 };
            var ex1 = Assert.Throws<ArgumentException>(() => PowerSet.UniqueSubsetSums(numberList, 3));
            Assert.That(ex1.Message, Is.EqualTo("Input list is too large. Maximum supported size is 15 to avoid excessive memory usage."));

            List<double> numberList_2 = new List<double> {
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
            var result = PowerSet.UniqueSubsetSums(numberList_2, -3);
            List<double> expectedList = new List<double> { 0 };
            Assert.That(result.Count == 1);
        }

        [Test]
        public void TestPowerSetGeneration_ThreeNumber()
        {
            List<double> numberList = new List<double> { 1.0, 2.0, 3.0 };
            var uniqueSumOfSubsets = PowerSet.UniqueSubsetSums(numberList, 4);
            List<double> expectedSums = new List<double> { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
            Assert.AreEqual(expectedSums.Count, uniqueSumOfSubsets.Count);
            Assert.IsTrue(expectedSums.SequenceEqual(uniqueSumOfSubsets));
        }

        [Test]
        public void TestPowerSetGeneration_DuplicateNumbers()
        {
            List<double> numberList = new List<double> { 1.0, 2.0, 2.0 };
            var uniqueSumOfSubsets_Max3 = PowerSet.UniqueSubsetSums(numberList, 3);
            List<double> expectedSums_Max3 = new List<double> { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
            Assert.AreEqual(expectedSums_Max3.Count, uniqueSumOfSubsets_Max3.Count);
            Assert.IsTrue(expectedSums_Max3.SequenceEqual(uniqueSumOfSubsets_Max3));

            var uniqueSumOfSubsets_Max2 = PowerSet.UniqueSubsetSums(numberList, 2);
            List<double> expectedSums_Max2 = new List<double> { 0.0, 1.0, 2.0, 3.0, 4.0 };
            Assert.AreEqual(expectedSums_Max2.Count, uniqueSumOfSubsets_Max2.Count);
            Assert.IsTrue(expectedSums_Max2.SequenceEqual(uniqueSumOfSubsets_Max2));

            var uniqueSumOfSubsets_Max1 = PowerSet.UniqueSubsetSums(numberList, 1);
            List<double> expectedSums_Max1 = new List<double> { 0.0, 1.0, 2.0 };
            Assert.AreEqual(expectedSums_Max1.Count, uniqueSumOfSubsets_Max1.Count);
            Assert.IsTrue(expectedSums_Max1.SequenceEqual(uniqueSumOfSubsets_Max1));

        }

        [Test]
        public void TestPowerSetGeneration_complexSet()
        {
            List<double> numberList = new List<double> { 1.0, 2.0, 3.0, 4.0 };
            var uniqueSumOfSubsets = PowerSet.UniqueSubsetSums(numberList, 4);
            List<double> expectedSums = new List<double> { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };
            Assert.AreEqual(expectedSums.Count, uniqueSumOfSubsets.Count);
            Assert.IsTrue(expectedSums.SequenceEqual(uniqueSumOfSubsets));
        }
    }
}
