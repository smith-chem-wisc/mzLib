using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using MzLibUtil;
using Readers;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestMzLibUtil
    {
        [Test]
        [TestCase(@"C:\Users\bubba\Documents\Projects\K562\K562_2\20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase("C:\\Users\bubba\\Documents\\Projects\\K562\\K562_2\\20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase("20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase("20100730_Velos1_TaGe_SA_K565_4", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase(@"C:\Users\bubba\Documents\Projects\K562\K562_2\20100730_Velos1_TaGe_SA_K565_4", "20100730_Velos1_TaGe_SA_K565_4")]
        //test extra period in folder name of path
        [TestCase(@"C:\Users\bubba\Documents.docs\Projects\K562\K562_2\20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        //test extra period in filename
        [TestCase(@"C:\Users\bubba\Documents\Projects\K562\K562_2\20100730_Velos1_.TaGe_SA_K565_4.raw", "20100730_Velos1_.TaGe_SA_K565_4")]
        [TestCase("/home/seth/Pictures/penguin.jpg","penguin")]
        [TestCase("/home/seth/Pictures/penguin", "penguin")]
        [TestCase("penguin.jpg", "penguin")]
        [TestCase("penguin", "penguin")]
        [TestCase("penguin.jpg.gz", "penguin")]
        [TestCase("penguin.jpg.zip", "penguin")]
        [TestCase("penguin.jpg.mzXML", "penguin.jpg")]
        public static void TestPeriodTolerantFilenameWithoutExtension(string filenameAndOrPath, string expectedResult)
        {
            string result = PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(filenameAndOrPath);
            string extensionResult = filenameAndOrPath.GetPeriodTolerantFilenameWithoutExtension();
            Assert.AreEqual(expectedResult, result);
            Assert.AreEqual(expectedResult, extensionResult);
        }

        [Test]
        public static void TestToEnum()
        {
            Assert.IsTrue(0.ToEnum<TimsTofMsMsType>(out var result));
            Assert.AreEqual(TimsTofMsMsType.MS, result);

            Assert.IsTrue(2.ToEnum<TimsTofMsMsType>(out result));
            Assert.AreEqual(TimsTofMsMsType.MSMSFragment, result);

            Assert.IsTrue(8.ToEnum<TimsTofMsMsType>(out result));
            Assert.AreEqual(TimsTofMsMsType.PASEF, result);

            Assert.IsTrue(9.ToEnum<TimsTofMsMsType>(out result));
            Assert.AreEqual(TimsTofMsMsType.DIA, result);

            Assert.IsTrue(10.ToEnum<TimsTofMsMsType>(out result));
            Assert.AreEqual(TimsTofMsMsType.PRM, result);

            Assert.IsTrue(0.ToEnum<TimsTofAcquisitionMode>(out var result2));
            Assert.AreEqual(TimsTofAcquisitionMode.MS, result2);

            Assert.IsFalse(1.ToEnum<TimsTofMsMsType>(out result));
            Assert.IsFalse(11.ToEnum<TimsTofMsMsType>(out result));
            Assert.IsFalse(7.ToEnum<TimsTofMsMsType>(out result));
            
        }

        [Test]
        public void IsNullOrEmpty_IEnumerable_Null_ReturnsTrue()
        {
            IEnumerable<int> nullEnumerable = null;
            Assert.IsTrue(nullEnumerable.IsNullOrEmpty());
        }

        [Test]
        public void IsNullOrEmpty_IEnumerable_Empty_ReturnsTrue()
        {
            IEnumerable<int> emptyEnumerable = new List<int>();
            Assert.IsTrue(emptyEnumerable.IsNullOrEmpty());
        }

        [Test]
        public void IsNullOrEmpty_IEnumerable_NotEmpty_ReturnsFalse()
        {
            IEnumerable<int> notEmptyEnumerable = new List<int> { 1 };
            Assert.IsFalse(notEmptyEnumerable.IsNullOrEmpty());
        }

        [Test]
        public void IsNullOrEmpty_IDictionary_Null_ReturnsTrue()
        {
            IDictionary<int, int> nullDictionary = null;
            Assert.IsTrue(nullDictionary.IsNullOrEmpty());
        }

        [Test]
        public void IsNullOrEmpty_IDictionary_Empty_ReturnsTrue()
        {
            IDictionary<int, int> emptyDictionary = new Dictionary<int, int>();
            Assert.IsTrue(emptyDictionary.IsNullOrEmpty());
        }

        [Test]
        public void IsNullOrEmpty_IDictionary_NotEmpty_ReturnsFalse()
        {
            IDictionary<int, int> notEmptyDictionary = new Dictionary<int, int> { { 1, 1 } };
            Assert.IsFalse(notEmptyDictionary.IsNullOrEmpty());
        }

        [Test]
        public void IsDefaultOrNull_WithNullReferenceType_ReturnsTrue()
        {
            string value = null;
            bool result = value.IsDefaultOrNull();
            Assert.IsTrue(result);
        }

        [Test]
        public void IsDefaultOrNull_WithNonNullReferenceType_ReturnsFalse()
        {
            string value = "test";
            bool result = value.IsDefaultOrNull();
            Assert.IsFalse(result);
        }

        [Test]
        public void IsDefaultOrNull_WithDefaultStruct_ReturnsTrue()
        {
            TestStruct value = default(TestStruct);
            bool result = value.IsDefaultOrNull();
            Assert.IsTrue(result);
        }

        [Test]
        public void IsDefaultOrNull_WithNonDefaultStruct_ReturnsFalse()
        {
            TestStruct value = new TestStruct { X = 1, Y = 2 };
            bool result = value.IsDefaultOrNull();
            Assert.IsFalse(result);
        }

        [Test]
        public void IsNotDefaultOrNull_WithNullReferenceType_ReturnsFalse()
        {
            string value = null;
            bool result = value.IsNotDefaultOrNull();
            Assert.IsFalse(result);
        }

        [Test]
        public void TestRemoveSpecialCharacters()
        {
            // Test default pipe removal
            string seqWithPipes = "PE|PTI|DE";
            string seqNoPipes = seqWithPipes.ToString();
            ClassExtensions.RemoveSpecialCharacters(ref seqNoPipes);
            Assert.AreEqual("PEPTIDE", seqNoPipes);


            // Test specified character replacement
            string seqWithHash = seqWithPipes.ToString();
            ClassExtensions.RemoveSpecialCharacters(ref seqWithHash, replacement: "#", specialCharacter: @"\|");
            Assert.AreEqual("PE#PTI#DE", seqWithHash);

            // Test specified character removal
            string cleanSeq = seqWithHash.ToString();
            ClassExtensions.RemoveSpecialCharacters(ref cleanSeq, specialCharacter: "#");
            Assert.AreEqual("PEPTIDE", cleanSeq);


        }

        public struct TestStruct
        {
            public int X { get; set; }
            public int Y { get; set; }
        }
    }
}
