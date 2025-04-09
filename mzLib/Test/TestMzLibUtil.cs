using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using MzLibUtil;
using Readers;
using System.Collections.Generic;
using System.Collections.Concurrent;
using System.Threading.Tasks;

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
        public void Increment_IncrementsExistingKey()
        {
            // Arrange
            var dictionary = new Dictionary<string, int>
            {
                { "key1", 1 }
            };

            // Act
            dictionary.Increment("key1");

            // Assert
            Assert.That(dictionary["key1"], Is.EqualTo(2));
        }

        [Test]
        public void Increment_AddsNewKeyWithInitialValue()
        {
            // Arrange
            var dictionary = new Dictionary<string, int>();

            // Act
            dictionary.Increment("key1");

            // Assert
            Assert.That(dictionary.ContainsKey("key1"));
            Assert.That(dictionary["key1"], Is.EqualTo(1));
        }

        [Test]
        public void Increment_IncrementsMultipleTimes()
        {
            // Arrange
            var dictionary = new Dictionary<string, int>
            {
                { "key1", 1 }
            };

            // Act
            dictionary.Increment("key1");
            dictionary.Increment("key1");
            dictionary.Increment("key1");

            // Assert
            Assert.That(dictionary["key1"], Is.EqualTo(4));
        }

        [Test]
        public void Increment_AddsAndIncrementsNewKey()
        {
            // Arrange
            var dictionary = new Dictionary<string, int>();

            // Act
            dictionary.Increment("key1");
            dictionary.Increment("key1");

            // Assert
            Assert.That(dictionary.ContainsKey("key1"));
            Assert.That(dictionary["key1"], Is.EqualTo(2));
        }

        [Test]
        public void Increment_ThreadSafeWithConcurrentDictionary()
        {
            // Arrange
            var dictionary = new ConcurrentDictionary<string, int>();
            var tasks = new List<Task>();

            // Act
            for (int i = 0; i < 1000; i++)
            {
                tasks.Add(Task.Run(() => dictionary.Increment("key1", 1)));
            }
            Task.WaitAll(tasks.ToArray());

            // Assert
            Assert.That(dictionary["key1"], Is.EqualTo(1000));
        }

        [Test]
        public void Increment_IncrementsBySpecifiedValue()
        {
            // Arrange
            var dictionary = new Dictionary<string, int>
            {
                { "key1", 1 }
            };

            // Act
            dictionary.Increment("key1", 5);

            // Assert
            Assert.That(dictionary["key1"], Is.EqualTo(6));
        }

        [Test]
        public void Increment_AddsNewKeyWithSpecifiedValue()
        {
            // Arrange
            var dictionary = new Dictionary<string, int>();

            // Act
            dictionary.Increment("key1", 5);

            // Assert
            Assert.That(dictionary.ContainsKey("key1"));
            Assert.That(dictionary["key1"], Is.EqualTo(5));
        }

        [Test]
        public void Increment_IncrementsBySpecifiedValueMultipleTimes()
        {
            // Arrange
            var dictionary = new Dictionary<string, int>
            {
                { "key1", 1 }
            };

            // Act
            dictionary.Increment("key1", 2);
            dictionary.Increment("key1", 3);
            dictionary.Increment("key1", 4);

            // Assert
            Assert.That(dictionary["key1"], Is.EqualTo(10));
        }

        [Test]
        public void Increment_AddsAndIncrementsNewKeyBySpecifiedValue()
        {
            // Arrange
            var dictionary = new Dictionary<string, int>();

            // Act
            dictionary.Increment("key1", 2);
            dictionary.Increment("key1", 3);

            // Assert
            Assert.That(dictionary.ContainsKey("key1"));
            Assert.That(dictionary["key1"], Is.EqualTo(5));
        }

        [Test]
        public void Increment_ThreadSafeWithConcurrentDictionaryBySpecifiedValue()
        {
            // Arrange
            var dictionary = new ConcurrentDictionary<string, int>();
            var tasks = new List<Task>();

            // Act
            for (int i = 0; i < 1000; i++)
            {
                int value = i % 10 + 1; // Increment by values from 1 to 10
                tasks.Add(Task.Run(() => dictionary.Increment("key1", value)));
            }
            Task.WaitAll(tasks.ToArray());

            // Assert
            Assert.That(dictionary["key1"], Is.EqualTo(5500));
        }
    }
}
