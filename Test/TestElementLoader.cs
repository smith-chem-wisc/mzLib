using Chemistry;
using NUnit.Framework;
using System.IO;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    [Parallelizable(ParallelScope.None)] // Element loader isn't thread safe
    public class TestElementLoader
    {
        [OneTimeSetUp]
        public void Setup()
        {
            var elementLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "SetUp.dat");
            Loaders.LoadElements(elementLocation);
        }

        [Test]
        public void TestUpdateElements()
        {
            var elementLocation = Path.Combine(TestContext.CurrentContext.TestDirectory, "lal.dat");
            Loaders.UpdateElements(elementLocation);
            Loaders.UpdateElements(elementLocation);
            Assert.IsTrue(PeriodicTable.ValidateAbundances(1e-15));
            Assert.IsTrue(PeriodicTable.ValidateAverageMasses(1e-2));
        }

        [Test]
        public static void ValidatePeriodicTable()
        {
            Assert.IsTrue(PeriodicTable.ValidateAverageMasses(1e-2));
            Assert.IsFalse(PeriodicTable.ValidateAverageMasses(1e-3));
            Assert.IsTrue(PeriodicTable.ValidateAbundances(1e-15));
            Assert.IsFalse(PeriodicTable.ValidateAbundances(0));
        }
    }
}