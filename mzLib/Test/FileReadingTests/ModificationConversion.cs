using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    public static class ModificationConversions
    {
        [Test]
        public static void FirstGo()
        {
            var allMods = ModificationExtensions.AllKnownMods;
            Assert.That(allMods.Count, Is.GreaterThanOrEqualTo(3110)); // number at time of creation
        }
    }
}
