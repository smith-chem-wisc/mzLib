using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    public static class ModificationConversion
    {
        [Test]
        public static void FirstGo()
        {
            var allMods = ModificationExtensions.AllKnownMods;
            Assert.That(allMods.Count, Is.GreaterThanOrEqualTo(3110)); // number at time of creation
        }
    }
}
