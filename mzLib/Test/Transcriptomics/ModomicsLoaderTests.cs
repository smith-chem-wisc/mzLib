using System;
using System.IO;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using UsefulProteomicsDatabases.Transcriptomics;
using Omics.Modifications;

namespace Test.Transcriptomics
{
    [TestFixture]
    public class ModomicsLoaderTests
    {
        [Test]
        public void LoadModomics_ParsesEntriesCorrectly()
        {
            // Adjust resource name to match your project/namespace
            var resourceName = "UsefulProteomicsDatabases.Transcriptomics.modomics.json";

            var mods = ModomicsLoader.LoadModomics();

            Assert.That(mods, Is.Not.Empty);
            Assert.That(mods.All(m => m is Modification), Is.True);
            Assert.That(mods.All(m => !string.IsNullOrWhiteSpace(m.IdWithMotif)), Is.True);

            // Check a known modification
            var m5Cm = mods.FirstOrDefault(m => m.IdWithMotif == "m5Cm");
            Assert.That(m5Cm, Is.Not.Null);
            Assert.That(m5Cm.ModificationType, Is.EqualTo("RNA"));
            Assert.That(m5Cm.FeatureType, Is.EqualTo("MODOMICS"));
            Assert.That(m5Cm.ChemicalFormula, Is.Not.Null);
            Assert.That(m5Cm.MonoisotopicMass, Is.GreaterThan(0));
        }

        [Test]
        public void LoadModomics_SkipsUnknownBases()
        {
            var resourceName = "UsefulProteomicsDatabases.Transcriptomics.modomics.json";

            var mods = ModomicsLoader.LoadModomics();

            // Should skip modifications with unknown bases (e.g. "X")
            Assert.That(mods.Any(m => m.IdWithMotif == "Xm" || m.IdWithMotif == "xX"), Is.False);
        }
    }
}
