using NUnit.Framework;
using Omics.Modifications;
using Omics.Modifications.Conversion;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using UsefulProteomicsDatabases;

namespace Test.UsefulProteomicsDatabases;

[TestFixture]
[ExcludeFromCodeCoverage]
public class ProteinDbWriterIntegrationTests
{
    [Test]
    public void ProteinDbWriter_DefaultConversionOptions_PreservesModifications()
    {
        var mod = Mods.GetModification("Carbamidomethyl on C");
        Assert.That(mod, Is.Not.Null, "Test requires Carbamidomethyl modification.");

        var proteinMods = new Dictionary<int, List<Modification>>
        {
            { 2, new List<Modification> { mod! } }
        };

        var protein = new Protein("ACDE", "P3", oneBasedModifications: proteinMods);
        var tempFile = Path.Combine(TestContext.CurrentContext.WorkDirectory, Guid.NewGuid() + ".xml");

        try
        {
            ProteinDbWriter.WriteXmlDatabase(
                new Dictionary<string, HashSet<Tuple<int, Modification>>>(),
                new List<Protein> { protein },
                tempFile,
                updateTimeStamp: false);

            var contents = File.ReadAllText(tempFile);
            Assert.That(contents, Does.Contain(mod!.IdWithMotif));
        }
        finally
        {
            if (File.Exists(tempFile))
            {
                File.Delete(tempFile);
            }
        }
    }

    [Test]
    public void ProteinDbWriter_KeepOriginalAnnotation_PreservesUnsupportedModifications()
    {
        ModificationMotif.TryGetMotif("M", out var motif);
        var customMod = new Modification(
            _originalId: "CustomWriter",
            _modificationType: "CustomWriter",
            _target: motif,
            _monoisotopicMass: null,
            _chemicalFormula: null);

        var extraMods = new Dictionary<string, HashSet<Tuple<int, Modification>>>
        {
            { "P4", new HashSet<Tuple<int, Modification>> { new Tuple<int, Modification>(2, customMod) } }
        };

        var protein = new Protein("MA", "P4");
        var options = new ProteinDbWriterConversionOptions
        {
            Enabled = true,
            TargetConvention = ModificationNamingConvention.Unimod,
            HandlingMode = SequenceConversionHandlingMode.KeepOriginalAnnotation
        };

        var tempFile = Path.Combine(TestContext.CurrentContext.WorkDirectory, Guid.NewGuid() + ".xml");

        try
        {
            ProteinDbWriter.WriteXmlDatabase(
                extraMods,
                new List<Protein> { protein },
                tempFile,
                updateTimeStamp: false,
                conversionOptions: options);

            var contents = File.ReadAllText(tempFile);
            Assert.That(contents, Does.Contain(customMod.IdWithMotif));
        }
        finally
        {
            if (File.Exists(tempFile))
            {
                File.Delete(tempFile);
            }
        }
    }
}
