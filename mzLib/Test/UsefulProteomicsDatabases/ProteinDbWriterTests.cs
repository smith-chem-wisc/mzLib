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
public class ProteinDbWriterTests
{
    [Test]
    public void ProteinDbWriter_ConversionOptions_ConvertsModifications()
    {
        var mod = Mods.GetModification("Carbamidomethyl on C");
        Assert.That(mod, Is.Not.Null, "Test requires Carbamidomethyl modification.");

        var proteinMods = new Dictionary<int, List<Modification>>
        {
            { 2, new List<Modification> { mod! } }
        };

        var protein = new Protein("ACDE", "P1", oneBasedModifications: proteinMods);
        var options = new ProteinDbWriterConversionOptions
        {
            Enabled = true,
            TargetConvention = ModificationNamingConvention.Unimod,
            HandlingMode = SequenceConversionHandlingMode.RemoveIncompatibleMods
        };

        var expected = SequenceConverter.Default.ConvertModificationDefinition(mod!, ModificationNamingConvention.Unimod);
        var tempFile = Path.Combine(TestContext.CurrentContext.WorkDirectory, Guid.NewGuid() + ".xml");

        try
        {
            ProteinDbWriter.WriteXmlDatabase(
                new Dictionary<string, HashSet<Tuple<int, Modification>>>(),
                new List<Protein> { protein },
                tempFile,
                updateTimeStamp: false,
                conversionOptions: options);

            var contents = File.ReadAllText(tempFile);
            Assert.That(contents, Does.Contain(expected.IdWithMotif));
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
    public void ProteinDbWriter_WithoutConversion_KeepsOriginalModifications()
    {
        var mod = Mods.GetModification("Carbamidomethyl on C");
        Assert.That(mod, Is.Not.Null, "Test requires Carbamidomethyl modification.");

        var proteinMods = new Dictionary<int, List<Modification>>
        {
            { 2, new List<Modification> { mod! } }
        };

        var protein = new Protein("ACDE", "P2", oneBasedModifications: proteinMods);
        var tempFile = Path.Combine(TestContext.CurrentContext.WorkDirectory, Guid.NewGuid() + ".xml");

        try
        {
            ProteinDbWriter.WriteXmlDatabase(
                new Dictionary<string, HashSet<Tuple<int, Modification>>>(),
                new List<Protein> { protein },
                tempFile,
                updateTimeStamp: false,
                conversionOptions: null);

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
}
