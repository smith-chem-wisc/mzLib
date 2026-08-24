using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Omics.SpectralMatch;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Quantification;
using Quantification.Strategies;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Test.Omics;

namespace Test.Quantification;

/// <summary>
/// Covers the delivery half of the quantification pipeline: writing results back onto the entities
/// that opt in via <see cref="IHasSampleIntensities"/>, populating <see cref="QuantificationResults"/>,
/// and the output files produced by <see cref="QuantificationWriter"/>.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class QuantificationDeliveryTests
{
    #region Fixture

    private class TestExperimentalDesign : IExperimentalDesign
    {
        public Dictionary<string, ISampleInfo[]> FileNameSampleInfoDictionary { get; init; }
    }

    private const string File1 = "file1.raw";
    private const string File2 = "file2.raw";
    private static readonly string[] Channels = { "126", "127N", "127C" };

    /// <summary>
    /// Two files x three TMT channels, one peptide per protein, with distinct intensities per cell so
    /// that a transposed or mis-ordered result is detectable.
    /// </summary>
    private static void BuildFixture(
        out IExperimentalDesign design,
        out List<ISpectralMatch> spectralMatches,
        out List<IBioPolymerWithSetMods> peptides,
        out List<IBioPolymerGroup> proteinGroups)
    {
        var dict = new Dictionary<string, ISampleInfo[]>();
        foreach (string file in new[] { File1, File2 })
        {
            var samples = new ISampleInfo[Channels.Length];
            for (int c = 0; c < Channels.Length; c++)
            {
                samples[c] = new IsobaricQuantSampleInfo(
                    file, "Cond" + c, 1, 1, 0, 0, Channels[c], 126.0 + c * 0.1, c == 0);
            }
            dict[file] = samples;
        }
        design = new TestExperimentalDesign { FileNameSampleInfoDictionary = dict };

        var protein1 = new Protein("PEPTIDEK", "P001");
        var protein2 = new Protein("SAMPLERK", "P002");
        var digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5);
        var noMods = new List<Modification>();

        var pep1 = protein1.Digest(digestionParams, noMods, noMods).First();
        var pep2 = protein2.Digest(digestionParams, noMods, noMods).First();
        peptides = new List<IBioPolymerWithSetMods> { pep1, pep2 };

        proteinGroups = new List<IBioPolymerGroup>
        {
            new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein1 },
                new HashSet<IBioPolymerWithSetMods> { pep1 },
                new HashSet<IBioPolymerWithSetMods> { pep1 }),
            new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein2 },
                new HashSet<IBioPolymerWithSetMods> { pep2 },
                new HashSet<IBioPolymerWithSetMods> { pep2 })
        };

        spectralMatches = new List<ISpectralMatch>
        {
            new MockSpectralMatch(File1, pep1.FullSequence, pep1.BaseSequence, 10.0, 1, new[] { pep1 })
                { Intensities = new double[] { 100, 200, 300 } },
            new MockSpectralMatch(File1, pep2.FullSequence, pep2.BaseSequence, 10.0, 2, new[] { pep2 })
                { Intensities = new double[] { 400, 500, 600 } },
            new MockSpectralMatch(File2, pep1.FullSequence, pep1.BaseSequence, 10.0, 3, new[] { pep1 })
                { Intensities = new double[] { 1000, 2000, 3000 } },
            new MockSpectralMatch(File2, pep2.FullSequence, pep2.BaseSequence, 10.0, 4, new[] { pep2 })
                { Intensities = new double[] { 4000, 5000, 6000 } }
        };
    }

    private static QuantificationParameters SimpleParameters(string outputDirectory = "")
    {
        return new QuantificationParameters
        {
            SpectralMatchNormalizationStrategy = new NoNormalization(),
            SpectralMatchToPeptideRollUpStrategy = new SumRollUp(),
            PeptideNormalizationStrategy = new NoNormalization(),
            CollapseStrategy = new NoCollapse(),
            PeptideToProteinRollUpStrategy = new SumRollUp(),
            ProteinNormalizationStrategy = new NoNormalization(),
            OutputDirectory = outputDirectory,
            UseSharedPeptidesForProteinQuant = false
        };
    }

    #endregion

    /// <summary>
    /// The three write flags must default to off. They previously defaulted to on while
    /// <see cref="QuantificationWriter"/> threw, which made a default-configured Run() throw.
    /// </summary>
    [Test]
    public void WriteFlags_DefaultToOff()
    {
        var fresh = new QuantificationParameters();
        Assert.Multiple(() =>
        {
            Assert.That(fresh.WriteRawInformation, Is.False);
            Assert.That(fresh.WritePeptideInformation, Is.False);
            Assert.That(fresh.WriteProteinInformation, Is.False);
        });
    }

    /// <summary>
    /// A default-configured engine must run to completion rather than throwing out of the writer.
    /// </summary>
    [Test]
    public void Run_WithDefaultParameters_DoesNotThrow()
    {
        BuildFixture(out var design, out var spectralMatches, out var peptides, out var proteinGroups);

        var parameters = SimpleParameters();
        var engine = new QuantificationEngine(parameters, design, spectralMatches, peptides, proteinGroups);

        QuantificationResults results = null;
        Assert.DoesNotThrow(() => results = engine.Run());
        Assert.That(results.Success, Is.True);
        Assert.That(results.WrittenFiles, Is.Empty);
    }

    /// <summary>
    /// The engine must write its final values onto each protein group's
    /// <see cref="IHasSampleIntensities.IntensitiesBySample"/> and set the column order.
    /// </summary>
    [Test]
    public void Run_PopulatesIntensitiesBySampleOnProteinGroups()
    {
        BuildFixture(out var design, out var spectralMatches, out var peptides, out var proteinGroups);

        foreach (var group in proteinGroups)
        {
            Assert.That(group.IntensitiesBySample, Is.Null, "fixture should start unquantified");
            Assert.That(group.SamplesForQuantification, Is.Null);
        }

        new QuantificationEngine(SimpleParameters(), design, spectralMatches, peptides, proteinGroups)
            .Run(out var proteinMatrix);

        foreach (var group in proteinGroups)
        {
            Assert.That(group.SamplesForQuantification, Is.Not.Null);
            Assert.That(group.SamplesForQuantification, Has.Count.EqualTo(6), "2 files x 3 channels");
            Assert.That(group.IntensitiesBySample, Is.Not.Null);
        }

        var p001 = proteinGroups.Single(g => g.BioPolymerGroupName.Contains("P001"));
        var row = proteinMatrix.GetRow(p001);
        var samples = p001.SamplesForQuantification;

        for (int i = 0; i < samples.Count; i++)
        {
            Assert.That(p001.IntensitiesBySample[samples[i]], Is.EqualTo(row[i]).Within(1e-9),
                $"column {i} ({samples[i]}) must match the matrix");
        }

        Assert.That(p001.IntensitiesBySample[samples[0]], Is.EqualTo(100.0).Within(1e-9));
        Assert.That(p001.IntensitiesBySample[samples[5]], Is.EqualTo(3000.0).Within(1e-9));
    }

    /// <summary>
    /// The matrix uses 0 to mean "not observed", so zero-valued cells must be absent from the
    /// dictionary rather than present with a zero — otherwise every sample looks quantified.
    /// </summary>
    [Test]
    public void Run_OmitsUnobservedSamplesRatherThanStoringZero()
    {
        BuildFixture(out var design, out var spectralMatches, out var peptides, out var proteinGroups);

        // P002's peptide is only seen in file1, so its file2 columns should be absent, not zero.
        string p002Peptide = peptides.Single(p => ((Protein)p.Parent).Accession == "P002").BaseSequence;
        spectralMatches = spectralMatches
            .Where(sm => !(sm.FullFilePath == File2 && sm.BaseSequence == p002Peptide))
            .ToList();
        Assert.That(spectralMatches, Has.Count.EqualTo(3), "one match should have been removed");

        new QuantificationEngine(SimpleParameters(), design, spectralMatches, peptides, proteinGroups).Run();

        var p002 = proteinGroups.Single(g => g.BioPolymerGroupName.Contains("P002"));
        var file2Samples = p002.SamplesForQuantification
            .Where(s => s.FullFilePathWithExtension == File2)
            .ToList();

        Assert.That(file2Samples, Has.Count.EqualTo(3));
        foreach (var sample in file2Samples)
            Assert.That(p002.IntensitiesBySample.ContainsKey(sample), Is.False,
                "an unobserved sample must be absent, not stored as zero");

        Assert.That(p002.IntensitiesBySample, Has.Count.EqualTo(3), "only the file1 channels are observed");
    }

    /// <summary>
    /// The point of the write-back: BioPolymerGroup's existing output machinery
    /// (PopulateSampleGroupResults / GetTabSeparatedHeader / ToString) can now be fed by the engine.
    /// </summary>
    [Test]
    public void Run_MakesProteinGroupOutputMachineryRenderIntensities()
    {
        BuildFixture(out var design, out var spectralMatches, out var peptides, out var proteinGroups);
        new QuantificationEngine(SimpleParameters(), design, spectralMatches, peptides, proteinGroups).Run();

        var p001 = (BioPolymerGroup)proteinGroups.Single(g => g.BioPolymerGroupName.Contains("P001"));
        p001.PopulateSampleGroupResults();

        Assert.That(p001.SampleGroupResults, Is.Not.Null.And.Not.Empty);
        Assert.That(p001.SampleGroupResults.Any(r => r.HasIntensityData), Is.True,
            "at least one sample group must carry intensity data after a run");
        Assert.That(p001.SampleGroupResults.Sum(r => r.Intensity), Is.GreaterThan(0));
    }

    /// <summary>
    /// Peptide values are delivered on the results object rather than on the peptide objects, because
    /// PeptideWithSetModifications is [Serializable] and serialized in bulk during indexing.
    /// </summary>
    [Test]
    public void Run_ReturnsPeptideAndProteinTables()
    {
        BuildFixture(out var design, out var spectralMatches, out var peptides, out var proteinGroups);
        var results = new QuantificationEngine(SimpleParameters(), design, spectralMatches, peptides, proteinGroups).Run();

        Assert.Multiple(() =>
        {
            Assert.That(results.Success, Is.True);
            Assert.That(results.Summary, Is.EqualTo("Quantification completed successfully."));
            Assert.That(results.Samples, Has.Count.EqualTo(6));
            Assert.That(results.PeptideIntensities, Has.Count.EqualTo(2));
            Assert.That(results.ProteinIntensities, Has.Count.EqualTo(2));
        });

        var pep1 = peptides.Single(p => p.BaseSequence == "PEPTIDEK");
        var pep1Values = results.PeptideIntensities[pep1];
        Assert.That(pep1Values.Values.Sum(), Is.EqualTo(100 + 200 + 300 + 1000 + 2000 + 3000).Within(1e-9));

        var p001 = proteinGroups.Single(g => g.BioPolymerGroupName.Contains("P001"));
        Assert.That(results.ProteinIntensities[p001], Is.EquivalentTo(p001.IntensitiesBySample),
            "the returned table and the write-back must agree");
    }

    /// <summary>
    /// Validation failures must surface as an unsuccessful result with empty tables, not as a throw.
    /// </summary>
    [Test]
    public void Run_OnInvalidInput_ReturnsFailureWithEmptyTables()
    {
        BuildFixture(out _, out var spectralMatches, out var peptides, out var proteinGroups);
        var results = new QuantificationEngine(SimpleParameters(), null, spectralMatches, peptides, proteinGroups).Run();

        Assert.Multiple(() =>
        {
            Assert.That(results.Success, Is.False);
            Assert.That(results.Summary, Does.Contain("Experimental design is null"));
            Assert.That(results.PeptideIntensities, Is.Empty);
            Assert.That(results.ProteinIntensities, Is.Empty);
            Assert.That(results.Samples, Is.Empty);
        });
    }

    /// <summary>
    /// With all three write flags on, the engine writes three files, reports their paths, and the
    /// protein file's contents match the matrix it was built from.
    /// </summary>
    [Test]
    public void Run_WithWritingEnabled_WritesAllThreeFiles()
    {
        BuildFixture(out var design, out var spectralMatches, out var peptides, out var proteinGroups);

        string outputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory,
            "QuantWriterOutput_" + TestContext.CurrentContext.Test.ID);
        try
        {
            var parameters = SimpleParameters(outputDirectory);
            parameters.WriteRawInformation = true;
            parameters.WritePeptideInformation = true;
            parameters.WriteProteinInformation = true;

            var results = new QuantificationEngine(parameters, design, spectralMatches, peptides, proteinGroups).Run();

            Assert.That(results.WrittenFiles, Has.Count.EqualTo(3));
            foreach (string file in results.WrittenFiles)
                Assert.That(File.Exists(file), Is.True, file + " should exist");

            string rawPath = Path.Combine(outputDirectory, QuantificationWriter.RawFileName);
            var rawLines = File.ReadAllLines(rawPath);
            Assert.That(rawLines[0].Split('\t'), Does.Contain("Reporter_3"));
            Assert.That(rawLines, Has.Length.EqualTo(1 + spectralMatches.Count));

            string proteinPath = Path.Combine(outputDirectory, QuantificationWriter.ProteinGroupFileName);
            var proteinLines = File.ReadAllLines(proteinPath);
            var header = proteinLines[0].Split('\t');
            Assert.That(header[0], Is.EqualTo("Protein Group"));
            Assert.That(header, Has.Length.EqualTo(7), "row label + 6 sample columns");
            Assert.That(header.Skip(1).Distinct().Count(), Is.EqualTo(6), "sample columns must be distinct");
            Assert.That(proteinLines, Has.Length.EqualTo(1 + proteinGroups.Count));

            var p001Line = proteinLines.Single(l => l.StartsWith("P001", StringComparison.Ordinal)).Split('\t');
            Assert.That(double.Parse(p001Line[1]), Is.EqualTo(100.0).Within(1e-9));
            Assert.That(double.Parse(p001Line[6]), Is.EqualTo(3000.0).Within(1e-9));

            string peptidePath = Path.Combine(outputDirectory, QuantificationWriter.PeptideFileName);
            Assert.That(File.ReadAllLines(peptidePath), Has.Length.EqualTo(1 + peptides.Count));
        }
        finally
        {
            if (Directory.Exists(outputDirectory)) Directory.Delete(outputDirectory, recursive: true);
        }
    }

    /// <summary>
    /// A label-free run has one value per match and gets a single Intensity column; an isobaric run
    /// gets one Reporter_n column per channel. Carried over from the writer in #1001.
    /// </summary>
    [TestCase(1, "Intensity")]
    [TestCase(3, "Reporter_1")]
    [TestCase(11, "Reporter_11")]
    public void RawHeader_NamesLabelFreeAndIsobaricColumnsDifferently(int channelCount, string expectedColumn)
    {
        var columns = QuantificationWriter.RawHeader(channelCount).Split('\t');

        Assert.That(columns, Does.Contain(expectedColumn));
        Assert.That(columns, Has.Length.EqualTo(7 + channelCount), "seven identity columns plus the values");
        Assert.That(columns[0], Is.EqualTo("File"));
    }

    /// <summary>A single-channel header must not use the plural Reporter naming.</summary>
    [Test]
    public void RawHeader_SingleChannel_HasNoReporterColumns()
    {
        var columns = QuantificationWriter.RawHeader(1).Split('\t');

        Assert.That(columns, Has.None.StartsWith("Reporter_"));
        Assert.That(columns[^1], Is.EqualTo("Intensity"));
    }

    /// <summary>
    /// Two channels from different files can share a label; the writer must still emit distinct headers.
    /// </summary>
    [Test]
    public void UniqueColumnLabels_DisambiguatesRepeatedLabels()
    {
        var duplicated = new ISampleInfo[]
        {
            new IsobaricQuantSampleInfo("a.raw", "C", 1, 1, 0, 0, "126", 126.0, false),
            new IsobaricQuantSampleInfo("a.raw", "C", 1, 1, 0, 0, "126", 126.0, false),
            new IsobaricQuantSampleInfo("a.raw", "C", 1, 1, 0, 0, "127N", 127.0, false)
        };

        var labels = QuantificationWriter.UniqueColumnLabels(duplicated);

        Assert.That(labels, Has.Count.EqualTo(3));
        Assert.That(labels.Distinct().Count(), Is.EqualTo(3));
        Assert.That(labels[0], Is.EqualTo("a_126"));
        Assert.That(labels[1], Is.EqualTo("a_126_2"));
    }

    /// <summary>
    /// A label-free design gives one column per file, and the writer labels it by file name.
    /// </summary>
    [Test]
    public void SampleColumnLabel_UsesFileNameForLabelFreeSamples()
    {
        var spectraFile = new SpectraFileInfo(@"C:\data\run7.mzML", "Control", 0, 0, 0);
        Assert.That(QuantificationWriter.SampleColumnLabel(spectraFile), Is.EqualTo("run7"));
    }

    /// <summary>
    /// Nothing quantified means nothing written, rather than an empty file.
    /// </summary>
    [Test]
    public void WriteRawData_WithNoQuantifiedMatches_WritesNothing()
    {
        var unquantified = new List<ISpectralMatch>
        {
            new MockSpectralMatch(File1, "PEPTIDEK", "PEPTIDEK", 10.0, 1, Array.Empty<IBioPolymerWithSetMods>())
        };
        Assert.That(QuantificationWriter.WriteRawData(unquantified, string.Empty), Is.Null);
        Assert.That(QuantificationWriter.WriteRawData(new List<ISpectralMatch>(), string.Empty), Is.Null);
    }
}
