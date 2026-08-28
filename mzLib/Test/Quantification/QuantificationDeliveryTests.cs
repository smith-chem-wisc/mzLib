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
    /// <param name="sourceDirectory">
    /// When given, the spectral matches carry absolute file paths inside this directory, which is what
    /// lets the engine derive a default output directory. The experimental design stays keyed by bare
    /// file name either way, as <see cref="QuantificationEngine.PivotByFile"/> requires.
    /// </param>
    private static void BuildFixture(
        out IExperimentalDesign design,
        out List<ISpectralMatch> spectralMatches,
        out List<IBioPolymerWithSetMods> peptides,
        out List<IBioPolymerGroup> proteinGroups,
        string sourceDirectory = null)
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

        string path1 = sourceDirectory == null ? File1 : Path.Combine(sourceDirectory, File1);
        string path2 = sourceDirectory == null ? File2 : Path.Combine(sourceDirectory, File2);

        spectralMatches = new List<ISpectralMatch>
        {
            new MockSpectralMatch(path1, pep1.FullSequence, pep1.BaseSequence, 10.0, 1, new[] { pep1 })
                { Intensities = new double[] { 100, 200, 300 } },
            new MockSpectralMatch(path1, pep2.FullSequence, pep2.BaseSequence, 10.0, 2, new[] { pep2 })
                { Intensities = new double[] { 400, 500, 600 } },
            new MockSpectralMatch(path2, pep1.FullSequence, pep1.BaseSequence, 10.0, 3, new[] { pep1 })
                { Intensities = new double[] { 1000, 2000, 3000 } },
            new MockSpectralMatch(path2, pep2.FullSequence, pep2.BaseSequence, 10.0, 4, new[] { pep2 })
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
            UseSharedPeptidesForProteinQuant = false,
            // Most tests here are about the values, not the files. Writing is on by default and
            // requires an output directory, so it is turned off unless a test asks for it.
            WriteRawInformation = false,
            WritePeptideInformation = false,
            WriteProteinInformation = false
        };
    }

    #endregion

    /// <summary>
    /// Writing stays on by default: the raw file is what makes re-processing a search under different
    /// strategies possible, so a caller has to opt out of it rather than opt in.
    /// </summary>
    [Test]
    public void WriteFlags_DefaultToOn()
    {
        var fresh = new QuantificationParameters();
        Assert.Multiple(() =>
        {
            Assert.That(fresh.WriteRawInformation, Is.True);
            Assert.That(fresh.WritePeptideInformation, Is.True);
            Assert.That(fresh.WriteProteinInformation, Is.True);
        });
    }

    /// <summary>
    /// Writing is on by default. With no OutputDirectory AND no usable source-file directory — this
    /// fixture's matches carry bare file names — the engine has nowhere defensible to write, so it must
    /// say so rather than scattering files into the working directory.
    /// </summary>
    [Test]
    public void Run_WithWritingOnAndNothingToDeriveFrom_ReturnsClearFailure()
    {
        BuildFixture(out var design, out var spectralMatches, out var peptides, out var proteinGroups);

        var parameters = SimpleParameters();
        parameters.WriteRawInformation = true;   // OutputDirectory is still empty

        QuantificationResults results = null;
        Assert.DoesNotThrow(() =>
            results = new QuantificationEngine(parameters, design, spectralMatches, peptides, proteinGroups).Run());

        Assert.Multiple(() =>
        {
            Assert.That(results.Success, Is.False);
            Assert.That(results.Summary, Does.Contain("OutputDirectory is not set"));
            Assert.That(results.WrittenFiles, Is.Empty);
            Assert.That(results.ProteinIntensities, Is.Empty);
        });

        foreach (var group in proteinGroups)
            Assert.That(group.IntensitiesBySample, Is.Null, "a rejected run must not touch the entities");
    }

    /// <summary>
    /// The default: writing is on, no OutputDirectory was given, and the matches name real files — so
    /// output lands beside the data rather than failing or landing in the working directory.
    /// </summary>
    [Test]
    public void Run_WithWritingOnAndNoOutputDirectory_WritesBesideTheSourceFiles()
    {
        string dataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory,
            "QuantDefaultOutput_" + TestContext.CurrentContext.Test.ID);
        Directory.CreateDirectory(dataDirectory);
        try
        {
            BuildFixture(out var design, out var spectralMatches, out var peptides, out var proteinGroups,
                dataDirectory);

            var parameters = SimpleParameters();   // OutputDirectory left empty
            parameters.WriteRawInformation = true;
            parameters.WritePeptideInformation = true;
            parameters.WriteProteinInformation = true;

            var results = new QuantificationEngine(parameters, design, spectralMatches, peptides, proteinGroups).Run();

            Assert.Multiple(() =>
            {
                Assert.That(results.Success, Is.True);
                Assert.That(results.OutputDirectory, Is.EqualTo(dataDirectory));
                Assert.That(results.WrittenFiles, Has.Count.EqualTo(3));
            });

            foreach (string file in results.WrittenFiles)
            {
                Assert.That(File.Exists(file), Is.True, file + " should exist");
                Assert.That(Path.GetDirectoryName(file), Is.EqualTo(dataDirectory));
            }
        }
        finally
        {
            if (Directory.Exists(dataDirectory)) Directory.Delete(dataDirectory, recursive: true);
        }
    }

    /// <summary>
    /// The derived directory is only a default. An explicit OutputDirectory wins, even when the source
    /// files would have supplied one.
    /// </summary>
    [Test]
    public void Run_WithOutputDirectorySet_IgnoresTheSourceFileDirectory()
    {
        string root = Path.Combine(TestContext.CurrentContext.TestDirectory,
            "QuantExplicitOutput_" + TestContext.CurrentContext.Test.ID);
        string dataDirectory = Path.Combine(root, "data");
        string chosenDirectory = Path.Combine(root, "chosen");
        Directory.CreateDirectory(dataDirectory);
        try
        {
            BuildFixture(out var design, out var spectralMatches, out var peptides, out var proteinGroups,
                dataDirectory);

            var parameters = SimpleParameters(chosenDirectory);
            parameters.WriteProteinInformation = true;

            var results = new QuantificationEngine(parameters, design, spectralMatches, peptides, proteinGroups).Run();

            Assert.That(results.OutputDirectory, Is.EqualTo(chosenDirectory));
            Assert.That(File.Exists(Path.Combine(chosenDirectory, QuantificationWriter.ProteinGroupFileName)), Is.True);
            Assert.That(Directory.GetFiles(dataDirectory), Is.Empty, "nothing should be written beside the data");
        }
        finally
        {
            if (Directory.Exists(root)) Directory.Delete(root, recursive: true);
        }
    }

    /// <summary>
    /// Bare file names and relative paths resolve against the process's working directory, which is the
    /// ambiguity the default exists to avoid, so neither may produce one.
    /// </summary>
    [Test]
    public void TryGetSourceFileDirectory_WithNoUsablePath_DerivesNothing()
    {
        Assert.Multiple(() =>
        {
            Assert.That(QuantificationEngine.TryGetSourceFileDirectory(null, out _), Is.False);
            Assert.That(QuantificationEngine.TryGetSourceFileDirectory(new List<ISpectralMatch>(), out _), Is.False);
            Assert.That(QuantificationEngine.TryGetSourceFileDirectory(MatchesFrom("file1.raw"), out _), Is.False,
                "a bare file name has no directory");
            Assert.That(QuantificationEngine.TryGetSourceFileDirectory(
                    MatchesFrom(Path.Combine("data", "file1.raw")), out _), Is.False,
                "a relative directory depends on the working directory");
        });
    }

    /// <summary>
    /// One directory is used as-is; several sibling directories collapse to their nearest common
    /// ancestor, so a fractionated search writes to the folder holding the fractions.
    /// </summary>
    [Test]
    public void TryGetSourceFileDirectory_DerivesTheDirectoryOrItsNearestCommonAncestor()
    {
        string root = Path.Combine(TestContext.CurrentContext.TestDirectory, "QuantDerive");
        string fraction1 = Path.Combine(root, "fraction1");
        string fraction2 = Path.Combine(root, "fraction2");

        Assert.Multiple(() =>
        {
            Assert.That(QuantificationEngine.TryGetSourceFileDirectory(
                MatchesFrom(Path.Combine(fraction1, "a.raw"), Path.Combine(fraction1, "b.raw")),
                out string single), Is.True);
            Assert.That(single, Is.EqualTo(fraction1));

            Assert.That(QuantificationEngine.TryGetSourceFileDirectory(
                MatchesFrom(Path.Combine(fraction1, "a.raw"), Path.Combine(fraction2, "b.raw")),
                out string ancestor), Is.True);
            Assert.That(ancestor, Is.EqualTo(root));
        });
    }

    /// <summary>
    /// Files whose only shared ancestor is the drive or filesystem root have no shared home. Writing to
    /// the root would be worse than asking, so nothing is derived.
    /// </summary>
    [Test]
    public void TryGetSourceFileDirectory_RejectsARootAsTheCommonAncestor()
    {
        string root = Path.GetPathRoot(TestContext.CurrentContext.TestDirectory);
        Assert.That(root, Is.Not.Null.And.Not.Empty);

        bool derived = QuantificationEngine.TryGetSourceFileDirectory(
            MatchesFrom(Path.Combine(root, "alpha", "a.raw"), Path.Combine(root, "beta", "b.raw")),
            out _);

        Assert.That(derived, Is.False);
    }

    /// <summary>Spectral matches that carry nothing but the given file paths.</summary>
    private static List<ISpectralMatch> MatchesFrom(params string[] filePaths) =>
        filePaths
            .Select((path, i) => (ISpectralMatch)new MockSpectralMatch(path, "SEQ", "SEQ", 10.0, i + 1))
            .ToList();

    /// <summary>
    /// Quantifying without writing anything is legitimate, and must not require an output directory.
    /// </summary>
    [Test]
    public void Run_WithWritingOff_NeedsNoOutputDirectory()
    {
        BuildFixture(out var design, out var spectralMatches, out var peptides, out var proteinGroups);

        QuantificationResults results = null;
        Assert.DoesNotThrow(() =>
            results = new QuantificationEngine(SimpleParameters(), design, spectralMatches, peptides, proteinGroups).Run());

        Assert.That(results.Success, Is.True);
        Assert.That(results.WrittenFiles, Is.Empty);
        Assert.That(results.OutputDirectory, Is.Null, "a run that writes nothing resolves no directory");
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

    /// <summary>
    /// The rendered table must be rectangular. Every row of a protein table is written by
    /// <see cref="BioPolymerGroup.ToString"/> while the header is written once, from the first group
    /// in the list -- so if two groups disagree about which columns exist, every value after the
    /// disagreement lands under the wrong header and the file is silently wrong rather than obviously
    /// broken.
    ///
    /// They disagreed because <c>HasIntensityData</c> gated the two intensity columns per sample
    /// group, and it is false for a group whose samples were all unobserved. P002 here is seen only in
    /// file1, so its three file2 groups carry no intensity: it emitted two columns per group where
    /// P001 emitted four.
    ///
    /// The gate is now whether the entity was quantified at all, not whether a particular group has a
    /// value, so an unobserved sample yields an empty cell rather than a missing column. That keeps
    /// absent distinguishable from zero -- see
    /// <see cref="Run_OmitsUnobservedSamplesRatherThanStoringZero"/>, which pins the same distinction
    /// on the data model.
    /// </summary>
    [Test]
    public void RenderedTable_IsRectangular_WhenAGroupHasUnobservedSamples()
    {
        BuildFixture(out var design, out var spectralMatches, out var peptides, out var proteinGroups);

        // P002 is observed in file1 only, so its three file2 samples have no intensity.
        string p002Peptide = peptides.Single(p => ((Protein)p.Parent).Accession == "P002").BaseSequence;
        spectralMatches = spectralMatches
            .Where(sm => !(sm.FullFilePath == File2 && sm.BaseSequence == p002Peptide))
            .ToList();

        new QuantificationEngine(SimpleParameters(), design, spectralMatches, peptides, proteinGroups).Run();

        // Render the table the way a writer does: one header, taken from the first group, then a row
        // per group from its own ToString().
        string header = proteinGroups.First().GetTabSeparatedHeader();
        var rows = proteinGroups.Select(g => g.ToString()).ToList();

        int headerFields = header.Split('\t').Length;
        Assert.Multiple(() =>
        {
            for (int i = 0; i < rows.Count; i++)
            {
                Assert.That(rows[i].Split('\t'), Has.Length.EqualTo(headerFields),
                    $"row {i} ({proteinGroups[i].BioPolymerGroupName}) does not line up with the header; "
                    + "every value after the first missing column is written under the wrong name");
            }
        });
    }

    /// <summary>
    /// The same invariant for a group the engine quantified but found nothing for -- a decoy, or a
    /// group whose peptides all fell below threshold. It is a matrix row like any other
    /// (<see cref="QuantificationEngine.GetUniquePeptideToProteinMap"/> gives every group an entry,
    /// empty or not), so it gets a full row of zeros, no intensities, and previously no intensity
    /// columns at all -- while its neighbours emitted them. That is the whole-file version of the
    /// same defect, because the header may be taken from either kind of group.
    /// </summary>
    [Test]
    public void RenderedTable_IsRectangular_WhenAGroupWasQuantifiedButHasNoValues()
    {
        BuildFixture(out var design, out var spectralMatches, out var peptides, out var proteinGroups);

        // A third group with no peptide in the matrix: quantified, but nothing to show for it.
        var orphan = new Protein("ORPHANK", "P003");
        var orphanPeptide = orphan
            .Digest(new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5), new List<Modification>(), new List<Modification>())
            .First();
        proteinGroups.Add(new BioPolymerGroup(
            new HashSet<IBioPolymer> { orphan },
            new HashSet<IBioPolymerWithSetMods> { orphanPeptide },
            new HashSet<IBioPolymerWithSetMods> { orphanPeptide }));

        new QuantificationEngine(SimpleParameters(), design, spectralMatches, peptides, proteinGroups).Run();

        var p003 = proteinGroups.Single(g => g.BioPolymerGroupName.Contains("P003"));
        Assert.That(p003.SamplesForQuantification, Is.Not.Null,
            "every protein group is a matrix row, so the engine should have written the sample list");
        Assert.That(p003.IntensitiesBySample, Is.Empty,
            "nothing was observed for it, so no sample carries a value");

        string header = proteinGroups.First().GetTabSeparatedHeader();
        int headerFields = header.Split('\t').Length;

        Assert.Multiple(() =>
        {
            foreach (var group in proteinGroups)
                Assert.That(group.ToString().Split('\t'), Has.Length.EqualTo(headerFields),
                    $"{group.BioPolymerGroupName} does not line up with the header");
        });

        // And the header must be the same whichever group happens to be first -- the writer takes it
        // from proteinGroups.First(), which is not necessarily a quantified one.
        Assert.That(p003.GetTabSeparatedHeader(), Is.EqualTo(header),
            "an unquantified group must describe the same columns as a quantified one");
    }

    /// <summary>
    /// An entity the engine never touched keeps the spectral-count-only shape. This is what
    /// distinguishes "quantification did not run" from "it ran and found nothing". The gate is both
    /// <see cref="IHasSampleIntensities.SamplesForQuantification"/> being non-empty AND
    /// <see cref="IHasSampleIntensities.IntensitiesBySample"/> being non-null, rather than a
    /// per-group flag: a run either quantified its groups or it did not, and every group of one
    /// quantified entity agrees. Both fields are needed because ConstructSubsetBioPolymerGroup can
    /// leave an empty sample list beside a non-null dictionary.
    /// </summary>
    [Test]
    public void RenderedTable_OmitsIntensityColumnsEntirely_WhenNothingWasQuantified()
    {
        BuildFixture(out _, out _, out _, out var proteinGroups);

        var unquantified = proteinGroups.First();
        Assert.That(unquantified.SamplesForQuantification, Is.Null, "fixture starts unquantified");

        string header = unquantified.GetTabSeparatedHeader();
        Assert.Multiple(() =>
        {
            Assert.That(header, Does.Not.Contain("Intensity_"),
                "no quantification ran, so there is nothing for an intensity column to hold");
            Assert.That(header, Does.Not.Contain("IntensityOccupancy_"));
            Assert.That(unquantified.ToString().Split('\t'), Has.Length.EqualTo(header.Split('\t').Length));
        });
    }

    /// <summary>
    /// The rectangular table must not have been bought by inventing zeros. An unobserved sample gets
    /// an empty cell; a sample observed at zero would get "0". Rectangularity is a field-count
    /// property and this is a field-content one, so the tests above cannot see the difference -- a
    /// change that wrote 0 into the absent cells would keep them all passing while claiming a
    /// measurement that was never made.
    /// </summary>
    [Test]
    public void RenderedRow_LeavesAnUnobservedSampleEmpty_RatherThanWritingZero()
    {
        BuildFixture(out var design, out var spectralMatches, out var peptides, out var proteinGroups);

        string p002Peptide = peptides.Single(p => ((Protein)p.Parent).Accession == "P002").BaseSequence;
        spectralMatches = spectralMatches
            .Where(sm => !(sm.FullFilePath == File2 && sm.BaseSequence == p002Peptide))
            .ToList();

        new QuantificationEngine(SimpleParameters(), design, spectralMatches, peptides, proteinGroups).Run();

        var p002 = proteinGroups.Single(g => g.BioPolymerGroupName.Contains("P002"));
        var header = p002.GetTabSeparatedHeader().Split('	');
        var row = p002.ToString().Split('	');

        // The three file2 channels are unobserved for P002; the three file1 channels are not.
        var absent = header
            .Select((name, i) => (name, i))
            .Where(x => x.name.StartsWith("Intensity_") && x.name.Contains(File2.Replace(".raw", "")))
            .ToList();
        var present = header
            .Select((name, i) => (name, i))
            .Where(x => x.name.StartsWith("Intensity_") && x.name.Contains(File1.Replace(".raw", "")))
            .ToList();

        Assert.That(absent, Is.Not.Empty, "the unobserved channels should still have columns");
        Assert.That(present, Is.Not.Empty, "the observed channels should too");

        Assert.Multiple(() =>
        {
            foreach (var (name, i) in absent)
                Assert.That(row[i], Is.Empty,
                    $"{name} was never observed, so its cell must be empty rather than a fabricated zero");
            foreach (var (name, i) in present)
                Assert.That(row[i], Is.Not.Empty, $"{name} was observed and must carry its value");
        });
    }

    /// <summary>
    /// A subset group built for a file that none of the parent's samples match gets an empty sample
    /// list beside an empty-but-non-null intensity dictionary
    /// (<see cref="BioPolymerGroup.ConstructSubsetBioPolymerGroup"/>). Gating the columns on the
    /// dictionary alone would then advertise intensity columns that
    /// <see cref="BioPolymerGroup.PopulateSampleGroupResults"/> can never fill, because with no
    /// samples it falls back to grouping by spectral-match file path.
    ///
    /// So the gate requires a non-empty sample list as well. This pins that: the subset describes no
    /// intensity columns rather than permanently empty ones.
    /// </summary>
    [Test]
    public void SubsetGroupForAnUnmatchedFile_DoesNotAdvertiseIntensityColumnsItCannotFill()
    {
        BuildFixture(out var design, out var spectralMatches, out var peptides, out var proteinGroups);
        new QuantificationEngine(SimpleParameters(), design, spectralMatches, peptides, proteinGroups).Run();

        var parent = proteinGroups.First();
        Assert.That(parent.SamplesForQuantification, Is.Not.Null.And.Not.Empty, "parent was quantified");

        // The fixture's samples carry bare file names, so an absolute path matches none of them --
        // the path/bare-name mismatch this repo's own fixtures produce.
        var subset = parent.ConstructSubsetBioPolymerGroup(@"C:\somewhereile1.raw");

        Assert.That(subset.SamplesForQuantification, Is.Empty, "no sample matched that path");
        Assert.That(subset.IntensitiesBySample, Is.Not.Null.And.Empty,
            "the subset still gets a dictionary, which is why the dictionary alone cannot be the gate");

        string header = subset.GetTabSeparatedHeader();
        Assert.Multiple(() =>
        {
            Assert.That(header, Does.Not.Contain("Intensity_"),
                "nothing can fill these, so they must not be advertised");
            Assert.That(subset.ToString().Split('	'), Has.Length.EqualTo(header.Split('	').Length));
        });
    }
}
