using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Readers;
using Quantification.Differential;
using StatisticalModels;

namespace Test.Quantification.Differential;

/// <summary>
/// The design layer built from a pipeline's <see cref="IExperimentalDesign"/> (<see cref="AnalysisSample.FromExperimentalDesign"/>,
/// <see cref="AnalysisSample.IdFor"/>), and <see cref="Contrast.ToWeights"/>, which hands a contrast to the model fits.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class AnalysisSampleTests
{
    private static SpectraFileInfo Run(string file, string condition, int biorep, int techrep = 0, int fraction = 0) =>
        new($@"C:\data\{file}.raw", condition, biorep, techrep, fraction);

    private static IsobaricQuantSampleInfo Channel(string file, int plex, string channel, string condition, string? name,
        bool reference = false, int fraction = 0) =>
        new($@"C:\data\{file}.raw", condition, 0, 0, fraction, plex, channel, 126.127726, reference) { SampleName = name };

    private static SampleExperimentalDesign Design(params ISampleInfo[][] files)
    {
        var design = new SampleExperimentalDesign();
        foreach (var samples in files)
            design.Add(samples[0].FullFilePathWithExtension, samples);
        return design;
    }

    [Test]
    public void ALabelFreeSamplesFractionsAndReplicatesAreOneSample()
    {
        var design = Design(
            new ISampleInfo[] { Run("a_f1", "A", 0, fraction: 1) },
            new ISampleInfo[] { Run("a_f2", "A", 0, fraction: 2) },
            new ISampleInfo[] { Run("a_t2", "A", 0, techrep: 1, fraction: 1) },
            new ISampleInfo[] { Run("a2", "A", 1) },
            new ISampleInfo[] { Run("b1", "B", 0) });

        var samples = AnalysisSample.FromExperimentalDesign(design);

        Assert.That(samples.Select(s => s.SampleId), Is.EqualTo(new[] { "A_1", "A_2", "B_1" }));
        Assert.That(samples.Select(s => s.Factors[AnalysisSample.ConditionFactor]), Is.EqualTo(new[] { "A", "A", "B" }));
        Assert.That(samples.All(s => s.Batch is null && s.Individual is null && s.Strata.Count == 0 && s.Covariates.Count == 0));
    }

    [Test]
    public void ALabelFreeSampleIsKeyedAsTheQuantifiedTablesLabelIt()
    {
        // The observation table names a run's sample with SampleGroupLabels.ForSample too, so the two layers agree.
        var runs = new[] { Run("x", "ctrl", 0), Run("y", "ctrl", 2), Run("z", "drug_2", 0) };
        foreach (var run in runs)
            Assert.That(AnalysisSample.IdFor(run), Is.EqualTo(SampleGroupLabels.ForSample(run)));
        Assert.That(runs.Select(AnalysisSample.IdFor), Is.EqualTo(new[] { "ctrl_1", "ctrl_3", "drug_2_1" }));
    }

    [Test]
    public void AnIsobaricChannelIsOneSampleInEveryFractionOfItsPlexAndThePlexIsItsBatch()
    {
        ISampleInfo[] Plex1(string file, int fraction) => new ISampleInfo[]
        {
            Channel(file, 1, "126", "A", "mouse1", fraction: fraction),
            Channel(file, 1, "127N", "B", "mouse2", fraction: fraction),
            Channel(file, 1, "128C", "pool", "pool", reference: true, fraction: fraction),
        };
        var design = Design(Plex1("p1_f1", 1), Plex1("p1_f2", 2), new ISampleInfo[]
        {
            Channel("p2", 2, "126", "A", "mouse3"),
            Channel("p2", 2, "127N", "B", null),
            Channel("p2", 2, "128C", "pool", "pool2", reference: true),
        });

        var samples = AnalysisSample.FromExperimentalDesign(design);

        Assert.That(samples.Select(s => s.SampleId), Is.EqualTo(new[] { "mouse1", "mouse2", "mouse3", "plex2_127N" }),
            "the pooled reference channels are not biological samples");
        Assert.That(samples.Select(s => s.Batch), Is.EqualTo(new[] { "plex1", "plex1", "plex2", "plex2" }));
        Assert.That(samples.Select(s => s.Factors[AnalysisSample.ConditionFactor]), Is.EqualTo(new[] { "A", "B", "A", "B" }));
    }

    /// <summary>
    /// Review round 4 (@pcruzparri): MetaMorpheus's TMT design writes every kit channel for every file, and a channel the
    /// plex does not use comes with no condition and no sample name. It was returned as a sample, so the design refused
    /// the whole experiment over a sample nobody declared.
    /// </summary>
    [Test]
    public void AChannelThePlexDoesNotUseIsNotASample()
    {
        ISampleInfo[] Plex(string file, int plex, string first, string second) => new ISampleInfo[]
        {
            Channel(file, plex, "126", "ctrl", first),
            Channel(file, plex, "127N", "drug", second),
            Channel(file, plex, "131C", "", null),
        };
        var design = Design(Plex("p1", 1, "s1", "s2"), Plex("p2", 2, "s3", "s4"));

        var samples = AnalysisSample.FromExperimentalDesign(design);

        Assert.That(samples.Select(s => s.SampleId), Is.EqualTo(new[] { "s1", "s2", "s3", "s4" }));
        var analysis = AnalysisDesign.Create(samples,
            new DesignSpec { Factors = new[] { new FactorSpec(AnalysisSample.ConditionFactor, "ctrl") } });
        Assert.That(analysis.Samples, Has.Count.EqualTo(4));
    }

    /// <summary>The same for mzLib's own SDRF reader: PXD008841's first plex read with the TMT11 kit leaves 131C empty in every file.</summary>
    [Test]
    public void AnSdrfPlexWithAnEmptyChannelGivesOnlyItsAnnotatedSamples()
    {
        Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(IsobaricMassTagType.TMT11, out var tmt11), Is.True);
        string sdrf = Path.Combine(TestContext.CurrentContext.TestDirectory, "FileReadingTests", "ExternalFileTypes",
            "SdrfDesign", "PXD008841.subset.sdrf.tsv");
        var read = SdrfIsobaricDesign.Read(sdrf,
            new SdrfIsobaricDesignOptions { Tag = tmt11!, PlexFileNamePattern = @"TMT_?pool(?<plex>\d+)" });
        Assert.That(read.Refusals, Is.Empty, read.Report());
        var plexOne = new SampleExperimentalDesign();
        foreach (var (file, channels) in read.ToExperimentalDesign().FileNameSampleInfoDictionary)
            if (channels.Cast<IsobaricQuantSampleInfo>().First().PlexId == 1)
                plexOne.Add(file, channels);
        Assert.That(plexOne.FileNameSampleInfoDictionary.Values.SelectMany(c => c).Count(c => c.Condition.Length == 0),
            Is.GreaterThan(0), "the fixture needs the empty channel");

        var samples = AnalysisSample.FromExperimentalDesign(plexOne);

        Assert.That(samples, Has.Count.EqualTo(10));
        Assert.That(samples.Select(s => s.SampleId), Has.None.EqualTo("plex1_131C"));
        Assert.That(samples.All(s => s.Factors[AnalysisSample.ConditionFactor].Length > 0));
    }

    /// <summary>A channel that is named but has no condition is a sample someone declared, so it is kept and the design refuses it by name.</summary>
    [Test]
    public void ANamedChannelWithoutAConditionIsKeptAndRefusedByName()
    {
        var design = Design(new ISampleInfo[]
        {
            Channel("p1", 1, "126", "ctrl", "s1"), Channel("p1", 1, "127N", "drug", "s2"), Channel("p1", 1, "127C", "", "s9"),
        });

        var samples = AnalysisSample.FromExperimentalDesign(design);

        Assert.That(samples.Select(s => s.SampleId), Is.EqualTo(new[] { "s1", "s2", "s9" }));
        var e = Assert.Throws<ArgumentException>(() => AnalysisDesign.Create(samples,
            new DesignSpec { Factors = new[] { new FactorSpec(AnalysisSample.ConditionFactor, "ctrl") } }));
        Assert.That(e!.Message, Does.Contain("'s9'"));
    }

    [Test]
    public void TheSamplesFeedTheDesign()
    {
        var design = Design(
            new ISampleInfo[] { Channel("p1", 1, "126", "ctrl", "s1"), Channel("p1", 1, "127N", "drug", "s2") },
            new ISampleInfo[] { Channel("p2", 2, "126", "ctrl", "s3"), Channel("p2", 2, "127N", "drug", "s4") });

        var analysis = AnalysisDesign.Create(AnalysisSample.FromExperimentalDesign(design),
            new DesignSpec { Factors = new[] { new FactorSpec(AnalysisSample.ConditionFactor, "ctrl") } });

        Assert.That(analysis.ColumnNames, Is.EqualTo(new[] { "intercept", "condition=drug", "batch=plex2" }));
        Assert.That(analysis.DefaultContrasts().Select(c => c.Label), Is.EqualTo(new[] { "condition=drug vs condition=ctrl" }));
    }

    [Test]
    public void ADesignMixingLabelFreeFilesAndChannelsIsRefused()
    {
        var design = Design(new ISampleInfo[] { Run("a", "A", 0) }, new ISampleInfo[] { Channel("p", 1, "126", "A", "s") });
        var e = Assert.Throws<ArgumentException>(() => AnalysisSample.FromExperimentalDesign(design));
        Assert.That(e!.Message, Does.Contain("mixes label-free files with isobaric channels"));
    }

    [Test]
    public void ASampleWhoseRunsDisagreeOnItsConditionIsRefusedAndNamed()
    {
        var design = Design(
            new ISampleInfo[] { Channel("p_f1", 1, "126", "A", "mouse1", fraction: 1) },
            new ISampleInfo[] { Channel("p_f2", 1, "126", "B", "mouse1", fraction: 2) });
        var e = Assert.Throws<ArgumentException>(() => AnalysisSample.FromExperimentalDesign(design));
        Assert.That(e!.Message, Does.Contain("Sample 'mouse1' is condition 'A' in 'p_f1.raw' but 'B' in 'p_f2.raw'"));
    }

    [Test]
    public void ANameOnTwoChannelsIsRefusedRatherThanMerged()
    {
        // A sample name is display text (IsobaricQuantSampleInfo.SampleName), so channels are never merged by it.
        var design = Design(
            new ISampleInfo[] { Channel("p1", 1, "126", "A", "mouse1") },
            new ISampleInfo[] { Channel("p2", 2, "126", "A", "mouse1") });
        var e = Assert.Throws<ArgumentException>(() => AnalysisSample.FromExperimentalDesign(design));
        Assert.That(e!.Message, Does.Contain("Sample name 'mouse1' is on channel 126 of plex 1 and channel 126 of plex 2"));
    }

    [Test]
    public void AChannelNamedTwoWaysIsRefused()
    {
        var design = Design(
            new ISampleInfo[] { Channel("p_f1", 1, "126", "A", "mouse1", fraction: 1) },
            new ISampleInfo[] { Channel("p_f2", 1, "126", "A", null, fraction: 2) });
        var e = Assert.Throws<ArgumentException>(() => AnalysisSample.FromExperimentalDesign(design));
        Assert.That(e!.Message, Does.Contain("Channel 126 of plex 1 is named 'mouse1' and 'plex1_126'"));
    }

    [Test]
    public void AContrastBecomesTheWeightsTheFitsTakeNamedByItsId()
    {
        var analysis = AnalysisDesign.Create(new[]
        {
            new AnalysisSample("s1") { Factors = new Dictionary<string, string> { ["condition"] = "ctrl" } },
            new AnalysisSample("s2") { Factors = new Dictionary<string, string> { ["condition"] = "drug" } },
            new AnalysisSample("s3") { Factors = new Dictionary<string, string> { ["condition"] = "ctrl" } },
            new AnalysisSample("s4") { Factors = new Dictionary<string, string> { ["condition"] = "drug" } },
        }, new DesignSpec { Factors = new[] { new FactorSpec("condition", "ctrl") } });
        var contrast = analysis.DefaultContrasts().Single();

        ContrastWeights weights = contrast.ToWeights();
        Assert.That(weights.Name, Is.EqualTo("c1"));
        Assert.That(weights.Weights, Is.EqualTo(contrast.Weights));

        // And the moderated result maps back to the contrast by that name.
        var responses = new double[3, 4];
        for (int f = 0; f < 3; f++)
            for (int s = 0; s < 4; s++)
                responses[f, s] = 20 + f + (s % 2 == 1 ? 1.0 : 0) + 0.1 * ((f + s) % 3);
        var fit = LinearModel.Fit(responses, analysis.Matrix, analysis.ColumnNames);
        var moderated = EmpiricalBayes.ModerateContrasts(fit, new[] { weights }, trend: false);
        Assert.That(moderated.Single().Name, Is.EqualTo(contrast.Id));
    }
}
