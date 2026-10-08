using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Quantification.Differential;
using Readers;
using DetectionType = FlashLFQ.DetectionType;
using FlashLfqEngine = FlashLFQ.FlashLfqEngine;
using FlashLfqIdentification = FlashLFQ.Identification;
using FlashLfqObservations = FlashLFQ.FlashLfqObservations;
using FlashLfqProteinGroup = FlashLFQ.ProteinGroup;
using FlashLfqResults = FlashLFQ.FlashLfqResults;

namespace Test.Quantification.Differential;

/// <summary>
/// GR-18: the in-memory <see cref="FlashLfqResults"/> and the peptide table FlashLFQ writes from them give the same
/// <see cref="ObservationTable"/>, under both bases and whatever the design makes of the files.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class FlashLfqObservationsTests
{
    private static string FlashLfqData => Path.Combine(TestContext.CurrentContext.TestDirectory, "FlashLFQ", "TestData");

    private static string ExternalFile(string name) => Path.Combine(TestContext.CurrentContext.TestDirectory,
        "FileReadingTests", "ExternalFileTypes", name);

    private string _outputDirectory = "";

    [SetUp]
    public void SetUp()
    {
        _outputDirectory = Path.Combine(Path.GetTempPath(), "FlashLfqObservationsTests_" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_outputDirectory);
    }

    [TearDown]
    public void TearDown()
    {
        if (Directory.Exists(_outputDirectory))
            Directory.Delete(_outputDirectory, true);
    }

    /// <summary>Writes the peptide table and reads it back, as a pipeline that stored it would.</summary>
    private QuantifiedPeptideFile WriteAndRead(FlashLfqResults results)
    {
        string path = Path.Combine(_outputDirectory, "AllQuantifiedPeptides.tsv");
        results.WriteResults(null, path, null, null, true);
        return new QuantifiedPeptideFile(path);
    }

    private static void AssertSameTable(ObservationTable expected, ObservationTable actual)
    {
        Assert.That(actual.Basis, Is.EqualTo(expected.Basis));
        Assert.That(actual.Runs, Is.EqualTo(expected.Runs));
        Assert.That(actual.Samples, Is.EqualTo(expected.Samples));
        Assert.That(actual.SamplesWithoutValues, Is.EqualTo(expected.SamplesWithoutValues));
        Assert.That(actual.Warnings, Is.EqualTo(expected.Warnings));
        Assert.That(actual.ProteinGroups.Values, Is.EquivalentTo(expected.ProteinGroups.Values));
        Assert.That(actual.Peptides.Select(p => (p.FullSequence, p.BaseSequence, string.Join(";", p.ProteinGroups))),
            Is.EqualTo(expected.Peptides.Select(p => (p.FullSequence, p.BaseSequence, string.Join(";", p.ProteinGroups)))));

        // Bit for bit: the stored table is the same numbers, not numbers within a tolerance.
        Assert.That(actual.Observations().Select(o => (BitConverter.DoubleToInt64Bits(o.Log2Intensity), o.State)),
            Is.EqualTo(expected.Observations().Select(o => (BitConverter.DoubleToInt64Bits(o.Log2Intensity), o.State))));
    }

    /// <summary>
    /// The K562 PSMs, each mapped to one group per accession it matches, so a PSM matching two proteins makes its peptide
    /// shared. One FlashLFQ group object per name, since FlashLFQ compares groups by reference. (Multi-accession group
    /// names are covered by <see cref="EveryDetectionTypeMapsToOneStateFromEitherSource"/>.)
    /// </summary>
    private static List<FlashLfqIdentification> K562Identifications(List<SpectraFileInfo> files)
    {
        var groups = new Dictionary<string, FlashLfqProteinGroup>();
        var identifications = new List<FlashLfqIdentification>();
        foreach (var psm in SpectrumMatchTsvReader.ReadPsmTsv(Path.Combine(FlashLfqData, "AllPSMs.psmtsv"), out _))
        {
            SpectraFileInfo? file = files.FirstOrDefault(f => f.FilenameWithoutExtension == psm.FileNameWithoutExtension);
            if (file == null) continue;

            var psmGroups = psm.ProteinAccession.Split('|').Distinct().Select(accession =>
            {
                if (!groups.TryGetValue(accession, out var group))
                    groups[accession] = group = new FlashLfqProteinGroup(accession, "", "");
                return group;
            }).ToList();

            identifications.Add(new FlashLfqIdentification(file, psm.BaseSeq, psm.FullSequence, (double)psm.MonoisotopicMass,
                (double)psm.RetentionTime, psm.PrecursorCharge, psmGroups, decoy: psm.DecoyContamTarget.Contains('D')));
        }
        return identifications;
    }

    private static IEnumerable<ProteinGroupInfo> GroupsOf(FlashLfqResults results) =>
        results.ProteinGroups.Keys.Select(name => new ProteinGroupInfo(name, 0.001, false, false));

    [Test]
    [NonParallelizable] // full FlashLFQ runs
    [TestCase(false, TestName = "BothSourcesGiveTheSameTable_WithoutMbr")]
    [TestCase(true, TestName = "BothSourcesGiveTheSameTable_WithMbr")]
    public void BothSourcesGiveTheSameTableOnOneSearch(bool matchBetweenRuns)
    {
        var files = new List<SpectraFileInfo>
        {
            new(Path.Combine(FlashLfqData, "20100614_Velos1_TaGe_SA_K562_3.mzML"), "a", 0, 0, 0),
            new(Path.Combine(FlashLfqData, "20100614_Velos1_TaGe_SA_K562_4.mzML"), "b", 0, 0, 0),
        };
        FlashLfqResults results = new FlashLfqEngine(K562Identifications(files), matchBetweenRuns: matchBetweenRuns,
            maxThreads: 1, silent: true).Run();
        QuantifiedPeptideFile stored = WriteAndRead(results);
        var groups = GroupsOf(results).ToList();
        Assert.That(results.PeptideModifiedSequences.Values.Any(p => p.ProteinGroups.Count > 1), "some peptides are shared");

        string first = files[0].FilenameWithoutExtension, second = files[1].FilenameWithoutExtension;
        var designs = new Dictionary<string, ObservationRun[]>
        {
            ["two samples"] = ObservationRun.FromSpectraFiles(files).ToArray(),
            ["two fractions"] = new[] { new ObservationRun(first, "S", Fraction: 1), new ObservationRun(second, "S", Fraction: 2) },
            ["two replicates"] = new[]
            {
                new ObservationRun(first, "S", TechnicalReplicate: 1), new ObservationRun(second, "S", TechnicalReplicate: 2),
            },
        };

        foreach (var (design, runs) in designs)
        foreach (QuantBasis basis in Enum.GetValues<QuantBasis>())
        {
            var inMemory = FlashLfqObservations.FromResults(results, runs, groups, basis);
            var fromFile = FlashLfqObservations.FromPeptideTable(stored, runs, groups, basis);

            Assert.That(inMemory.Peptides, Has.Count.EqualTo(results.PeptideModifiedSequences.Count), design);
            Assert.That(inMemory.Observations().Count(o => !double.IsNaN(o.Log2Intensity)), Is.GreaterThan(0), design);
            AssertSameTable(inMemory, fromFile);

            if (matchBetweenRuns && basis == QuantBasis.MbrKept && design == "two samples")
                Assert.That(inMemory.Observations().Count(o => o.State == ObservationState.MbrTransferred && !double.IsNaN(o.Log2Intensity)),
                    Is.GreaterThan(0), "the run transfers something, so parity covers transfers");
        }
    }

    /// <summary>
    /// Hand-set results covering every detection type MetaMorpheus writes, a calibrated file name, a file with no value
    /// and an intensity that needs all 17 digits to round-trip.
    /// </summary>
    private static (FlashLfqResults Results, ObservationRun[] Runs) EveryDetectionType()
    {
        var a = new SpectraFileInfo(@"C:\data\A-calib.mzML", "c", 0, 0, 0);
        var b = new SpectraFileInfo(@"C:\data\B.mzML", "c", 1, 0, 0);
        var c = new SpectraFileInfo(@"C:\data\C.mzML", "d", 0, 0, 0);
        var g1 = new FlashLfqProteinGroup("P1|P2", "", "");
        var g2 = new FlashLfqProteinGroup("Q1", "", "");
        var ids = new List<FlashLfqIdentification>
        {
            new(a, "PEPMSMS", "PEPMSMS", 800, 10, 2, new List<FlashLfqProteinGroup> { g1 }),
            new(a, "PEPZERO", "PEPZERO", 800, 10, 2, new List<FlashLfqProteinGroup> { g1 }),
            new(b, "PEPZERO", "PEPZERO", 800, 10, 2, new List<FlashLfqProteinGroup> { g2 }),
            new(a, "PEPAMB", "PEP[Common Variable:Oxidation on M]AMB", 800, 10, 2, new List<FlashLfqProteinGroup> { g2 }),
        };
        var results = new FlashLfqResults(new List<SpectraFileInfo> { a, b, c }, ids);
        results.CalculatePeptideResults(quantifyAmbiguousPeptides: false); // every file NotDetected, intensity 0

        void Set(string sequence, SpectraFileInfo file, double intensity, DetectionType type)
        {
            results.PeptideModifiedSequences[sequence].SetIntensity(file, intensity);
            results.PeptideModifiedSequences[sequence].SetDetectionType(file, type);
        }

        Set("PEPMSMS", a, 12345.678901234567, DetectionType.MSMS);
        Set("PEPMSMS", b, 2000.25, DetectionType.MBR);
        Set("PEPZERO", a, 0, DetectionType.MSMS); // zeroed because the sample's strongest fraction was ambiguous
        Set("PEPZERO", b, 0, DetectionType.MSMSIdentifiedButNotQuantified);
        Set("PEP[Common Variable:Oxidation on M]AMB", a, 3000, DetectionType.MSMSAmbiguousPeakfinding);
        Set("PEP[Common Variable:Oxidation on M]AMB", b, 0, DetectionType.MBR);

        return (results, ObservationRun.FromSpectraFiles(new[]
        {
            new SpectraFileInfo(@"C:\data\A.mzML", "c", 0, 0, 0), b, c,
        }).ToArray());
    }

    [Test]
    public void EveryDetectionTypeMapsToOneStateFromEitherSource()
    {
        var (results, runs) = EveryDetectionType();
        QuantifiedPeptideFile stored = WriteAndRead(results);
        var groups = new[] { new ProteinGroupInfo("P1|P2", 0.001, false, false), new ProteinGroupInfo("Q1", 0.02, false, false) };

        foreach (QuantBasis basis in Enum.GetValues<QuantBasis>())
            AssertSameTable(FlashLfqObservations.FromResults(results, runs, groups, basis),
                FlashLfqObservations.FromPeptideTable(stored, runs, groups, basis));

        var table = FlashLfqObservations.FromPeptideTable(stored, runs, groups, QuantBasis.MbrKept);
        Assert.That(table.Samples, Is.EqualTo(new[] { "c_1", "c_2", "d_1" }));
        Assert.That(table.SamplesWithoutValues, Is.EqualTo(new[] { "d_1" }), "the file with no value is kept");

        var cells = table.Observations().ToDictionary(o => (o.Peptide.BaseSequence, o.SampleId));
        Assert.That(cells[("PEPMSMS", "c_1")], Has.Property("State").EqualTo(ObservationState.Quantified)
            .And.Property("Log2Intensity").EqualTo(Math.Log2(12345.678901234567)));
        Assert.That(cells[("PEPMSMS", "c_2")], Has.Property("State").EqualTo(ObservationState.MbrTransferred)
            .And.Property("Log2Intensity").EqualTo(Math.Log2(2000.25)));
        Assert.That(cells[("PEPMSMS", "d_1")].State, Is.EqualTo(ObservationState.NotDetected));
        Assert.That(cells[("PEPZERO", "c_1")].State, Is.EqualTo(ObservationState.AmbiguousPeak), "MSMS written with 0");
        Assert.That(cells[("PEPZERO", "c_2")].State, Is.EqualTo(ObservationState.IdentifiedNotQuantified));
        Assert.That(cells[("PEPAMB", "c_1")].State, Is.EqualTo(ObservationState.AmbiguousPeak), "ambiguous, though written with 3000");
        Assert.That(cells[("PEPAMB", "c_2")].State, Is.EqualTo(ObservationState.AmbiguousPeak), "MBR written with 0");
        Assert.That(cells.Values.Where(o => o.State != ObservationState.Quantified && o.State != ObservationState.MbrTransferred)
            .Select(o => o.Log2Intensity), Is.All.NaN);

        Assert.That(table.Peptides.Single(p => p.BaseSequence == "PEPZERO").UniqueProteinGroup, Is.Null, "two groups");
        Assert.That(table.Peptides.Single(p => p.BaseSequence == "PEPMSMS").UniqueProteinGroup, Is.EqualTo("P1|P2"));
    }

    [Test]
    public void IsoTrackerResultsAreRefused()
    {
        var (results, runs) = EveryDetectionType();
        results.IsoTracker = true;
        Assert.Throws<ArgumentException>(() =>
            FlashLfqObservations.FromResults(results, runs, Array.Empty<ProteinGroupInfo>(), QuantBasis.MsmsOnly));
    }

    private QuantifiedPeptideFile Table(params string[] lines)
    {
        string path = Path.Combine(_outputDirectory, "QuantifiedPeptides.tsv");
        File.WriteAllLines(path, lines);
        return new QuantifiedPeptideFile(path);
    }

    [TestCase("IsoTrack_MSMS")]
    [TestCase("Foo")]
    [TestCase("3")]
    [TestCase("")]
    public void AnUnreadableDetectionTypeIsRefused(string detectionType)
    {
        var table = Table("Sequence\tBase Sequence\tProtein Groups\tGene Names\tOrganism\tIntensity_f1\tDetection Type_f1",
            $"PEPA\tPEPA\tP1\t\t\t100\t{detectionType}");

        var ex = Assert.Throws<ArgumentException>(() => FlashLfqObservations.FromPeptideTable(table,
            new[] { new ObservationRun("f1", "S1") }, Array.Empty<ProteinGroupInfo>(), QuantBasis.MsmsOnly));
        Assert.That(ex!.Message, Does.Contain("PEPA").And.Contain("f1"));
    }

    [Test]
    public void ABlankIntensityIsRefused()
    {
        var table = Table("Sequence\tBase Sequence\tProtein Groups\tGene Names\tOrganism\tIntensity_f1\tDetection Type_f1",
            "PEPA\tPEPA\tP1\t\t\t\tMSMS");

        var ex = Assert.Throws<ArgumentException>(() => FlashLfqObservations.FromPeptideTable(table,
            new[] { new ObservationRun("f1", "S1") }, Array.Empty<ProteinGroupInfo>(), QuantBasis.MsmsOnly));
        Assert.That(ex!.Message, Does.Contain("PEPA").And.Contain("f1"));
    }

    [Test]
    public void AnIsoTrackerPeptideTableIsRefused()
    {
        var table = Table(
            "Sequence\tBase Sequence\tPeak Order\tProtein Groups\tGene Names\tOrganism\tIntensity_f1\tRetentionTime (min)_f1\tDetection Type_f1",
            "PEPA\tPEPA\t1\tP1\t\t\t100\t12.5\tMSMS");

        Assert.Throws<ArgumentException>(() => FlashLfqObservations.FromPeptideTable(table,
            new[] { new ObservationRun("f1", "S1") }, Array.Empty<ProteinGroupInfo>(), QuantBasis.MsmsOnly));
    }

    [Test]
    public void AnEmptyPeptideTableGivesEverySampleWithoutValues()
    {
        var table = Table("Sequence\tBase Sequence\tProtein Groups\tGene Names\tOrganism\tIntensity_f1\tDetection Type_f1");

        var result = FlashLfqObservations.FromPeptideTable(table, new[] { new ObservationRun("f1", "S1") },
            Array.Empty<ProteinGroupInfo>(), QuantBasis.MsmsOnly);

        Assert.That(result.Peptides, Is.Empty);
        Assert.That(result.SamplesWithoutValues, Is.EqualTo(new[] { "S1" }));
    }

    /// <summary>The 18-file MetaMorpheus 1.1.11 peptide table in the reader tests: calibrated column names, one transfer.</summary>
    [Test]
    public void ReadsAMetaMorpheusPeptideTable()
    {
        var stored = new QuantifiedPeptideFile(ExternalFile("MetaMorpheus_1.1.11_AllQuantifiedPeptides.tsv"));
        var runs = stored.First().Samples.Keys
            .Select(label => label.EndsWith("-calib", StringComparison.Ordinal) ? label[..^"-calib".Length] : label)
            .Select(name => new ObservationRun(name, name)).ToArray();

        var mbrKept = FlashLfqObservations.FromPeptideTable(stored, runs, Array.Empty<ProteinGroupInfo>(), QuantBasis.MbrKept);
        var msmsOnly = FlashLfqObservations.FromPeptideTable(stored, runs, Array.Empty<ProteinGroupInfo>(), QuantBasis.MsmsOnly);

        Assert.That(mbrKept.Samples, Has.Count.EqualTo(18));
        Assert.That(mbrKept.Peptides, Has.Count.EqualTo(3));
        int peptide = mbrKept.Peptides.Select(p => p.BaseSequence).ToList().IndexOf("CACASHVAK");
        int sample = mbrKept.Samples.ToList().IndexOf("QE-002118_GM6_a");
        Assert.That(mbrKept.Log2Intensity(peptide, sample), Is.EqualTo(Math.Log2(143944.34375)));
        Assert.That(mbrKept.State(peptide, sample), Is.EqualTo(ObservationState.MbrTransferred));
        Assert.That(msmsOnly.Log2Intensity(peptide, sample), Is.NaN);
        Assert.That(msmsOnly.State(peptide, sample), Is.EqualTo(ObservationState.MbrTransferred));
        Assert.That(mbrKept.Warnings, Has.Count.EqualTo(1), "no protein table given, so its groups are undescribed");
    }

    /// <summary>The six-row MetaMorpheus 1.1.11 protein-group table in the reader tests: a two-protein group, a contaminant and a decoy.</summary>
    [Test]
    public void ReadsAMetaMorpheusProteinGroupTable()
    {
        var file = new ProteinGroupFromTsvFile(ExternalFile("MetaMorpheus_1.1.11_AllQuantifiedProteinGroups.tsv"));

        var groups = FlashLfqObservations.ReadProteinGroups(file);

        Assert.That(groups.Select(g => g.Name), Is.EqualTo(file.Select(r => r.ProteinGroupName)));
        Assert.That(groups.Select(g => g.QValue), Is.EqualTo(file.Select(r => r.QValue)));
        Assert.That(groups.Single(g => g.Name == "P0C0S5|Q71UI9").Accessions, Is.EqualTo(new[] { "P0C0S5", "Q71UI9" }));
        Assert.That(groups.Where(g => g.IsDecoy).Select(g => g.Name), Is.EqualTo(new[] { "DECOY_P62750" }));
        Assert.That(groups.Where(g => g.IsContaminant).Select(g => g.Name), Is.EqualTo(new[] { "P02769" }));
        Assert.That(groups.Select(g => g.Genes), Is.EqualTo(file.Select(r => r.Gene)));
    }
}
