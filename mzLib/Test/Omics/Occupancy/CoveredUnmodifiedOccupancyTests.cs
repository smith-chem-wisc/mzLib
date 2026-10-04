using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Omics.SpectralMatch;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

namespace Test.Omics.Occupancy;

/// <summary>
/// A site covered without being modified reports 0/N rather than nothing, once the caller says which sites
/// were seen modified anywhere in the search. Protein P00001 is ACDEFGHIK; phosphorylation on D sits at
/// AllModsOneIsNterminus position 4 (residue 3), which ACDEF (residues 1-5) covers and GHIK (6-9) does not.
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class CoveredUnmodifiedOccupancyTests
{
    private MockBioPolymer _protein;
    private Modification _phospho;
    private MockBioPolymerWithSetMods _modified;
    private MockBioPolymerWithSetMods _unmodified;
    private MockBioPolymerWithSetMods _otherRegion;

    [SetUp]
    public void SetUp()
    {
        _protein = new MockBioPolymer("ACDEFGHIK", "P00001");
        ModificationMotif.TryGetMotif("D", out var motif);
        _phospho = new Modification("Phosphorylation", null, "Biological", null, motif, "Anywhere.", null, 79.966);
        _modified = new MockBioPolymerWithSetMods("ACDEF", "ACD[Phosphorylation]EF", _protein, 1, 5,
            new Dictionary<int, Modification> { { 4, _phospho } });
        _unmodified = new MockBioPolymerWithSetMods("ACDEF", "ACDEF", _protein, 1, 5);
        _otherRegion = new MockBioPolymerWithSetMods("GHIK", "GHIK", _protein, 6, 9);
    }

    private MockSpectralMatch Modified(string file, int scan, double? intensity = null) =>
        WithIntensity(new MockSpectralMatch(file, "ACD[Phosphorylation]EF", "ACDEF", 1.0, scan, [_modified]), intensity);

    private MockSpectralMatch Unmodified(string file, int scan, double? intensity = null) =>
        WithIntensity(new MockSpectralMatch(file, "ACDEF", "ACDEF", 1.0, scan, [_unmodified]), intensity);

    private MockSpectralMatch OtherRegion(string file, int scan) =>
        new MockSpectralMatch(file, "GHIK", "GHIK", 1.0, scan, [_otherRegion]);

    private static MockSpectralMatch WithIntensity(MockSpectralMatch psm, double? intensity)
    {
        if (intensity.HasValue)
            psm.Intensities = [intensity.Value];
        return psm;
    }

    private (int, string)[] PhosphoAtFour => [(4, _phospho.IdWithMotif)];

    #region Calculator

    [Test]
    public void ACoveredSiteNotModifiedHereReportsZeroOfItsTotal()
    {
        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            _protein, [Unmodified("b.raw", 1, 1_000_000), Unmodified("b.raw", 2, 3_000_000)], PhosphoAtFour);

        var site = result[4].Single();
        Assert.That(site.ModificationIdWithMotif, Is.EqualTo(_phospho.IdWithMotif));
        Assert.That(site.ModifiedCount, Is.EqualTo(0));
        Assert.That(site.TotalCount, Is.EqualTo(2));
        Assert.That(site.CountBasedOccupancy, Is.EqualTo(0));
        Assert.That(site.ModifiedIntensity, Is.EqualTo(0));
        Assert.That(site.TotalIntensity, Is.EqualTo(4_000_000));
        Assert.That(site.IntensityBasedStoichiometry, Is.EqualTo(0));
        Assert.That(site.ToModInfoString(), Is.EqualTo("pos3[Phosphorylation on D,info:fraction=0.00(0/2)]"));
    }

    [Test]
    public void ASiteNoPsmCoversIsStillOmitted()
    {
        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            _protein, [OtherRegion("b.raw", 1)], PhosphoAtFour);

        Assert.That(result, Is.Empty);
    }

    [Test]
    public void ASiteModifiedByEveryCoveringPsmReportsAllOfThem()
    {
        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            _protein, [Modified("a.raw", 1), Modified("a.raw", 2)], PhosphoAtFour);

        var site = result[4].Single();
        Assert.That(site.ModifiedCount, Is.EqualTo(2));
        Assert.That(site.TotalCount, Is.EqualTo(2));
        Assert.That(site.CountBasedOccupancy, Is.EqualTo(1.0));
    }

    [Test]
    public void AModifiedSiteCountsTheSameWhetherOrNotItIsAlsoListed()
    {
        var psms = new List<ISpectralMatch> { Modified("a.raw", 1, 1_000_000), Unmodified("a.raw", 2, 3_000_000) };

        var listed = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(_protein, psms, PhosphoAtFour)[4].Single();
        var unlisted = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(_protein, psms)[4].Single();

        Assert.That(listed.ModifiedCount, Is.EqualTo(unlisted.ModifiedCount).And.EqualTo(1));
        Assert.That(listed.TotalCount, Is.EqualTo(unlisted.TotalCount).And.EqualTo(2));
        Assert.That(listed.ModifiedIntensity, Is.EqualTo(unlisted.ModifiedIntensity));
        Assert.That(listed.TotalIntensity, Is.EqualTo(unlisted.TotalIntensity));
    }

    [Test]
    public void WithoutSitesToReportACoveredUnmodifiedSiteIsOmittedAsBefore()
    {
        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            _protein, [Unmodified("b.raw", 1), Unmodified("b.raw", 2)]);

        Assert.That(result, Is.Empty);
    }

    [Test]
    public void EntriesAtOnePositionAreOrderedByModificationWhenSitesAreListed()
    {
        ModificationMotif.TryGetMotif("D", out var motif);
        var aaa = new Modification("Aaa", null, "Biological", null, motif, "Anywhere.", null, 1.0);
        var zzzForm = new MockBioPolymerWithSetMods("ACDEF", "ACD[Zzz]EF", _protein, 1, 5,
            new Dictionary<int, Modification> { { 4, new Modification("Zzz", null, "Biological", null, motif, "Anywhere.", null, 2.0) } });
        var zzz = new MockSpectralMatch("a.raw", "ACD[Zzz]EF", "ACDEF", 1.0, 1, [zzzForm]);

        // Zzz is seen here; Aaa only elsewhere. Listed, Aaa still comes first.
        var result = ModificationOccupancyCalculator.CalculateParentLevelOccupancy(
            _protein, [zzz], [(4, "Zzz on D"), (4, aaa.IdWithMotif)]);

        Assert.That(result[4].Select(o => o.ModificationIdWithMotif), Is.EqualTo(new[] { "Aaa on D", "Zzz on D" }));
        Assert.That(result[4].Select(o => o.ModifiedCount), Is.EqualTo(new[] { 0, 1 }));
    }

    [Test]
    public void SitesSeenModifiedAreOnlyThePairsSomePsmModifies()
    {
        ModificationMotif.TryGetMotif("D", out var motif);
        var common = new Modification("Oxidation", null, "Common Variable", null, motif, "Anywhere.", null, 15.995);
        var commonForm = new MockBioPolymerWithSetMods("ACDEF", "ACD[Oxidation]EF", _protein, 1, 5,
            new Dictionary<int, Modification> { { 4, common } });

        var sites = ModificationOccupancyCalculator.GetSitesSeenModified(_protein,
        [
            Modified("a.raw", 1),
            Unmodified("b.raw", 2),
            OtherRegion("b.raw", 3),
            new MockSpectralMatch("c.raw", "ACD[Oxidation]EF", "ACDEF", 1.0, 4, [commonForm])
        ]);

        Assert.That(sites, Is.EquivalentTo(new[] { (4, _phospho.IdWithMotif) }));
    }

    #endregion

    #region BioPolymerGroup

    private BioPolymerGroup GroupOf(params ISpectralMatch[] psms)
    {
        var forms = new HashSet<IBioPolymerWithSetMods> { _modified, _unmodified, _otherRegion };
        return new BioPolymerGroup(new HashSet<IBioPolymer> { _protein }, forms, forms)
        {
            AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch>(psms)
        };
    }

    private static SampleGroupResult ResultFor(BioPolymerGroup group, string file) =>
        group.SampleGroupResults!.Single(r => r.Identity == file);

    [Test]
    public void AFileThatCoveredASiteModifiedInAnotherFileReportsZeroOfN()
    {
        // a.raw sees the site modified; b.raw covers it twice, unmodified; c.raw never covers it.
        var group = GroupOf(Modified("a.raw", 1), Unmodified("a.raw", 2),
            Unmodified("b.raw", 3), Unmodified("b.raw", 4), OtherRegion("c.raw", 5));

        group.ReportSitesSeenModifiedInAnySampleGroup();
        group.PopulateSampleGroupResults();

        Assert.That(ResultFor(group, "a.raw").FormatOccupancy(["P00001"]),
            Is.EqualTo("pos3[Phosphorylation on D,info:fraction=0.50(1/2)]"));
        Assert.That(ResultFor(group, "b.raw").FormatOccupancy(["P00001"]),
            Is.EqualTo("pos3[Phosphorylation on D,info:fraction=0.00(0/2)]"));
        Assert.That(ResultFor(group, "c.raw").FormatOccupancy(["P00001"]), Is.Empty, "not covered stays absent");
    }

    [Test]
    public void ASiteModifiedInNoFileIsReportedInNone()
    {
        var group = GroupOf(Unmodified("a.raw", 1), Unmodified("b.raw", 2));

        group.ReportSitesSeenModifiedInAnySampleGroup();
        group.PopulateSampleGroupResults();

        Assert.That(group.SampleGroupResults!.Select(r => r.FormatOccupancy(["P00001"])), Is.All.Empty);
    }

    [Test]
    public void WithoutOptingInTheGroupOmitsACoveredUnmodifiedSiteAsBefore()
    {
        var group = GroupOf(Modified("a.raw", 1), Unmodified("b.raw", 2));

        group.PopulateSampleGroupResults();

        Assert.That(group.OccupancySitesToReport, Is.Null);
        Assert.That(ResultFor(group, "b.raw").FormatOccupancy(["P00001"]), Is.Empty);
    }

    [Test]
    public void ZeroOfNWithNoMeasuredIntensityIsCountedButNeverWrittenAsAnIntensityFraction()
    {
        // b.raw's PSMs carry no intensity: its count cell reports 0/2, and its intensity cell stays empty
        // rather than printing 0/0 as if a zero had been measured.
        var group = GroupOf(Modified("a.raw", 1, 1_000_000), Unmodified("b.raw", 2), Unmodified("b.raw", 3));

        group.ReportSitesSeenModifiedInAnySampleGroup();
        group.PopulateSampleGroupResults();

        var b = ResultFor(group, "b.raw");
        Assert.That(b.FormatOccupancy(["P00001"]), Is.EqualTo("pos3[Phosphorylation on D,info:fraction=0.00(0/2)]"));
        Assert.That(b.FormatOccupancy(["P00001"], intensityBased: true), Is.Empty);
    }

    [Test]
    public void ZeroOfNWithMeasuredIntensityIsWrittenAsAZeroIntensityFraction()
    {
        var group = GroupOf(Modified("a.raw", 1, 1_000_000), Unmodified("b.raw", 2, 2_000_000));

        group.ReportSitesSeenModifiedInAnySampleGroup();
        group.PopulateSampleGroupResults();

        Assert.That(ResultFor(group, "b.raw").FormatOccupancy(["P00001"], intensityBased: true),
            Is.EqualTo("pos3[Phosphorylation on D,info:fraction=0.0000(0/2E+06)]"));
    }

    [Test]
    public void AFileSubsetKeepsTheSitesSeenModifiedInEveryFile()
    {
        var group = GroupOf(Modified("a.raw", 1), Unmodified("b.raw", 2));
        group.ReportSitesSeenModifiedInAnySampleGroup();

        var subset = (BioPolymerGroup)group.ConstructSubsetBioPolymerGroup("b.raw");
        subset.PopulateSampleGroupResults();

        Assert.That(subset.OccupancySitesToReport, Is.SameAs(group.OccupancySitesToReport));
        Assert.That(ResultFor(subset, "b.raw").FormatOccupancy(["P00001"]),
            Is.EqualTo("pos3[Phosphorylation on D,info:fraction=0.00(0/1)]"));
    }

    [Test]
    public void SettingTheSitesInvalidatesTheResults()
    {
        var group = GroupOf(Modified("a.raw", 1), Unmodified("b.raw", 2));
        group.PopulateSampleGroupResults();
        Assert.That(group.SampleGroupResults, Is.Not.Null);

        group.ReportSitesSeenModifiedInAnySampleGroup();

        Assert.That(group.SampleGroupResults, Is.Null);
    }

    #endregion
}
