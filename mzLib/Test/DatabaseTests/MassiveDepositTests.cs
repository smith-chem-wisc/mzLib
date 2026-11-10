using MzLibUtil;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests;

/// <summary>
/// M1b: what a MassIVE deposit's listing and params.xml say about it: its runs (<see cref="MassiveRuns"/>), its
/// deposited SDRF and design-table candidates (<see cref="MassiveFileExtensions"/>), and its submission form
/// (<see cref="MassiveParameters"/>). The listings are MassiveClient's own, captured 2026-11-10 (MassiveFixtures/PROVENANCE.md).
/// </summary>
[TestFixture]
[ExcludeFromCodeCoverage]
public class MassiveDepositTests
{
    private static readonly Dictionary<string, string> Roots = new()
    {
        ["MSV000093798"] = "/v06/MSV000093798", ["MSV000080763"] = "/x01/MSV000080763", ["MSV000098930"] = "/v10/MSV000098930",
        ["MSV000090777"] = "/v05/MSV000090777", ["MSV000096226"] = "/v07/MSV000096226", ["MSV000092529"] = "/v06/MSV000092529",
        ["MSV000096205"] = "/v07/MSV000096205", ["MSV000100606"] = "/v12/MSV000100606", ["MSV000102340"] = "/v13/MSV000102340",
        ["MSV000103504"] = "/v15/MSV000103504", ["MSV000102494"] = "/v13/MSV000102494", ["MSV000094031"] = "/v07/MSV000094031",
        ["MSV000097296"] = "/v09/MSV000097296",
    };

    private static string FixturePath(string name) =>
        Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "MassiveFixtures", name);

    /// <summary>A captured listing, as MassiveClient returned it.</summary>
    private static List<MassiveFile> Listing(string msv) =>
        File.ReadAllLines(FixturePath(Path.Combine("listings", msv + ".listing.tsv")), Encoding.UTF8)
            .Skip(1)
            .Where(line => line.Length > 0)
            .Select(line => line.Split('\t'))
            .Select(c => new MassiveFile(msv, Roots[msv], c[0],
                c[1].Length == 0 ? null : long.Parse(c[1], CultureInfo.InvariantCulture),
                c[3].Length == 0 ? null : DateTime.ParseExact(c[3], "yyyy-MM-dd HH:mm", CultureInfo.InvariantCulture),
                isLink: c[2] == "1"))
            .ToList();

    private static MassiveFile File_(string path, long size = 1) => new("MSV000000001", "/v01/MSV000000001", path, size);

    // ---- the run rule ------------------------------------------------------------------------------------------

    [Test]
    public void Runs_ThermoDeposit_EveryRawFileIsARun_TheLinkIsNot()
    {
        IReadOnlyList<MassiveRun> runs = MassiveRuns.Of(Listing("MSV000093798"));

        Assert.Multiple(() =>
        {
            Assert.That(runs, Has.Count.EqualTo(192));
            Assert.That(runs.All(r => r.Kind == MassiveRunKind.Vendor && !r.IsFolder && r.FileName.EndsWith(".raw")), Is.True);
            Assert.That(runs.All(r => r.Collection == "raw" && r.Files.Count == 1), Is.True);
            Assert.That(runs.Any(r => r.Path.StartsWith("ccms_peak")), Is.False);
        });
    }

    [Test]
    public void Runs_DepositWhoseSpectraAreOnlyInCcmsPeak_HasNone()
    {
        // ccms_peak is a link MassIVE's FTP server will not open; the originals are packed in search/SEARCH/Saliva_txt.rar.
        Assert.That(MassiveRuns.Of(Listing("MSV000080763")), Is.Empty);
    }

    [Test]
    public void Runs_RawCollectionZips_AreTheRuns_ASearchZipIsNot()
    {
        IReadOnlyList<MassiveRun> runs = MassiveRuns.Of(Listing("MSV000098930"));

        Assert.Multiple(() =>
        {
            Assert.That(runs, Has.Count.EqualTo(508));
            Assert.That(runs.All(r => r.Kind == MassiveRunKind.Archive && r.Path.StartsWith("raw/RAW/") && r.FileName.EndsWith(".zip")), Is.True);
            Assert.That(runs.Any(r => r.Path == "search/20181028_021159_SEC_HeLa_CCL2_Kyoto.zip"), Is.False);
        });
    }

    [Test]
    public void Runs_SciexDeposit_OneRunPerStem_TheWiffPrimary_ThreeMembers()
    {
        List<MassiveFile> listing = Listing("MSV000090777");
        IReadOnlyList<MassiveRun> runs = MassiveRuns.Of(listing);

        MassiveRun first = runs.First(r => r.FileName == "220314_0008_SSL2_30_ND_INJ_4_L_iRT_1000ng_120min_01.wiff");
        Assert.Multiple(() =>
        {
            Assert.That(runs, Has.Count.EqualTo(50));
            Assert.That(runs.All(r => r.Kind == MassiveRunKind.Vendor && r.Primary!.FileName.EndsWith(".wiff")), Is.True);
            Assert.That(runs.All(r => r.Files.Count == 4), Is.True, "the .wiff, its .wiff.scan, its .wiff2 and its .timeseries.data");
            Assert.That(first.Members.Select(m => Path.GetExtension(m.FileName)), Is.EquivalentTo(new[] { ".scan", ".wiff2", ".data" }));
            Assert.That(runs.SelectMany(r => r.Files).Any(f => f.FileName == ".DS_Store"), Is.False);
            Assert.That(runs.Select(r => r.Path).Distinct().Count(), Is.EqualTo(50));
        });
    }

    [Test]
    public void Runs_BrukerFolders_AreOneRunEach_WithEveryFileBelow()
    {
        List<MassiveFile> listing = Listing("MSV000092529");
        IReadOnlyList<MassiveRun> runs = MassiveRuns.Of(listing);

        Assert.Multiple(() =>
        {
            Assert.That(runs, Has.Count.EqualTo(6));
            Assert.That(runs.All(r => r.IsFolder && r.Primary == null && r.FileName.EndsWith(".d") && r.Kind == MassiveRunKind.Vendor), Is.True);
            Assert.That(runs.Sum(r => r.Files.Count), Is.EqualTo(listing.Count(f => f.Path.Contains(".d/"))));
            Assert.That(runs[0].Path, Is.EqualTo("raw/230515_untreatedBMVEC_rep1.2_A5_1_5516.d"));
            Assert.That(runs[0].SizeBytes, Is.EqualTo(runs[0].Files.Sum(f => f.SizeBytes)));
        });
    }

    [Test]
    public void Runs_WatersRawFolder_IsOneRun()
    {
        // Synthetic: no Waters .raw folder is among sdrf's 121 deposits; this is the shape MassLynx writes.
        var listing = new[] { File_("raw/Sample A.raw/_FUNC001.DAT"), File_("raw/Sample A.raw/_FUNC001.IDX"), File_("raw/Sample A.raw/_HEADER.TXT") };

        IReadOnlyList<MassiveRun> runs = MassiveRuns.Of(listing);

        Assert.That(runs.Single().Path, Is.EqualTo("raw/Sample A.raw"));
        Assert.That(runs.Single().Files, Has.Count.EqualTo(3));
    }

    [Test]
    public void Runs_PeakListsOnly_OneRunPerAcquisition_MzmlFirst()
    {
        // Synthetic, in the shape of a peak-list-only deposit.
        var listing = new[] { File_("peak/a.mgf"), File_("peak/a.mzML"), File_("peak/b.mzXML.gz"), File_("search/results.mzid"), File_("peak/notes.txt") };

        IReadOnlyList<MassiveRun> runs = MassiveRuns.Of(listing);

        Assert.Multiple(() =>
        {
            Assert.That(runs.Select(r => r.Path), Is.EqualTo(new[] { "peak/a.mzML", "peak/b.mzXML.gz" }));
            Assert.That(runs[0].Members.Single().FileName, Is.EqualTo("a.mgf"));
            Assert.That(runs.All(r => r.Kind == MassiveRunKind.PeakList), Is.True);
        });
    }

    [Test]
    public void Runs_VendorFilesPresent_PeakListsAndZipsAreNotRuns()
    {
        var listing = new[] { File_("raw/x.raw"), File_("raw/y.zip"), File_("peak/z.mzML"), File_("raw/x.mzML") };

        IReadOnlyList<MassiveRun> runs = MassiveRuns.Of(listing);

        Assert.That(runs.Single().Path, Is.EqualTo("raw/x.raw"));
        Assert.That(runs.Single().Members.Single().FileName, Is.EqualTo("x.mzML"), "another format of the same acquisition is a member");
    }

    [Test]
    public void Runs_SearchOnlyDeposit_HasNone()
    {
        Assert.That(MassiveRuns.Of(new[] { File_("search/out.mzid"), File_("ccms_parameters/params.xml") }), Is.Empty);
    }

    [Test]
    public void Runs_NullListingOrEntry_IsRefused()
    {
        Assert.Throws<ArgumentNullException>(() => MassiveRuns.Of(null!));
        Assert.Throws<ArgumentNullException>(() => MassiveRuns.Of(new MassiveFile[] { null! }));
    }

    [TestCase("X.timeseries.data", "X")]
    [TestCase("X.wiff.scan", "X")]
    [TestCase("X.wiff2", "X")]
    [TestCase("X.raw.zip", "X")]
    [TestCase("X.mzML.gz", "X")]
    public void RunStem_SharesTheDraftersStem_PlusTheTimeseriesSidecar(string fileName, string stem)
    {
        Assert.That(MassiveRuns.StemOf(fileName), Is.EqualTo(stem));
    }

    [Test]
    public void RunStem_IsTheDraftersStem_ForEveryDataExtensionItKnows()
    {
        // The run rule must name runs exactly as sdrf's drafter names files, so it uses the drafter's own list.
        var extensions = (string[])typeof(SdrfFileNamePattern)
            .GetField("DataExtensions", BindingFlags.NonPublic | BindingFlags.Static)!.GetValue(null)!;

        Assert.That(extensions, Is.Not.Empty);
        foreach (string extension in extensions)
            Assert.That(MassiveRuns.StemOf("Sample_1" + extension), Is.EqualTo(SdrfFileNamePattern.Stem("Sample_1" + extension)), extension);
    }

    // ---- deposited SDRF and design-table candidates ------------------------------------------------------------

    [TestCase("MSV000096205", "metadata/SBF912MAG_samplesheet.sdrf.tsv")]
    [TestCase("MSV000100606", "metadata/PATL2_KO.sdrf.tsv")]
    [TestCase("MSV000102340", "metadata/sdrf.tsv")]
    [TestCase("MSV000103504", "metadata/sdrf.tsv")]
    [TestCase("MSV000102494", "metadata/SDRF_Metadata.txt")]
    [TestCase("MSV000094031", "metadata/decellularization manuscript.sdrf.xlsx")]
    [TestCase("MSV000097296", "metadata/SDRF_submission.xlsx")]
    public void DepositedSdrf_FindsEachKnownOne_AndNothingElse(string msv, string path)
    {
        Assert.That(Listing(msv).WhereDepositedSdrf().Select(f => f.Path), Is.EqualTo(new[] { path }));
    }

    [TestCase("MSV000090777", 27)]
    [TestCase("MSV000096226", 2)]
    [TestCase("MSV000098930", 1)]
    [TestCase("MSV000093798", 0)]
    public void DesignTableCandidates_InEveryCollection_NoSdrfAmongThem(string msv, int expected)
    {
        List<MassiveFile> listing = Listing(msv);
        List<MassiveFile> candidates = listing.WhereDesignTableCandidate().ToList();

        Assert.Multiple(() =>
        {
            Assert.That(candidates, Has.Count.EqualTo(expected));
            Assert.That(listing.WhereDepositedSdrf(), Is.Empty);
        });
    }

    [Test]
    public void DesignTableCandidate_AnnotationBesideAnSdrf_IsACandidate_NotAnSdrf()
    {
        MassiveFile annotation = Listing("MSV000103504").Single(f => f.Path == "metadata/experiment_annotation.tsv");

        Assert.That(annotation.IsDesignTableCandidate(), Is.True);
        Assert.That(annotation.IsDepositedSdrf(), Is.False);
    }

    [TestCase("metadata/run.sdrf.tsv", true, false)]
    [TestCase("metadata/SDRF.TSV", true, false)]
    [TestCase("metadata/x.sdrf.csv", false, false)]
    [TestCase("metadata/annotation.sdrf.tsv", true, false)]
    [TestCase("metadata/experimentalDesign.txt", false, true)]
    [TestCase("search/combined/txt/experimentalDesignTemplate.txt", false, true)]
    [TestCase("metadata/experimental.design.txt", false, true)]
    [TestCase("quant/a_ConditionSetup.tsv", false, true)]
    [TestCase("quant/a_ConditionSetup.xlsx", false, false)]
    [TestCase("metadata/annotation.xlsx", false, false)]
    [TestCase("raw/x.raw", false, false)]
    public void NameRules_ByName(string path, bool sdrf, bool candidate)
    {
        MassiveFile file = File_(path);

        Assert.That(file.IsDepositedSdrf(), Is.EqualTo(sdrf));
        Assert.That(file.IsDesignTableCandidate(), Is.EqualTo(candidate));
    }

    [TestCase("MSV000090777", 27, 1)]
    [TestCase("MSV000092529", 1, 0)]
    [TestCase("MSV000096226", 2, 0)]
    [TestCase("MSV000093798", 0, 0)]
    [TestCase("MSV000098930", 0, 0)]
    [TestCase("MSV000080763", 0, 0)]
    [TestCase("MSV000103504", 0, 0)]
    public void SpectronautOutput_InRealListings(string msv, int analyses, int libraries)
    {
        List<MassiveFile> listing = Listing(msv);
        List<MassiveFile> output = listing.WhereSpectronautOutput().ToList();

        Assert.Multiple(() =>
        {
            Assert.That(output.Count(f => f.FileName.EndsWith(".sne")), Is.EqualTo(analyses));
            Assert.That(output.Count(f => f.FileName.EndsWith(".kit")), Is.EqualTo(libraries));
            Assert.That(output, Has.Count.EqualTo(analyses + libraries));
            Assert.That(listing.HasSpectronautOutput(), Is.EqualTo(analyses + libraries > 0));
        });
    }

    [Test]
    public void SpectronautOutput_InAnyCollection_AsAging037Found()
    {
        // MSV000090777's library sits under library/, as aging's four datasets keep theirs; one of MSV000096226's
        // analyses sits under raw/, so the rule does not look at the collection.
        Assert.Multiple(() =>
        {
            Assert.That(Listing("MSV000090777").WhereSpectronautOutput().Single(f => f.FileName.EndsWith(".kit")).Path,
                Is.EqualTo("library/Spectral_Library/2022_0513_SSL2_DDA_Zeno_120min_spectral_library.kit"));
            Assert.That(Listing("MSV000096226").WhereSpectronautOutput().Select(f => f.Collection), Is.EquivalentTo(new[] { "raw", "search" }));
        });
    }

    [TestCase("search/x.sne", true)]
    [TestCase("library/X.KIT", true)]
    [TestCase("search/x.sne.zip", false)]
    [TestCase("search/x_Report.tsv", false)]
    [TestCase("raw/x.raw", false)]
    public void SpectronautOutput_ByName(string path, bool expected)
    {
        Assert.That(File_(path).IsSpectronautOutput(), Is.EqualTo(expected));
    }

    [Test]
    public void SpectronautOutput_ALinkIsNot_AndNullsAreRefused()
    {
        Assert.That(new MassiveFile("MSV000000001", "/v01/MSV000000001", "x.sne", null, isLink: true).IsSpectronautOutput(), Is.False);
        Assert.Throws<ArgumentNullException>(() => ((MassiveFile)null!).IsSpectronautOutput());
        Assert.Throws<ArgumentNullException>(() => ((IEnumerable<MassiveFile>)null!).HasSpectronautOutput());
    }

    [Test]
    public void NameRules_ALinkIsNeither()
    {
        var link = new MassiveFile("MSV000000001", "/v01/MSV000000001", "metadata_sdrf.tsv", null, isLink: true);

        Assert.That(link.IsDepositedSdrf(), Is.False);
        Assert.That(link.IsDesignTableCandidate(), Is.False);
        Assert.Throws<ArgumentNullException>(() => ((MassiveFile)null!).IsDepositedSdrf());
        Assert.Throws<ArgumentNullException>(() => ((MassiveFile)null!).IsDesignTableCandidate());
        Assert.Throws<ArgumentNullException>(() => ((IEnumerable<MassiveFile>)null!).WhereDepositedSdrf().ToList());
        Assert.Throws<ArgumentNullException>(() => ((IEnumerable<MassiveFile>)null!).WhereDesignTableCandidate().ToList());
    }

    // ---- params.xml --------------------------------------------------------------------------------------------

    [Test]
    public void Params_Dataset_ReadsTheMappingsSpeciesAndInstrument()
    {
        MassiveParameters p = MassiveParameters.Read(FixturePath("MSV000080763.params.xml"));

        Assert.Multiple(() =>
        {
            Assert.That(p.Workflow, Is.EqualTo("MASSIVE-COMPLETE"));
            Assert.That(p.UploadFileMappings, Has.Count.EqualTo(31), "consecutive mappings must not be skipped");
            Assert.That(p.UploadFileMappings![0], Is.EqualTo(new MassiveUploadMapping("RAW-00000.raw",
                "ccms/dataset/pride_imports/PXD004319/RAW/20150722_QE5_UPLC8_RJC_COLLAB_DB_4075_01.raw")));
            Assert.That(p.UploadFileMappings[30], Is.EqualTo(new MassiveUploadMapping("SEARCH-00000.rar",
                "ccms/dataset/pride_imports/PXD004319/SEARCH/Saliva_txt.rar")));
            Assert.That(p.Species, Is.EqualTo(new[] { "NCBITaxon:9606" }));
            Assert.That(p.Instruments, Is.EqualTo(new[] { "MS:1001911" }));
        });
    }

    [Test]
    public void Params_LargerDataset_MapsEveryKindOfFile()
    {
        using FileStream stream = File.OpenRead(FixturePath("MSV000096226.params.xml"));
        MassiveParameters p = MassiveParameters.Read(stream);

        Assert.Multiple(() =>
        {
            Assert.That(p.UploadFileMappings, Has.Count.EqualTo(85));
            Assert.That(p.UploadFileMappings!.Single(m => m.StoredName == "metadata-00000.tsv").UploadPath, Does.EndWith("_ConditionSetup.tsv"));
            Assert.That(p.UploadFileMappings.Single(m => m.StoredName == "methods-00000.docx").UploadPath, Does.EndWith("/M&M_Ub.docx"));
            Assert.That(p.Species, Is.EqualTo(new[] { "NCBITaxon:10090" }));
            Assert.That(p.Instruments, Is.EqualTo(new[] { "MS:1003028" }));
        });
    }

    [Test]
    public void Params_Update_HasItsOwnWorkflow_AndStatesNoSpecies()
    {
        MassiveParameters p = MassiveParameters.Read(FixturePath("MSV000090777_update.params.xml"));

        Assert.Multiple(() =>
        {
            Assert.That(p.Workflow, Is.EqualTo("UPDATE-MASSIVE-DATASET"));
            Assert.That(p.UploadFileMappings, Has.Count.EqualTo(6));
            Assert.That(p.Species, Is.Null, "an update's form does not restate the species");
            Assert.That(p.Instruments, Is.Null);
        });
    }

    [TestCase("MSV000080763.params.xml")]
    [TestCase("MSV000096226.params.xml")]
    [TestCase("MSV000090777_update.params.xml")]
    public void Params_NoEmailNameOrPassword_LeavesTheParser(string fixture)
    {
        string raw = File.ReadAllText(FixturePath(fixture));
        Assert.That(raw, Does.Contain("@example.org"), "the fixture must hold the fakes, or this test proves nothing");

        MassiveParameters p = MassiveParameters.Read(FixturePath(fixture));
        IEnumerable<string> everything = new[] { p.Workflow }
            .Concat(p.UploadFileMappings?.SelectMany(m => new[] { m.StoredName, m.UploadPath }) ?? Enumerable.Empty<string>())
            .Concat(p.Species ?? Enumerable.Empty<string>())
            .Concat(p.Instruments ?? Enumerable.Empty<string>())
            .Where(v => v != null)!;

        Assert.That(everything.Any(v => v.Contains('@') || v.Contains("Redacted") || v.Contains("redacted-password")), Is.False);
        Assert.That(typeof(MassiveParameters).GetProperties().Select(pr => pr.Name),
            Is.EquivalentTo(new[] { "Workflow", "UploadFileMappings", "Species", "Instruments" }), "no property may hold anything else");
    }

    [TestCase("<foo/>")]
    [TestCase("not xml")]
    [TestCase("<?xml version=\"1.0\"?><!DOCTYPE parameters [<!ENTITY x SYSTEM \"file:///c:/windows/win.ini\">]><parameters><parameter name=\"workflow\">&x;</parameter></parameters>")]
    [TestCase("<parameters><parameter name=\"upload_file_mapping\">no bar here</parameter></parameters>")]
    public void Params_NotAParametersFile_IsAnMzLibException(string xml)
    {
        Assert.Throws<MzLibException>(() => MassiveParameters.Read(new MemoryStream(Encoding.UTF8.GetBytes(xml))));
    }

    [Test]
    public void Params_ParametersWithNothingBetweenThem_AreEachRead()
    {
        // The real files put a line break between elements; a reader that moved past an element and then read again
        // would skip every second one only when there is none. This pins it.
        const string xml = "<parameters><parameter name=\"upload_file_mapping\">A-0.raw|a.raw</parameter>"
                         + "<parameter name=\"upload_file_mapping\">B-0.raw|b.raw</parameter>"
                         + "<parameter name=\"upload_file_mapping\">C-0.raw|c.raw</parameter>"
                         + "<parameter name=\"workflow\">MASSIVE-COMPLETE</parameter></parameters>";

        MassiveParameters p = MassiveParameters.Read(new MemoryStream(Encoding.UTF8.GetBytes(xml)));

        Assert.That(p.UploadFileMappings!.Select(m => m.StoredName), Is.EqualTo(new[] { "A-0.raw", "B-0.raw", "C-0.raw" }));
        Assert.That(p.Workflow, Is.EqualTo("MASSIVE-COMPLETE"));
    }

    [Test]
    public void Params_Arguments_AreChecked()
    {
        Assert.Throws<ArgumentNullException>(() => MassiveParameters.Read((Stream)null!));
        Assert.Throws<ArgumentException>(() => MassiveParameters.Read(" "));
    }
}
