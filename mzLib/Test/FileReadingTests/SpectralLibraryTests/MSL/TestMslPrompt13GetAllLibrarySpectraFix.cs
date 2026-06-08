using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using SpectralLibraryClass = Readers.SpectralLibrary.SpectralLibrary;

namespace Test.MslSpectralLibrary;

/// <summary>
/// NUnit 4 regression tests for the fix to
/// <see cref="SpectralLibrary.GetAllLibrarySpectra"/> that caused .msl inputs
/// to be silently excluded from bulk enumeration.
///
/// Coverage matrix:
///   • GetAllLibrarySpectra with .msl-only input  → must yield all spectra
///   • LoadResults           with .msl-only input  → Results must be non-empty
///   • GetAllLibrarySpectra with mixed .msl + .msp → all spectra from both sources
///   • LoadResults           with mixed input       → Results count = msl + msp totals
///   • Decoys included by default in GetAllLibrarySpectra via .msl path
///   • Sequence/Name round-trip for .msl-sourced spectra
/// </summary>
[TestFixture]
public sealed class TestMslPrompt13GetAllLibrarySpectraFix
{
    // ── Fixture paths ─────────────────────────────────────────────────────────

    private static readonly string OutputDir =
        Path.Combine(Path.GetTempPath(), "MslPrompt13Tests");

    [OneTimeSetUp]
    public void OneTimeSetUp() => Directory.CreateDirectory(OutputDir);

    [OneTimeTearDown]
    public void OneTimeTearDown()
    {
        if (Directory.Exists(OutputDir))
            Directory.Delete(OutputDir, recursive: true);
    }

    // ── Helpers ───────────────────────────────────────────────────────────────

    /// <summary>
    /// Builds a list of MslLibraryEntry values with distinct sequences and optional decoys.
    /// </summary>
    private static List<MslLibraryEntry> BuildEntries(
        IEnumerable<(string seq, int charge, bool isDecoy)> specs)
    {
        var entries = new List<MslLibraryEntry>();
        foreach (var (seq, charge, isDecoy) in specs)
        {
            entries.Add(new MslLibraryEntry
            {
                FullSequence = seq,
                BaseSequence = seq,
                PrecursorMz = 449.74,
                ChargeState = charge,
                IsDecoy = isDecoy,
                MatchedFragmentIons = new List<MslFragmentIon>
                {
                    new MslFragmentIon
                    {
                        Mz             = 175.119f,
                        Intensity      = 1.0f,
                        ProductType    = ProductType.y,
                        FragmentNumber = 1,
                        Charge         = 1
                    }
                }
            });
        }
        return entries;
    }

    /// <summary>Saves entries to a .msl file and returns the path.</summary>
    private static string SaveMsl(string stem, List<MslLibraryEntry> entries)
    {
        string path = Path.Combine(OutputDir, stem + ".msl");
        MslLibrary.Save(path, entries);
        return path;
    }

    /// <summary>Writes a minimal MSP file with one entry and returns the path.</summary>
    private static string WriteMsp(string stem, string sequence = "MSPEPTIDE", int charge = 2)
    {
        string path = Path.Combine(OutputDir, stem + ".msp");
        File.WriteAllText(path, string.Join("\n",
            $"Name: {sequence}/{charge}",
            "MW: 530.26",
            "Comment: Parent=530.26 iRT=40.0",
            "Num peaks: 1",
            "120.081\t1.0\t\"b2^1/0ppm\""
        ));
        return path;
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Core regression: .msl-only input
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// REGRESSION: Before the fix, GetAllLibrarySpectra() returned 0 items for an
    /// .msl-only SpectralLibrary.  After the fix it must return all entries.
    /// </summary>
    [Test]
    public void GetAllLibrarySpectra_MslOnly_YieldsAllSpectra()
    {
        var entries = BuildEntries(new[]
        {
            ("PEPTIDE",  2, false),
            ("ACDEFG",   3, false),
        });
        string mslPath = SaveMsl("get_all_msl_only", entries);

        using var lib = new SpectralLibraryClass(new List<string> { mslPath });
        var all = lib.GetAllLibrarySpectra().ToList();

        Assert.That(all, Has.Count.EqualTo(2),
            "GetAllLibrarySpectra must yield 2 spectra from an .msl-only library.");

        lib.CloseConnections();
    }

    /// <summary>
    /// REGRESSION: Before the fix, LoadResults() left Results empty for .msl-only input.
    /// After the fix Results.Count must equal the number of entries in the .msl file.
    /// </summary>
    [Test]
    public void LoadResults_MslOnly_ResultsNonEmpty()
    {
        var entries = BuildEntries(new[]
        {
            ("LOADPEP",  2, false),
            ("LOADPEP2", 3, false),
        });
        string mslPath = SaveMsl("load_results_msl_only", entries);

        var lib = new SpectralLibraryClass(new List<string> { mslPath });
        lib.LoadResults();

        Assert.That(lib.Results, Has.Count.EqualTo(2),
            "LoadResults must populate Results for an .msl-only SpectralLibrary.");

        lib.CloseConnections();
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Mixed input: .msl + .msp
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// A SpectralLibrary constructed with one .msl (2 entries) and one .msp (1 entry)
    /// must yield 3 spectra from GetAllLibrarySpectra.
    /// </summary>
    [Test]
    public void GetAllLibrarySpectra_MixedMslAndMsp_YieldsAll()
    {
        var mslEntries = BuildEntries(new[]
        {
            ("MIXPEP1", 2, false),
            ("MIXPEP2", 3, false),
        });
        string mslPath = SaveMsl("mixed_msl", mslEntries);
        string mspPath = WriteMsp("mixed_msp", "MSPEPTIDE", 2);

        using var lib = new SpectralLibraryClass(new List<string> { mslPath, mspPath });
        var all = lib.GetAllLibrarySpectra().ToList();

        Assert.That(all, Has.Count.EqualTo(3),
            "GetAllLibrarySpectra must yield spectra from both .msl and .msp sources.");

        lib.CloseConnections();
    }

    /// <summary>
    /// LoadResults on a mixed .msl + .msp SpectralLibrary must populate Results with
    /// spectra from both sources (msl count + msp count).
    /// </summary>
    [Test]
    public void LoadResults_MixedMslAndMsp_ResultsCountIsSum()
    {
        var mslEntries = BuildEntries(new[]
        {
            ("LRPEP1", 2, false),
            ("LRPEP2", 3, false),
        });
        string mslPath = SaveMsl("lr_mixed_msl", mslEntries);
        string mspPath = WriteMsp("lr_mixed_msp", "LRMSP", 2);

        var lib = new SpectralLibraryClass(new List<string> { mslPath, mspPath });
        lib.LoadResults();

        Assert.That(lib.Results, Has.Count.EqualTo(3),
            "LoadResults must include spectra from both .msl and .msp inputs.");

        lib.CloseConnections();
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Sequence / Name round-trip for .msl-sourced spectra
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// Spectra yielded by GetAllLibrarySpectra from an .msl source must have the correct
    /// Sequence and Name (Sequence/Charge) values — confirming ToLibrarySpectrum() is
    /// called correctly and fields are not lost.
    /// </summary>
    [Test]
    public void GetAllLibrarySpectra_MslEntry_NameRoundTripCorrect()
    {
        var entries = BuildEntries(new[] { ("ROUNDTRIP", 2, false) });
        string mslPath = SaveMsl("name_roundtrip", entries);

        using var lib = new SpectralLibraryClass(new List<string> { mslPath });
        var spectrum = lib.GetAllLibrarySpectra().Single();

        Assert.That(spectrum.Sequence, Is.EqualTo("ROUNDTRIP"),
            "Sequence from .msl-sourced LibrarySpectrum must match the original FullSequence.");
        Assert.That(spectrum.Name, Is.EqualTo("ROUNDTRIP/2"),
            "Name must be 'Sequence/Charge' for .msl-sourced LibrarySpectrum.");

        lib.CloseConnections();
    }

    /// <summary>
    /// Fragment ions must be present on spectra yielded via GetAllLibrarySpectra from .msl.
    /// </summary>
    [Test]
    public void GetAllLibrarySpectra_MslEntry_FragmentIonsPreserved()
    {
        var entries = BuildEntries(new[] { ("FRAGMENTS", 2, false) });
        string mslPath = SaveMsl("fragments_present", entries);

        using var lib = new SpectralLibraryClass(new List<string> { mslPath });
        var spectrum = lib.GetAllLibrarySpectra().Single();

        Assert.That(spectrum.MatchedFragmentIons, Is.Not.Empty,
            "MatchedFragmentIons must not be empty for .msl-sourced spectra.");

        lib.CloseConnections();
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Decoy handling
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// GetAllLibrarySpectra must include decoy entries from .msl libraries.
    /// MslLibrary.GetAllEntries(includeDecoys: true) is the code path under test.
    /// </summary>
    [Test]
    public void GetAllLibrarySpectra_MslOnly_IncludesDecoys()
    {
        var entries = BuildEntries(new[]
        {
            ("TARGET", 2, false),
            ("DECOY",  2, true),
        });
        string mslPath = SaveMsl("decoys_included", entries);

        using var lib = new SpectralLibraryClass(new List<string> { mslPath });
        var all = lib.GetAllLibrarySpectra().ToList();

        Assert.That(all, Has.Count.EqualTo(2),
            "GetAllLibrarySpectra must include decoy entries from .msl sources.");

        lib.CloseConnections();
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Regression guard: MSP-only behaviour is unchanged
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// Confirms that the fix does not break the existing MSP-only path.
    /// GetAllLibrarySpectra on an MSP-only library must still yield all MSP entries.
    /// </summary>
    [Test]
    public void GetAllLibrarySpectra_MspOnly_StillWorks()
    {
        string mspPath = WriteMsp("msp_only_regression", "MSPONLY", 2);

        using var lib = new SpectralLibraryClass(new List<string> { mspPath });
        var all = lib.GetAllLibrarySpectra().ToList();

        Assert.That(all, Has.Count.EqualTo(1),
            "GetAllLibrarySpectra must still yield MSP entries after the fix.");

        lib.CloseConnections();
    }

    /// <summary>
    /// LoadResults on an MSP-only library must still populate Results correctly.
    /// </summary>
    [Test]
    public void LoadResults_MspOnly_StillPopulatesResults()
    {
        string mspPath = WriteMsp("load_results_msp_regression", "MSPREG", 2);

        var lib = new SpectralLibraryClass(new List<string> { mspPath });
        lib.LoadResults();

        Assert.That(lib.Results, Has.Count.EqualTo(1),
            "LoadResults MSP-only regression: Results must contain 1 entry.");

        lib.CloseConnections();
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Edge case: empty .msl library
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// An .msl file with zero entries must not throw and must yield an empty sequence.
    /// </summary>
    [Test]
    public void GetAllLibrarySpectra_EmptyMsl_YieldsNothing()
    {
        string mslPath = SaveMsl("empty_msl", new List<MslLibraryEntry>());

        using var lib = new SpectralLibraryClass(new List<string> { mslPath });
        var all = lib.GetAllLibrarySpectra().ToList();

        Assert.That(all, Is.Empty,
            "An empty .msl library must yield zero spectra from GetAllLibrarySpectra.");

        lib.CloseConnections();
    }

    // ══════════════════════════════════════════════════════════════════════════
    // Multiple .msl files
    // ══════════════════════════════════════════════════════════════════════════

    /// <summary>
    /// Two .msl files with 2 entries each must yield 4 spectra in total.
    /// </summary>
    [Test]
    public void GetAllLibrarySpectra_TwoMslFiles_YieldsAllFromBoth()
    {
        var entries1 = BuildEntries(new[] { ("PEPFILE1A", 2, false), ("PEPFILE1B", 3, false) });
        var entries2 = BuildEntries(new[] { ("PEPFILE2A", 2, false), ("PEPFILE2B", 3, false) });
        string mslPath1 = SaveMsl("two_msl_file1", entries1);
        string mslPath2 = SaveMsl("two_msl_file2", entries2);

        using var lib = new SpectralLibraryClass(new List<string> { mslPath1, mslPath2 });
        var all = lib.GetAllLibrarySpectra().ToList();

        Assert.That(all, Has.Count.EqualTo(4),
            "GetAllLibrarySpectra must yield entries from all .msl files provided.");

        lib.CloseConnections();
    }
}