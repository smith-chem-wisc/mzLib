using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;

namespace Test.MslSpectralLibrary;

/// <summary>
/// NUnit 4 test suite for <see cref="MslLibrary.UpdateAndSave"/> and
/// <see cref="MslUpdateResult"/>.
///
/// <para>
/// The fixture uses a five-entry original library with known iRT values. Incoming
/// <see cref="LibrarySpectrum"/> objects are constructed with controlled observed-RT and
/// fragment-ion counts so every code path can be exercised deterministically.
/// </para>
///
/// Coverage groups:
/// <list type="number">
///   <item>Return value — <see cref="MslUpdateResult"/> fields are correct</item>
///   <item>Retained entries — original kept when no incoming spectrum matches</item>
///   <item>Replaced entries — incoming wins when it has more fragment ions</item>
///   <item>Original kept — tie goes to the original library</item>
///   <item>Added entries — novel incoming spectra appended</item>
///   <item>RT normalisation — regression applied when enough anchors exist</item>
///   <item>RT normalisation skipped — below anchor threshold</item>
///   <item>RT identity fallback — degenerate anchor data (all same RT)</item>
///   <item>Duplicate incoming — best spectrum kept per sequence/charge key</item>
///   <item>Empty incoming — output is a verbatim copy of the original library</item>
///   <item>Metadata preservation — IonMobility / ProteinAccession / IsDecoy copied</item>
///   <item>Output file is readable — <see cref="MslLibrary.Load"/> succeeds</item>
///   <item>MslUpdateResult helpers — TotalCount and ToString</item>
///   <item>Guard clauses — null arguments and disposed instance</item>
///   <item>minRegressionAnchors clamping — values below 2 clamped to 2</item>
/// </list>
/// </summary>
[TestFixture]
public sealed class TestMslUpdateAndSave
{
	// ── Fixture paths ─────────────────────────────────────────────────────────

	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "MslUpdateAndSaveTests");

	private static string SharedLibraryPath =>
		Path.Combine(OutputDirectory, "original_fixture.msl");

	private static string TempPath(string testName) =>
		Path.Combine(OutputDirectory, testName + ".msl");

	// ── Fixture constants ─────────────────────────────────────────────────────

	// The five original entries have these iRT values (used to verify regression output).
	// Observed RT values in incoming spectra are chosen so that the regression
	// iRT = 2.0 × observedRT + 1.0 holds exactly for all five anchor spectra.
	private const double IrtA = 11.0;   // AAAA/2  — observedRT = 5.0
	private const double IrtB = 21.0;   // BBBB/2  — observedRT = 10.0
	private const double IrtC = 31.0;   // CCCC/2  — observedRT = 15.0
	private const double IrtD = 41.0;   // DDDD/2  — observedRT = 20.0
	private const double IrtE = 51.0;   // EEEE/2  — observedRT = 25.0

	private const double RtSlope = 2.0;
	private const double RtIntercept = 1.0;

	// ── One-time setup / teardown ─────────────────────────────────────────────

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		Directory.CreateDirectory(OutputDirectory);
		MslLibrary.Save(SharedLibraryPath, BuildOriginalEntries());
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDirectory))
			Directory.Delete(OutputDirectory, recursive: true);
	}

	// ── Fixture data factories ────────────────────────────────────────────────

	/// <summary>
	/// Five-entry original library. Each entry has a known iRT, 2 fragment ions,
	/// and metadata fields (IonMobility, ProteinAccession, IsDecoy) that are checked
	/// by the metadata-preservation tests.
	/// </summary>
	private static List<MslLibraryEntry> BuildOriginalEntries() =>
		new()
		{
			MakeEntry("AAAA", 2, mz: 200.0, irt: IrtA, fragCount: 2,
				ionMobility: 0.85, protein: "P00001", isDecoy: false),
			MakeEntry("BBBB", 2, mz: 250.0, irt: IrtB, fragCount: 2,
				ionMobility: 0.90, protein: "P00002", isDecoy: false),
			MakeEntry("CCCC", 2, mz: 300.0, irt: IrtC, fragCount: 2,
				ionMobility: 0.95, protein: "P00003", isDecoy: false),
			MakeEntry("DDDD", 2, mz: 350.0, irt: IrtD, fragCount: 2,
				ionMobility: 1.00, protein: "P00004", isDecoy: false),
			MakeEntry("EEEE", 2, mz: 400.0, irt: IrtE, fragCount: 2,
				ionMobility: 1.05, protein: "P00005", isDecoy: true),
		};

	/// <summary>
	/// Builds a <see cref="LibrarySpectrum"/> for an incoming search result.
	/// <paramref name="observedRt"/> is the raw chromatographic RT from the search run.
	/// <paramref name="fragCount"/> controls how many fragment ions the spectrum carries.
	/// </summary>
	private static LibrarySpectrum MakeSpectrum(string seq, int charge, double precMz,
		double observedRt, int fragCount)
	{
		var ions = MakeMatchedIons(fragCount);
		return new LibrarySpectrum(seq, precMz, charge, ions, observedRt);
	}

	/// <summary>
	/// Builds a list of <paramref name="count"/> y-ions with dummy m/z values.
	/// </summary>
	private static List<MatchedFragmentIon> MakeMatchedIons(int count)
	{
		var ions = new List<MatchedFragmentIon>(count);
		for (int i = 1; i <= count; i++)
		{
			var product = new Product(
				ProductType.y,
				FragmentationTerminus.C,
				neutralMass: 0.0,
				fragmentNumber: i,
				residuePosition: i,
				neutralLoss: 0.0);
			ions.Add(new MatchedFragmentIon(product, experMz: 100.0 + i * 10, experIntensity: 1.0f, charge: 1));
		}
		return ions;
	}

	/// <summary>
	/// Creates a minimal <see cref="MslLibraryEntry"/> with <paramref name="fragCount"/>
	/// fragment ions and the supplied metadata.
	/// </summary>
	private static MslLibraryEntry MakeEntry(string seq, int charge, double mz, double irt,
		int fragCount, double ionMobility = 0.0, string protein = "", bool isDecoy = false)
	{
		var frags = new List<MslFragmentIon>(fragCount);
		for (int i = 1; i <= fragCount; i++)
		{
			frags.Add(new MslFragmentIon
			{
				Mz = 100.0f + i * 10,
				Intensity = 1.0f,
				ProductType = ProductType.y,
				FragmentNumber = i,
				ResiduePosition = i,
				Charge = 1
			});
		}

		return new MslLibraryEntry
		{
			ModifiedSequence = seq,
			StrippedSequence = seq,
			PrecursorMz = mz,
			Charge = charge,
			Irt = irt,
			IonMobility = ionMobility,
			ProteinAccession = protein,
			IsDecoy = isDecoy,
			IsProteotypic = !isDecoy,
			QValue = float.NaN,
			Source = MslFormat.SourceType.Predicted,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Fragments = frags
		};
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group 1 — Return value correctness
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// The returned <see cref="MslUpdateResult.OutputPath"/> must equal the path passed in.
	/// </summary>
	[Test]
	public void UpdateAndSave_OutputPath_MatchesArgument()
	{
		string output = TempPath(nameof(UpdateAndSave_OutputPath_MatchesArgument));
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		MslUpdateResult result = lib.UpdateAndSave(Array.Empty<LibrarySpectrum>(), output);

		Assert.That(result.OutputPath, Is.EqualTo(output));
	}

	/// <summary>
	/// When all five incoming spectra are anchors and the regression can be computed,
	/// <see cref="MslUpdateResult.RtNormalisationApplied"/> must be <see langword="true"/>.
	/// </summary>
	[Test]
	public void UpdateAndSave_RtNormalisationApplied_TrueWhenEnoughAnchors()
	{
		string output = TempPath(nameof(UpdateAndSave_RtNormalisationApplied_TrueWhenEnoughAnchors));
		var incoming = BuildFiveAnchorSpectra(fragCountOverride: 2); // same count → no replacement

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = lib.UpdateAndSave(incoming, output, minRegressionAnchors: 5);

		Assert.That(result.RtNormalisationApplied, Is.True,
			"Normalisation must be applied when anchor count >= minRegressionAnchors.");
	}

	/// <summary>
	/// The <see cref="MslUpdateResult.AnchorCount"/> must equal the number of incoming
	/// spectra whose sequence/charge key is present in the original library.
	/// </summary>
	[Test]
	public void UpdateAndSave_AnchorCount_EqualsOverlapCount()
	{
		string output = TempPath(nameof(UpdateAndSave_AnchorCount_EqualsOverlapCount));

		// Three of the five incoming spectra overlap with the original library.
		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 5.0,  fragCount: 2),
			MakeSpectrum("BBBB", 2, 250.0, observedRt: 10.0, fragCount: 2),
			MakeSpectrum("CCCC", 2, 300.0, observedRt: 15.0, fragCount: 2),
			MakeSpectrum("XXXX", 2, 999.0, observedRt: 50.0, fragCount: 2), // novel
        };

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = lib.UpdateAndSave(incoming, output, minRegressionAnchors: 2);

		Assert.That(result.AnchorCount, Is.EqualTo(3),
			"AnchorCount must equal the number of sequence/charge keys shared between the " +
			"incoming list and the original library.");
	}

	/// <summary>
	/// <see cref="MslUpdateResult.TotalCount"/> must equal
	/// Retained + Replaced + Added.
	/// </summary>
	[Test]
	public void UpdateAndSave_TotalCount_EqualsRetainedPlusReplacedPlusAdded()
	{
		string output = TempPath(nameof(UpdateAndSave_TotalCount_EqualsRetainedPlusReplacedPlusAdded));

		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 5.0,  fragCount: 5), // replaces (5 > 2)
            MakeSpectrum("ZZZZ", 2, 999.0, observedRt: 99.0, fragCount: 2), // novel
        };

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = lib.UpdateAndSave(incoming, output, minRegressionAnchors: 2);

		Assert.That(result.TotalCount,
			Is.EqualTo(result.RetainedCount + result.ReplacedCount + result.AddedCount),
			"TotalCount must equal the sum of the three component counts.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group 2 — Retained entries
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// An entry in the original library with no matching incoming spectrum must appear
	/// in the output with its original iRT unchanged.
	/// </summary>
	[Test]
	public void UpdateAndSave_RetainedEntry_IrtUnchanged()
	{
		string output = TempPath(nameof(UpdateAndSave_RetainedEntry_IrtUnchanged));

		// Only AAAA is in the incoming list; BBBB–EEEE have no match.
		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 5.0, fragCount: 2),
		};

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		lib.UpdateAndSave(incoming, output, minRegressionAnchors: 10); // anchors too few → no normalisation

		using MslLibrary updated = MslLibrary.Load(output);
		updated.TryGetEntry("BBBB", 2, out MslLibraryEntry? entry);

		Assert.That(entry, Is.Not.Null, "BBBB/2 must be present in the updated library.");
		Assert.That(entry!.Irt, Is.EqualTo(IrtB).Within(0.01),
			"Retained entry iRT must be unchanged.");
	}

	/// <summary>
	/// <see cref="MslUpdateResult.RetainedCount"/> must count every original entry that
	/// was not replaced.
	/// </summary>
	[Test]
	public void UpdateAndSave_RetainedCount_CorrectWhenNoIncoming()
	{
		string output = TempPath(nameof(UpdateAndSave_RetainedCount_CorrectWhenNoIncoming));
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);

		MslUpdateResult result = lib.UpdateAndSave(Array.Empty<LibrarySpectrum>(), output);

		Assert.That(result.RetainedCount, Is.EqualTo(5),
			"All five original entries must be retained when incoming is empty.");
		Assert.That(result.ReplacedCount, Is.EqualTo(0));
		Assert.That(result.AddedCount, Is.EqualTo(0));
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group 3 — Replaced entries
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// When an incoming spectrum has more fragment ions than the original entry,
	/// <see cref="MslUpdateResult.ReplacedCount"/> must be incremented and the
	/// output entry must carry the new fragment ions.
	/// </summary>
	[Test]
	public void UpdateAndSave_IncomingHasMoreFragments_EntryReplaced()
	{
		string output = TempPath(nameof(UpdateAndSave_IncomingHasMoreFragments_EntryReplaced));

		// AAAA original has 2 fragments; incoming has 5.
		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 5.0, fragCount: 5),
		};

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = lib.UpdateAndSave(incoming, output, minRegressionAnchors: 10);

		Assert.That(result.ReplacedCount, Is.EqualTo(1), "One entry must be reported replaced.");

		using MslLibrary updated = MslLibrary.Load(output);
		updated.TryGetEntry("AAAA", 2, out MslLibraryEntry? entry);

		Assert.That(entry!.Fragments.Count, Is.EqualTo(5),
			"Replaced entry must carry the incoming fragment count.");
	}

	/// <summary>
	/// When the incoming spectrum has fewer fragments than the original, the original
	/// is kept and <see cref="MslUpdateResult.ReplacedCount"/> stays at zero.
	/// </summary>
	[Test]
	public void UpdateAndSave_IncomingHasFewerFragments_OriginalKept()
	{
		string output = TempPath(nameof(UpdateAndSave_IncomingHasFewerFragments_OriginalKept));

		// AAAA original has 2 fragments; incoming has only 1.
		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 5.0, fragCount: 1),
		};

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = lib.UpdateAndSave(incoming, output, minRegressionAnchors: 10);

		Assert.That(result.ReplacedCount, Is.EqualTo(0));
		Assert.That(result.RetainedCount, Is.EqualTo(5));
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group 4 — Ties go to the original
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// When the incoming spectrum has the same fragment count as the original entry,
	/// the original is kept (no replacement).
	/// </summary>
	[Test]
	public void UpdateAndSave_EqualFragmentCount_OriginalKept()
	{
		string output = TempPath(nameof(UpdateAndSave_EqualFragmentCount_OriginalKept));

		// Both original and incoming have 2 fragments — tie → keep original.
		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 5.0, fragCount: 2),
		};

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = lib.UpdateAndSave(incoming, output, minRegressionAnchors: 10);

		Assert.That(result.ReplacedCount, Is.EqualTo(0),
			"Tie in fragment count must not cause a replacement.");
		Assert.That(result.RetainedCount, Is.EqualTo(5));
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group 5 — Added entries
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// An incoming spectrum whose sequence/charge key is not present in the original
	/// library must be appended to the output and counted in
	/// <see cref="MslUpdateResult.AddedCount"/>.
	/// </summary>
	[Test]
	public void UpdateAndSave_NovelSpectrum_Appended()
	{
		string output = TempPath(nameof(UpdateAndSave_NovelSpectrum_Appended));

		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("NOVEL", 2, 555.0, observedRt: 30.0, fragCount: 3),
		};

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = lib.UpdateAndSave(incoming, output, minRegressionAnchors: 10);

		Assert.That(result.AddedCount, Is.EqualTo(1), "Novel entry must increment AddedCount.");
		Assert.That(result.RetainedCount, Is.EqualTo(5), "All original entries must be retained.");

		using MslLibrary updated = MslLibrary.Load(output);
		Assert.That(updated.PrecursorCount, Is.EqualTo(6),
			"Updated library must have one more entry than the original.");
		Assert.That(updated.TryGetEntry("NOVEL", 2, out _), Is.True,
			"Novel entry must be findable by sequence/charge in the updated library.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group 6 — RT normalisation applied
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// When five anchor spectra exist and the regression iRT = 2×observedRT + 1 is exact,
	/// a replaced entry's iRT in the output must equal slope×newRT + intercept.
	/// </summary>
	[Test]
	public void UpdateAndSave_RtNormalisation_ReplacedEntryIrtTransformed()
	{
		string output = TempPath(nameof(UpdateAndSave_RtNormalisation_ReplacedEntryIrtTransformed));

		// Five anchors with exact linear relationship: iRT = 2×observedRT + 1
		// AAAA gets replaced (incoming has 5 fragments > original 2).
		var incoming = BuildFiveAnchorSpectra(fragCountOverride: 5);

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = lib.UpdateAndSave(incoming, output, minRegressionAnchors: 5);

		Assert.That(result.RtNormalisationApplied, Is.True);
		Assert.That(result.RtSlope, Is.EqualTo(RtSlope).Within(1e-6));
		Assert.That(result.RtIntercept, Is.EqualTo(RtIntercept).Within(1e-4));

		// AAAA observed RT was 5.0 → expected iRT = 2.0×5.0 + 1.0 = 11.0
		using MslLibrary updated = MslLibrary.Load(output);
		updated.TryGetEntry("AAAA", 2, out MslLibraryEntry? aaaa);
		double expectedIrt = RtSlope * 5.0 + RtIntercept;
		Assert.That(aaaa!.Irt, Is.EqualTo(expectedIrt).Within(0.01),
			"Replaced entry iRT must equal slope × observedRT + intercept.");
	}

	/// <summary>
	/// When RT normalisation is applied, a novel added entry's iRT must also be transformed.
	/// </summary>
	[Test]
	public void UpdateAndSave_RtNormalisation_AddedEntryIrtTransformed()
	{
		string output = TempPath(nameof(UpdateAndSave_RtNormalisation_AddedEntryIrtTransformed));

		double novelObservedRt = 30.0;
		double expectedIrt = RtSlope * novelObservedRt + RtIntercept; // = 61.0

		// Five anchors to satisfy regression threshold + one novel entry.
		var incoming = BuildFiveAnchorSpectra(fragCountOverride: 2); // same count → retain originals
		incoming.Add(MakeSpectrum("NOVEL", 2, 555.0, novelObservedRt, fragCount: 2));

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = lib.UpdateAndSave(incoming, output, minRegressionAnchors: 5);

		Assert.That(result.RtNormalisationApplied, Is.True);

		using MslLibrary updated = MslLibrary.Load(output);
		updated.TryGetEntry("NOVEL", 2, out MslLibraryEntry? novel);
		Assert.That(novel!.Irt, Is.EqualTo(expectedIrt).Within(0.01),
			"Added entry iRT must be normalised via the regression.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group 7 — RT normalisation skipped
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// When the anchor count is below <c>minRegressionAnchors</c>,
	/// <see cref="MslUpdateResult.RtNormalisationApplied"/> must be <see langword="false"/>
	/// and regression parameters must be identity (slope=1, intercept=0).
	/// </summary>
	[Test]
	public void UpdateAndSave_TooFewAnchors_NormalisationSkipped()
	{
		string output = TempPath(nameof(UpdateAndSave_TooFewAnchors_NormalisationSkipped));

		// Only one overlap — below the default threshold of 5.
		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 5.0, fragCount: 5),
		};

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = lib.UpdateAndSave(incoming, output, minRegressionAnchors: 5);

		Assert.That(result.RtNormalisationApplied, Is.False);
		Assert.That(result.RtSlope, Is.EqualTo(1.0).Within(1e-10));
		Assert.That(result.RtIntercept, Is.EqualTo(0.0).Within(1e-10));
	}

	/// <summary>
	/// When RT normalisation is skipped, the added entry's iRT must equal the raw observed RT.
	/// </summary>
	[Test]
	public void UpdateAndSave_NormalisationSkipped_AddedEntryIrtIsRawRt()
	{
		string output = TempPath(nameof(UpdateAndSave_NormalisationSkipped_AddedEntryIrtIsRawRt));

		double rawRt = 77.0;
		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("NOVEL", 2, 555.0, observedRt: rawRt, fragCount: 2),
		};

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		lib.UpdateAndSave(incoming, output, minRegressionAnchors: 5);

		using MslLibrary updated = MslLibrary.Load(output);
		updated.TryGetEntry("NOVEL", 2, out MslLibraryEntry? novel);
		Assert.That(novel!.Irt, Is.EqualTo(rawRt).Within(0.01),
			"When normalisation is skipped, raw observed RT must be stored as-is.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group 8 — RT identity fallback (degenerate regression)
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// When all anchor observed-RT values are identical (zero variance in x), the OLS
	/// denominator is zero. The method must return the identity transform (slope=1,
	/// intercept=0) rather than producing NaN or throwing.
	/// </summary>
	[Test]
	public void UpdateAndSave_AllAnchorsHaveSameRt_IdentityFallback()
	{
		string output = TempPath(nameof(UpdateAndSave_AllAnchorsHaveSameRt_IdentityFallback));

		// All five anchors share the same observed RT = 10.0 → degenerate regression.
		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 10.0, fragCount: 2),
			MakeSpectrum("BBBB", 2, 250.0, observedRt: 10.0, fragCount: 2),
			MakeSpectrum("CCCC", 2, 300.0, observedRt: 10.0, fragCount: 2),
			MakeSpectrum("DDDD", 2, 350.0, observedRt: 10.0, fragCount: 2),
			MakeSpectrum("EEEE", 2, 400.0, observedRt: 10.0, fragCount: 2),
		};

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = null!;
		Assert.That(() =>
		{
			result = lib.UpdateAndSave(incoming, output, minRegressionAnchors: 5);
		}, Throws.Nothing, "Degenerate regression must not throw.");

		// The OLS guard returns identity.
		Assert.That(result.RtSlope, Is.EqualTo(1.0).Within(1e-10));
		Assert.That(result.RtIntercept, Is.EqualTo(0.0).Within(1e-10));
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group 9 — Duplicate incoming spectra
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// When the same sequence/charge appears multiple times in the incoming list, only the
	/// one with the most fragment ions is used (mirrors MetaMorpheus deduplication).
	/// </summary>
	[Test]
	public void UpdateAndSave_DuplicateIncoming_BestFragmentCountWins()
	{
		string output = TempPath(nameof(UpdateAndSave_DuplicateIncoming_BestFragmentCountWins));

		// AAAA appears twice; the second copy has more fragments and must win.
		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 5.0, fragCount: 3),
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 5.0, fragCount: 7), // best
        };

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		lib.UpdateAndSave(incoming, output, minRegressionAnchors: 10);

		using MslLibrary updated = MslLibrary.Load(output);
		updated.TryGetEntry("AAAA", 2, out MslLibraryEntry? entry);
		Assert.That(entry!.Fragments.Count, Is.EqualTo(7),
			"The duplicate with the most fragment ions must be selected.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group 10 — Empty incoming list
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// When the incoming list is empty, the output must be a verbatim copy of the
	/// original library (same precursor count and same iRT values).
	/// </summary>
	[Test]
	public void UpdateAndSave_EmptyIncoming_OutputIsVerbatimCopy()
	{
		string output = TempPath(nameof(UpdateAndSave_EmptyIncoming_OutputIsVerbatimCopy));

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = lib.UpdateAndSave(Array.Empty<LibrarySpectrum>(), output);

		Assert.That(result.RetainedCount, Is.EqualTo(5));
		Assert.That(result.ReplacedCount, Is.EqualTo(0));
		Assert.That(result.AddedCount, Is.EqualTo(0));

		using MslLibrary updated = MslLibrary.Load(output);
		Assert.That(updated.PrecursorCount, Is.EqualTo(5));

		// Spot-check iRTs are preserved.
		updated.TryGetEntry("AAAA", 2, out MslLibraryEntry? a);
		Assert.That(a!.Irt, Is.EqualTo(IrtA).Within(0.01));
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group 11 — Metadata preservation
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// When an entry is replaced, the original IonMobility must be copied onto the
	/// replacement because <see cref="LibrarySpectrum"/> has no ion-mobility field.
	/// </summary>
	[Test]
	public void UpdateAndSave_ReplacedEntry_IonMobilityPreserved()
	{
		string output = TempPath(nameof(UpdateAndSave_ReplacedEntry_IonMobilityPreserved));

		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 5.0, fragCount: 5), // replaces AAAA
        };

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		lib.UpdateAndSave(incoming, output, minRegressionAnchors: 10);

		using MslLibrary updated = MslLibrary.Load(output);
		updated.TryGetEntry("AAAA", 2, out MslLibraryEntry? entry);
		Assert.That(entry!.IonMobility, Is.EqualTo(0.85).Within(1e-4),
			"IonMobility must be copied from the original entry to the replacement.");
	}

	/// <summary>
	/// When an entry is replaced, the original ProteinAccession must be preserved.
	/// </summary>
	[Test]
	public void UpdateAndSave_ReplacedEntry_ProteinAccessionPreserved()
	{
		string output = TempPath(nameof(UpdateAndSave_ReplacedEntry_ProteinAccessionPreserved));

		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 5.0, fragCount: 5),
		};

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		lib.UpdateAndSave(incoming, output, minRegressionAnchors: 10);

		using MslLibrary updated = MslLibrary.Load(output);
		updated.TryGetEntry("AAAA", 2, out MslLibraryEntry? entry);
		Assert.That(entry!.ProteinAccession, Is.EqualTo("P00001"),
			"ProteinAccession must be copied from the original entry to the replacement.");
	}

	/// <summary>
	/// When an entry is replaced, the original IsDecoy flag must be preserved.
	/// </summary>
	[Test]
	public void UpdateAndSave_ReplacedEntry_IsDecoyPreserved()
	{
		string output = TempPath(nameof(UpdateAndSave_ReplacedEntry_IsDecoyPreserved));

		// EEEE is the only decoy (IsDecoy = true). Replace it with a higher-count spectrum.
		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("EEEE", 2, 400.0, observedRt: 25.0, fragCount: 5),
		};

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		lib.UpdateAndSave(incoming, output, minRegressionAnchors: 10);

		using MslLibrary updated = MslLibrary.Load(output);
		updated.TryGetEntry("EEEE", 2, out MslLibraryEntry? entry);
		Assert.That(entry!.IsDecoy, Is.True,
			"IsDecoy must be copied from the original entry to the replacement.");
	}

	/// <summary>
	/// When an entry is replaced, its <see cref="MslFormat.SourceType"/> must be set to
	/// <see cref="MslFormat.SourceType.Empirical"/>.
	/// </summary>
	[Test]
	public void UpdateAndSave_ReplacedEntry_SourceTypeIsEmpirical()
	{
		string output = TempPath(nameof(UpdateAndSave_ReplacedEntry_SourceTypeIsEmpirical));

		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 5.0, fragCount: 5),
		};

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		lib.UpdateAndSave(incoming, output, minRegressionAnchors: 10);

		using MslLibrary updated = MslLibrary.Load(output);
		updated.TryGetEntry("AAAA", 2, out MslLibraryEntry? entry);
		Assert.That(entry!.Source, Is.EqualTo(MslFormat.SourceType.Empirical),
			"Replaced entries must be tagged as Empirical.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group 12 — Output file is readable
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// The file written by <see cref="MslLibrary.UpdateAndSave"/> must be a valid .msl
	/// file that can be opened by <see cref="MslLibrary.Load"/> without error.
	/// </summary>
	[Test]
	public void UpdateAndSave_OutputFile_IsReadableByLoad()
	{
		string output = TempPath(nameof(UpdateAndSave_OutputFile_IsReadableByLoad));
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		lib.UpdateAndSave(Array.Empty<LibrarySpectrum>(), output);

		Assert.That(() =>
		{
			using MslLibrary reloaded = MslLibrary.Load(output);
			_ = reloaded.PrecursorCount;
		}, Throws.Nothing, "Output file must be openable by MslLibrary.Load.");
	}

	/// <summary>
	/// The output file must be a valid .msl file openable by
	/// <see cref="MslLibrary.LoadIndexOnly"/> without error.
	/// </summary>
	[Test]
	public void UpdateAndSave_OutputFile_IsReadableByLoadIndexOnly()
	{
		string output = TempPath(nameof(UpdateAndSave_OutputFile_IsReadableByLoadIndexOnly));
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		lib.UpdateAndSave(Array.Empty<LibrarySpectrum>(), output);

		Assert.That(() =>
		{
			using MslLibrary reloaded = MslLibrary.LoadIndexOnly(output);
			_ = reloaded.PrecursorCount;
		}, Throws.Nothing, "Output file must be openable by MslLibrary.LoadIndexOnly.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group 13 — MslUpdateResult helpers
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// <see cref="MslUpdateResult.ToString"/> must contain the output path, the three
	/// count fields, and the normalisation summary.
	/// </summary>
	[Test]
	public void MslUpdateResult_ToString_ContainsKeyFields()
	{
		string output = TempPath(nameof(MslUpdateResult_ToString_ContainsKeyFields));
		var incoming = BuildFiveAnchorSpectra(fragCountOverride: 2);

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = lib.UpdateAndSave(incoming, output, minRegressionAnchors: 5);

		string s = result.ToString();
		Assert.That(s, Does.Contain(output), "ToString must contain the output path.");
		Assert.That(s, Does.Contain("retained="), "ToString must contain retained count.");
		Assert.That(s, Does.Contain("replaced="), "ToString must contain replaced count.");
		Assert.That(s, Does.Contain("added="), "ToString must contain added count.");
		Assert.That(s, Does.Contain("normalised"), "ToString must mention RT normalisation.");
	}

	/// <summary>
	/// <see cref="MslUpdateResult.TotalCount"/> must equal the precursor count of the
	/// reloaded output library.
	/// </summary>
	[Test]
	public void MslUpdateResult_TotalCount_MatchesOutputLibraryPrecursorCount()
	{
		string output = TempPath(nameof(MslUpdateResult_TotalCount_MatchesOutputLibraryPrecursorCount));

		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("AAAA",  2, 200.0, observedRt: 5.0,  fragCount: 5), // replace
            MakeSpectrum("NOVEL", 2, 555.0, observedRt: 99.0, fragCount: 2), // add
        };

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = lib.UpdateAndSave(incoming, output, minRegressionAnchors: 10);

		using MslLibrary updated = MslLibrary.Load(output);
		Assert.That(result.TotalCount, Is.EqualTo(updated.PrecursorCount),
			"TotalCount must equal the precursor count of the written file.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group 14 — Guard clauses
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Passing <see langword="null"/> for <c>newSpectra</c> must throw
	/// <see cref="ArgumentNullException"/>.
	/// </summary>
	[Test]
	public void UpdateAndSave_NullNewSpectra_ThrowsArgumentNullException()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		Assert.That(
			() => lib.UpdateAndSave(null!, TempPath("null_spectra")),
			Throws.ArgumentNullException,
			"null newSpectra must throw ArgumentNullException.");
	}

	/// <summary>
	/// Passing <see langword="null"/> for <c>outputPath</c> must throw
	/// <see cref="ArgumentNullException"/>.
	/// </summary>
	[Test]
	public void UpdateAndSave_NullOutputPath_ThrowsArgumentNullException()
	{
		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		Assert.That(
			() => lib.UpdateAndSave(Array.Empty<LibrarySpectrum>(), null!),
			Throws.ArgumentNullException,
			"null outputPath must throw ArgumentNullException.");
	}

	/// <summary>
	/// Calling <see cref="MslLibrary.UpdateAndSave"/> on a disposed library must throw
	/// <see cref="ObjectDisposedException"/>.
	/// </summary>
	[Test]
	public void UpdateAndSave_DisposedLibrary_ThrowsObjectDisposedException()
	{
		MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		lib.Dispose();

		Assert.That(
			() => lib.UpdateAndSave(Array.Empty<LibrarySpectrum>(), TempPath("disposed")),
			Throws.InstanceOf<ObjectDisposedException>(),
			"UpdateAndSave on a disposed library must throw ObjectDisposedException.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Group 15 — minRegressionAnchors clamping
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Passing <c>minRegressionAnchors = 0</c> must be silently clamped to 2; the method
	/// must not throw and normalisation must be applied when at least 2 anchors are present.
	/// </summary>
	[Test]
	public void UpdateAndSave_MinAnchorsZero_ClampedToTwo()
	{
		string output = TempPath(nameof(UpdateAndSave_MinAnchorsZero_ClampedToTwo));

		// Two overlapping spectra — enough to satisfy the clamped threshold of 2.
		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 5.0,  fragCount: 2),
			MakeSpectrum("BBBB", 2, 250.0, observedRt: 10.0, fragCount: 2),
		};

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = null!;
		Assert.That(
			() => result = lib.UpdateAndSave(incoming, output, minRegressionAnchors: 0),
			Throws.Nothing,
			"minRegressionAnchors = 0 must be accepted without throwing.");

		Assert.That(result.RtNormalisationApplied, Is.True,
			"With minRegressionAnchors clamped to 2 and 2 anchors present, " +
			"normalisation must be applied.");
	}

	/// <summary>
	/// Passing <c>minRegressionAnchors = 1</c> must also be clamped to 2.
	/// </summary>
	[Test]
	public void UpdateAndSave_MinAnchorsOne_ClampedToTwo()
	{
		string output = TempPath(nameof(UpdateAndSave_MinAnchorsOne_ClampedToTwo));

		var incoming = new List<LibrarySpectrum>
		{
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 5.0,  fragCount: 2),
			MakeSpectrum("BBBB", 2, 250.0, observedRt: 10.0, fragCount: 2),
		};

		using MslLibrary lib = MslLibrary.Load(SharedLibraryPath);
		MslUpdateResult result = lib.UpdateAndSave(incoming, output, minRegressionAnchors: 1);

		Assert.That(result.RtNormalisationApplied, Is.True,
			"With minRegressionAnchors clamped to 2 and 2 anchors present, " +
			"normalisation must be applied.");
	}

	// ────────────────────────────────────────────────────────────────────────
	// Private helpers
	// ────────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Returns five <see cref="LibrarySpectrum"/> objects that each match one original entry
	/// and satisfy iRT = 2×observedRT + 1 exactly. The <paramref name="fragCountOverride"/>
	/// controls how many fragment ions each spectrum carries.
	/// </summary>
	private static List<LibrarySpectrum> BuildFiveAnchorSpectra(int fragCountOverride)
		=> new()
		{
			MakeSpectrum("AAAA", 2, 200.0, observedRt: 5.0,  fragCount: fragCountOverride),
			MakeSpectrum("BBBB", 2, 250.0, observedRt: 10.0, fragCount: fragCountOverride),
			MakeSpectrum("CCCC", 2, 300.0, observedRt: 15.0, fragCount: fragCountOverride),
			MakeSpectrum("DDDD", 2, 350.0, observedRt: 20.0, fragCount: fragCountOverride),
			MakeSpectrum("EEEE", 2, 400.0, observedRt: 25.0, fragCount: fragCountOverride),
		};
}