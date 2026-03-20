using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MspReader = Readers.SpectralLibrary.SpectralLibrary;

namespace Test.MslSpectralLibrary;

/// <summary>
/// NUnit 4 tests for <see cref="MslConverter"/>.
///
/// <para>
/// All tests that require a real MSP file use <c>syntheticTest2_PredictedLibrary.msp</c>
/// (42 precursors, 968 fragments — the same baseline as <c>TestMslRegressionMsp</c>).
/// The file is located by <see cref="FindMspFile"/>, which walks up the directory tree
/// from the NUnit output directory. If the file is absent the fixture is skipped with a
/// clear message rather than failing.
/// </para>
///
/// Coverage groups:
/// <list type="number">
///   <item>Guard conditions — null/empty/missing path arguments</item>
///   <item><see cref="MslConverter.FromMspFile"/> — happy path against real MSP</item>
///   <item><see cref="MslConverter.FromMspLibrarySpectra"/> — happy path with synthetic spectra</item>
///   <item>Query API — verify the returned <see cref="MslLibrary"/> behaves like a file-loaded one</item>
///   <item>Empty input — zero spectra produces a valid empty library</item>
/// </list>
/// </summary>
[TestFixture]
public sealed class TestMslConverter
{
	// ── Baseline constants (must match syntheticTest2_PredictedLibrary.msp) ──

	/// <summary>Number of "Name:" entries in <c>syntheticTest2_PredictedLibrary.msp</c>.</summary>
	private const int ExpectedPrecursorCount = 42;

	/// <summary>Sum of all "Num peaks:" values in the MSP file.</summary>
	private const int ExpectedTotalFragments = 968;

	/// <summary>Float tolerance for precursor and fragment m/z comparisons (float32 rounding).</summary>
	private const double MzTolerance = 1e-4;

	// ── File paths ────────────────────────────────────────────────────────────

	private const string MspFileName = "syntheticTest2_PredictedLibrary.msp";
	private static readonly string MspFilePath = FindMspFile();

	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "MslConverterTests");

	// ── Shared state populated in OneTimeSetUp ────────────────────────────────

	/// <summary>
	/// Library created by <see cref="MslConverter.FromMspFile"/> in
	/// <see cref="OneTimeSetUp"/>. Null when the MSP file is absent.
	/// </summary>
	private static MslLibrary? _fromFile;

	/// <summary>
	/// Spectra read directly from the MSP file; used as ground truth in comparison tests.
	/// </summary>
	private static List<LibrarySpectrum> _mspSpectra = new();

	// ── One-time setup / teardown ─────────────────────────────────────────────

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		if (!File.Exists(MspFilePath))
		{
			Assert.Ignore(
				$"MSP test file not found at '{MspFilePath}'. " +
				$"Ensure {MspFileName} is present in TestData.");
			return;
		}

		Directory.CreateDirectory(OutputDirectory);

		// Read ground-truth spectra directly from MSP for later comparison
		var mspLib = new MspReader(new List<string> { MspFilePath });
		_mspSpectra = mspLib.GetAllLibrarySpectra().ToList();
		mspLib.CloseConnections();

		// Build the shared library under test via MslConverter
		_fromFile = MslConverter.FromMspFile(MspFilePath);
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		_fromFile?.Dispose();

		if (Directory.Exists(OutputDirectory))
			Directory.Delete(OutputDirectory, recursive: true);
	}

	// ── Group 1: Guard conditions ─────────────────────────────────────────────

	/// <summary><see cref="MslConverter.FromMspFile"/> must throw on null path.</summary>
	[Test]
	public void FromMspFile_NullPath_ThrowsArgumentNullException()
	{
		Assert.That(() => MslConverter.FromMspFile(null!),
			Throws.InstanceOf<ArgumentNullException>());
	}

	/// <summary><see cref="MslConverter.FromMspFile"/> must throw on empty string.</summary>
	[Test]
	public void FromMspFile_EmptyPath_ThrowsArgumentNullException()
	{
		Assert.That(() => MslConverter.FromMspFile(string.Empty),
			Throws.InstanceOf<ArgumentNullException>());
	}

	/// <summary><see cref="MslConverter.FromMspFile"/> must throw when the file does not exist.</summary>
	[Test]
	public void FromMspFile_MissingFile_ThrowsFileNotFoundException()
	{
		Assert.That(() => MslConverter.FromMspFile(@"C:\does\not\exist.msp"),
			Throws.InstanceOf<FileNotFoundException>());
	}

	/// <summary><see cref="MslConverter.FromMspLibrarySpectra"/> must throw on null.</summary>
	[Test]
	public void FromMspLibrarySpectra_NullList_ThrowsArgumentNullException()
	{
		Assert.That(() => MslConverter.FromMspLibrarySpectra(null!),
			Throws.InstanceOf<ArgumentNullException>());
	}

	// ── Group 2: FromMspFile happy path ───────────────────────────────────────

	/// <summary>
	/// <see cref="MslConverter.FromMspFile"/> must return a non-null <see cref="MslLibrary"/>.
	/// </summary>
	[Test]
	public void FromMspFile_ReturnsNonNullLibrary()
	{
		Assert.That(_fromFile, Is.Not.Null);
	}

	/// <summary>
	/// The returned library must be in full-load mode (<see cref="MslLibrary.IsIndexOnly"/>
	/// must be <see langword="false"/>).
	/// </summary>
	[Test]
	public void FromMspFile_IsIndexOnly_IsFalse()
	{
		Assert.That(_fromFile!.IsIndexOnly, Is.False);
	}

	/// <summary>
	/// Precursor count must match the number of "Name:" entries in the MSP file.
	/// </summary>
	[Test]
	public void FromMspFile_PrecursorCount_MatchesMsp()
	{
		Assert.That(_fromFile!.PrecursorCount, Is.EqualTo(ExpectedPrecursorCount));
	}

	/// <summary>
	/// Total fragment count across all entries must match the MSP baseline.
	/// </summary>
	[Test]
	public void FromMspFile_TotalFragmentCount_MatchesMsp()
	{
		int totalFragments = _fromFile!.GetAllEntries()
			.Sum(e => e.Fragments.Count);

		Assert.That(totalFragments, Is.EqualTo(ExpectedTotalFragments));
	}

	/// <summary>
	/// For every spectrum in the MSP, <see cref="MslLibrary.TryGetLibrarySpectrum"/> must
	/// find the entry and its precursor m/z must match within float tolerance.
	/// </summary>
	[Test]
	public void FromMspFile_AllPrecursorMz_MatchMsp()
	{
		foreach (LibrarySpectrum mspSpectrum in _mspSpectra)
		{
			(string seq, int charge) = ParseName(mspSpectrum.Name);

			bool found = _fromFile!.TryGetLibrarySpectrum(seq, charge, out LibrarySpectrum? lib);

			Assert.That(found, Is.True,
				$"Precursor '{mspSpectrum.Name}' not found in converted library.");
			Assert.That(lib!.PrecursorMz,
				Is.EqualTo(mspSpectrum.PrecursorMz).Within(MzTolerance),
				$"PrecursorMz mismatch for '{mspSpectrum.Name}'.");
		}
	}

	/// <summary>
	/// For every spectrum in the MSP, fragment m/z values (sorted ascending) must match
	/// those in the converted library within float tolerance.
	/// </summary>
	[Test]
	public void FromMspFile_AllFragmentMz_MatchMsp()
	{
		foreach (LibrarySpectrum mspSpectrum in _mspSpectra)
		{
			(string seq, int charge) = ParseName(mspSpectrum.Name);
			_fromFile!.TryGetLibrarySpectrum(seq, charge, out LibrarySpectrum? lib);
			if (lib is null) continue; // caught by AllPrecursorMz test

			List<double> mspMzs = mspSpectrum.MatchedFragmentIons
				.Select(f => f.Mz).OrderBy(x => x).ToList();
			List<double> libMzs = lib.MatchedFragmentIons
				.Select(f => f.Mz).OrderBy(x => x).ToList();

			Assert.That(libMzs.Count, Is.EqualTo(mspMzs.Count),
				$"Fragment count mismatch for '{mspSpectrum.Name}'.");

			for (int i = 0; i < mspMzs.Count; i++)
				Assert.That(libMzs[i], Is.EqualTo(mspMzs[i]).Within(MzTolerance),
					$"Fragment [{i}] m/z mismatch for '{mspSpectrum.Name}'.");
		}
	}

	// ── Group 3: FromMspLibrarySpectra happy path ─────────────────────────────

	/// <summary>
	/// <see cref="MslConverter.FromMspLibrarySpectra"/> with an empty list must return a
	/// valid <see cref="MslLibrary"/> with zero precursors rather than throwing.
	/// </summary>
	[Test]
	public void FromMspLibrarySpectra_EmptyList_ReturnsEmptyLibrary()
	{
		MslLibrary lib = MslConverter.FromMspLibrarySpectra(new List<LibrarySpectrum>());

		Assert.That(lib, Is.Not.Null);
		Assert.That(lib.PrecursorCount, Is.EqualTo(0));
	}

	/// <summary>
	/// A single synthetic <see cref="LibrarySpectrum"/> must round-trip through
	/// <see cref="MslConverter.FromMspLibrarySpectra"/> and be retrievable by
	/// <see cref="MslLibrary.TryGetLibrarySpectrum"/>.
	/// </summary>
	[Test]
	public void FromMspLibrarySpectra_SingleSpectrum_IsRetrievable()
	{
		var ion = MakePeptideIon(ProductType.b, 2, mz: 200.0, intensity: 1.0f);
		var spectrum = new LibrarySpectrum("PEPTIDE", 400.0, 2,
			new List<MatchedFragmentIon> { ion }, 10.0);

		MslLibrary lib = MslConverter.FromMspLibrarySpectra(
			new List<LibrarySpectrum> { spectrum });

		bool found = lib.TryGetLibrarySpectrum("PEPTIDE", 2, out LibrarySpectrum? retrieved);

		Assert.That(found, Is.True);
		Assert.That(retrieved!.PrecursorMz, Is.EqualTo(400.0).Within(MzTolerance));
		Assert.That(retrieved.MatchedFragmentIons, Has.Count.EqualTo(1));
		Assert.That(retrieved.MatchedFragmentIons[0].Mz,
			Is.EqualTo(200.0).Within(MzTolerance));
	}

	/// <summary>
	/// Multiple spectra passed to <see cref="MslConverter.FromMspLibrarySpectra"/> must all
	/// be present and queryable in the returned library.
	/// </summary>
	[Test]
	public void FromMspLibrarySpectra_MultipleSpectra_AllRetrievable()
	{
		var spectra = new List<LibrarySpectrum>
		{
			new LibrarySpectrum("PEPTIDE",    400.0, 2,
				new List<MatchedFragmentIon> { MakePeptideIon(ProductType.b, 2, 200.0, 1.0f) }, 10.0),
			new LibrarySpectrum("ACDEFGHIK",  530.0, 2,
				new List<MatchedFragmentIon> { MakePeptideIon(ProductType.y, 3, 364.2, 1.0f) }, 25.0),
			new LibrarySpectrum("SAMPLER",    450.0, 3,
				new List<MatchedFragmentIon> { MakePeptideIon(ProductType.b, 4, 410.0, 0.8f) }, 18.0),
		};

		MslLibrary lib = MslConverter.FromMspLibrarySpectra(spectra);

		Assert.That(lib.PrecursorCount, Is.EqualTo(3));
		Assert.That(lib.TryGetLibrarySpectrum("PEPTIDE", 2, out _), Is.True);
		Assert.That(lib.TryGetLibrarySpectrum("ACDEFGHIK", 2, out _), Is.True);
		Assert.That(lib.TryGetLibrarySpectrum("SAMPLER", 3, out _), Is.True);
	}

	// ── Group 4: Query API behaves like a file-loaded library ─────────────────

	/// <summary>
	/// <see cref="MslLibrary.QueryMzWindow"/> must return entries whose
	/// <c>PrecursorMz</c> falls within the requested window.
	/// </summary>
	[Test]
	public void FromMspLibrarySpectra_QueryMzWindow_ReturnsCorrectEntries()
	{
		// Build two precursors at known m/z values
		var spectra = new List<LibrarySpectrum>
		{
			new LibrarySpectrum("AAAA", 500.0, 2,
				new List<MatchedFragmentIon> { MakePeptideIon(ProductType.b, 2, 200.0, 1.0f) }, 5.0),
			new LibrarySpectrum("BBBB", 800.0, 2,
				new List<MatchedFragmentIon> { MakePeptideIon(ProductType.y, 2, 300.0, 1.0f) }, 5.0),
		};
		MslLibrary lib = MslConverter.FromMspLibrarySpectra(spectra);

		// Window that covers only the first precursor
		var hits = lib.QueryMzWindow(490f, 510f);

		Assert.That(hits.Length, Is.EqualTo(1));
		Assert.That(hits[0].PrecursorMz, Is.EqualTo(500.0f).Within(0.01f));
	}

	/// <summary>
	/// The <see cref="MslLibrary"/> produced by <see cref="MslConverter.FromMspFile"/>
	/// must produce identical results to loading the same data via the write-then-read
	/// round-trip used by <see cref="TestMslRegressionMsp"/>.
	/// </summary>
	[Test]
	public void FromMspFile_QueryResults_MatchSaveAndLoadRoundTrip()
	{
		// Build a file-based library for comparison
		string mslPath = Path.Combine(OutputDirectory, "roundtrip_compare.msl");
		MslLibrary.SaveFromLibrarySpectra(mslPath, _mspSpectra);
		using MslLibrary fileLib = MslLibrary.Load(mslPath);

		Assert.That(_fromFile!.PrecursorCount, Is.EqualTo(fileLib.PrecursorCount),
			"In-memory and file-loaded libraries must have the same precursor count.");

		// Spot-check: every entry in the file-loaded library must be findable in the
		// in-memory library with a matching precursor m/z
		foreach (MslLibraryEntry fileEntry in fileLib.GetAllEntries())
		{
			bool found = _fromFile.TryGetEntry(
				fileEntry.ModifiedSequence, fileEntry.Charge, out MslLibraryEntry? memEntry);

			Assert.That(found, Is.True,
				$"Entry '{fileEntry.ModifiedSequence}/{fileEntry.Charge}' missing from in-memory library.");
			Assert.That(memEntry!.PrecursorMz,
				Is.EqualTo(fileEntry.PrecursorMz).Within(MzTolerance),
				$"PrecursorMz mismatch for '{fileEntry.ModifiedSequence}/{fileEntry.Charge}'.");
		}
	}

	// ── Helpers ───────────────────────────────────────────────────────────────

	/// <summary>
	/// Parses a "SEQUENCE/CHARGE" library spectrum name into its two components.
	/// </summary>
	private static (string sequence, int charge) ParseName(string name)
	{
		int slash = name.LastIndexOf('/');
		return (name[..slash], int.Parse(name[(slash + 1)..]));
	}

	/// <summary>
	/// Builds a minimal peptide <see cref="MatchedFragmentIon"/> for use in synthetic
	/// test spectra.
	/// </summary>
	private static MatchedFragmentIon MakePeptideIon(
		ProductType type, int number, double mz, float intensity)
	{
		var product = new Product(
			type,
			type is ProductType.b or ProductType.a or ProductType.c
				? FragmentationTerminus.N
				: FragmentationTerminus.C,
			neutralMass: mz.ToMass(1),
			fragmentNumber: number,
			residuePosition: number,
			neutralLoss: 0.0);

		return new MatchedFragmentIon(product, mz, intensity, charge: 1);
	}

	/// <summary>
	/// Locates <c>syntheticTest2_PredictedLibrary.msp</c> by walking up the directory tree
	/// from the NUnit output directory, checking TestData and SpectralLibraryData subfolders
	/// at each level. Returns a sentinel path (that will not exist) if not found, causing
	/// <see cref="OneTimeSetUp"/> to skip gracefully.
	/// </summary>
	private static string FindMspFile()
	{
		string outputDir = TestContext.CurrentContext.TestDirectory;

		string[] subfolders =
		{
			"TestData",
			Path.Combine("SpectralLibrary", "SpectralLibraryData"),
			"SpectralLibraryData"
		};

		foreach (string sub in subfolders)
		{
			string candidate = Path.Combine(outputDir, sub, MspFileName);
			if (File.Exists(candidate)) return candidate;
		}

		string? dir = outputDir;
		for (int depth = 0; depth < 8 && dir != null; depth++)
		{
			string direct = Path.Combine(dir, MspFileName);
			if (File.Exists(direct)) return direct;

			foreach (string sub in subfolders)
			{
				string sibling = Path.Combine(dir, sub, MspFileName);
				if (File.Exists(sibling)) return sibling;
			}

			dir = Path.GetDirectoryName(dir);
		}

		return Path.Combine(outputDir, "TestData", MspFileName);
	}
}