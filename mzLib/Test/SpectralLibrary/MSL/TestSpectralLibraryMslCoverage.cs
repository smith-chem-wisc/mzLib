using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Fragmentation.Peptide;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;

namespace Test.MslSpectralLibrary;

/// <summary>
/// Coverage tests for the MSL-specific code paths added to SpectralLibrary.cs and
/// MslFileTypeHandler.cs by PR #1036.
///
/// Following the convention of the existing test suite, this class tests the MSL
/// infrastructure directly via MslFileTypeHandler and MslLibrary rather than
/// instantiating the SpectralLibrary facade class (which carries a namespace/type
/// name collision that the existing tests intentionally avoid).
///
/// Surfaces covered:
///   1. MslFileTypeHandler.IsMslFile   — all branches
///   2. MslFileTypeHandler.Open        — full-load path, index-only threshold, FileNotFoundException
///   3. SpectralLibrary.ReadFragmentIon — the new internal-ion path added by this PR:
///        a. charge-absent  (Groups[5].Success == false -> default 1)
///        b. charge-explicit (Groups[5].Success == true  -> parsed)
///        c. distinct primary / secondary product types
///        d. mutual exclusion with the terminal-ion path
///        e. null peptide sequence tolerance
/// </summary>
[TestFixture]
public sealed class TestSpectralLibraryMslCoverage
{
	// ── Fixture setup ─────────────────────────────────────────────────────────

	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "SpectralLibraryMslCoverageTests");

	[OneTimeSetUp]
	public void OneTimeSetUp() => Directory.CreateDirectory(OutputDirectory);

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDirectory))
			Directory.Delete(OutputDirectory, recursive: true);
	}

	// ── Helpers ───────────────────────────────────────────────────────────────

	private static string TempPath(string name) =>
		Path.Combine(OutputDirectory, name + ".msl");

	/// <summary>
	/// Split chars that exactly match what ReadLibrarySpectrum passes to ReadFragmentIon.
	/// </summary>
	private static readonly char[] FragmentSplit = { '\t', '"', ')', '/' };
	private static readonly char[] NeutralLossSplit = { '-' };

	private static MslLibraryEntry MakeEntry(string sequence = "PEPTIDE", int charge = 2,
		double precursorMz = 449.74, double irt = 35.4) =>
		new MslLibraryEntry
		{
			ModifiedSequence = sequence,
			StrippedSequence = sequence,
			PrecursorMz = precursorMz,
			Charge = charge,
			Irt = irt,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Fragments = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 175.119f, Intensity = 1.0f,
					ProductType = ProductType.y, FragmentNumber = 1,
					ResiduePosition = 6, Charge = 1
				}
			}
		};

	// ── MslFileTypeHandler.IsMslFile ──────────────────────────────────────────

	[Test]
	public void IsMslFile_MslExtension_ReturnsTrue()
	{
		Assert.That(MslFileTypeHandler.IsMslFile("library.msl"), Is.True);
	}

	[Test]
	public void IsMslFile_MslExtensionUpperCase_ReturnsTrue()
	{
		Assert.That(MslFileTypeHandler.IsMslFile("LIBRARY.MSL"), Is.True,
			"Extension check must be case-insensitive");
	}

	[Test]
	public void IsMslFile_MspExtension_ReturnsFalse()
	{
		Assert.That(MslFileTypeHandler.IsMslFile("library.msp"), Is.False);
	}

	[Test]
	public void IsMslFile_NullOrEmpty_ReturnsFalse()
	{
		Assert.That(MslFileTypeHandler.IsMslFile(null), Is.False);
		Assert.That(MslFileTypeHandler.IsMslFile(""), Is.False);
	}

	[Test]
	public void IsMslFile_NoExtension_ReturnsFalse()
	{
		Assert.That(MslFileTypeHandler.IsMslFile("library"), Is.False);
	}

	// ── MslFileTypeHandler.Open — FileNotFoundException ───────────────────────

	[Test]
	public void Open_NonExistentFile_ThrowsFileNotFoundException()
	{
		Assert.Throws<FileNotFoundException>(
			() => MslFileTypeHandler.Open(Path.Combine(OutputDirectory, "does_not_exist.msl")),
			"Open must throw FileNotFoundException for a missing file");
	}

	// ── MslFileTypeHandler.Open — full-load path (file < threshold) ───────────

	[Test]
	public void Open_SmallFile_ReturnsFullyLoadedLibrary()
	{
		string path = TempPath("open_full");
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry() });

		// Default threshold is 1 GiB; a tiny test file is always below it
		using MslLibrary lib = MslFileTypeHandler.Open(path);

		Assert.That(lib.PrecursorCount, Is.EqualTo(1));
		Assert.That(lib.IsIndexOnly, Is.False,
			"A file below the size threshold must be opened with full Load, not LoadIndexOnly");
	}

	[Test]
	public void Open_SmallFile_EntryIsReachable()
	{
		string path = TempPath("open_entry");
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry() });

		using MslLibrary lib = MslFileTypeHandler.Open(path);

		Assert.That(lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? entry), Is.True);
		Assert.That(entry!.PrecursorMz, Is.EqualTo(449.74f).Within(0.01f));
	}

	// ── MslFileTypeHandler.Open — index-only path (threshold override) ────────

	[Test]
	public void Open_ThresholdSetToZero_ReturnsIndexOnlyLibrary()
	{
		string path = TempPath("open_indexonly");
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry() });

		// Passing threshold = 0 forces index-only regardless of file size
		using MslLibrary lib = MslFileTypeHandler.Open(path, indexOnlyThresholdBytes: 0);

		Assert.That(lib.IsIndexOnly, Is.True,
			"A threshold of 0 must route every file through LoadIndexOnly");
		Assert.That(lib.PrecursorCount, Is.EqualTo(1));
	}

	[Test]
	public void Open_IndexOnlyMode_EntryIsReachableOnDemand()
	{
		string path = TempPath("open_indexonly_entry");
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry() });

		using MslLibrary lib = MslFileTypeHandler.Open(path, indexOnlyThresholdBytes: 0);

		Assert.That(lib.TryGetEntry("PEPTIDE", 2, out MslLibraryEntry? entry), Is.True,
			"Demand-loaded entry must be reachable in index-only mode");
		Assert.That(entry!.Fragments, Is.Not.Empty,
			"Fragments must be fetched on demand in index-only mode");
	}

	// ── SpectralLibrary.ReadFragmentIon — internal ion path ───────────────────
	//
	// ReadFragmentIon is public static. The fully qualified call below sidesteps
	// the CS0118 ambiguity between the Readers.SpectralLibrary namespace and the
	// Readers.SpectralLibrary.SpectralLibrary class — the same reason the existing
	// test suite never instantiates SpectralLibrary directly.

	/// <summary>
	/// bIb[3-6] — charge group absent -> must default to 1.
	/// Covers the Groups[5].Success == false branch.
	/// </summary>
	[Test]
	public void ReadFragmentIon_InternalIon_ChargeAbsent_DefaultsToOne()
	{
		string line = "350.175\t0.85\tbIb[3-6]";

		MatchedFragmentIon ion = Readers.SpectralLibrary.SpectralLibrary.ReadFragmentIon(
			line, FragmentSplit, NeutralLossSplit, "PEPTIDER");

		Assert.That(ion.Charge, Is.EqualTo(1),
			"Charge must default to 1 when absent from the internal ion annotation");
		Assert.That(ion.Mz, Is.EqualTo(350.175).Within(1e-3));
		Assert.That(ion.Intensity, Is.EqualTo(0.85).Within(1e-6));
		Assert.That(ion.NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.None),
			"Internal ions must carry FragmentationTerminus.None");
		Assert.That(ion.NeutralTheoreticalProduct.FragmentNumber, Is.EqualTo(3),
			"FragmentNumber must equal the start residue");
		Assert.That(ion.NeutralTheoreticalProduct.SecondaryFragmentNumber, Is.EqualTo(6),
			"SecondaryFragmentNumber must equal the end residue");
	}

	/// <summary>
	/// aIb[2-5]^2 — charge group present -> must parse the explicit value.
	/// Covers the Groups[5].Success == true branch.
	/// </summary>
	[Test]
	public void ReadFragmentIon_InternalIon_ChargeExplicit_ParsedCorrectly()
	{
		string line = "263.140\t0.60\taIb[2-5]^2";

		MatchedFragmentIon ion = Readers.SpectralLibrary.SpectralLibrary.ReadFragmentIon(
			line, FragmentSplit, NeutralLossSplit, "PEPTIDER");

		Assert.That(ion.Charge, Is.EqualTo(2),
			"Explicit ^2 suffix must be parsed from the internal ion annotation");
		Assert.That(ion.NeutralTheoreticalProduct.Terminus, Is.EqualTo(FragmentationTerminus.None));
		Assert.That(ion.NeutralTheoreticalProduct.FragmentNumber, Is.EqualTo(2));
		Assert.That(ion.NeutralTheoreticalProduct.SecondaryFragmentNumber, Is.EqualTo(5));
	}

	/// <summary>
	/// aIy[1-4] — primary=a, secondary=y. Both Enum.Parse calls must succeed with
	/// distinct type strings.
	/// </summary>
	[Test]
	public void ReadFragmentIon_InternalIon_DistinctPrimaryAndSecondaryTypes()
	{
		string line = "215.105\t1.00\taIy[1-4]";

		MatchedFragmentIon ion = Readers.SpectralLibrary.SpectralLibrary.ReadFragmentIon(
			line, FragmentSplit, NeutralLossSplit, "PEPTIDER");

		Assert.That(ion.NeutralTheoreticalProduct.ProductType,
			Is.EqualTo(ProductType.a), "Primary product type must be 'a'");
		Assert.That(ion.NeutralTheoreticalProduct.SecondaryProductType,
			Is.EqualTo(ProductType.y), "Secondary product type must be 'y'");
	}

	/// <summary>
	/// y5 — a terminal ion must not match the internal-ion regex and must fall
	/// through to the original terminal-ion parser, producing FragmentationTerminus.C.
	/// Confirms the two paths are mutually exclusive.
	/// </summary>
	[Test]
	public void ReadFragmentIon_TerminalIon_DoesNotMatchInternalIonPath()
	{
		string line = "689.378\t1.00\ty5";

		MatchedFragmentIon ion = Readers.SpectralLibrary.SpectralLibrary.ReadFragmentIon(
			line, FragmentSplit, NeutralLossSplit, "PEPTIDER");

		Assert.That(ion.NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.C),
			"Terminal y ion must resolve to C terminus, not None");
		Assert.That(ion.NeutralTheoreticalProduct.ProductType, Is.EqualTo(ProductType.y));
	}

	/// <summary>
	/// Internal ions must parse correctly even when peptideSequence is null —
	/// they never reach the peptideLength branch that uses that parameter.
	/// </summary>
	[Test]
	public void ReadFragmentIon_InternalIon_NullPeptideSequence_DoesNotThrow()
	{
		string line = "350.175\t0.90\tbIb[3-6]";

		Assert.DoesNotThrow(() =>
		{
			MatchedFragmentIon ion = Readers.SpectralLibrary.SpectralLibrary.ReadFragmentIon(
				line, FragmentSplit, NeutralLossSplit, peptideSequence: null);
			Assert.That(ion.NeutralTheoreticalProduct.Terminus,
				Is.EqualTo(FragmentationTerminus.None));
		}, "Internal ion parsing must not throw when peptideSequence is null");
	}
}