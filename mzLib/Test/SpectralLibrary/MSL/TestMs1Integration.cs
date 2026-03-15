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
/// NUnit 4 integration tests for Prompt 6: file-type routing, internal ion MSP parsing,
/// MSP → MSL round-trip, <see cref="ToLibrarySpectrum"/> terminus coverage,
/// <see cref="FromLibrarySpectrum"/> completeness, and <see cref="SpectralLibrary"/>
/// MSL routing.
///
/// All file-system writes go to a dedicated temp directory cleaned up in
/// <see cref="OneTimeTearDown"/>.  All assertions use <c>Assert.That</c> exclusively.
/// </summary>
[TestFixture]
public sealed class TestMslIntegration
{
	// ── Fixture paths ─────────────────────────────────────────────────────────

	/// <summary>
	/// Root directory under which all test output files are written.
	/// Isolates this suite from other test suites that also write to temp.
	/// </summary>
	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "MslIntegrationTests");

	/// <summary>
	/// Shared .msl library written once in <see cref="OneTimeSetUp"/> and used by
	/// <see cref="SpectralLibrary"/> routing tests.
	/// </summary>
	private static string SharedMslPath =>
		Path.Combine(OutputDirectory, "shared_integration.msl");

	/// <summary>
	/// Returns a unique per-test .msl path that does not yet exist on disk.
	/// </summary>
	/// <param name="name">Stem used as the filename.</param>
	private static string TempMsl(string name) =>
		Path.Combine(OutputDirectory, name + ".msl");

	/// <summary>
	/// Returns a unique per-test .msp path that does not yet exist on disk.
	/// </summary>
	/// <param name="name">Stem used as the filename.</param>
	private static string TempMsp(string name) =>
		Path.Combine(OutputDirectory, name + ".msp");

	// ── One-time setup / teardown ─────────────────────────────────────────────

	/// <summary>
	/// Creates the output directory and writes the shared fixture .msl file used by
	/// <see cref="SpectralLibrary"/> routing tests.
	/// </summary>
	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		Directory.CreateDirectory(OutputDirectory);
		MslLibrary.Save(SharedMslPath, BuildFixtureEntries());
	}

	/// <summary>Deletes all test output files created during this run.</summary>
	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDirectory))
			Directory.Delete(OutputDirectory, recursive: true);
	}

	// ── Fixture helpers ───────────────────────────────────────────────────────

	/// <summary>
	/// Builds a small two-entry library used by routing and round-trip tests.
	/// Entry 0: PEPTIDE/2, target, b2 and y1 fragments.
	/// Entry 1: ACDEFGHIK/2, decoy, y3 fragment.
	/// </summary>
	private static List<MslLibraryEntry> BuildFixtureEntries()
	{
		var entry0 = new MslLibraryEntry
		{
			ModifiedSequence = "PEPTIDE",
			StrippedSequence = "PEPTIDE",
			PrecursorMz = 449.74,
			Charge = 2,
			Irt = 35.4,
			IsDecoy = false,
			IsProteotypic = true,
			Source = MslFormat.SourceType.Predicted,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Fragments = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 98.060f, Intensity = 1.0f,
					ProductType = ProductType.b, FragmentNumber = 2,
					ResiduePosition = 1, Charge = 1
				},
				new MslFragmentIon
				{
					Mz = 175.119f, Intensity = 0.8f,
					ProductType = ProductType.y, FragmentNumber = 1,
					ResiduePosition = 6, Charge = 1
				}
			}
		};

		var entry1 = new MslLibraryEntry
		{
			ModifiedSequence = "ACDEFGHIK",
			StrippedSequence = "ACDEFGHIK",
			PrecursorMz = 529.76,
			Charge = 2,
			Irt = 42.1,
			IsDecoy = true,
			Source = MslFormat.SourceType.Predicted,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Fragments = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz = 364.197f, Intensity = 1.0f,
					ProductType = ProductType.y, FragmentNumber = 3,
					ResiduePosition = 6, Charge = 1
				}
			}
		};

		return new List<MslLibraryEntry> { entry0, entry1 };
	}

	/// <summary>
	/// Builds a <see cref="LibrarySpectrum"/> containing a mix of terminal and internal
	/// ions.  Used by MSP round-trip tests and <see cref="FromLibrarySpectrum"/> tests.
	/// </summary>
	private static LibrarySpectrum BuildMixedSpectrum()
	{
		// b2 terminal ion
		var b2Product = new Product(ProductType.b, FragmentationTerminus.N,
			neutralMass: 0, fragmentNumber: 2, residuePosition: 2, neutralLoss: 0);
		var b2Ion = new MatchedFragmentIon(b2Product, experMz: 227.1, experIntensity: 1.0, charge: 1);

		// y3 terminal ion
		var y3Product = new Product(ProductType.y, FragmentationTerminus.C,
			neutralMass: 0, fragmentNumber: 3, residuePosition: 5, neutralLoss: 0);
		var y3Ion = new MatchedFragmentIon(y3Product, experMz: 391.2, experIntensity: 0.8, charge: 1);

		// bIb[2-4] internal ion — terminus None
		var bibProduct = new Product(ProductType.b, FragmentationTerminus.None,
			neutralMass: 0, fragmentNumber: 2, residuePosition: 2, neutralLoss: 0);
		var bibIon = new MatchedFragmentIon(bibProduct, experMz: 312.1, experIntensity: 0.3, charge: 1);

		return new LibrarySpectrum(
			sequence: "PEPTIDE",
			precursorMz: 449.74,
			chargeState: 2,
			peaks: new List<MatchedFragmentIon> { b2Ion, y3Ion, bibIon },
			rt: 35.4);
	}

	/// <summary>
	/// Writes a minimal MSP file containing one entry with a mix of terminal and internal
	/// ion annotations.  Returns the path of the written file.
	/// </summary>
	/// <param name="name">Stem used for the file name.</param>
	private static string WriteMspWithInternalIons(string name)
	{
		string path = TempMsp(name);

		// Internal ion annotation format: {PType}I{SType}[{start}-{end}]
		// Terminal annotation format: {type}{number}^{charge}
		string msp = string.Join("\n",
			"Name: PEPTIDE/2",
			"MW: 449.74",
			"Comment: Parent=449.74 iRT=35.4",
			"Num peaks: 3",
			"227.1\t1.0\t\"b2^1/0ppm\"",
			"391.2\t0.8\t\"y3^1/0ppm\"",
			"312.1\t0.3\t\"bIb[2-4]^1/0ppm\"",
			"");

		File.WriteAllText(path, msp);
		return path;
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Group 1 — File-type detection
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>Verifies that <c>.msl</c> extension is recognised.</summary>
	[Test]
	public void IsMslFile_ReturnsTrue_ForMslExtension()
	{
		Assert.That(MslFileTypeHandler.IsMslFile("library.msl"), Is.True);
	}

	/// <summary>Verifies that <c>.msp</c> extension is not recognised as .msl.</summary>
	[Test]
	public void IsMslFile_ReturnsFalse_ForMspExtension()
	{
		Assert.That(MslFileTypeHandler.IsMslFile("library.msp"), Is.False);
	}

	/// <summary>Verifies case-insensitive matching for <c>.MSL</c> and <c>.Msl</c>.</summary>
	[Test]
	public void IsMslFile_CaseInsensitive()
	{
		// Upper case
		Assert.That(MslFileTypeHandler.IsMslFile("Library.MSL"), Is.True);
		// Mixed case
		Assert.That(MslFileTypeHandler.IsMslFile("Library.Msl"), Is.True);
	}

	/// <summary>Verifies that null input returns false rather than throwing.</summary>
	[Test]
	public void IsMslFile_ReturnsFalse_ForNullInput()
	{
		Assert.That(MslFileTypeHandler.IsMslFile(null), Is.False);
	}

	/// <summary>Verifies that empty string returns false rather than throwing.</summary>
	[Test]
	public void IsMslFile_ReturnsFalse_ForEmptyString()
	{
		Assert.That(MslFileTypeHandler.IsMslFile(string.Empty), Is.False);
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Group 2 — SpectralLibrary routing
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// A <see cref="SpectralLibrary"/> constructed with a .msl path must be able to find an
	/// entry via <see cref="SpectralLibrary.TryGetSpectrum"/>.
	/// </summary>
	[Test]
	public void SpectralLibrary_WithMslPath_TryGetSpectrum_FindsEntry()
	{
		var lib = new SpectralLibraryClass(new List<string> { SharedMslPath });

		bool found = lib.TryGetSpectrum("PEPTIDE", 2, out LibrarySpectrum spectrum);
		lib.CloseConnections();

		Assert.That(found, Is.True);
		Assert.That(spectrum, Is.Not.Null);
		Assert.That(spectrum.Name, Is.EqualTo("PEPTIDE/2"));
	}

	/// <summary>
	/// <see cref="SpectralLibrary.ContainsSpectrum"/> must return true for a sequence present
	/// in an .msl library and false for one that is absent.
	/// </summary>
	[Test]
	public void SpectralLibrary_WithMslPath_ContainsSpectrum_CorrectResults()
	{
		var lib = new SpectralLibraryClass(new List<string> { SharedMslPath });

		Assert.That(lib.ContainsSpectrum("PEPTIDE", 2), Is.True);
		Assert.That(lib.ContainsSpectrum("NOTINLIB", 2), Is.False);
		lib.CloseConnections();
	}

	/// <summary>
	/// When a mix of .msl and .msp paths are supplied, spectra from both sources must be
	/// reachable via <see cref="SpectralLibrary.TryGetSpectrum"/>.
	/// </summary>
	[Test]
	public void SpectralLibrary_WithMixedPaths_MspAndMsl_BothSearched()
	{
		// Write a minimal MSP with a different sequence
		string mspPath = TempMsp("mixed_test");
		File.WriteAllText(mspPath, string.Join("\n",
			"Name: ACDEFGHIK/2",
			"MW: 529.76",
			"Comment: Parent=529.76 iRT=42.1",
			"Num peaks: 1",
			"364.197\t1.0\t\"y3^1/0ppm\"",
			""));

		var lib = new SpectralLibraryClass(new List<string> { SharedMslPath, mspPath });

		// MSL entry
		bool foundMsl = lib.TryGetSpectrum("PEPTIDE", 2, out LibrarySpectrum mslSpectrum);
		// MSP entry
		bool foundMsp = lib.TryGetSpectrum("ACDEFGHIK", 2, out LibrarySpectrum mspSpectrum);
		lib.CloseConnections();

		Assert.That(foundMsl, Is.True, "MSL entry not found");
		Assert.That(foundMsp, Is.True, "MSP entry not found");
		Assert.That(mslSpectrum.Name, Is.EqualTo("PEPTIDE/2"));
		Assert.That(mspSpectrum.Name, Is.EqualTo("ACDEFGHIK/2"));
	}

	/// <summary>
	/// <see cref="SpectralLibrary.CloseConnections"/> must not throw, and the underlying
	/// .msl file must be deletable afterwards (confirming the file handle was released).
	/// </summary>
	[Test]
	public void SpectralLibrary_CloseConnections_DisposesAllMslLibraries()
	{
		// Write a dedicated .msl for this test so we can attempt deletion
		string path = TempMsl("close_connections_test");
		MslLibrary.Save(path, BuildFixtureEntries());

		var lib = new SpectralLibraryClass(new List<string> { path });

		Assert.That(() => lib.CloseConnections(), Throws.Nothing);

		// On Windows, if the file handle is still open the delete will throw
		Assert.That(() => File.Delete(path), Throws.Nothing);
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Group 3 — Internal ion MSP parsing via ReadFragmentIon
	// ═════════════════════════════════════════════════════════════════════════

	// Shared split chars that match those used by ReadLibrarySpectrum
	private static readonly char[] FragSplit = { '\t', '"', ')', '/' };
	private static readonly char[] NeutralLossSplit = { '-' };

	/// <summary>
	/// <c>bIb[3-6]^1</c> must be parsed without throwing and produce a non-null ion.
	/// </summary>
	[Test]
	public void ReadFragmentIon_InternalBIb_ParsesCorrectly()
	{
		// Tab-separated: mz \t intensity \t annotation
		string line = "312.1\t0.3\t\"bIb[3-6]^1/0ppm\"";

		MatchedFragmentIon ion = SpectralLibraryClass.ReadFragmentIon(
			line, FragSplit, NeutralLossSplit, peptideSequence: "PEPTIDE");

		Assert.That(ion, Is.Not.Null);
	}

	/// <summary>
	/// <c>aIb[2-5]^1</c> must parse to a product with ProductType.a.
	/// </summary>
	[Test]
	public void ReadFragmentIon_InternalAIb_ParsesCorrectly()
	{
		string line = "256.1\t0.2\t\"aIb[2-5]^1/0ppm\"";

		MatchedFragmentIon ion = SpectralLibraryClass.ReadFragmentIon(
			line, FragSplit, NeutralLossSplit, peptideSequence: "PEPTIDE");

		Assert.That(ion.NeutralTheoreticalProduct.ProductType, Is.EqualTo(ProductType.a));
	}

	/// <summary>
	/// Internal ions must always have <see cref="FragmentationTerminus.None"/> as their
	/// terminus regardless of which product type letters are used.
	/// </summary>
	[Test]
	public void ReadFragmentIon_InternalIon_HasTerminusNone()
	{
		string line = "312.1\t0.3\t\"bIb[3-6]^1/0ppm\"";

		MatchedFragmentIon ion = SpectralLibraryClass.ReadFragmentIon(
			line, FragSplit, NeutralLossSplit, peptideSequence: "PEPTIDE");

		Assert.That(ion.NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.None));
	}

	/// <summary>
	/// The start residue parsed from <c>bIb[3-6]</c> must become
	/// <c>Product.FragmentNumber</c> == 3.
	/// </summary>
	[Test]
	public void ReadFragmentIon_InternalIon_StartAndEndResidue_Correct()
	{
		string line = "312.1\t0.3\t\"bIb[3-6]^1/0ppm\"";

		MatchedFragmentIon ion = SpectralLibraryClass.ReadFragmentIon(
			line, FragSplit, NeutralLossSplit, peptideSequence: "PEPTIDE");

		// FragmentNumber carries the start residue for internal ions
		Assert.That(ion.NeutralTheoreticalProduct.FragmentNumber, Is.EqualTo(3));
	}

	/// <summary>
	/// A standard terminal ion annotation (<c>b5^1</c>) must still parse correctly after
	/// the internal-ion regex was added — verifying no regression.
	/// </summary>
	[Test]
	public void ReadFragmentIon_TerminalIon_Unchanged_AfterUpdate()
	{
		string line = "567.3\t0.9\t\"b5^1/0ppm\"";

		MatchedFragmentIon ion = SpectralLibraryClass.ReadFragmentIon(
			line, FragSplit, NeutralLossSplit, peptideSequence: "PEPTIDE");

		Assert.That(ion.NeutralTheoreticalProduct.ProductType, Is.EqualTo(ProductType.b));
		Assert.That(ion.NeutralTheoreticalProduct.FragmentNumber, Is.EqualTo(5));
		Assert.That(ion.NeutralTheoreticalProduct.Terminus, Is.EqualTo(FragmentationTerminus.N));
	}

	/// <summary>
	/// An internal ion annotated with charge 2 (<c>bIb[3-6]^2</c>) must parse with
	/// <see cref="MatchedFragmentIon.Charge"/> == 2.
	/// </summary>
	[Test]
	public void ReadFragmentIon_InternalIon_WithCharge2_ParsesCharge()
	{
		string line = "156.6\t0.1\t\"bIb[3-6]^2/0ppm\"";

		MatchedFragmentIon ion = SpectralLibraryClass.ReadFragmentIon(
			line, FragSplit, NeutralLossSplit, peptideSequence: "PEPTIDE");

		Assert.That(ion.Charge, Is.EqualTo(2));
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Group 4 — MSP → MSL round-trip (internal ions)
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// A <see cref="LibrarySpectrum"/> read from an MSP containing internal ions can be
	/// saved to .msl and read back without error.
	/// </summary>
	[Test]
	public void MspLibraryWithInternalIons_WritesToMsl_AndReadsBack()
	{
		string mspPath = WriteMspWithInternalIons("roundtrip_basic");
		string mslPath = TempMsl("roundtrip_basic");

		// Read the MSP
		var mspLib = new SpectralLibraryClass(new List<string> { mspPath });
		bool found = mspLib.TryGetSpectrum("PEPTIDE", 2, out LibrarySpectrum spectrum);
		mspLib.CloseConnections();

		Assert.That(found, Is.True, "MSP spectrum not found");

		// Convert and save to MSL
		MslLibrary.SaveFromLibrarySpectra(mslPath, new List<LibrarySpectrum> { spectrum });

		// Read back
		using MslLibrary mslLib = MslLibrary.Load(mslPath);
		bool foundBack = mslLib.TryGetLibrarySpectrum("PEPTIDE", 2, out LibrarySpectrum back);

		Assert.That(foundBack, Is.True);
		Assert.That(back, Is.Not.Null);
		Assert.That(back.MatchedFragmentIons.Count, Is.GreaterThan(0));
	}

	/// <summary>
	/// After an MSP → MSL round-trip the number of fragment ions is preserved.
	/// </summary>
	[Test]
	public void MspToMsl_FragmentCount_Preserved()
	{
		string mspPath = WriteMspWithInternalIons("roundtrip_fragcount");
		string mslPath = TempMsl("roundtrip_fragcount");

		var mspLib = new SpectralLibraryClass(new List<string> { mspPath });
		mspLib.TryGetSpectrum("PEPTIDE", 2, out LibrarySpectrum spectrum);
		mspLib.CloseConnections();

		int originalCount = spectrum.MatchedFragmentIons.Count;

		MslLibrary.SaveFromLibrarySpectra(mslPath, new List<LibrarySpectrum> { spectrum });

		using MslLibrary mslLib = MslLibrary.Load(mslPath);
		mslLib.TryGetLibrarySpectrum("PEPTIDE", 2, out LibrarySpectrum back);

		Assert.That(back.MatchedFragmentIons.Count, Is.EqualTo(originalCount));
	}

	/// <summary>
	/// After an MSP → MSL round-trip terminal ions retain their original
	/// <see cref="FragmentationTerminus"/> values.
	/// </summary>
	[Test]
	public void MspToMsl_TerminalIon_Terminus_Preserved()
	{
		string mspPath = WriteMspWithInternalIons("roundtrip_terminus");
		string mslPath = TempMsl("roundtrip_terminus");

		var mspLib = new SpectralLibraryClass(new List<string> { mspPath });
		mspLib.TryGetSpectrum("PEPTIDE", 2, out LibrarySpectrum spectrum);
		mspLib.CloseConnections();

		MslLibrary.SaveFromLibrarySpectra(mslPath, new List<LibrarySpectrum> { spectrum });

		using MslLibrary mslLib = MslLibrary.Load(mslPath);
		mslLib.TryGetLibrarySpectrum("PEPTIDE", 2, out LibrarySpectrum back);

		// b2 ion must still be N-terminal
		var bIon = back.MatchedFragmentIons
			.FirstOrDefault(i => i.NeutralTheoreticalProduct.ProductType == ProductType.b
							   && i.NeutralTheoreticalProduct.FragmentNumber == 2);
		Assert.That(bIon, Is.Not.Null, "b2 ion not found after round-trip");
		Assert.That(bIon.NeutralTheoreticalProduct.Terminus, Is.EqualTo(FragmentationTerminus.N));
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Group 5 — ToLibrarySpectrum terminus coverage
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Helper: builds a minimal <see cref="MslLibraryEntry"/> with one fragment of the
	/// specified type and calls <see cref="MslLibraryEntry.ToLibrarySpectrum"/>.
	/// </summary>
	private static LibrarySpectrum MakeSpectrumForTerminusTest(
		ProductType productType,
		MslFormat.MoleculeType molType,
		bool isInternal = false)
	{
		var entry = new MslLibraryEntry
		{
			ModifiedSequence = "PEPTIDE",
			StrippedSequence = "PEPTIDE",
			PrecursorMz = 449.74,
			Charge = 2,
			MoleculeType = molType,
			Fragments = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz                   = 300.0f,
					Intensity            = 1.0f,
					ProductType          = productType,
					SecondaryProductType = isInternal ? ProductType.b : (ProductType?)null,
					FragmentNumber       = 3,
					SecondaryFragmentNumber = isInternal ? 6 : 0,
					ResiduePosition      = 3,
					Charge               = 1
				}
			}
		};
		return entry.ToLibrarySpectrum();
	}

	/// <summary>Peptide b ions must have <see cref="FragmentationTerminus.N"/>.</summary>
	[Test]
	public void ToLibrarySpectrum_PeptideBIon_HasNTerminus()
	{
		var spectrum = MakeSpectrumForTerminusTest(ProductType.b, MslFormat.MoleculeType.Peptide);
		Assert.That(spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.N));
	}

	/// <summary>Peptide y ions must have <see cref="FragmentationTerminus.C"/>.</summary>
	[Test]
	public void ToLibrarySpectrum_PeptideYIon_HasCTerminus()
	{
		var spectrum = MakeSpectrumForTerminusTest(ProductType.y, MslFormat.MoleculeType.Peptide);
		Assert.That(spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.C));
	}

	/// <summary>Peptide c ions must have <see cref="FragmentationTerminus.N"/>.</summary>
	[Test]
	public void ToLibrarySpectrum_PeptideCIon_HasNTerminus()
	{
		var spectrum = MakeSpectrumForTerminusTest(ProductType.c, MslFormat.MoleculeType.Peptide);
		Assert.That(spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.N));
	}

	/// <summary>Peptide a ions must have <see cref="FragmentationTerminus.N"/>.</summary>
	[Test]
	public void ToLibrarySpectrum_PeptideAIon_HasNTerminus()
	{
		var spectrum = MakeSpectrumForTerminusTest(ProductType.a, MslFormat.MoleculeType.Peptide);
		Assert.That(spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.N));
	}

	/// <summary>Peptide z ions must have <see cref="FragmentationTerminus.C"/>.</summary>
	[Test]
	public void ToLibrarySpectrum_PeptideZIon_HasCTerminus()
	{
		var spectrum = MakeSpectrumForTerminusTest(ProductType.z, MslFormat.MoleculeType.Peptide);
		Assert.That(spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.C));
	}

	/// <summary>
	/// zDot is not in the terminus dictionary so it falls through to
	/// <see cref="FragmentationTerminus.None"/>.
	/// </summary>
	[Test]
	public void ToLibrarySpectrum_PeptideZDotIon_HasTerminusNone()
	{
		var spectrum = MakeSpectrumForTerminusTest(ProductType.zDot, MslFormat.MoleculeType.Peptide);
		Assert.That(spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.None));
	}

	/// <summary>Oligo a ions must have <see cref="FragmentationTerminus.FivePrime"/>.</summary>
	[Test]
	public void ToLibrarySpectrum_OligoAIon_HasFivePrimeTerminus()
	{
		var spectrum = MakeSpectrumForTerminusTest(ProductType.a, MslFormat.MoleculeType.Oligonucleotide);
		Assert.That(spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.FivePrime));
	}

	/// <summary>Oligo w ions must have <see cref="FragmentationTerminus.ThreePrime"/>.</summary>
	[Test]
	public void ToLibrarySpectrum_OligoWIon_HasThreePrimeTerminus()
	{
		var spectrum = MakeSpectrumForTerminusTest(ProductType.w, MslFormat.MoleculeType.Oligonucleotide);
		Assert.That(spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.ThreePrime));
	}

	/// <summary>D (diagnostic) ions must have <see cref="FragmentationTerminus.None"/>.</summary>
	[Test]
	public void ToLibrarySpectrum_DiagnosticDIon_HasTerminusNone()
	{
		var spectrum = MakeSpectrumForTerminusTest(ProductType.D, MslFormat.MoleculeType.Peptide);
		Assert.That(spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.None));
	}

	/// <summary>M (precursor) ions must have <see cref="FragmentationTerminus.None"/>.</summary>
	[Test]
	public void ToLibrarySpectrum_MIon_HasTerminusNone()
	{
		var spectrum = MakeSpectrumForTerminusTest(ProductType.M, MslFormat.MoleculeType.Peptide);
		Assert.That(spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.None));
	}

	/// <summary>
	/// Internal fragment ions (SecondaryProductType != null) must always have
	/// <see cref="FragmentationTerminus.None"/>.
	/// </summary>
	[Test]
	public void ToLibrarySpectrum_InternalIon_HasTerminusNone()
	{
		var spectrum = MakeSpectrumForTerminusTest(
			ProductType.b, MslFormat.MoleculeType.Peptide, isInternal: true);
		Assert.That(spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus,
			Is.EqualTo(FragmentationTerminus.None));
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Group 6 — FromLibrarySpectrum completeness
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// <see cref="MslLibraryEntry.FromLibrarySpectrum"/> must preserve
	/// <see cref="ProductType"/> for all fragment ions.
	/// </summary>
	[Test]
	public void FromLibrarySpectrum_AllFragmentTypes_ProductTypePreserved()
	{
		LibrarySpectrum spectrum = BuildMixedSpectrum();
		MslLibraryEntry entry = MslLibraryEntry.FromLibrarySpectrum(spectrum);

		var originalTypes = spectrum.MatchedFragmentIons
			.Select(i => i.NeutralTheoreticalProduct.ProductType)
			.ToList();
		var roundTrippedTypes = entry.Fragments
			.Select(f => f.ProductType)
			.ToList();

		Assert.That(roundTrippedTypes, Is.EquivalentTo(originalTypes));
	}

	/// <summary>
	/// <see cref="MslLibraryEntry.FromLibrarySpectrum"/> must preserve the neutral-loss
	/// mass stored in each <see cref="Product"/>.
	/// </summary>
	[Test]
	public void FromLibrarySpectrum_NeutralLoss_Preserved()
	{
		// Build a spectrum with an H2O loss on the y ion
		double h2oLoss = -18.010565;
		var yProduct = new Product(ProductType.y, FragmentationTerminus.C,
			neutralMass: 0, fragmentNumber: 3, residuePosition: 5, neutralLoss: h2oLoss);
		var yIon = new MatchedFragmentIon(yProduct, experMz: 373.2, experIntensity: 0.6, charge: 1);

		var spectrum = new LibrarySpectrum("PEPTIDE", 449.74, 2,
			new List<MatchedFragmentIon> { yIon }, rt: 35.4);

		MslLibraryEntry entry = MslLibraryEntry.FromLibrarySpectrum(spectrum);

		Assert.That(entry.Fragments[0].NeutralLoss,
			Is.EqualTo(h2oLoss).Within(1e-4));
	}

	/// <summary>
	/// <see cref="MslLibraryEntry.FromLibrarySpectrum"/> must preserve fragment charge.
	/// </summary>
	[Test]
	public void FromLibrarySpectrum_Charge_Preserved()
	{
		var product = new Product(ProductType.b, FragmentationTerminus.N,
			neutralMass: 0, fragmentNumber: 2, residuePosition: 2, neutralLoss: 0);
		var ion = new MatchedFragmentIon(product, experMz: 114.1, experIntensity: 1.0, charge: 2);

		var spectrum = new LibrarySpectrum("PEPTIDE", 449.74, 2,
			new List<MatchedFragmentIon> { ion }, rt: 35.4);

		MslLibraryEntry entry = MslLibraryEntry.FromLibrarySpectrum(spectrum);

		Assert.That(entry.Fragments[0].Charge, Is.EqualTo(2));
	}

	/// <summary>
	/// <see cref="MslLibraryEntry.FromLibrarySpectrum"/> must preserve the
	/// <see cref="LibrarySpectrum.Sequence"/> as <see cref="MslLibraryEntry.ModifiedSequence"/>.
	/// </summary>
	[Test]
	public void FromLibrarySpectrum_Sequence_Preserved()
	{
		LibrarySpectrum spectrum = BuildMixedSpectrum();
		MslLibraryEntry entry = MslLibraryEntry.FromLibrarySpectrum(spectrum);

		Assert.That(entry.ModifiedSequence, Is.EqualTo(spectrum.Sequence));
	}

	/// <summary>
	/// <see cref="MslLibraryEntry.FromLibrarySpectrum"/> must throw
	/// <see cref="ArgumentNullException"/> when passed null.
	/// </summary>
	[Test]
	public void FromLibrarySpectrum_NullInput_Throws()
	{
		Assert.That(() => MslLibraryEntry.FromLibrarySpectrum(null),
			Throws.ArgumentNullException);
	}

	// ═════════════════════════════════════════════════════════════════════════
	// Group 7 — MslFileTypeHandler.Open smoke test
	// ═════════════════════════════════════════════════════════════════════════

	/// <summary>
	/// <see cref="MslFileTypeHandler.Open"/> on an existing file must return a non-null
	/// <see cref="MslLibrary"/> with the correct precursor count.
	/// </summary>
	[Test]
	public void MslFileTypeHandler_Open_ReturnsLibraryWithCorrectCount()
	{
		using MslLibrary lib = MslFileTypeHandler.Open(SharedMslPath);

		Assert.That(lib, Is.Not.Null);
		Assert.That(lib.PrecursorCount, Is.EqualTo(2));
	}

	/// <summary>
	/// <see cref="MslFileTypeHandler.Open"/> on a missing file must throw
	/// <see cref="FileNotFoundException"/>.
	/// </summary>
	[Test]
	public void MslFileTypeHandler_Open_MissingFile_Throws()
	{
		Assert.That(
			() => MslFileTypeHandler.Open(Path.Combine(OutputDirectory, "does_not_exist.msl")),
			Throws.TypeOf<FileNotFoundException>());
	}

	/// <summary>
	/// A very small threshold forces index-only mode; the library must still find entries.
	/// </summary>
	[Test]
	public void MslFileTypeHandler_Open_SmallThreshold_UsesIndexOnlyMode()
	{
		// Threshold of 1 byte forces any real file to use index-only mode
		using MslLibrary lib = MslFileTypeHandler.Open(SharedMslPath, indexOnlyThresholdBytes: 1);

		Assert.That(lib.IsIndexOnly, Is.True);

		bool found = lib.TryGetLibrarySpectrum("PEPTIDE", 2, out LibrarySpectrum spectrum);
		Assert.That(found, Is.True);
		Assert.That(spectrum, Is.Not.Null);
	}
}