using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using Readers.SpectralLibrary;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.MslSpectralLibrary;

/// <summary>
/// NUnit 4 tests for Prompt 12: verifies that <see cref="MslLibraryEntry.ToLibrarySpectrum"/>
/// and <see cref="MslLibraryEntry.FromLibrarySpectrum"/> correctly propagate
/// <see cref="Product.SecondaryProductType"/> and <see cref="Product.SecondaryFragmentNumber"/>
/// for internal fragment ions.
///
/// Prior to Prompt 12 both conversion methods silently dropped these fields, causing internal-ion
/// boundary information (e.g. bIy[3-6]) to be lost on every <c>LibrarySpectrum</c> round-trip.
/// These tests confirm the fix is in place and regression-safe.
///
/// File-system writes (binary round-trip test) go to a temp directory cleaned up in
/// <see cref="OneTimeTearDown"/>.
/// </summary>
[TestFixture]
public sealed class TestMslInternalIonRoundTrip
{
	// ── Temp directory ────────────────────────────────────────────────────────

	private static readonly string OutputDirectory =
		Path.Combine(Path.GetTempPath(), "MslInternalIonRoundTripTests");

	private static string TempMsl(string name) =>
		Path.Combine(OutputDirectory, name + ".msl");

	[OneTimeSetUp]
	public void OneTimeSetUp() => Directory.CreateDirectory(OutputDirectory);

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDirectory))
			Directory.Delete(OutputDirectory, recursive: true);
	}

	// ── Shared factory ────────────────────────────────────────────────────────

	/// <summary>
	/// Canonical test entry containing one terminal b3 ion and one bIy[3-6] internal ion.
	/// Used across tests that need a predictable, minimal internal-ion entry.
	/// </summary>
	private static MslLibraryEntry MakeInternalIonEntry() =>
		new MslLibraryEntry
		{
			ModifiedSequence = "PEPTIDE",
			StrippedSequence = "PEPTIDE",
			PrecursorMz = 400.2,
			Charge = 2,
			Irt = 25.0,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Fragments = new List<MslFragmentIon>
			{
                // Terminal b ion
                new MslFragmentIon
				{
					ProductType             = ProductType.b,
					SecondaryProductType    = null,
					FragmentNumber          = 3,
					SecondaryFragmentNumber = 0,
					ResiduePosition         = 3,
					Charge                  = 1,
					Mz                      = 312.15f,
					Intensity               = 0.9f,
					NeutralLoss             = 0.0
				},
                // Internal bIy ion spanning residues 3–6
                new MslFragmentIon
				{
					ProductType             = ProductType.b,
					SecondaryProductType    = ProductType.y,
					FragmentNumber          = 3,
					SecondaryFragmentNumber = 6,
					ResiduePosition         = 0,
					Charge                  = 1,
					Mz                      = 415.22f,
					Intensity               = 0.5f,
					NeutralLoss             = 0.0
				}
			}
		};

	/// <summary>
	/// Returns the single internal-ion <see cref="MatchedFragmentIon"/> from a
	/// <see cref="LibrarySpectrum"/> produced by <see cref="MakeInternalIonEntry"/>.
	/// Throws if none is found.
	/// </summary>
	private static MatchedFragmentIon GetInternalIon(LibrarySpectrum spectrum) =>
		spectrum.MatchedFragmentIons.Single(i => i.NeutralTheoreticalProduct.IsInternalFragment);

	/// <summary>
	/// Returns the single terminal <see cref="MatchedFragmentIon"/> from a
	/// <see cref="LibrarySpectrum"/> produced by <see cref="MakeInternalIonEntry"/>.
	/// Throws if none is found.
	/// </summary>
	private static MatchedFragmentIon GetTerminalIon(LibrarySpectrum spectrum) =>
		spectrum.MatchedFragmentIons.Single(i => !i.NeutralTheoreticalProduct.IsInternalFragment);

	// ── ToLibrarySpectrum tests ───────────────────────────────────────────────

	/// <summary>
	/// A bIy[3-6] internal ion must produce a Product whose SecondaryProductType is y.
	/// </summary>
	[Test]
	public void InternalIon_ToLibrarySpectrum_SecondaryProductType_Preserved()
	{
		var spectrum = MakeInternalIonEntry().ToLibrarySpectrum();
		var internalIon = GetInternalIon(spectrum);

		Assert.That(internalIon.NeutralTheoreticalProduct.SecondaryProductType,
			Is.EqualTo(ProductType.y));
	}

	/// <summary>
	/// A bIy[3-6] internal ion must produce a Product whose SecondaryFragmentNumber is 6.
	/// </summary>
	[Test]
	public void InternalIon_ToLibrarySpectrum_SecondaryFragmentNumber_Preserved()
	{
		var spectrum = MakeInternalIonEntry().ToLibrarySpectrum();
		var internalIon = GetInternalIon(spectrum);

		Assert.That(internalIon.NeutralTheoreticalProduct.SecondaryFragmentNumber, Is.EqualTo(6));
	}

	/// <summary>
	/// Product.IsInternalFragment must return true for the converted bIy ion.
	/// </summary>
	[Test]
	public void InternalIon_ToLibrarySpectrum_IsInternalFragment_True()
	{
		var spectrum = MakeInternalIonEntry().ToLibrarySpectrum();
		var internalIon = GetInternalIon(spectrum);

		Assert.That(internalIon.NeutralTheoreticalProduct.IsInternalFragment, Is.True);
	}

	/// <summary>
	/// The annotation on the converted Product must equal "bIy[3-6]", matching the
	/// mzLib Product.Annotation format for internal ions.
	/// </summary>
	[Test]
	public void InternalIon_ToLibrarySpectrum_Annotation_Correct()
	{
		var spectrum = MakeInternalIonEntry().ToLibrarySpectrum();
		var internalIon = GetInternalIon(spectrum);

		Assert.That(internalIon.NeutralTheoreticalProduct.Annotation, Is.EqualTo("bIy[3-6]"));
	}

	/// <summary>
	/// A terminal b3 ion must produce a Product with SecondaryProductType == null after
	/// conversion — confirming that the secondary-field pass-through does not corrupt
	/// terminal ions.
	/// </summary>
	[Test]
	public void TerminalIon_ToLibrarySpectrum_SecondaryProductType_Null()
	{
		var spectrum = MakeInternalIonEntry().ToLibrarySpectrum();
		var terminalIon = GetTerminalIon(spectrum);

		Assert.That(terminalIon.NeutralTheoreticalProduct.SecondaryProductType, Is.Null);
	}

	/// <summary>
	/// Product.IsInternalFragment must return false for a terminal b3 ion.
	/// </summary>
	[Test]
	public void TerminalIon_ToLibrarySpectrum_IsInternalFragment_False()
	{
		var spectrum = MakeInternalIonEntry().ToLibrarySpectrum();
		var terminalIon = GetTerminalIon(spectrum);

		Assert.That(terminalIon.NeutralTheoreticalProduct.IsInternalFragment, Is.False);
	}

	// ── FromLibrarySpectrum tests ─────────────────────────────────────────────

	/// <summary>
	/// FromLibrarySpectrum must recover SecondaryProductType from a Product that carries it.
	/// </summary>
	[Test]
	public void InternalIon_FromLibrarySpectrum_SecondaryProductType_Preserved()
	{
		// Build a LibrarySpectrum with a Product that has SecondaryProductType set
		var internalProduct = new Product(
			ProductType.b,
			FragmentationTerminus.None,
			neutralMass: 0.0,
			fragmentNumber: 3,
			residuePosition: 0,
			neutralLoss: 0.0,
			secondaryProductType: ProductType.y,
			secondaryFragmentNumber: 6);

		var spectrum = new LibrarySpectrum(
			sequence: "PEPTIDE",
			precursorMz: 400.2,
			chargeState: 2,
			peaks: new List<MatchedFragmentIon>
			{
				new MatchedFragmentIon(internalProduct, experMz: 415.22, experIntensity: 0.5, charge: 1)
			},
			rt: 25.0);

		var entry = MslLibraryEntry.FromLibrarySpectrum(spectrum);
		var mslInternal = entry.Fragments.Single(f => f.IsInternalFragment);

		Assert.That(mslInternal.SecondaryProductType, Is.EqualTo(ProductType.y));
	}

	/// <summary>
	/// FromLibrarySpectrum must recover SecondaryFragmentNumber from a Product that carries it.
	/// </summary>
	[Test]
	public void InternalIon_FromLibrarySpectrum_SecondaryFragmentNumber_Preserved()
	{
		var internalProduct = new Product(
			ProductType.b,
			FragmentationTerminus.None,
			neutralMass: 0.0,
			fragmentNumber: 3,
			residuePosition: 0,
			neutralLoss: 0.0,
			secondaryProductType: ProductType.y,
			secondaryFragmentNumber: 6);

		var spectrum = new LibrarySpectrum(
			sequence: "PEPTIDE",
			precursorMz: 400.2,
			chargeState: 2,
			peaks: new List<MatchedFragmentIon>
			{
				new MatchedFragmentIon(internalProduct, experMz: 415.22, experIntensity: 0.5, charge: 1)
			},
			rt: 25.0);

		var entry = MslLibraryEntry.FromLibrarySpectrum(spectrum);
		var mslInternal = entry.Fragments.Single(f => f.IsInternalFragment);

		Assert.That(mslInternal.SecondaryFragmentNumber, Is.EqualTo(6));
	}

	// ── Full round-trip tests ─────────────────────────────────────────────────

	/// <summary>
	/// MslLibraryEntry → ToLibrarySpectrum() → FromLibrarySpectrum() must preserve both
	/// secondary fields end-to-end through the LibrarySpectrum representation.
	/// </summary>
	[Test]
	public void InternalIon_FullRoundTrip_MslEntry_LibrarySpectrum_MslEntry()
	{
		var original = MakeInternalIonEntry();

		// Round-trip through LibrarySpectrum
		var spectrum = original.ToLibrarySpectrum();
		var recovered = MslLibraryEntry.FromLibrarySpectrum(spectrum);

		// Find the internal ion in the recovered entry
		var recoveredInternal = recovered.Fragments.Single(f => f.IsInternalFragment);

		Assert.That(recoveredInternal.SecondaryProductType, Is.EqualTo(ProductType.y),
			"SecondaryProductType must survive MslEntry → LibrarySpectrum → MslEntry round-trip");
		Assert.That(recoveredInternal.SecondaryFragmentNumber, Is.EqualTo(6),
			"SecondaryFragmentNumber must survive MslEntry → LibrarySpectrum → MslEntry round-trip");
	}

	/// <summary>
	/// Write an entry containing internal ions to a .msl file, read it back,
	/// call ToLibrarySpectrum(), and assert that both secondary fields are correct.
	/// This exercises the full path: MslWriter → MslReader → ToLibrarySpectrum.
	/// </summary>
	[Test]
	public void InternalIon_BinaryRoundTrip_WriteThenRead_SecondaryFieldsIntact()
	{
		string path = TempMsl("internal_ion_binary_roundtrip");
		var original = MakeInternalIonEntry();

		// Write then read back
		MslLibrary.Save(path, new List<MslLibraryEntry> { original });
		var lib = MslLibrary.Load(path);
		var readBack = lib.GetAllEntries().Single();

		// Convert to LibrarySpectrum and locate the internal ion
		var spectrum = readBack.ToLibrarySpectrum();
		var internalIon = GetInternalIon(spectrum);

		Assert.That(internalIon.NeutralTheoreticalProduct.SecondaryProductType,
			Is.EqualTo(ProductType.y),
			"SecondaryProductType must survive binary write → read → ToLibrarySpectrum");
		Assert.That(internalIon.NeutralTheoreticalProduct.SecondaryFragmentNumber,
			Is.EqualTo(6),
			"SecondaryFragmentNumber must survive binary write → read → ToLibrarySpectrum");
	}

	/// <summary>
	/// An entry with 3 terminal ions and 2 internal ions must convert all 5 correctly
	/// in a single ToLibrarySpectrum() call — no cross-contamination between ion types.
	/// </summary>
	[Test]
	public void MixedIons_TerminalAndInternal_AllConvertCorrectly()
	{
		var entry = new MslLibraryEntry
		{
			ModifiedSequence = "ACDEFGHIK",
			StrippedSequence = "ACDEFGHIK",
			PrecursorMz = 529.76,
			Charge = 2,
			Irt = 42.0,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			DissociationType = DissociationType.HCD,
			Fragments = new List<MslFragmentIon>
			{
                // Terminal ions
                new MslFragmentIon
				{
					ProductType = ProductType.b, FragmentNumber = 2,
					SecondaryProductType = null, SecondaryFragmentNumber = 0,
					ResiduePosition = 2, Charge = 1, Mz = 188.07f, Intensity = 1.0f
				},
				new MslFragmentIon
				{
					ProductType = ProductType.y, FragmentNumber = 3,
					SecondaryProductType = null, SecondaryFragmentNumber = 0,
					ResiduePosition = 6, Charge = 1, Mz = 364.20f, Intensity = 0.9f
				},
				new MslFragmentIon
				{
					ProductType = ProductType.b, FragmentNumber = 4,
					SecondaryProductType = null, SecondaryFragmentNumber = 0,
					ResiduePosition = 4, Charge = 1, Mz = 450.18f, Intensity = 0.7f
				},
                // Internal ions
                new MslFragmentIon
				{
					ProductType = ProductType.b, SecondaryProductType = ProductType.y,
					FragmentNumber = 2, SecondaryFragmentNumber = 5,
					ResiduePosition = 0, Charge = 1, Mz = 390.16f, Intensity = 0.4f
				},
				new MslFragmentIon
				{
					ProductType = ProductType.a, SecondaryProductType = ProductType.b,
					FragmentNumber = 3, SecondaryFragmentNumber = 7,
					ResiduePosition = 0, Charge = 1, Mz = 512.24f, Intensity = 0.3f
				}
			}
		};

		var spectrum = entry.ToLibrarySpectrum();

		Assert.That(spectrum.MatchedFragmentIons, Has.Count.EqualTo(5),
			"All 5 ions must be present after conversion");

		// Terminal ions: SecondaryProductType must be null, IsInternalFragment false
		var terminalIons = spectrum.MatchedFragmentIons
			.Where(i => !i.NeutralTheoreticalProduct.IsInternalFragment)
			.ToList();

		Assert.That(terminalIons, Has.Count.EqualTo(3),
			"Exactly 3 terminal ions expected");
		Assert.That(terminalIons.All(i => i.NeutralTheoreticalProduct.SecondaryProductType == null),
			Is.True, "All terminal ions must have SecondaryProductType == null");

		// Internal ions: IsInternalFragment true, SecondaryProductType set
		var internalIons = spectrum.MatchedFragmentIons
			.Where(i => i.NeutralTheoreticalProduct.IsInternalFragment)
			.ToList();

		Assert.That(internalIons, Has.Count.EqualTo(2),
			"Exactly 2 internal ions expected");

		var byIon = internalIons.Single(i =>
			i.NeutralTheoreticalProduct.ProductType == ProductType.b &&
			i.NeutralTheoreticalProduct.SecondaryProductType == ProductType.y);
		Assert.That(byIon.NeutralTheoreticalProduct.SecondaryFragmentNumber, Is.EqualTo(5));

		var abIon = internalIons.Single(i =>
			i.NeutralTheoreticalProduct.ProductType == ProductType.a &&
			i.NeutralTheoreticalProduct.SecondaryProductType == ProductType.b);
		Assert.That(abIon.NeutralTheoreticalProduct.SecondaryFragmentNumber, Is.EqualTo(7));
	}
}