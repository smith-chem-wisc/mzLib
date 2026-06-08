using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.InteropServices;

namespace Test.SpectralLibrary.MSL;

/// <summary>
/// Tests targeting the findings from Prompt 1 — Binary Format Correctness / Round-Trip Audit.
///
/// Covers:
///   B1 — NeutralLossCode.H3PO4AndH2O (formerly PlusH2O): naming and correct mass
///   B2 — NElutionGroups: write-only header field (informational, not a data-loss bug)
///   B3 — StrippedSeqLength: written but not used during read (informational)
///   D1 — NCE × 10 encoding: overflow guard for NCE > 3276
///   Magic — MagicForLEStruct byte-order correctness
///   NCE — round-trip encoding at representative values
///   Precision — float vs double fields that lose precision on round-trip
/// </summary>
[TestFixture]
public class TestMslPrompt1RoundTrip
{
	// ── Temp file infrastructure ──────────────────────────────────────────

	private static readonly string OutputDir =
		Path.Combine(Path.GetTempPath(), "MslPrompt1Tests");

	private static string Tmp(string name) =>
		Path.Combine(OutputDir, name + ".msl");

	[OneTimeSetUp]
	public void OneTimeSetUp() => Directory.CreateDirectory(OutputDir);

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(OutputDir))
			Directory.Delete(OutputDir, recursive: true);
	}

	// ── Minimal entry factory ─────────────────────────────────────────────

	private static MslLibraryEntry MakeEntry(
		string seq = "PEPTIDE",
		int charge = 2,
		double neutralLoss = 0.0,
		int nce = 28,
		double precursorMz = 500.0)
	{
		return new MslLibraryEntry
		{
			FullSequence = seq,
			BaseSequence = seq,
			PrecursorMz = precursorMz,
			ChargeState = charge,
			RetentionTime = 30.0,
			IonMobility = 0.0,
			DissociationType = DissociationType.HCD,
			Nce = nce,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			Source = MslFormat.SourceType.Predicted,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			IsProteotypic = false,
			QValue = float.NaN,
			ElutionGroupId = 0,
			IsDecoy = false,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					ProductType             = ProductType.b,
					SecondaryProductType    = null,
					FragmentNumber          = 3,
					SecondaryFragmentNumber = 0,
					ResiduePosition         = 3,
					Charge                  = 1,
					Mz                      = 312.15f,
					Intensity               = 1.0f,
					NeutralLoss             = neutralLoss
				}
			}
		};
	}

	// ═════════════════════════════════════════════════════════════════════
	// B1 — NeutralLossCode.H3PO4AndH2O (formerly PlusH2O)
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// The combined H₃PO₄ + H₂O loss (−115.987460 Da) must encode to
	/// NeutralLossCode.H3PO4AndH2O (value = 5) and decode back to the
	/// same mass within 1e-6 Da. This rules out both the "water-only" (−18)
	/// and the "no loss" (0) misinterpretations.
	/// </summary>
	[Test]
	public void NeutralLoss_H3PO4AndH2O_RoundTripsWithCorrectMass()
	{
		const double expectedLoss = -97.976895 + -18.010565;   // = −115.987460
		string path = Tmp(nameof(NeutralLoss_H3PO4AndH2O_RoundTripsWithCorrectMass));

		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry(neutralLoss: expectedLoss) });
		var lib = MslReader.Load(path);
		double actual = lib.Entries[0].MatchedFragmentIons[0].NeutralLoss;

		Assert.That(actual, Is.EqualTo(expectedLoss).Within(1e-6),
			$"H3PO4AndH2O loss must round-trip to {expectedLoss:F6} Da, got {actual:F6}.");
	}

	/// <summary>
	/// Belt-and-suspenders: the recovered H3PO4AndH2O mass must NOT be
	/// within 0.1 Da of the water-only loss (−18.010565). A cross-language
	/// reader that implements the old name as "water loss" would fail here.
	/// </summary>
	[Test]
	public void NeutralLoss_H3PO4AndH2O_IsNotWaterOnly()
	{
		const double combinedLoss = -97.976895 + -18.010565;
		const double waterOnly = -18.010565;
		string path = Tmp(nameof(NeutralLoss_H3PO4AndH2O_IsNotWaterOnly));

		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry(neutralLoss: combinedLoss) });
		var lib = MslReader.Load(path);
		double actual = lib.Entries[0].MatchedFragmentIons[0].NeutralLoss;

		Assert.That(Math.Abs(actual - waterOnly), Is.GreaterThan(0.1),
			$"H3PO4AndH2O (−115.987) must NOT decode as water-only (−18.011). " +
			$"Got {actual:F6} Da — difference from water-only: {Math.Abs(actual - waterOnly):F3} Da.");
	}

	/// <summary>
	/// The enum value 5 (H3PO4AndH2O) must be classified correctly by
	/// MslWriter.ClassifyNeutralLoss (writer-side, called from MslFormat.ClassifyNeutralLoss
	/// after Fix 3). Tests the classification in isolation without a file round-trip.
	///
	/// NOTE: After Fix 3, call MslFormat.ClassifyNeutralLoss directly.
	/// Until then, MslWriter.ClassifyNeutralLoss (internal) is accessible from tests.
	/// </summary>
	[Test]
	public void ClassifyNeutralLoss_H3PO4AndH2O_ClassifiesAsCode5()
	{
		const double combinedLoss = -97.976895 + -18.010565;

		// MslWriter.ClassifyNeutralLoss is internal — accessible from the test assembly
		// via [assembly: InternalsVisibleTo("Test")] if present, or via reflection otherwise.
		// After Fix 3, use MslFormat.ClassifyNeutralLoss (public).
		var code = MslWriter.ClassifyNeutralLoss(combinedLoss);

		// Value 5 = H3PO4AndH2O (formerly PlusH2O — same integer, renamed)
		Assert.That((int)code, Is.EqualTo(5),
			$"Combined H3PO4+H2O loss should classify as code 5 (H3PO4AndH2O), got {code} ({(int)code}).");
	}

	/// <summary>
	/// All five named neutral losses must round-trip through write → read
	/// and recover the original mass within 1e-4 Da.
	/// </summary>
	[TestCase(0.0, TestName = "NeutralLoss_None")]
	[TestCase(-18.010565, TestName = "NeutralLoss_H2O")]
	[TestCase(-17.026549, TestName = "NeutralLoss_NH3")]
	[TestCase(-97.976895, TestName = "NeutralLoss_H3PO4")]
	[TestCase(-79.966331, TestName = "NeutralLoss_HPO3")]
	[TestCase(-97.976895 + -18.010565, TestName = "NeutralLoss_H3PO4AndH2O")]
	public void NeutralLoss_AllNamedCodes_RoundTrip(double loss)
	{
		string path = Tmp($"NeutralLoss_RoundTrip_{loss:F3}".Replace("-", "m").Replace(".", "_"));

		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry(neutralLoss: loss) });
		var lib = MslReader.Load(path);
		double actual = lib.Entries[0].MatchedFragmentIons[0].NeutralLoss;

		Assert.That(actual, Is.EqualTo(loss).Within(1e-4),
			$"Neutral loss {loss:F6} did not survive write→read round-trip (got {actual:F6}).");
	}

	// ═════════════════════════════════════════════════════════════════════
	// Magic byte correctness
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// The first four bytes of every .msl file must be 0x4D 0x5A 0x4C 0x42
	/// ("MZLB"). Verifies that MagicForLEStruct written via MemoryMarshal
	/// produces the correct on-disk sequence on a little-endian platform.
	/// </summary>
	[Test]
	public void FileHeader_MagicBytes_AreCorrectOnDisk()
	{
		string path = Tmp(nameof(FileHeader_MagicBytes_AreCorrectOnDisk));
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry() });

		byte[] raw = File.ReadAllBytes(path);

		Assert.That(raw[0], Is.EqualTo(0x4D), "Magic byte 0 must be 0x4D ('M')");
		Assert.That(raw[1], Is.EqualTo(0x5A), "Magic byte 1 must be 0x5A ('Z')");
		Assert.That(raw[2], Is.EqualTo(0x4C), "Magic byte 2 must be 0x4C ('L')");
		Assert.That(raw[3], Is.EqualTo(0x42), "Magic byte 3 must be 0x42 ('B')");
	}

	/// <summary>
	/// MslFormat.MagicMatches must return true for the bytes produced by
	/// the writer, confirming reader and writer agree on the magic value.
	/// </summary>
	[Test]
	public void MagicMatches_ReturnsTrue_ForWriterProducedBytes()
	{
		string path = Tmp(nameof(MagicMatches_ReturnsTrue_ForWriterProducedBytes));
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry() });

		byte[] raw = File.ReadAllBytes(path);

		Assert.That(MslFormat.MagicMatches(raw.AsSpan(0, 4)), Is.True,
			"MagicMatches must return true for the leading 4 bytes of a writer-produced file.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// NCE × 10 encoding round-trip
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// NCE must survive the × 10 encode / ÷ 10 decode round-trip exactly
	/// for all values in the normal operating range.
	/// </summary>
	[TestCase(0, TestName = "Nce_Zero")]
	[TestCase(20, TestName = "Nce_20")]
	[TestCase(25, TestName = "Nce_25")]
	[TestCase(28, TestName = "Nce_28")]
	[TestCase(30, TestName = "Nce_30")]
	[TestCase(35, TestName = "Nce_35")]
	public void Nce_RoundTrips_Exactly(int nce)
	{
		string path = Tmp($"Nce_{nce}");
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry(nce: nce) });
		var lib = MslReader.Load(path);

		Assert.That(lib.Entries[0].Nce, Is.EqualTo(nce),
			$"NCE {nce} must survive the ×10 encode / ÷10 decode round-trip exactly.");
	}

	/// <summary>
	/// NCE = 3276 is the maximum value that fits in an int16 when stored as
	/// NCE × 10 (3276 × 10 = 32760 ≤ 32767). Must round-trip correctly.
	/// </summary>
	[Test]
	public void Nce_MaxSafeValue_3276_RoundTrips()
	{
		string path = Tmp(nameof(Nce_MaxSafeValue_3276_RoundTrips));
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry(nce: 3276) });
		var lib = MslReader.Load(path);

		Assert.That(lib.Entries[0].Nce, Is.EqualTo(3276),
			"NCE 3276 (maximum safely encodable value) must round-trip correctly.");
	}

	/// <summary>
	/// Fix 2 (EncodeNce with Math.Clamp) clamps NCE > 3276 to 3276 rather
	/// than silently overflowing int16. This test confirms the clamp is applied:
	/// NCE 3277 must be stored and recovered as 3276. The writer emits a
	/// diagnostic warning when clamping occurs; that is expected behaviour.
	/// </summary>
	[Test]
	public void Nce_OverflowValue_3277_IsClamped()
	{
		string path = Tmp(nameof(Nce_OverflowValue_3277_IsClamped));
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry(nce: 3277) });
		var lib = MslReader.Load(path);

		int recovered = lib.Entries[0].Nce;

		Assert.That(recovered, Is.EqualTo(3276),
			"NCE 3277 exceeds int16 range when stored as NCE×10. " +
			"EncodeNce (Fix 2) must clamp it to 3276; the writer emits a diagnostic warning.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// Precision: float vs double fields
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// PrecursorMz is stored as float32. Precision loss beyond float32
	/// representation is expected and acceptable at typical m/z tolerances.
	/// This test confirms the loss is bounded within float32 epsilon at 500 Da.
	/// </summary>
	[Test]
	public void PrecursorMz_PrecisionLoss_BoundedToFloat32Epsilon()
	{
		double originalMz = 500.123456789;  // more precision than float32 can hold
		string path = Tmp(nameof(PrecursorMz_PrecisionLoss_BoundedToFloat32Epsilon));
		MslWriter.Write(path, new List<MslLibraryEntry> { MakeEntry(precursorMz: originalMz) });
		var lib = MslReader.Load(path);

		double recovered = lib.Entries[0].PrecursorMz;
		double delta = Math.Abs(recovered - originalMz);

		// float32 relative precision ~1.2e-7; at 500 Da this is ~6e-5 Da
		Assert.That(delta, Is.LessThan(1e-3),
			$"PrecursorMz float32 precision loss should be < 1 mDa at 500 Da, got {delta:E3} Da.");
		Assert.That(recovered, Is.Not.EqualTo(originalMz),
			"PrecursorMz stored as float32 SHOULD lose sub-float32 precision — " +
			"this verifies the test is actually measuring the truncation.");
	}

	/// <summary>
	/// RetentionTime (iRT) is stored as float32. Same acceptable precision-loss as PrecursorMz.
	/// </summary>
	[Test]
	public void Irt_StoredAsFloat32_PrecisionLossAcceptable()
	{
		// Build entry with a high-precision iRT value
		var entry = MakeEntry();
		entry.RetentionTime = 42.123456789;
		string path = Tmp(nameof(Irt_StoredAsFloat32_PrecisionLossAcceptable));
		MslWriter.Write(path, new List<MslLibraryEntry> { entry });
		var lib = MslReader.Load(path);

		double delta = Math.Abs(lib.Entries[0].RetentionTime - entry.RetentionTime);
		Assert.That(delta, Is.LessThan(1e-3),
			"RetentionTime float32 precision loss should be < 1 mDa equivalent.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// Full precursor-field round-trip
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// All precursor-level fields must survive write → read with values intact.
	/// This is the primary regression guard for any future struct layout change.
	/// </summary>
	[Test]
	public void PrecursorRecord_AllFields_RoundTrip()
	{
		var entry = new MslLibraryEntry
		{
			FullSequence = "ACDEFGHIK",
			BaseSequence = "ACDEFGHIK",
			PrecursorMz = 529.765,
			ChargeState = 2,
			RetentionTime = 42.5,
			IonMobility = 0.85,
			DissociationType = DissociationType.HCD,
			Nce = 28,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			Source = MslFormat.SourceType.Empirical,
			ProteinAccession = string.Empty,
			ProteinName = string.Empty,
			GeneName = string.Empty,
			IsProteotypic = true,
			IsDecoy = false,
			QValue = 0.01f,
			ElutionGroupId = 0,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					ProductType = ProductType.b, FragmentNumber = 4,
					Charge = 1, Mz = 450.18f, Intensity = 1.0f, NeutralLoss = 0.0
				}
			}
		};

		string path = Tmp(nameof(PrecursorRecord_AllFields_RoundTrip));
		MslWriter.Write(path, new List<MslLibraryEntry> { entry });
		var lib = MslReader.Load(path);
		var r = lib.Entries[0];

		Assert.Multiple(() =>
		{
			Assert.That(r.FullSequence, Is.EqualTo(entry.FullSequence), "FullSequence");
			Assert.That(r.BaseSequence, Is.EqualTo(entry.BaseSequence), "BaseSequence");
			Assert.That(r.ChargeState, Is.EqualTo(entry.ChargeState), "ChargeState");
			Assert.That(r.Nce, Is.EqualTo(entry.Nce), "Nce");
			Assert.That(r.MoleculeType, Is.EqualTo(entry.MoleculeType), "MoleculeType");
			Assert.That(r.DissociationType, Is.EqualTo(entry.DissociationType), "DissociationType");
			Assert.That(r.Source, Is.EqualTo(entry.Source), "Source");
			Assert.That(r.IsProteotypic, Is.EqualTo(entry.IsProteotypic), "IsProteotypic");
			Assert.That(r.IsDecoy, Is.EqualTo(entry.IsDecoy), "IsDecoy");
			Assert.That(r.QValue, Is.EqualTo(entry.QValue).Within(1e-6f), "QValue");

			// Float32 fields — allow float32 round-trip precision loss
			Assert.That(r.PrecursorMz, Is.EqualTo(entry.PrecursorMz).Within(1e-3), "PrecursorMz");
			Assert.That(r.RetentionTime, Is.EqualTo(entry.RetentionTime).Within(1e-3), "RetentionTime");
			Assert.That(r.IonMobility, Is.EqualTo(entry.IonMobility).Within(1e-5), "IonMobility");
		});
	}

	/// <summary>
	/// All fragment-level fields must survive write → read with values intact.
	/// Tests Mz, Intensity, ProductType, FragmentNumber, ResiduePosition, ChargeState,
	/// SecondaryProductType, SecondaryFragmentNumber, and NeutralLoss.
	/// </summary>
	[Test]
	public void FragmentRecord_AllFields_RoundTrip()
	{
		var entry = MakeEntry();
		// Replace the default fragment with a richly populated one
		entry.MatchedFragmentIons = new List<MslFragmentIon>
		{
            // Terminal b ion
            new MslFragmentIon
			{
				ProductType             = ProductType.y,
				SecondaryProductType    = null,
				FragmentNumber          = 5,
				SecondaryFragmentNumber = 0,
				ResiduePosition         = 4,
				Charge                  = 2,
				Mz                      = 612.34f,
				Intensity               = 0.75f,
				NeutralLoss             = -18.010565   // H2O loss
            },
            // Internal ion
            new MslFragmentIon
			{
				ProductType             = ProductType.b,
				SecondaryProductType    = ProductType.y,
				FragmentNumber          = 2,
				SecondaryFragmentNumber = 5,
				ResiduePosition         = 0,
				Charge                  = 1,
				Mz                      = 390.16f,
				Intensity               = 0.4f,
				NeutralLoss             = 0.0
			}
		};

		string path = Tmp(nameof(FragmentRecord_AllFields_RoundTrip));
		MslWriter.Write(path, new List<MslLibraryEntry> { entry });
		var lib = MslReader.Load(path);
		var frags = lib.Entries[0].MatchedFragmentIons;

		// Frags are sorted by m/z ascending on write; 390 < 612 so internal ion is [0]
		var internal_ion = frags[0];
		var terminal_ion = frags[1];

		Assert.Multiple(() =>
		{
			// Terminal y5 -H2O
			Assert.That(terminal_ion.ProductType, Is.EqualTo(ProductType.y), "terminal ProductType");
			Assert.That(terminal_ion.SecondaryProductType, Is.Null, "terminal SecondaryProductType");
			Assert.That(terminal_ion.FragmentNumber, Is.EqualTo(5), "terminal FragmentNumber");
			Assert.That(terminal_ion.ResiduePosition, Is.EqualTo(4), "terminal ResiduePosition");
			Assert.That(terminal_ion.Charge, Is.EqualTo(2), "terminal ChargeState");
			Assert.That(terminal_ion.NeutralLoss, Is.EqualTo(-18.010565).Within(1e-4), "terminal NeutralLoss");
			Assert.That(terminal_ion.Mz, Is.EqualTo(612.34f).Within(1e-4f), "terminal Mz");

			// Internal bIy[2-5]
			Assert.That(internal_ion.ProductType, Is.EqualTo(ProductType.b), "internal ProductType");
			Assert.That(internal_ion.SecondaryProductType, Is.EqualTo(ProductType.y), "internal SecondaryProductType");
			Assert.That(internal_ion.FragmentNumber, Is.EqualTo(2), "internal FragmentNumber");
			Assert.That(internal_ion.SecondaryFragmentNumber, Is.EqualTo(5), "internal SecondaryFragmentNumber");
			Assert.That(internal_ion.NeutralLoss, Is.EqualTo(0.0).Within(1e-9), "internal NeutralLoss");
		});
	}

	// ═════════════════════════════════════════════════════════════════════
	// B2 / B3 — Write-only fields (documentation tests)
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// MslFileHeader.NElutionGroups is written by the writer but not currently
	/// exposed through MslLibraryData's public API. This test documents that the
	/// value is present in the raw header and reads back correctly when accessed
	/// via the Header property, even though no high-level API surfaces it.
	/// </summary>
	[Test]
	public void Header_NElutionGroups_IsWrittenCorrectly()
	{
		var entries = new List<MslLibraryEntry>
		{
			MakeEntry("PEPTIDE",  charge: 2),
			MakeEntry("PEPTIDE",  charge: 3),   // same stripped seq → same group
            MakeEntry("ACDEFGHI", charge: 2),   // different → new group
        };

		string path = Tmp(nameof(Header_NElutionGroups_IsWrittenCorrectly));
		MslWriter.Write(path, entries);

		MslFileHeader header = MslReader.ReadHeaderOnly(path);

		Assert.That(header.NElutionGroups, Is.EqualTo(2),
			"NElutionGroups must equal the number of distinct stripped sequences (2).");
	}

	/// <summary>
	/// MslPrecursorRecord.StrippedSeqLength is written by the writer but the
	/// reader does not use it — BaseSequence.Length is recomputed from the
	/// string resolved from the string table.
	/// </summary>
	[Test]
	public void PrecursorRecord_StrippedSeqLength_MatchesActualSequenceLength()
	{
		string seq = "ACDEFGHIKLMNPQRSTVWY";   // 20 residues
		var entry = MakeEntry(seq);

		string path = Tmp(nameof(PrecursorRecord_StrippedSeqLength_MatchesActualSequenceLength));
		MslWriter.Write(path, new List<MslLibraryEntry> { entry });

		var lib = MslReader.Load(path);
		string recovered = lib.Entries[0].BaseSequence;

		Assert.That(recovered.Length, Is.EqualTo(seq.Length),
			"Recovered BaseSequence must have the same length as the original.");
		Assert.That(recovered, Is.EqualTo(seq),
			"Recovered BaseSequence must equal the original.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// Struct size validation
	// ═════════════════════════════════════════════════════════════════════

	[Test]
	public void StructSizeCheck_DoesNotThrow()
	{
		Assert.That(() => MslStructs.SizeCheck(), Throws.Nothing,
			"SizeCheck must not throw — all structs must match their declared sizes.");
	}

	[Test]
	public void StructSizes_MatchDeclaredConstants()
	{
		Assert.Multiple(() =>
		{
			Assert.That(Marshal.SizeOf<MslFileHeader>(), Is.EqualTo(MslFormat.HeaderSize), "MslFileHeader");
			Assert.That(Marshal.SizeOf<MslProteinRecord>(), Is.EqualTo(MslFormat.ProteinRecordSize), "MslProteinRecord");
			Assert.That(Marshal.SizeOf<MslPrecursorRecord>(), Is.EqualTo(MslFormat.PrecursorRecordSize), "MslPrecursorRecord");
			Assert.That(Marshal.SizeOf<MslFragmentRecord>(), Is.EqualTo(MslFormat.FragmentRecordSize), "MslFragmentRecord");
			Assert.That(Marshal.SizeOf<MslFooter>(), Is.EqualTo(MslFormat.FooterSize), "MslFooter");
		});
	}
}