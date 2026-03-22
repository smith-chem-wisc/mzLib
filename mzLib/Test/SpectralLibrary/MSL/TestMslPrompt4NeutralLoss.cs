using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.SpectralLibrary.MSL;

/// <summary>
/// Tests targeting the findings from Prompt 4 — Neutral Loss Encoding:
/// Writer/Reader Consistency and Custom Loss Handling.
///
/// Covers:
///   NL1 — All five named losses round-trip with correct mass (belt-and-suspenders)
///   NL2 — NeutralLossCode.H3PO4AndH2O is correctly renamed from PlusH2O (Fix 1)
///   NL3 — WriteStreaming throws for custom losses (documents current behaviour);
///          after Fix 7 this becomes a round-trip success test
///   NL4 — Write() supports custom losses (non-streaming path)
///   NL5 — Flags byte neutral loss field: encoding/decoding consistency
///   NL6 — No two named losses are within 0.01 Da of each other (no mis-classification)
///   NL7 — 3-bit field is nearly full: adding a 7th named loss would exhaust headroom
/// </summary>
[TestFixture]
public class TestMslPrompt4NeutralLoss
{
	private string _tempDir = null!;

	[OneTimeSetUp]
	public void OneTimeSetUp()
	{
		_tempDir = Path.Combine(Path.GetTempPath(),
			$"TestMslPrompt4_{Guid.NewGuid():N}");
		Directory.CreateDirectory(_tempDir);
	}

	[OneTimeTearDown]
	public void OneTimeTearDown()
	{
		if (Directory.Exists(_tempDir))
			Directory.Delete(_tempDir, recursive: true);
	}

	private string TempPath(string name) =>
		Path.Combine(_tempDir, name + ".msl");

	// ── Shared entry builder ──────────────────────────────────────────────

	private static MslLibraryEntry EntryWithLoss(double neutralLoss, string seq = "PEPTIDE") =>
		new MslLibraryEntry
		{
			FullSequence = seq,
			BaseSequence = seq,
			PrecursorMz = 500.0,
			ChargeState = 2,
			RetentionTime = 30.0,
			DissociationType = DissociationType.HCD,
			Nce = 28,
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
					ProductType     = ProductType.b,
					FragmentNumber  = 3,
					Charge          = 1,
					Mz              = 312.15f,
					Intensity       = 1.0f,
					NeutralLoss     = neutralLoss,
					ResiduePosition = 3
				}
			}
		};

	// ═════════════════════════════════════════════════════════════════════
	// NL1 — All five named losses round-trip with correct mass
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Each of the five named neutral losses must survive write → read and
	/// recover to within 1e-6 Da. Belt-and-suspenders complement to the
	/// all-codes parameterized test in TestMslPrompt1RoundTrip.
	/// </summary>
	[TestCase(0.0, "None", TestName = "NamedLoss_None")]
	[TestCase(-18.010565, "H2O", TestName = "NamedLoss_H2O")]
	[TestCase(-17.026549, "NH3", TestName = "NamedLoss_NH3")]
	[TestCase(-97.976895, "H3PO4", TestName = "NamedLoss_H3PO4")]
	[TestCase(-79.966331, "HPO3", TestName = "NamedLoss_HPO3")]
	[TestCase(-97.976895 + -18.010565, "H3PO4AndH2O", TestName = "NamedLoss_H3PO4AndH2O")]
	public void NamedLoss_RoundTrips_WithCorrectMass(double loss, string codeName)
	{
		string path = TempPath($"named_{codeName}");
		MslWriter.Write(path, new[] { EntryWithLoss(loss) });
		var lib = MslReader.Load(path);

		double actual = lib.Entries[0].MatchedFragmentIons[0].NeutralLoss;

		Assert.That(actual, Is.EqualTo(loss).Within(1e-6),
			$"NeutralLossCode.{codeName} ({loss:F6} Da) must survive write→read round-trip.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// NL2 — H3PO4AndH2O is the renamed PlusH2O (Fix 1 verification)
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// NeutralLossCode.H3PO4AndH2O (value = 5, renamed from PlusH2O in Fix 1)
	/// must represent −115.987460 Da — NOT −18.010565 Da (water-only).
	/// </summary>
	[Test]
	public void H3PO4AndH2O_IsNotWaterOnly_MassIsCorrect()
	{
		const double combinedLoss = -97.976895 + -18.010565;   // = −115.987460
		const double waterOnly = -18.010565;

		string path = TempPath("H3PO4AndH2O_not_water");
		MslWriter.Write(path, new[] { EntryWithLoss(combinedLoss) });
		var lib = MslReader.Load(path);

		double actual = lib.Entries[0].MatchedFragmentIons[0].NeutralLoss;

		Assert.That(actual, Is.EqualTo(combinedLoss).Within(1e-6),
			$"H3PO4AndH2O must decode to {combinedLoss:F6} Da.");
		Assert.That(Math.Abs(actual - waterOnly), Is.GreaterThan(0.1),
			$"H3PO4AndH2O must NOT be within 0.1 Da of water-only loss ({waterOnly:F6} Da).");
	}

	/// <summary>
	/// The integer value of NeutralLossCode.H3PO4AndH2O must be 5 — the same as
	/// the old PlusH2O — to maintain binary compatibility with existing .msl files.
	/// </summary>
	[Test]
	public void H3PO4AndH2O_EnumValue_IsStillFive()
	{
		Assert.That((int)MslFormat.NeutralLossCode.H3PO4AndH2O, Is.EqualTo(5),
			"H3PO4AndH2O must retain integer value 5 for binary compatibility with " +
			"files written before the rename.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// NL3 — WriteStreaming custom loss behaviour (Fix 7 applied)
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// After Fix 7: WriteStreaming supports custom neutral losses without throwing.
	/// The round-trip mass must exactly equal the written value.
	/// </summary>
	[Test]
	public void WriteStreaming_CustomLoss_RoundTrips_Correctly()
	{
		const double customLoss = -203.0794;   // glycan loss (HexNAc)
		string path = TempPath("streaming_custom_roundtrip");

		var entries = new List<MslLibraryEntry> { EntryWithLoss(customLoss) };

		Assert.That(() => MslWriter.WriteStreaming(path, entries),
			Throws.Nothing,
			"After Fix 7, WriteStreaming must not throw for custom neutral losses.");

		var lib = MslReader.Load(path);
		double actual = lib.Entries[0].MatchedFragmentIons[0].NeutralLoss;

		Assert.That(actual, Is.EqualTo(customLoss).Within(1e-9),
			$"Custom loss {customLoss:F6} Da must survive WriteStreaming→Read round-trip exactly.");
	}

	/// <summary>
	/// Write() (non-streaming path) supports custom neutral losses via the
	/// extended annotation table. Both paths must produce the same result.
	/// </summary>
	[Test]
	public void Write_CustomLoss_RoundTrips_Correctly()
	{
		const double customLoss = -203.0794;
		string path = TempPath("write_custom_roundtrip");

		Assert.That(() => MslWriter.Write(path, new[] { EntryWithLoss(customLoss) }),
			Throws.Nothing,
			"Write() must support custom neutral losses without throwing.");

		var lib = MslReader.Load(path);
		double actual = lib.Entries[0].MatchedFragmentIons[0].NeutralLoss;

		Assert.That(actual, Is.EqualTo(customLoss).Within(1e-9),
			$"Custom loss {customLoss:F6} Da must survive Write→Read round-trip exactly.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// NL4 — Flags byte: correct encoding and decoding for all codes
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// EncodeFragmentFlags / DecodeFragmentFlags must round-trip every defined
	/// NeutralLossCode value through the 3-bit field in the flags byte.
	/// This directly tests the bit-packing layer independently of file I/O.
	/// </summary>
	[TestCase(MslFormat.NeutralLossCode.None, TestName = "FlagsByte_None")]
	[TestCase(MslFormat.NeutralLossCode.H2O, TestName = "FlagsByte_H2O")]
	[TestCase(MslFormat.NeutralLossCode.NH3, TestName = "FlagsByte_NH3")]
	[TestCase(MslFormat.NeutralLossCode.H3PO4, TestName = "FlagsByte_H3PO4")]
	[TestCase(MslFormat.NeutralLossCode.HPO3, TestName = "FlagsByte_HPO3")]
	[TestCase(MslFormat.NeutralLossCode.H3PO4AndH2O, TestName = "FlagsByte_H3PO4AndH2O")]
	[TestCase(MslFormat.NeutralLossCode.Custom, TestName = "FlagsByte_Custom")]
	public void FlagsByte_NeutralLossCode_RoundTrips(MslFormat.NeutralLossCode code)
	{
		byte encoded = MslFormat.EncodeFragmentFlags(
			isInternal: false, isDiagnostic: false,
			lossCode: code, excludeFromQuant: false);

		var (_, _, decodedCode, _) = MslFormat.DecodeFragmentFlags(encoded);

		Assert.That(decodedCode, Is.EqualTo(code),
			$"NeutralLossCode.{code} must survive EncodeFragmentFlags → DecodeFragmentFlags.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// NL5 — No mis-classification risk between named losses
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Every pair of named losses must be more than 0.01 Da apart.
	/// This confirms the 0.01 Da tolerance cannot cause one named loss to
	/// be mis-classified as another.
	/// </summary>
	[Test]
	public void NamedLosses_AllPairwiseDistances_ExceedTolerance()
	{
		var losses = new (string name, double mass)[]
		{
			("None",        0.0),
			("H2O",        -18.010565),
			("NH3",        -17.026549),
			("H3PO4",      -97.976895),
			("HPO3",       -79.966331),
			("H3PO4AndH2O",-115.987460),
		};

		const double tolerance = 0.01;

		for (int i = 0; i < losses.Length; i++)
			for (int j = i + 1; j < losses.Length; j++)
			{
				double dist = Math.Abs(losses[i].mass - losses[j].mass);
				Assert.That(dist, Is.GreaterThan(tolerance),
					$"Named losses {losses[i].name} ({losses[i].mass:F6}) and " +
					$"{losses[j].name} ({losses[j].mass:F6}) are only {dist:F6} Da apart — " +
					$"within the 0.01 Da tolerance. Mis-classification is possible.");
			}
	}

	// ═════════════════════════════════════════════════════════════════════
	// NL6 — 3-bit field headroom
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// The 3-bit neutral loss field supports values 0–7. Six values are currently
	/// assigned (None=0 through Custom=6); value 7 is the only remaining slot.
	/// This test documents the constraint so that any future addition of a named
	/// loss is visible as a failing assertion.
	/// </summary>
	[Test]
	public void NeutralLossCode_ThreeBitField_HasOneRemainingSlot()
	{
		var defined = Enum.GetValues<MslFormat.NeutralLossCode>();

		// 3-bit field: max representable value = 7; values 0-7 = 8 slots
		const int maxSlots = 8;   // 2^3

		Assert.That(defined.Length, Is.EqualTo(7),
			$"7 NeutralLossCode values are defined (None through Custom). " +
			$"The 3-bit field supports {maxSlots} values total, leaving 1 reserved slot (value 7). " +
			"Adding another named code requires widening the field to 4 bits — a breaking format change.");

		// Confirm no value exceeds 6 (value 7 is reserved)
		foreach (var code in defined)
		{
			Assert.That((int)code, Is.LessThanOrEqualTo(6),
				$"NeutralLossCode.{code} has value {(int)code}, which must not exceed 6 " +
				"(value 7 is the only reserved slot in the 3-bit field).");
		}
	}

	// ═════════════════════════════════════════════════════════════════════
	// NL7 — Extended annotation table: FileFlagHasExtAnnotations set correctly
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Writing a custom loss must set FileFlagHasExtAnnotations in the file header.
	/// Writing only named losses must NOT set it.
	/// </summary>
	[Test]
	public void CustomLoss_FileFlagHasExtAnnotations_SetWhenCustomPresent()
	{
		string path = TempPath("flag_set_custom");
		MslWriter.Write(path, new[] { EntryWithLoss(-203.0794) });

		MslFileHeader header = MslReader.ReadHeaderOnly(path);

		Assert.That((header.FileFlags & MslFormat.FileFlagHasExtAnnotations), Is.Not.EqualTo(0),
			"FileFlagHasExtAnnotations must be set when a custom neutral loss is written.");
	}

	[Test]
	public void NamedLoss_FileFlagHasExtAnnotations_ClearWhenNoCustom()
	{
		string path = TempPath("flag_clear_named");
		MslWriter.Write(path, new[] { EntryWithLoss(-18.010565) });   // H2O — named

		MslFileHeader header = MslReader.ReadHeaderOnly(path);

		Assert.That((header.FileFlags & MslFormat.FileFlagHasExtAnnotations), Is.EqualTo(0),
			"FileFlagHasExtAnnotations must NOT be set when only named neutral losses are written.");
	}

	// ═════════════════════════════════════════════════════════════════════
	// NL8 — ResiduePosition is preserved for non-custom losses
	// ═════════════════════════════════════════════════════════════════════

	/// <summary>
	/// For named neutral losses, ResiduePosition must round-trip normally.
	/// Only custom-loss fragments repurpose ResiduePosition as ExtAnnotationIdx.
	/// </summary>
	[Test]
	public void NamedLoss_ResiduePosition_PreservedInRoundTrip()
	{
		var entry = EntryWithLoss(-18.010565);   // H2O — named
		entry.MatchedFragmentIons[0].ResiduePosition = 7;

		string path = TempPath("residue_pos_named");
		MslWriter.Write(path, new[] { entry });
		var lib = MslReader.Load(path);

		Assert.That(lib.Entries[0].MatchedFragmentIons[0].ResiduePosition, Is.EqualTo(7),
			"ResiduePosition must be preserved for named neutral losses.");
	}

	/// <summary>
	/// For custom neutral losses, ResiduePosition is repurposed as ExtAnnotationIdx
	/// in the binary format and is recovered as 0 on read (documented trade-off).
	/// </summary>
	[Test]
	public void CustomLoss_ResiduePosition_IsZeroAfterRoundTrip()
	{
		var entry = EntryWithLoss(-203.0794);    // custom
		entry.MatchedFragmentIons[0].ResiduePosition = 7;

		string path = TempPath("residue_pos_custom");
		MslWriter.Write(path, new[] { entry });
		var lib = MslReader.Load(path);

		Assert.That(lib.Entries[0].MatchedFragmentIons[0].ResiduePosition, Is.EqualTo(0),
			"ResiduePosition is repurposed as ExtAnnotationIdx for custom-loss fragments " +
			"and is recovered as 0 on read. This is a documented trade-off.");
	}
}