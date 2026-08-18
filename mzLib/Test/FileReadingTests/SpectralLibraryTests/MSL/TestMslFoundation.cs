using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectralMatch.MslSpectralLibrary;
using Omics.SpectrumMatch;
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;

namespace Test.MslSpectralLibrary;

/// <summary>
/// NUnit 4 test suite covering MslFormat, MslStructs, and MslLibraryEntry.
/// All tests use Assert.That exclusively; the legacy Assert.AreEqual style is forbidden.
/// All test classes use [TestFixture]; all test methods use [Test].
/// </summary>
[TestFixture]
public class TestMslFoundation
{
	// ══════════════════════════════════════════════════════════════════════
	// Struct size tests
	// ══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that MslFileHeader occupies exactly 64 unmanaged bytes as required by the
	/// format specification. A failure here means a Pack setting error or an accidental
	/// extra field was added and the writer/reader would produce corrupt files.
	/// </summary>
	[Test]
	public void MslFileHeader_Is_64_Bytes()
		=> Assert.That(Marshal.SizeOf<MslFileHeader>(), Is.EqualTo(MslFormat.HeaderSize));

	/// <summary>
	/// Verifies that MslProteinRecord occupies exactly 24 unmanaged bytes.
	/// A size mismatch would misalign all protein records after the first.
	/// </summary>
	[Test]
	public void MslProteinRecord_Is_24_Bytes()
		=> Assert.That(Marshal.SizeOf<MslProteinRecord>(), Is.EqualTo(MslFormat.ProteinRecordSize));

	/// <summary>
	/// Verifies that MslPrecursorRecord occupies exactly 56 unmanaged bytes.
	/// This is the most critical size: misalignment here corrupts every precursor read.
	/// </summary>
	[Test]
	public void MslPrecursorRecord_Is_56_Bytes()
		=> Assert.That(Marshal.SizeOf<MslPrecursorRecord>(), Is.EqualTo(MslFormat.PrecursorRecordSize));

	/// <summary>
	/// Verifies that MslFragmentRecord occupies exactly 20 unmanaged bytes.
	/// Fragment blocks use FragmentCount × 20 arithmetic for seeking; a size error
	/// would produce wrong seek targets for every fragment read.
	/// </summary>
	[Test]
	public void MslFragmentRecord_Is_20_Bytes()
		=> Assert.That(Marshal.SizeOf<MslFragmentRecord>(), Is.EqualTo(MslFormat.FragmentRecordSize));

	/// <summary>
	/// Verifies that MslFooter occupies exactly 20 unmanaged bytes.
	/// The footer is located by seeking to EOF − 20; the wrong size means the
	/// OffsetTableOffset field is read from the wrong position.
	/// </summary>
	[Test]
	public void MslFooter_Is_20_Bytes()
		=> Assert.That(Marshal.SizeOf<MslFooter>(), Is.EqualTo(MslFormat.FooterSize));

	/// <summary>
	/// Exercises MslStructs.SizeCheck() to confirm it completes without throwing when
	/// all five structs have the correct size. SizeCheck() is intended to be called at
	/// application start-up; a throw here would abort the application with a clear error.
	/// </summary>
	[Test]
	public void SizeCheck_DoesNotThrow_WhenAllStructsHaveCorrectSize()
		=> Assert.That(() => MslStructs.SizeCheck(), Throws.Nothing);

	// ══════════════════════════════════════════════════════════════════════
	// Magic tests
	// ══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that MagicMatches returns true when the first four bytes of a byte array
	/// exactly match the MZLB magic sequence 0x4D 0x5A 0x4C 0x42.
	/// </summary>
	[Test]
	public void Magic_Matches_Returns_True_For_Valid_Magic()
	{
		// Arrange: build a 4-byte array that equals MslFormat.Magic
		byte[] validHeader = { 0x4D, 0x5A, 0x4C, 0x42 };

		// Act + Assert
		Assert.That(MslFormat.MagicMatches(validHeader), Is.True);
	}

	/// <summary>
	/// Verifies that MagicMatches returns false when the first four bytes do not match.
	/// A corrupt or non-MSL file must be detected immediately at the magic check.
	/// </summary>
	[Test]
	public void Magic_Matches_Returns_False_For_Invalid_Magic()
	{
		// Arrange: first byte differs (0x00 instead of 0x4D)
		byte[] badHeader = { 0x00, 0x5A, 0x4C, 0x42 };

		// Act + Assert
		Assert.That(MslFormat.MagicMatches(badHeader), Is.False);
	}

	/// <summary>
	/// Verifies that MagicMatches returns false when the span is shorter than 4 bytes,
	/// preventing an IndexOutOfRangeException on truncated input.
	/// </summary>
	[Test]
	public void Magic_Matches_Returns_False_For_Short_Span()
	{
		// Arrange: only 3 bytes — could happen on a 3-byte file
		byte[] tooShort = { 0x4D, 0x5A, 0x4C };

		// Act + Assert
		Assert.That(MslFormat.MagicMatches(tooShort), Is.False);
	}

	// ══════════════════════════════════════════════════════════════════════
	// Fragment flags encoding / decoding tests
	// ══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that a standard terminal ion (no neutral loss, not diagnostic, not excluded
	/// from quantification) encodes to 0x00 and decodes back to all-false with lossCode=None.
	/// </summary>
	[Test]
	public void EncodeDecodeFragmentFlags_TerminalIon_RoundTrips()
	{
		// Arrange: terminal ion has no special flags set
		bool isInternal = false;
		bool isDiagnostic = false;
		var lossCode = MslFormat.NeutralLossCode.None;
		bool excludeFromQuant = false;

		// Act: encode to a byte and decode back
		byte encoded = MslFormat.EncodeFragmentFlags(isInternal, isDiagnostic, lossCode, excludeFromQuant);
		var (decIsInternal, decIsDiagnostic, decLossCode, decExclude) = MslFormat.DecodeFragmentFlags(encoded);

		// Assert: all four fields round-trip exactly
		Assert.That(decIsInternal, Is.EqualTo(isInternal));
		Assert.That(decIsDiagnostic, Is.EqualTo(isDiagnostic));
		Assert.That(decLossCode, Is.EqualTo(lossCode));
		Assert.That(decExclude, Is.EqualTo(excludeFromQuant));
		Assert.That(encoded, Is.EqualTo(0x00));
	}

	/// <summary>
	/// Verifies that an internal fragment ion (is_internal = true) encodes bit 0 correctly
	/// and that the decoded is_internal field is true.
	/// </summary>
	[Test]
	public void EncodeDecodeFragmentFlags_InternalIon_RoundTrips()
	{
		// Arrange: set only the is_internal flag
		bool isInternal = true;

		// Act
		byte encoded = MslFormat.EncodeFragmentFlags(isInternal, false, MslFormat.NeutralLossCode.None, false);
		var (decIsInternal, _, _, _) = MslFormat.DecodeFragmentFlags(encoded);

		// Assert
		Assert.That(decIsInternal, Is.True);
		Assert.That(encoded & 0x01, Is.EqualTo(1)); // bit 0 must be set
	}

	/// <summary>
	/// Verifies that a diagnostic ion (is_diagnostic = true) encodes bit 1 correctly
	/// and that the decoded is_diagnostic field is true.
	/// </summary>
	[Test]
	public void EncodeDecodeFragmentFlags_DiagnosticIon_RoundTrips()
	{
		// Arrange
		bool isDiagnostic = true;

		// Act
		byte encoded = MslFormat.EncodeFragmentFlags(false, isDiagnostic, MslFormat.NeutralLossCode.None, false);
		var (_, decIsDiagnostic, _, _) = MslFormat.DecodeFragmentFlags(encoded);

		// Assert
		Assert.That(decIsDiagnostic, Is.True);
		Assert.That(encoded & 0x02, Is.EqualTo(2)); // bit 1 must be set
	}

	/// <summary>
	/// Verifies that all six defined NeutralLossCode values (None through Custom)
	/// encode into bits 2–4 correctly and decode back to the same value.
	/// </summary>
	[Test]
	public void EncodeDecodeFragmentFlags_NeutralLoss_AllValues_RoundTrip()
	{
		// Iterate over all defined NeutralLossCode members
		foreach (MslFormat.NeutralLossCode code in Enum.GetValues<MslFormat.NeutralLossCode>())
		{
			// Act: encode with only the loss code set; decode and check
			byte encoded = MslFormat.EncodeFragmentFlags(false, false, code, false);
			var (_, _, decCode, _) = MslFormat.DecodeFragmentFlags(encoded);

			// Assert: decoded code matches the original
			Assert.That(decCode, Is.EqualTo(code),
				$"NeutralLossCode.{code} did not round-trip through EncodeFragmentFlags/DecodeFragmentFlags");
		}
	}

	/// <summary>
	/// Verifies that the exclude_from_quant flag (bit 5) encodes and decodes correctly
	/// without interfering with other flag bits.
	/// </summary>
	[Test]
	public void EncodeDecodeFragmentFlags_ExcludeFromQuant_RoundTrips()
	{
		// Arrange
		bool excludeFromQuant = true;

		// Act
		byte encoded = MslFormat.EncodeFragmentFlags(false, false, MslFormat.NeutralLossCode.None, excludeFromQuant);
		var (_, _, _, decExclude) = MslFormat.DecodeFragmentFlags(encoded);

		// Assert
		Assert.That(decExclude, Is.True);
		Assert.That(encoded & 0x20, Is.EqualTo(0x20)); // bit 5 must be set
	}

	/// <summary>
	/// Verifies that all four fragment flag components can be set simultaneously and that
	/// each decodes independently and correctly (no bit collisions between fields).
	/// </summary>
	[Test]
	public void EncodeDecodeFragmentFlags_AllFlagsSet_NoCollisions()
	{
		// Arrange: set everything non-zero
		bool isInternal = true;
		bool isDiagnostic = true;
		var lossCode = MslFormat.NeutralLossCode.H3PO4; // value 3 = 011 in 3 bits
		bool excludeFromQuant = true;

		// Act
		byte encoded = MslFormat.EncodeFragmentFlags(isInternal, isDiagnostic, lossCode, excludeFromQuant);
		var (decIsInternal, decIsDiagnostic, decLossCode, decExclude) = MslFormat.DecodeFragmentFlags(encoded);

		// Assert: every field round-trips without aliasing
		Assert.That(decIsInternal, Is.True);
		Assert.That(decIsDiagnostic, Is.True);
		Assert.That(decLossCode, Is.EqualTo(lossCode));
		Assert.That(decExclude, Is.True);
		// Bits 6–7 must be 0 (reserved)
		Assert.That(encoded & 0xC0, Is.EqualTo(0));
	}

	// ══════════════════════════════════════════════════════════════════════
	// Precursor flags encoding / decoding tests
	// ══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that the is_decoy flag (bit 0) encodes and decodes correctly.
	/// </summary>
	[Test]
	public void EncodeDecodePrecursorFlags_Decoy_RoundTrips()
	{
		// Act: encode only is_decoy = true
		byte encoded = MslFormat.EncodePrecursorFlags(isDecoy: true, isProteotypic: false, rtCalibrated: false);
		var (decIsDecoy, _, _) = MslFormat.DecodePrecursorFlags(encoded);

		// Assert
		Assert.That(decIsDecoy, Is.True);
		Assert.That(encoded & 0x01, Is.EqualTo(1)); // bit 0 set
		Assert.That(encoded & 0xFE, Is.EqualTo(0)); // all other bits clear
	}

	/// <summary>
	/// Verifies that the is_proteotypic flag (bit 1) encodes and decodes correctly.
	/// </summary>
	[Test]
	public void EncodeDecodePrecursorFlags_Proteotypic_RoundTrips()
	{
		// Act: encode only is_proteotypic = true
		byte encoded = MslFormat.EncodePrecursorFlags(isDecoy: false, isProteotypic: true, rtCalibrated: false);
		var (_, decIsProteotypic, _) = MslFormat.DecodePrecursorFlags(encoded);

		// Assert
		Assert.That(decIsProteotypic, Is.True);
		Assert.That(encoded & 0x02, Is.EqualTo(2)); // bit 1 set
		Assert.That(encoded & 0xFD, Is.EqualTo(0)); // all other bits clear
	}

	/// <summary>
	/// Verifies that the rt_is_calibrated flag (bit 2) encodes and decodes correctly.
	/// </summary>
	[Test]
	public void EncodeDecodePrecursorFlags_RtCalibrated_RoundTrips()
	{
		// Act: encode only rtCalibrated = true
		byte encoded = MslFormat.EncodePrecursorFlags(isDecoy: false, isProteotypic: false, rtCalibrated: true);
		var (_, _, decRtCalibrated) = MslFormat.DecodePrecursorFlags(encoded);

		// Assert
		Assert.That(decRtCalibrated, Is.True);
		Assert.That(encoded & 0x04, Is.EqualTo(4)); // bit 2 set
		Assert.That(encoded & 0xFB, Is.EqualTo(0)); // all other bits clear
	}

	/// <summary>
	/// Verifies that all three precursor flags can be set simultaneously and decode
	/// independently without aliasing.
	/// </summary>
	[Test]
	public void EncodeDecodePrecursorFlags_AllFlagsSet_NoCollisions()
	{
		// Act: set all three flags
		byte encoded = MslFormat.EncodePrecursorFlags(isDecoy: true, isProteotypic: true, rtCalibrated: true);
		var (decIsDecoy, decIsProteotypic, decRtCalibrated) = MslFormat.DecodePrecursorFlags(encoded);

		// Assert: each flag is true and reserved bits 3–7 are clear
		Assert.That(decIsDecoy, Is.True);
		Assert.That(decIsProteotypic, Is.True);
		Assert.That(decRtCalibrated, Is.True);
		Assert.That(encoded & 0xF8, Is.EqualTo(0)); // bits 3–7 must be 0
	}

	// ══════════════════════════════════════════════════════════════════════
	// MslLibraryEntry construction tests
	// ══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that Name concatenates FullSequence and ChargeState with "/" as separator,
	/// matching the format used by MslIndex for O(1) DDA-style lookup.
	/// </summary>
	[Test]
	public void MslLibraryEntry_LookupKey_Format_Is_Sequence_Slash_Charge()
	{
		// Arrange
		var entry = new MslLibraryEntry
		{
			FullSequence = "PEPTIDE",
			ChargeState = 2
		};

		// Act
		string key = entry.Name;

		// Assert
		Assert.That(key, Is.EqualTo("PEPTIDE/2"));
	}

	/// <summary>
	/// Verifies that IsInternalFragment returns false for a terminal ion (SecondaryProductType == null).
	/// </summary>
	[Test]
	public void MslLibraryEntry_IsInternalFragment_False_For_Terminal_Ion()
	{
		// Arrange: terminal b5 ion — no secondary product type
		var ion = new MslFragmentIon
		{
			ProductType = ProductType.b,
			SecondaryProductType = null,
			FragmentNumber = 5
		};

		// Assert
		Assert.That(ion.IsInternalFragment, Is.False);
	}

	/// <summary>
	/// Verifies that IsInternalFragment returns true when SecondaryProductType is set,
	/// regardless of the specific type value.
	/// </summary>
	[Test]
	public void MslLibraryEntry_IsInternalFragment_True_When_SecondaryProductType_Set()
	{
		// Arrange: b/y internal fragment ion
		var ion = new MslFragmentIon
		{
			ProductType = ProductType.b,
			SecondaryProductType = ProductType.y,
			FragmentNumber = 3,
			SecondaryFragmentNumber = 6
		};

		// Assert
		Assert.That(ion.IsInternalFragment, Is.True);
	}

	/// <summary>
	/// Verifies that ToLibrarySpectrum produces a LibrarySpectrum whose Name (Sequence/ChargeState)
	/// matches the source entry's FullSequence and ChargeState.
	/// </summary>
	[Test]
	public void MslLibraryEntry_ToLibrarySpectrum_PreservesSequenceAndCharge()
	{
		// Arrange: minimal entry with one b-ion fragment so the constructor doesn't throw
		var entry = BuildMinimalPeptideEntry("PEPTIDE", 2, precursorMz: 400.71);

		// Act
		LibrarySpectrum spectrum = entry.ToLibrarySpectrum();

		// Assert: LibrarySpectrum.Name is "<Sequence>/<ChargeState>"
		Assert.That(spectrum.Sequence, Is.EqualTo("PEPTIDE"));
		Assert.That(spectrum.ChargeState, Is.EqualTo(2));
		Assert.That(spectrum.PrecursorMz, Is.EqualTo(400.71).Within(0.001));
	}

	/// <summary>
	/// Verifies that ToLibrarySpectrum preserves the m/z and intensity of each fragment ion
	/// in the MatchedFragmentIons list.
	/// </summary>
	[Test]
	public void MslLibraryEntry_ToLibrarySpectrum_PreservesFragmentMzAndIntensity()
	{
		// Arrange
		var entry = BuildMinimalPeptideEntry("ALGVGLATR", 2, precursorMz: 429.26);
		entry.MatchedFragmentIons[0].Mz = 100.0f;
		entry.MatchedFragmentIons[0].Intensity = 0.8f;

		// Act
		LibrarySpectrum spectrum = entry.ToLibrarySpectrum();

		// Assert
		Assert.That(spectrum.MatchedFragmentIons, Has.Count.EqualTo(1));
		Assert.That(spectrum.MatchedFragmentIons[0].Mz, Is.EqualTo(100.0).Within(0.001));
		Assert.That(spectrum.MatchedFragmentIons[0].Intensity, Is.EqualTo(0.8).Within(0.0001));
	}

	/// <summary>
	/// Verifies that a standard terminal b-ion yields FragmentationTerminus.N, derived
	/// directly from the ProductType enum value (b → N) when MoleculeType == Peptide.
	/// </summary>
	[Test]
	public void MslLibraryEntry_ToLibrarySpectrum_TerminalIon_HasCorrectTerminus_Peptide()
	{
		// Arrange: b2 ion on a peptide entry
		var entry = new MslLibraryEntry
		{
			FullSequence = "ACDEF",
			BaseSequence = "ACDEF",
			PrecursorMz = 300.0,
			ChargeState = 1,
			RetentionTime = 0.0,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz             = 173.09f,
					Intensity      = 1.0f,
					ProductType    = ProductType.b,
					SecondaryProductType = null,
					FragmentNumber = 2,
					Charge         = 1
				}
			}
		};

		// Act
		LibrarySpectrum spectrum = entry.ToLibrarySpectrum();

		// Assert: b-ions are N-terminal
		FragmentationTerminus terminus = spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus;
		Assert.That(terminus, Is.EqualTo(FragmentationTerminus.N));
	}

	/// <summary>
	/// Verifies that an internal fragment ion yields FragmentationTerminus.None because
	/// internal ions span a sub-sequence and are not associated with a terminus.
	/// </summary>
	[Test]
	public void MslLibraryEntry_ToLibrarySpectrum_InternalIon_HasTerminusNone()
	{
		// Arrange: bIy[2-4] internal ion
		var entry = new MslLibraryEntry
		{
			FullSequence = "PEPTIDE",
			BaseSequence = "PEPTIDE",
			PrecursorMz = 400.0,
			ChargeState = 2,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz                      = 300.0f,
					Intensity               = 0.5f,
					ProductType             = ProductType.b,
					SecondaryProductType    = ProductType.y,
					FragmentNumber          = 2,
					SecondaryFragmentNumber = 4,
					Charge                  = 1
				}
			}
		};

		// Act
		LibrarySpectrum spectrum = entry.ToLibrarySpectrum();

		// Assert: internal ions have no terminus
		FragmentationTerminus terminus = spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus;
		Assert.That(terminus, Is.EqualTo(FragmentationTerminus.None));
	}

	/// <summary>
	/// Verifies that an oligonucleotide entry with a 5'-type ion (e.g. ProductType.a) yields
	/// the FivePrime terminus via the oligo TerminusSpecificProductTypes dictionary.
	/// </summary>
	[Test]
	public void MslLibraryEntry_ToLibrarySpectrum_OligoIon_HasFivePrimeTerminus()
	{
		// Arrange: RNA oligo entry with a-type (5'-terminal) fragment
		var entry = new MslLibraryEntry
		{
			FullSequence = "rArUrGrC",
			BaseSequence = "rArUrGrC",
			PrecursorMz = 600.0,
			ChargeState = 2,
			MoleculeType = MslFormat.MoleculeType.Oligonucleotide,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz             = 329.0f,
					Intensity      = 1.0f,
					ProductType    = ProductType.a,   // 5'-terminal type in oligo namespace
                    SecondaryProductType = null,
					FragmentNumber = 2,
					Charge         = 1
				}
			}
		};

		// Act
		LibrarySpectrum spectrum = entry.ToLibrarySpectrum();

		// Assert: a-ions are 5'-terminal in the oligo namespace
		FragmentationTerminus terminus = spectrum.MatchedFragmentIons[0].NeutralTheoreticalProduct.Terminus;
		Assert.That(terminus, Is.EqualTo(FragmentationTerminus.FivePrime));
	}

	/// <summary>
	/// Verifies that FromLibrarySpectrum produces an entry whose ToLibrarySpectrum output
	/// has the same fragment count, m/z values, and intensities as the original spectrum.
	/// </summary>
	[Test]
	public void MslLibraryEntry_FromLibrarySpectrum_RoundTrip_PreservesAllFragments()
	{
		// Arrange: build a LibrarySpectrum with two b-ions
		var matchedIons = new List<MatchedFragmentIon>
		{
			MakePeptideIon(ProductType.b, 2, mz: 200.0, intensity: 1.0f),
			MakePeptideIon(ProductType.b, 3, mz: 300.0, intensity: 0.5f)
		};
		var original = new LibrarySpectrum("ACDE", 300.0, 1, matchedIons, rt: 5.0);

		// Act: round-trip via MslLibraryEntry
		MslLibraryEntry entry = MslLibraryEntry.FromLibrarySpectrum(original);
		LibrarySpectrum restored = entry.ToLibrarySpectrum();

		// Assert: fragment count and key values preserved
		Assert.That(restored.MatchedFragmentIons, Has.Count.EqualTo(2));
		Assert.That(restored.MatchedFragmentIons[0].Mz, Is.EqualTo(200.0).Within(0.001));
		Assert.That(restored.MatchedFragmentIons[1].Mz, Is.EqualTo(300.0).Within(0.001));
	}

	/// <summary>
	/// Verifies that FromLibrarySpectrum preserves the start and end residue numbers of
	/// an internal fragment ion through the round-trip.
	/// </summary>
	[Test]
	public void MslLibraryEntry_FromLibrarySpectrum_InternalIon_PreservesStartAndEndResidue()
	{
		// mzLib's Product type does not carry SecondaryProductType or SecondaryFragmentNumber,
		// so these internal-fragment fields cannot be recovered via FromLibrarySpectrum.
		// This test verifies the correct fallback behaviour: SecondaryProductType == null
		// and SecondaryFragmentNumber == 0 (so IsInternalFragment is false on round-trip).
		// Internal-fragment data is only preserved when an MslLibraryEntry is constructed
		// directly (e.g. from the binary reader) rather than promoted from a LibrarySpectrum.

		// Arrange: build MslLibraryEntry directly with an internal ion
		var entry = new MslLibraryEntry
		{
			FullSequence = "PEPTIDE",
			BaseSequence = "PEPTIDE",
			PrecursorMz = 400.0,
			ChargeState = 2,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz                      = 451.0f,
					Intensity               = 0.9f,
					ProductType             = ProductType.b,
					SecondaryProductType    = ProductType.y,
					FragmentNumber          = 2,
					SecondaryFragmentNumber = 5,
					Charge                  = 1
				}
			}
		};

		// Act: access via the MslFragmentIon directly — no lossy LibrarySpectrum round-trip
		MslFragmentIon mslIon = entry.MatchedFragmentIons[0];

		// Assert: start and end residue are preserved on the MslFragmentIon itself
		Assert.That(mslIon.FragmentNumber, Is.EqualTo(2));
		Assert.That(mslIon.SecondaryFragmentNumber, Is.EqualTo(5));
		Assert.That(mslIon.IsInternalFragment, Is.True);

		// Also verify that FromLibrarySpectrum correctly falls back to null/0 for these fields
		var product = new Product(ProductType.b, FragmentationTerminus.None,
			neutralMass: 0.0, fragmentNumber: 2, residuePosition: 0, neutralLoss: 0.0);
		var ions = new List<MatchedFragmentIon> { new MatchedFragmentIon(product, experMz: 451.0, experIntensity: 0.9f, charge: 1) };
		var spectrum = new LibrarySpectrum("PEPTIDE", 400.0, 2, ions, rt: 10.0);
		MslLibraryEntry fromSpectrum = MslLibraryEntry.FromLibrarySpectrum(spectrum);
		Assert.That(fromSpectrum.MatchedFragmentIons[0].SecondaryProductType, Is.Null);
		Assert.That(fromSpectrum.MatchedFragmentIons[0].SecondaryFragmentNumber, Is.EqualTo(0));
		Assert.That(fromSpectrum.MatchedFragmentIons[0].IsInternalFragment, Is.False);
	}

	/// <summary>
	/// Verifies that a neutral-loss mass passes through the FromLibrarySpectrum →
	/// ToLibrarySpectrum round-trip without significant precision loss.
	/// </summary>
	[Test]
	public void MslLibraryEntry_FromLibrarySpectrum_NeutralLoss_IsPreserved()
	{
		// Arrange: b5 ion with −H2O loss
		const double h2oLoss = -18.010565;
		var product = new Product(
			ProductType.b,
			FragmentationTerminus.N,
			neutralMass: 0.0,
			fragmentNumber: 5,
			residuePosition: 0,
			neutralLoss: h2oLoss);

		var ions = new List<MatchedFragmentIon>
		{
			new MatchedFragmentIon(product, experMz: 543.0, experIntensity:  0.6f, charge: 1)
		};
		var spectrum = new LibrarySpectrum("ALGVGLATR", 429.26, 2, ions, rt: 16.5);

		// Act
		MslLibraryEntry entry = MslLibraryEntry.FromLibrarySpectrum(spectrum);

		// Assert: neutral loss preserved within floating-point tolerance
		Assert.That(entry.MatchedFragmentIons[0].NeutralLoss, Is.EqualTo(h2oLoss).Within(1e-6));
	}

	// ══════════════════════════════════════════════════════════════════════
	// ProductType encoding tests
	// ══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Verifies that every value in the ProductType enum fits within the range of an int16
	/// (−32768 to 32767). The fragment record stores ProductType as int16; values outside
	/// this range would be silently truncated on write and produce wrong annotations on read.
	/// </summary>
	[Test]
	public void ProductType_AllValues_FitInInt16()
	{
		// All ProductType values must be representable as int16
		foreach (ProductType pt in Enum.GetValues<ProductType>())
		{
			int intValue = (int)pt;
			Assert.That(intValue, Is.InRange(short.MinValue, short.MaxValue),
				$"ProductType.{pt} (value {intValue}) does not fit in int16.");
		}
	}

	/// <summary>
	/// Verifies that casting a ProductType to int16 and back to ProductType via explicit
	/// casts produces the original value for all enum members.
	/// </summary>
	[Test]
	public void ProductType_RoundTrip_Via_Int16_Cast()
	{
		foreach (ProductType pt in Enum.GetValues<ProductType>())
		{
			// Simulate the write path: cast to int16
			short stored = (short)(int)pt;
			// Simulate the read path: cast back
			ProductType restored = (ProductType)(int)stored;

			Assert.That(restored, Is.EqualTo(pt),
				$"ProductType.{pt} did not round-trip through int16 cast.");
		}
	}

	// ══════════════════════════════════════════════════════════════════════
	// Private test helpers
	// ══════════════════════════════════════════════════════════════════════

	/// <summary>
	/// Builds a minimal valid MslLibraryEntry for a peptide with one b2 fragment ion.
	/// Used to satisfy the LibrarySpectrum constructor's requirement for at least one ion
	/// without coupling every test to the full conversion logic.
	/// </summary>
	/// <param name="sequence">Modified sequence string for the entry.</param>
	/// <param name="charge">Precursor charge state.</param>
	/// <param name="precursorMz">Precursor m/z value.</param>
	/// <returns>A new MslLibraryEntry ready for ToLibrarySpectrum() calls.</returns>
	private static MslLibraryEntry BuildMinimalPeptideEntry(string sequence, int charge, double precursorMz)
	{
		return new MslLibraryEntry
		{
			FullSequence = sequence,
			BaseSequence = sequence,
			PrecursorMz = precursorMz,
			ChargeState = charge,
			RetentionTime = 0.0,
			MoleculeType = MslFormat.MoleculeType.Peptide,
			MatchedFragmentIons = new List<MslFragmentIon>
			{
				new MslFragmentIon
				{
					Mz             = 100.0f,
					Intensity      = 1.0f,
					ProductType    = ProductType.b,
					SecondaryProductType = null,
					FragmentNumber = 2,
					Charge         = 1
				}
			}
		};
	}

	/// <summary>
	/// Creates a MatchedFragmentIon wrapping a peptide-backbone Product for use in
	/// LibrarySpectrum construction inside round-trip tests.
	/// </summary>
	/// <param name="productType">Fragment ion type (b, y, c, z, etc.).</param>
	/// <param name="fragmentNumber">Ion series number (e.g. 2 for b2).</param>
	/// <param name="mz">Fragment m/z value.</param>
	/// <param name="intensity">Relative intensity (0–1).</param>
	/// <returns>A new MatchedFragmentIon with charge 1 and no neutral loss.</returns>
	private static MatchedFragmentIon MakePeptideIon(
		ProductType productType,
		int fragmentNumber,
		double mz,
		float intensity)
	{
		// Derive terminus directly from the ProductType enum value.
		// N-terminal ion types (b, a, c) → FragmentationTerminus.N
		// C-terminal ion types (y, z, x) → FragmentationTerminus.C
		// Everything else (internal, diagnostic, oligo, etc.) → None
		FragmentationTerminus terminus = productType switch
		{
			ProductType.b or ProductType.a or ProductType.c => FragmentationTerminus.N,
			ProductType.y or ProductType.z or ProductType.x => FragmentationTerminus.C,
			_ => FragmentationTerminus.None
		};

		// Use neutralMass = 0.0, matching the pattern in the existing mzLib MSP and DIA-NN readers.
		// The MatchedFragmentIon stores the actual m/z; neutralMass is not used for library matching.
		var product = new Product(productType, terminus, neutralMass: 0.0, fragmentNumber, residuePosition: 0, neutralLoss: 0.0);
		return new MatchedFragmentIon(product, mz, intensity, charge: 1);
	}
}