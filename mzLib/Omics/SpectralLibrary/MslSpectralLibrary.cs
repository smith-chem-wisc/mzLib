using System.Runtime.CompilerServices;

namespace Omics.SpectralMatch.MslSpectralLibrary;

/// <summary>
/// All format-level constants, enumerations, and bitfield encoding/decoding helpers for the
/// .msl (mzLib Spectral Library) binary format. This class is the single authoritative source
/// for magic bytes, version numbers, record sizes, flag layouts, and molecule/source
/// classifications. No I/O logic lives here; every field and method is pure data definition
/// or pure bit manipulation.
/// </summary>
public static class MslFormat
{
	// ──────────────────────────────────────────────────────────────────────
	// Magic and versioning
	// ──────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Four-byte file magic: ASCII "MZLB" (0x4D 0x5A 0x4C 0x42).
	/// Every .msl file must begin with exactly these bytes; a mismatch indicates a
	/// corrupt or non-MSL file and must be rejected before any further parsing.
	/// Use <see cref="MagicMatches"/> to check raw bytes, or <see cref="MagicAsUInt32"/>
	/// to populate the uint fields in <c>MslFileHeader</c> and <c>MslFooter</c>.
	/// </summary>
	public static readonly byte[] Magic = { 0x4D, 0x5A, 0x4C, 0x42 };

	/// <summary>
	/// The file magic packed as a uint32 for writing into the <c>MslFileHeader.Magic</c>
	/// and <c>MslFooter.TrailingMagic</c> struct fields (which are uint rather than
	/// <c>unsafe fixed byte[4]</c> to avoid requiring the /unsafe compiler flag).
	/// Value is 0x4D5A4C42u; when written by <c>BinaryWriter.Write(uint)</c> on a
	/// little-endian system the on-disk bytes are 0x42 0x4C 0x5A 0x4D — use
	/// <c>BinaryPrimitives.WriteUInt32BigEndian</c> if byte order must match the
	/// canonical Magic array exactly.
	/// </summary>
	public const uint MagicAsUInt32 = 0x4D5A_4C42u;

	/// <summary>
	/// Format version stored at bytes 4–7 of the file header.
	/// Increment this constant whenever any struct layout, field offset, or semantic
	/// meaning changes in an incompatible way. Readers must reject files whose stored
	/// version differs from this value unless they explicitly handle the older version.
	/// <list type="bullet">
	///   <item>Version 1 — original format; <c>ExtAnnotationTableOffset</c> field is <c>Reserved</c> (must be 0).</item>
	///   <item>Version 2 — extended annotation table support (<see cref="FileFlagHasExtAnnotations"/>); <c>Reserved</c> repurposed as <c>ExtAnnotationTableOffset</c>.</item>
	///   <item>Version 3 — optional zstd block compression (<see cref="FileFlagIsCompressed"/>); 16-byte compression descriptor inserted at offset 64 when compressed.</item>
	/// </list>
	/// </summary>
	public const int CurrentVersion = 3;

	// ──────────────────────────────────────────────────────────────────────
	// Fixed record sizes (bytes)
	// ──────────────────────────────────────────────────────────────────────
	// All sizes are validated at startup by MslStructs.SizeCheck().
	// Do NOT change any size constant without also incrementing CurrentVersion.

	/// <summary>Size of the file header section in bytes (MslFileHeader struct).</summary>
	public const int HeaderSize = 64;

	/// <summary>Size of one protein table entry in bytes (MslProteinRecord struct).</summary>
	public const int ProteinRecordSize = 24;

	/// <summary>Size of one precursor entry in bytes (MslPrecursorRecord struct).</summary>
	public const int PrecursorRecordSize = 56;

	/// <summary>Size of one fragment ion entry in bytes (MslFragmentRecord struct).</summary>
	public const int FragmentRecordSize = 20;

	/// <summary>Size of the file footer in bytes (MslFooter struct).</summary>
	public const int FooterSize = 20;

	// ──────────────────────────────────────────────────────────────────────
	// Enumerations
	// ──────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Classifies what a precursor entry represents at the molecular level.
	/// Stored as int16 in the PrecursorRecord at offset 48.
	/// Determines which fragmentation namespace to use when reconstructing
	/// Product.Terminus and Product.NeutralMass at read time.
	/// </summary>
	public enum MoleculeType : short
	{
		/// <summary>Standard tryptic or non-tryptic peptide. Uses peptide fragmentation namespace.</summary>
		Peptide = 0,

		/// <summary>
		/// Top-down intact protein or large proteoform. Same backbone chemistry as Peptide
		/// but potentially very large StrippedSeqLength values.
		/// </summary>
		Proteoform = 1,

		/// <summary>RNA or DNA oligonucleotide. Uses Omics.Fragmentation.Oligo namespace.</summary>
		Oligonucleotide = 2,

		/// <summary>
		/// Peptide with intact glycan. Uses peptide backbone fragmentation plus
		/// glycan-specific Ycore / Y ions.
		/// </summary>
		Glycopeptide = 3
	}

	/// <summary>
	/// Describes how a precursor's fragment intensities were obtained.
	/// Stored as a byte in the PrecursorRecord at offset 55 (the SourceType field).
	/// </summary>
	public enum SourceType : byte
	{
		/// <summary>Fragment intensities predicted by an in-silico model (Prosit, internal ONNX, etc.).</summary>
		Predicted = 0,

		/// <summary>Fragment intensities measured directly from a real DIA or DDA experiment.</summary>
		Empirical = 1,

		/// <summary>
		/// Started as a predicted entry but the retention time and/or ion mobility have been
		/// overwritten with empirically measured values. Fragment intensities remain predicted.
		/// </summary>
		EmpiricalRefined = 2
	}

	/// <summary>
	/// Named neutral-loss codes for the 3-bit neutral_loss_code field occupying bits 2–4 of
	/// the fragment flags byte. Only the seven values below are assigned; all others are reserved.
	/// For the Custom value the actual mass shift is stored in the MslFragmentIon.NeutralLoss
	/// double field and is looked up via the extended_annotation_idx mechanism.
	/// </summary>
	public enum NeutralLossCode : byte
	{
		/// <summary>No neutral loss. NeutralLoss double is 0.</summary>
		None = 0,

		/// <summary>Loss of water: −18.010565 Da.</summary>
		H2O = 1,

		/// <summary>Loss of ammonia: −17.026549 Da.</summary>
		NH3 = 2,

		/// <summary>Loss of phosphoric acid: −97.976895 Da (phospho serine/threonine/tyrosine).</summary>
		H3PO4 = 3,

		/// <summary>Loss of metaphosphoric acid: −79.966331 Da (phospho, alternative loss).</summary>
		HPO3 = 4,

		/// <summary>Combined H3PO4 + H2O loss.</summary>
		PlusH2O = 5,

		/// <summary>
		/// Custom mass loss not covered by the named codes above.
		/// The exact mass shift (negative = loss, positive = gain) is stored in
		/// MslFragmentIon.NeutralLoss and must be written to / read from the binary
		/// extended annotation table.
		/// </summary>
		Custom = 6
	}

	// ──────────────────────────────────────────────────────────────────────
	// Fragment flags byte layout (byte 19 of each fragment record)
	// ──────────────────────────────────────────────────────────────────────
	// bit 0     : is_internal          – 1 if SecondaryProductType != null
	// bit 1     : is_diagnostic        – 1 if ProductType == D
	// bits 2–4  : neutral_loss_code    – NeutralLossCode (3 bits, values 0–6)
	// bit 5     : exclude_from_quant   – mirrors DIA-NN ExcludeFromAssay
	// bits 6–7  : reserved             – must be written as 0; ignored on read

	private const byte FragFlagIsInternal = 0b_0000_0001;
	private const byte FragFlagIsDiagnostic = 0b_0000_0010;
	private const byte FragFlagNeutralLossMask = 0b_0001_1100; // bits 2–4
	private const int FragFlagNeutralLossShift = 2;
	private const byte FragFlagExcludeFromQuant = 0b_0010_0000;

	// ──────────────────────────────────────────────────────────────────────
	// Precursor flags byte layout (byte 54 of each precursor record)
	// ──────────────────────────────────────────────────────────────────────
	// bit 0     : is_decoy
	// bit 1     : is_proteotypic
	// bit 2     : rt_is_calibrated – 0 = stored as iRT; 1 = run-specific RT in minutes
	// bits 3–7  : reserved

	private const byte PrecFlagIsDecoy = 0b_0000_0001;
	private const byte PrecFlagIsProteotypic = 0b_0000_0010;
	private const byte PrecFlagRtCalibrated = 0b_0000_0100;

	// ──────────────────────────────────────────────────────────────────────
	// File-level flags int32 layout (bytes 8–11 of file header)
	// ──────────────────────────────────────────────────────────────────────
	// bit 0     : has_ion_mobility    – at least one precursor has a meaningful 1/K0
	// bit 1     : has_protein_data    – protein table is populated
	// bit 2     : has_gene_data       – gene names are present in the string table
	// bit 3     : is_predicted        – all entries generated by an in-silico model
	// bit 4     : has_ext_annotations – extended annotation table section present (version 2+)
	// bit 5     : is_compressed       – fragment section is zstd-compressed (version 3+);
	//                                   FragmentBlockOffset values are offsets into the
	//                                   decompressed fragment buffer, not absolute file offsets;
	//                                   a 16-byte compression descriptor is present at offset 64
	// bits 6–31 : reserved            – must be written as 0; ignored on read

	/// <summary>File-level flag: at least one precursor has a meaningful ion-mobility (1/K0) value.</summary>
	public const int FileFlagHasIonMobility = 1 << 0;

	/// <summary>File-level flag: the protein table section is populated (NProteins > 0).</summary>
	public const int FileFlagHasProteinData = 1 << 1;

	/// <summary>File-level flag: gene names are present in the string table.</summary>
	public const int FileFlagHasGeneData = 1 << 2;

	/// <summary>File-level flag: the library was entirely generated by an in-silico prediction model.</summary>
	public const int FileFlagIsPredicted = 1 << 3;

	/// <summary>
	/// File-level flag: the extended annotation table section is present (version 2+).
	/// When set, <c>MslFileHeader.ExtAnnotationTableOffset</c> contains the absolute byte
	/// offset of the extended annotation table; the table holds custom neutral-loss masses
	/// indexed by <c>MslFragmentRecord.ResiduePosition</c> when <c>neutral_loss_code == Custom</c>.
	/// </summary>
	public const int FileFlagHasExtAnnotations = 1 << 4;

	/// <summary>
	/// File-level flag: the fragment section (plus extended annotation table if present) is
	/// compressed as a single zstd frame (version 3+).
	/// <para>
	/// When set, <c>MslPrecursorRecord.FragmentBlockOffset</c> values are byte offsets into
	/// the <em>decompressed</em> fragment buffer rather than absolute file positions. A 16-byte
	/// compression descriptor is inserted at file offset 64 (immediately after the file header)
	/// and contains the compressed and uncompressed sizes of the fragment data:
	/// </para>
	/// <code>
	/// [int64 CompressedFragmentSize]    byte count of the compressed zstd frame in the file
	/// [int64 UncompressedFragmentSize]  byte count after full decompression
	/// </code>
	/// <para>
	/// Index-only load (<see cref="MslLibrary.LoadIndexOnly"/>) is <b>not available</b> for
	/// compressed files; the reader always performs full decompression regardless of which
	/// load method was called, and <c>MslLibrary.IsIndexOnly</c> returns <c>false</c>.
	/// </para>
	/// </summary>
	public const int FileFlagIsCompressed = 1 << 5;

	// ──────────────────────────────────────────────────────────────────────
	// Static helpers
	// ──────────────────────────────────────────────────────────────────────

	/// <summary>
	/// Returns true if the first four bytes of a file match the MZLB magic sequence.
	/// Must be called before any further parsing; callers should reject the file immediately
	/// if this returns false.
	/// </summary>
	/// <param name="firstFourBytes">
	///   A span containing at least the first four bytes of the file being inspected.
	///   If fewer than four bytes are provided the method returns false.
	/// </param>
	/// <returns>True if bytes 0–3 are 0x4D, 0x5A, 0x4C, 0x42; otherwise false.</returns>
	[MethodImpl(MethodImplOptions.AggressiveInlining)]
	public static bool MagicMatches(ReadOnlySpan<byte> firstFourBytes)
	{
		// Guard: need at least 4 bytes to compare
		if (firstFourBytes.Length < 4)
			return false;

		// Compare each byte of the magic against the stored constant array
		return firstFourBytes[0] == Magic[0]
			&& firstFourBytes[1] == Magic[1]
			&& firstFourBytes[2] == Magic[2]
			&& firstFourBytes[3] == Magic[3];
	}

	/// <summary>
	/// Encodes the fragment flags byte (byte 19 of an MslFragmentRecord) from its four
	/// independent component values. The resulting byte is stored verbatim in the binary file.
	/// </summary>
	/// <param name="isInternal">
	///   True when the fragment is an internal ion (SecondaryProductType is not null).
	///   Sets bit 0 of the returned byte.
	/// </param>
	/// <param name="isDiagnostic">
	///   True when the fragment is a diagnostic ion (ProductType == D).
	///   Sets bit 1 of the returned byte.
	/// </param>
	/// <param name="lossCode">
	///   The neutral-loss category from the NeutralLossCode enumeration.
	///   Packed into bits 2–4 (maximum representable value is 7; values above 6 are reserved).
	/// </param>
	/// <param name="excludeFromQuant">
	///   True when the fragment should be excluded from quantification (mirrors DIA-NN
	///   ExcludeFromAssay). Sets bit 5 of the returned byte.
	/// </param>
	/// <returns>A single byte with bits set according to the parameters above.</returns>
	public static byte EncodeFragmentFlags(
		bool isInternal,
		bool isDiagnostic,
		NeutralLossCode lossCode,
		bool excludeFromQuant)
	{
		// Start with zero; OR in each flag bit individually for clarity
		byte result = 0;

		if (isInternal)
			result |= FragFlagIsInternal;

		if (isDiagnostic)
			result |= FragFlagIsDiagnostic;

		// Shift the 3-bit loss code into bits 2–4
		result |= (byte)(((byte)lossCode << FragFlagNeutralLossShift) & FragFlagNeutralLossMask);

		if (excludeFromQuant)
			result |= FragFlagExcludeFromQuant;

		return result;
	}

	/// <summary>
	/// Decodes the fragment flags byte back into its four independent component values.
	/// The inverse of <see cref="EncodeFragmentFlags"/>; used by the reader to populate
	/// MslFragmentIon fields from the stored binary byte.
	/// </summary>
	/// <param name="flags">
	///   The raw flags byte read from offset 19 of an MslFragmentRecord.
	/// </param>
	/// <returns>
	///   A value tuple with four named fields:
	///   <list type="bullet">
	///     <item><c>isInternal</c>       — true if bit 0 is set</item>
	///     <item><c>isDiagnostic</c>     — true if bit 1 is set</item>
	///     <item><c>lossCode</c>         — NeutralLossCode extracted from bits 2–4</item>
	///     <item><c>excludeFromQuant</c> — true if bit 5 is set</item>
	///   </list>
	/// </returns>
	public static (bool isInternal, bool isDiagnostic, NeutralLossCode lossCode, bool excludeFromQuant)
		DecodeFragmentFlags(byte flags)
	{
		// Extract each bit or bit-field independently
		bool isInternal = (flags & FragFlagIsInternal) != 0;
		bool isDiagnostic = (flags & FragFlagIsDiagnostic) != 0;
		bool excludeFromQuant = (flags & FragFlagExcludeFromQuant) != 0;

		// Shift bits 2–4 back down to the low bits, then cast to enum
		NeutralLossCode lossCode =
			(NeutralLossCode)((flags & FragFlagNeutralLossMask) >> FragFlagNeutralLossShift);

		return (isInternal, isDiagnostic, lossCode, excludeFromQuant);
	}

	/// <summary>
	/// Encodes the precursor flags byte (byte 54 of an MslPrecursorRecord) from its three
	/// independent Boolean properties.
	/// </summary>
	/// <param name="isDecoy">
	///   True for decoy precursors. Sets bit 0 of the returned byte.
	/// </param>
	/// <param name="isProteotypic">
	///   True when the peptide is expected to be uniquely detectable (proteotypic).
	///   Sets bit 1 of the returned byte.
	/// </param>
	/// <param name="rtCalibrated">
	///   True when the stored RT value is a run-specific calibrated retention time in minutes
	///   rather than an iRT prediction. Sets bit 2 of the returned byte.
	/// </param>
	/// <returns>A single byte with bits set according to the parameters above.</returns>
	public static byte EncodePrecursorFlags(bool isDecoy, bool isProteotypic, bool rtCalibrated)
	{
		// Build the flag byte by ORing individual bit constants
		byte result = 0;

		if (isDecoy)
			result |= PrecFlagIsDecoy;

		if (isProteotypic)
			result |= PrecFlagIsProteotypic;

		if (rtCalibrated)
			result |= PrecFlagRtCalibrated;

		return result;
	}

	/// <summary>
	/// Decodes the precursor flags byte back into its three independent Boolean properties.
	/// The inverse of <see cref="EncodePrecursorFlags"/>; used by the reader to populate
	/// MslLibraryEntry boolean fields from the stored binary byte.
	/// </summary>
	/// <param name="flags">
	///   The raw precursor flags byte read from offset 54 of an MslPrecursorRecord.
	/// </param>
	/// <returns>
	///   A value tuple with three named fields:
	///   <list type="bullet">
	///     <item><c>isDecoy</c>        — true if bit 0 is set</item>
	///     <item><c>isProteotypic</c>  — true if bit 1 is set</item>
	///     <item><c>rtCalibrated</c>   — true if bit 2 is set</item>
	///   </list>
	/// </returns>
	public static (bool isDecoy, bool isProteotypic, bool rtCalibrated) DecodePrecursorFlags(byte flags)
	{
		// Inspect each bit in isolation
		bool isDecoy = (flags & PrecFlagIsDecoy) != 0;
		bool isProteotypic = (flags & PrecFlagIsProteotypic) != 0;
		bool rtCalibrated = (flags & PrecFlagRtCalibrated) != 0;

		return (isDecoy, isProteotypic, rtCalibrated);
	}
}