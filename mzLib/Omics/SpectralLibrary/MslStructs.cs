using System.Runtime.InteropServices;

namespace Omics.SpectralMatch.MslSpectralLibrary;

/// <summary>
/// Fixed-layout binary structs that map directly to the on-disk record formats of the
/// .msl (mzLib Spectral Library) binary format. All structs use Pack = 1 to prevent
/// the CLR from inserting any alignment padding; every field therefore sits at exactly
/// the byte offset documented in its XML comment.
///
/// Consumers must never create struct instances by hand and then write them to disk
/// without first calling <see cref="SizeCheck"/> at application start-up, which
/// verifies that <c>Marshal.SizeOf&lt;T&gt;()</c> matches every documented size constant
/// in <see cref="MslFormat"/>. Any mismatch indicates a Pack setting error and is treated
/// as a programming bug rather than a recoverable runtime condition.
/// </summary>
public static class MslStructs
{
	/// <summary>
	/// Validates that all five binary structs have exactly the sizes declared in
	/// <see cref="MslFormat"/>. Should be called once at application start-up (e.g. from a
	/// static constructor on the writer or reader). Throws <see cref="InvalidOperationException"/>
	/// if any size does not match; a mismatch almost always indicates an incorrect Pack
	/// setting or an accidental extra field insertion.
	/// </summary>
	/// <exception cref="InvalidOperationException">
	///   Thrown when at least one struct's Marshal.SizeOf does not equal its declared constant.
	///   The message names the offending struct and shows both the expected and actual sizes.
	/// </exception>
	public static void SizeCheck()
	{
		// Verify each struct independently so the error message is unambiguous
		AssertSize<MslFileHeader>(MslFormat.HeaderSize, nameof(MslFileHeader));
		AssertSize<MslProteinRecord>(MslFormat.ProteinRecordSize, nameof(MslProteinRecord));
		AssertSize<MslPrecursorRecord>(MslFormat.PrecursorRecordSize, nameof(MslPrecursorRecord));
		AssertSize<MslFragmentRecord>(MslFormat.FragmentRecordSize, nameof(MslFragmentRecord));
		AssertSize<MslFooter>(MslFormat.FooterSize, nameof(MslFooter));
	}

	// Helper: compare Marshal.SizeOf<T> against an expected value and throw if they differ
	private static void AssertSize<T>(int expected, string typeName) where T : struct
	{
		// Compute the unmanaged size; this is the value the OS / BinaryReader will see
		int actual = Marshal.SizeOf<T>();

		if (actual != expected)
			throw new InvalidOperationException(
				$"Struct size mismatch for {typeName}: expected {expected} bytes, " +
				$"got {actual} bytes. Verify [StructLayout(LayoutKind.Sequential, Pack = 1)] " +
				$"is applied and no unintended fields were added.");
	}
}

// ──────────────────────────────────────────────────────────────────────────────
// MslFileHeader — 64 bytes
// ──────────────────────────────────────────────────────────────────────────────

/// <summary>
/// The 64-byte file header occupying the first bytes of every .msl file.
/// After reading the header the caller must validate that FormatVersion falls within
/// the supported range [MslFormat.MinSupportedVersion, MslFormat.CurrentVersion]
/// before proceeding. Version 1 is the original layout. Version 2 repurposes the
/// formerly-Reserved field at offset 28 as ExtAnnotationTableOffset for custom
/// neutral-loss masses. Version 3 adds optional zstd block compression of the
/// fragment section. An unsupported version should produce a clear error rather
/// than attempting to read a misaligned layout.
/// </summary>
[StructLayout(LayoutKind.Sequential, Pack = 1)]
public struct MslFileHeader
{
	/// <summary>
	/// Offset 0, 4 bytes. File magic stored as a uint for safe (non-unsafe) layout.
	/// The canonical value is 0x4D5A4C42 (little-endian on-disk: 0x42 0x4C 0x5A 0x4D).
	/// Writers must call <c>BinaryPrimitives.WriteUInt32BigEndian</c> or write the raw
	/// MslFormat.Magic bytes directly; readers should pass the first four bytes to
	/// MslFormat.MagicMatches() rather than comparing this uint directly.
	/// </summary>
	public uint Magic;

	/// <summary>
	/// Offset 4, 4 bytes. Format version. Valid values are MslFormat.MinSupportedVersion
	/// through MslFormat.CurrentVersion (inclusive). Readers must reject files outside
	/// this range with a clear error. Increment MslFormat.CurrentVersion — and add a
	/// corresponding entry in the version history in MslFormat.cs — whenever a layout
	/// change is made.
	/// Version history:
	///   1 — original layout; Reserved field at offset 28 is always 0.
	///   2 — offset 28 repurposed as ExtAnnotationTableOffset for custom neutral losses
	///       (FileFlagHasExtAnnotations). All v1 files remain readable by v2+ readers.
	///   3 — optional zstd block compression (FileFlagIsCompressed); 16-byte compression
	///       descriptor inserted at file offset 64 when compressed; FragmentBlockOffset
	///       values become decompressed-buffer-relative rather than absolute file offsets.
	/// </summary>
	public int FormatVersion;

	/// <summary>
	/// Offset 8, 4 bytes. Bitfield of file-level flags; see MslFormat.FileFlag* constants.
	/// Readers use this field to decide which optional sections are present without
	/// scanning ahead; writers must set this accurately before flushing the header.
	/// </summary>
	public int FileFlags;

	/// <summary>
	/// Offset 12, 4 bytes. Total number of precursor records in the file.
	/// Duplicated in the footer (MslFooter.NPrecursors) for validation; a mismatch
	/// between the two indicates a truncated or corrupt write.
	/// </summary>
	public int NPrecursors;

	/// <summary>
	/// Offset 16, 4 bytes. Number of protein records in the protein table.
	/// Zero means the protein table is absent (FileFlags bit 1 should also be 0).
	/// </summary>
	public int NProteins;

	/// <summary>
	/// Offset 20, 4 bytes. Count of unique elution group IDs across all precursors.
	/// Precursors that share the same stripped sequence share an elution group ID.
	/// </summary>
	public int NElutionGroups;

	/// <summary>
	/// Offset 24, 4 bytes. Total number of entries in the string table.
	/// The string table stores all modified sequences, stripped sequences, protein
	/// accessions, names, and gene names in a flat array indexed by int32 index.
	/// </summary>
	public int NStrings;

	/// <summary>
	/// Offset 28, 4 bytes. Absolute byte offset of the extended annotation table section.
	/// 0 when <c>MslFormat.FileFlagHasExtAnnotations</c> is not set (section absent).
	/// When the flag is set, this field holds the offset at which the extended annotation
	/// table begins; readers must seek here to read custom neutral-loss masses.
	///
	/// Formerly the <c>Reserved</c> field in format version 1 (always written as 0).
	/// Repurposed in format version 2. Version-1 readers see this field as Reserved=0
	/// for all version-1 files and will reject version-2 files on the version check
	/// before reaching this field, so backward compatibility is fully maintained.
	/// </summary>
	public int ExtAnnotationTableOffset;

	/// <summary>
	/// Offset 32, 8 bytes. Absolute file byte offset of the first MslProteinRecord.
	/// Zero when NProteins == 0. Readers seeking to the protein table jump directly
	/// to this offset rather than scanning past the header.
	/// </summary>
	public long ProteinTableOffset;

	/// <summary>
	/// Offset 40, 8 bytes. Absolute file byte offset of the first entry in the string table.
	/// All string-index fields in protein and precursor records resolve against this section.
	/// </summary>
	public long StringTableOffset;

	/// <summary>
	/// Offset 48, 8 bytes. Absolute file byte offset of the first MslPrecursorRecord.
	/// The precursor section immediately follows the string table in the streaming writer;
	/// readers skip directly to this offset for indexed or full-load access.
	/// </summary>
	public long PrecursorSectionOffset;

	/// <summary>
	/// Offset 56, 8 bytes. Absolute file byte offset of the first MslFragmentRecord.
	/// Each precursor record embeds its own FragmentBlockOffset for O(1) random access;
	/// this field allows readers to locate the start of the combined fragment section
	/// without iterating precursors.
	/// </summary>
	public long FragmentSectionOffset;
}

// ──────────────────────────────────────────────────────────────────────────────
// MslProteinRecord — 24 bytes
// ──────────────────────────────────────────────────────────────────────────────

/// <summary>
/// A 24-byte record describing one protein in the protein table section of the .msl file.
/// Precursor records reference proteins by zero-based index into this table (ProteinIdx field).
/// All string fields are stored as int32 indices into the file's string table rather than
/// as inline character data; this avoids variable-length records and allows multiple precursors
/// to reference the same protein strings without duplication.
/// </summary>
[StructLayout(LayoutKind.Sequential, Pack = 1)]
public struct MslProteinRecord
{
	/// <summary>
	/// Offset 0, 4 bytes. Index into the file string table for the UniProt accession string
	/// (e.g. "P12345" or "Q9UBT6-2" for isoforms). 0 is a valid index; -1 signals absent.
	/// </summary>
	public int AccessionStringIdx;

	/// <summary>
	/// Offset 4, 4 bytes. Index into the file string table for the human-readable protein
	/// name (e.g. "Serum albumin"). 0 is a valid index; -1 signals absent.
	/// </summary>
	public int NameStringIdx;

	/// <summary>
	/// Offset 8, 4 bytes. Index into the file string table for the gene symbol
	/// (e.g. "ALB"). 0 = absent (gene data not available for this entry; file-level
	/// FileFlag has_gene_data may still be set for other proteins).
	/// </summary>
	public int GeneStringIdx;

	/// <summary>
	/// Offset 12, 4 bytes. Protein group assignment ID. Proteins that belong to the same
	/// indistinguishable group share a ProteinGroupId value; 0 = ungrouped or unknown.
	/// </summary>
	public int ProteinGroupId;

	/// <summary>
	/// Offset 16, 4 bytes. Number of precursor records that reference this protein via their
	/// ProteinIdx field. Writers populate this during the pre-pass; readers may use it to
	/// pre-size per-protein collections.
	/// </summary>
	public int NPrecursors;

	/// <summary>
	/// Offset 20, 4 bytes. Bitfield of protein-level flags.
	/// bit 0: is_reviewed — 1 if the accession is a UniProt Swiss-Prot (reviewed) entry.
	/// bits 1–31: reserved, must be 0.
	/// </summary>
	public int ProteinFlags;
}

// ──────────────────────────────────────────────────────────────────────────────
// MslPrecursorRecord — 56 bytes
// ──────────────────────────────────────────────────────────────────────────────

/// <summary>
/// A 56-byte record that stores all per-precursor data needed to reconstruct a
/// LibrarySpectrum (plus mzLib-extended metadata) without reading any fragment data.
/// The PrecursorMz field is the primary sort key for the in-memory query engine; the
/// FragmentBlockOffset field enables O(1) random access to the associated fragment ions
/// without scanning the entire fragment section.
/// </summary>
[StructLayout(LayoutKind.Sequential, Pack = 1)]
public struct MslPrecursorRecord
{
	/// <summary>
	/// Offset 0, 4 bytes. Precursor m/z (mass-to-charge ratio).
	/// Stored as float32 to match Prosit/DIA-NN precision; precision loss vs double is
	/// negligible for library matching at standard tolerances (≤ 10 ppm).
	/// Maps to LibrarySpectrum.PrecursorMz.
	/// </summary>
	public float PrecursorMz;

	/// <summary>
	/// Offset 4, 4 bytes. Indexed retention time (iRT) or calibrated run-specific RT in minutes.
	/// Whether this is iRT or calibrated RT is indicated by bit 2 of PrecursorFlags.
	/// Maps to LibrarySpectrum.RetentionTime.
	/// </summary>
	public float Irt;

	/// <summary>
	/// Offset 8, 4 bytes. Collisional cross section expressed as 1/K0 (ion mobility).
	/// 0.0f means the value is not available; the file-level has_ion_mobility flag indicates
	/// whether any precursor in the file has a meaningful non-zero value.
	/// </summary>
	public float IonMobility;

	/// <summary>
	/// Offset 12, 2 bytes. Precursor charge state (e.g. 1, 2, 3, 4).
	/// Maps to LibrarySpectrum.ChargeState.
	/// </summary>
	public short Charge;

	/// <summary>
	/// Offset 14, 2 bytes. Number of MslFragmentRecord entries in the fragment block
	/// pointed to by FragmentBlockOffset. Readers use this count to know exactly how
	/// many 20-byte records to read after seeking.
	/// </summary>
	public short FragmentCount;

	/// <summary>
	/// Offset 16, 4 bytes. Elution group identifier. All precursors with the same
	/// stripped sequence share an ElutionGroupId, enabling the DIA scorer to co-elute
	/// charge states and modifications together.
	/// </summary>
	public int ElutionGroupId;

	/// <summary>
	/// Offset 20, 4 bytes. Zero-based index into the protein table for the source protein.
	/// −1 means no protein assignment (e.g. synthetic library entries with no FASTA origin).
	/// Maps to MslProteinRecord at index ProteinIdx.
	/// </summary>
	public int ProteinIdx;

	/// <summary>
	/// Offset 24, 4 bytes. String table index for the modified sequence string
	/// (mzLib notation, e.g. "PEPTM[Common Variable:Oxidation on M]IDE").
	/// Maps to LibrarySpectrum.Sequence.
	/// </summary>
	public int ModifiedSeqStringIdx;

	/// <summary>
	/// Offset 28, 4 bytes. String table index for the stripped (unmodified) amino-acid sequence.
	/// Used for elution group construction and DDA-style lookup by sequence.
	/// </summary>
	public int StrippedSeqStringIdx;

	/// <summary>
	/// Offset 32, 8 bytes. Absolute file byte offset of the first MslFragmentRecord belonging
	/// to this precursor. Combined with FragmentCount this fully specifies the fragment block
	/// extent, enabling O(1) seeking without scanning other precursors' fragments.
	/// For compressed files (FileFlagIsCompressed), this is an offset into the decompressed
	/// fragment buffer rather than an absolute file position.
	/// </summary>
	public long FragmentBlockOffset;

	/// <summary>
	/// Offset 40, 4 bytes. Library q-value (confidence score, lower = better).
	/// float.NaN when no q-value is available (e.g. purely predicted library).
	/// </summary>
	public float QValue;

	/// <summary>
	/// Offset 44, 4 bytes. Residue count of the stripped sequence.
	/// Stored explicitly because computing it from the string requires an additional
	/// string-table lookup; top-down proteoforms can have values in the millions.
	/// </summary>
	public int StrippedSeqLength;

	/// <summary>
	/// Offset 48, 2 bytes. Molecule classification; cast from MslFormat.MoleculeType (int16).
	/// Controls which fragmentation namespace is used at read time to reconstruct
	/// Product.Terminus and Product.NeutralMass.
	/// </summary>
	public short MoleculeType;

	/// <summary>
	/// Offset 50, 2 bytes. Dissociation type used to generate the spectrum
	/// (e.g. HCD, CID, EThcD); cast from MassSpectrometry.DissociationType to int16.
	/// Needed by fragment mass re-computation at read time.
	/// </summary>
	public short DissociationType;

	/// <summary>
	/// Offset 52, 2 bytes. Nominal collision energy × 10, stored as int16.
	/// E.g. NCE 28 → stored as 280. 0 = unknown or not applicable.
	/// </summary>
	public short Nce;

	/// <summary>
	/// Offset 54, 1 byte. Bitfield of precursor-level Boolean properties.
	/// bit 0: is_decoy; bit 1: is_proteotypic; bit 2: rt_is_calibrated; bits 3–7: reserved.
	/// Use MslFormat.EncodePrecursorFlags / DecodePrecursorFlags to read or write this field.
	/// </summary>
	public byte PrecursorFlags;

	/// <summary>
	/// Offset 55, 1 byte. Origin of fragment intensities; cast from MslFormat.SourceType.
	/// Predicted = 0, Empirical = 1, EmpiricalRefined = 2.
	/// </summary>
	public byte SourceType;
}

// ──────────────────────────────────────────────────────────────────────────────
// MslFragmentRecord — 20 bytes
// ──────────────────────────────────────────────────────────────────────────────

/// <summary>
/// A 20-byte record describing one fragment ion.  Each MslPrecursorRecord's fragment
/// block consists of exactly FragmentCount contiguous MslFragmentRecord entries starting
/// at FragmentBlockOffset. Fields that are derivable at read time (Product.NeutralMass,
/// Product.Terminus, and named neutral-loss masses) are NOT stored; they are recomputed
/// from the stored ProductType, SecondaryProductType, and the precursor's MoleculeType
/// to keep the record size minimal.
/// </summary>
[StructLayout(LayoutKind.Sequential, Pack = 1)]
public struct MslFragmentRecord
{
	/// <summary>
	/// Offset 0, 4 bytes. Observed or predicted m/z of the fragment ion.
	/// Maps to MatchedFragmentIon.Mz.
	/// </summary>
	public float Mz;

	/// <summary>
	/// Offset 4, 4 bytes. Relative intensity of the fragment, normalized 0–1 within a precursor
	/// (1.0 = most abundant fragment). Maps to MatchedFragmentIon.Intensity.
	/// </summary>
	public float Intensity;

	/// <summary>
	/// Offset 8, 2 bytes. N-terminal fragment ion type (cast from Omics.Fragmentation.ProductType).
	/// For terminal ions this is the sole product type; for internal ions it is the N-terminal
	/// terminus type of the internal series. Cast: (int16)(int)product.ProductType.
	/// </summary>
	public short ProductType;

	/// <summary>
	/// Offset 10, 2 bytes. C-terminal fragment ion type for internal ions; −1 for terminal ions.
	/// Cast: (int16)(int)product.SecondaryProductType if set, otherwise −1.
	/// At read time: cast back to ProductType? — null when the stored value is −1.
	/// </summary>
	public short SecondaryProductType;

	/// <summary>
	/// Offset 12, 2 bytes. Ion series number for terminal ions (e.g. 5 for b5 or y5);
	/// or the start residue position (0-based) for internal fragment ions.
	/// Maps to Product.FragmentNumber.
	/// </summary>
	public short FragmentNumber;

	/// <summary>
	/// Offset 14, 2 bytes. End residue position for internal ions (0-based, exclusive upper bound);
	/// 0 for terminal ions. Maps to Product.SecondaryFragmentNumber.
	/// </summary>
	public short SecondaryFragmentNumber;

	/// <summary>
	/// Offset 16, 2 bytes. Residue position within the peptide sequence from which the
	/// fragment was assigned. Maps to Product.ResiduePosition.
	/// </summary>
	public short ResiduePosition;

	/// <summary>
	/// Offset 18, 1 byte. ChargeState state of the fragment ion (typically 1 or 2).
	/// Maps to MatchedFragmentIon.ChargeState. Values above 127 are extremely rare for fragments;
	/// a byte is sufficient.
	/// </summary>
	public byte Charge;

	/// <summary>
	/// Offset 19, 1 byte. Packed bitfield for fragment-level Boolean and enum properties.
	/// Bit 0: is_internal; bit 1: is_diagnostic; bits 2–4: neutral_loss_code (NeutralLossCode);
	/// bit 5: exclude_from_quant; bits 6–7: reserved.
	/// Use MslFormat.EncodeFragmentFlags / DecodeFragmentFlags to read or write this field.
	/// </summary>
	public byte Flags;
}

// ──────────────────────────────────────────────────────────────────────────────
// MslFooter — 20 bytes
// ──────────────────────────────────────────────────────────────────────────────

/// <summary>
/// The 20-byte file footer, occupying the last bytes of every .msl file.
/// The footer exists to support streaming writes: the offset table is written last and its
/// position within the file is recorded here. Readers should verify that NPrecursors matches
/// the value in the file header and that TrailingMagic matches MslFormat.Magic before trusting
/// any data from the file.
/// </summary>
[StructLayout(LayoutKind.Sequential, Pack = 1)]
public struct MslFooter
{
	/// <summary>
	/// Offset from EOF: −20 (first byte of footer), 8 bytes.
	/// Absolute file byte offset of the per-precursor offset table. The offset table is a
	/// contiguous array of int64 values, one per precursor, giving the absolute byte offset
	/// of each precursor's MslPrecursorRecord. This enables O(log n) binary search on m/z
	/// without loading the entire precursor section into memory.
	/// </summary>
	public long OffsetTableOffset;

	/// <summary>
	/// Offset from EOF: −12, 4 bytes. Total precursor count, duplicated from the file header.
	/// A mismatch between this value and MslFileHeader.NPrecursors indicates a truncated write
	/// and the file should be treated as corrupt.
	/// </summary>
	public int NPrecursors;

	/// <summary>
	/// Offset from EOF: −8, 4 bytes. CRC-32 checksum of all file bytes from offset 0 to
	/// (OffsetTableOffset − 1) inclusive. Covers the header, protein table, string table,
	/// precursor section, and fragment section but not the offset table or the footer itself.
	/// </summary>
	public uint DataCrc32;

	/// <summary>
	/// Offset from EOF: −4, 4 bytes. Trailing magic stored as a uint, identical in value to
	/// MslFileHeader.Magic ("MZLB"). Presence of the correct trailing magic confirms the file
	/// was fully written and not truncated mid-stream. Writers populate this with the same
	/// uint written to MslFileHeader.Magic; readers verify it matches before trusting the footer.
	/// </summary>
	public uint TrailingMagic;
}