namespace Readers.SpectralLibrary.DiaNNSpectralLibrary
{
	/// <summary>
	/// Constants and enumerations for the DIA-NN 2.3.2 .speclib binary format (version −10).
	///
	/// The format is NOT a struct array. Every precursor is located by scanning for one of
	/// six 4-byte marker values embedded in the file. All layout knowledge here is the result
	/// of empirical binary analysis documented in Session3–Session7 findings documents.
	///
	/// Key format characteristics:
	///   - Version int32 = −10 at file offset 0
	///   - Anchor-based precursor records (no fixed-size array)
	///   - 12 fixed fragment slots per precursor (1 null placeholder + up to 11 real)
	///   - Fragment bytes 0–2 are opaque DIA-NN internal indices — ignore them
	///   - Protein 0 has 2 opaque prefix int32 fields; proteins 1..N−1 have 4
	/// </summary>
	public static class DiaNNBinaryStructs
	{
		// ─── Version ────────────────────────────────────────────────────────────────

		/// <summary>
		/// The only confirmed .speclib version produced by DIA-NN 2.3.2.
		/// The first int32 in the file is always −10.
		/// </summary>
		public const int ExpectedVersion = -10;

		// ─── Marker byte patterns ────────────────────────────────────────────────────

		/// <summary>
		/// The six 4-byte marker patterns. Each marks a precursor record anchor point.
		/// The marker appears at file position (nlp − 8), where nlp is the byte offset
		/// of the name-length int32.
		///
		/// IMPORTANT: Do NOT write these as float literals. The bit patterns do not
		/// correspond to clean decimal values (e.g. Internal = 2.9375016689300537f,
		/// not 2.9375f). Using approximate float literals produces wrong bytes and
		/// the scanner will never find them. Always write all six markers as raw bytes.
		///
		/// Patterns ending in 0x40 are normal floats; patterns ending in 0x00 are
		/// subnormal floats (near zero). A scanner that only looks for 0x40 in the
		/// fourth byte will silently miss all shared-peptide precursors.
		/// </summary>
		public static readonly byte[][] MarkerPatterns = new byte[][]
		{
			new byte[] { 0x07, 0x00, 0x3C, 0x40 },   // Internal        (~2.9375)
            new byte[] { 0x07, 0x00, 0x39, 0x40 },   // N-terminal      (~2.8906)
            new byte[] { 0x07, 0x00, 0x36, 0x40 },   // C-terminal      (~2.8438)
            new byte[] { 0x07, 0x00, 0x33, 0x40 },   // N+C terminal    (~2.7969)
            new byte[] { 0x07, 0x00, 0x3D, 0x00 },   // Shared conflict (subnormal)
            new byte[] { 0x07, 0x00, 0x39, 0x00 },   // Shared same     (subnormal)
        };

		// ─── Terminus type ───────────────────────────────────────────────────────────

		/// <summary>
		/// Encodes the terminus position of a peptide within its source protein(s),
		/// as determined by the marker float value in the binary file.
		/// </summary>
		public enum TerminusType
		{
			/// <summary>Peptide is neither at the N- nor C-terminus of any source protein.</summary>
			Internal,

			/// <summary>Peptide is at the N-terminus of all source proteins.</summary>
			NTerminal,

			/// <summary>Peptide is at the C-terminus of all source proteins.</summary>
			CTerminal,

			/// <summary>Peptide is simultaneously the N- and C-terminus (single-peptide protein).</summary>
			NAndCTerminal,

			/// <summary>
			/// Shared peptide whose terminus type conflicts across source proteins
			/// (e.g., N-terminal in protein A but internal in protein B).
			/// Marker bytes: 07 00 3D 00 (subnormal float — do NOT write as a float literal).
			/// </summary>
			SharedConflict,

			/// <summary>
			/// Shared peptide with the same terminus type in all source proteins
			/// (e.g., N-terminal in every protein that contains it).
			/// Marker bytes: 07 00 39 00 (subnormal float — do NOT write as a float literal).
			/// </summary>
			SharedSameType,
		}

		/// <summary>
		/// Returns the <see cref="TerminusType"/> encoded by four marker bytes.
		/// Returns null if the bytes do not match any known marker pattern.
		/// </summary>
		public static TerminusType? ClassifyMarker(byte b0, byte b1, byte b2, byte b3)
		{
			if (b0 != 0x07 || b1 != 0x00) return null;

			return (b2, b3) switch
			{
				(0x3C, 0x40) => TerminusType.Internal,
				(0x39, 0x40) => TerminusType.NTerminal,
				(0x36, 0x40) => TerminusType.CTerminal,
				(0x33, 0x40) => TerminusType.NAndCTerminal,
				(0x3D, 0x00) => TerminusType.SharedConflict,
				(0x39, 0x00) => TerminusType.SharedSameType,
				_ => null
			};
		}

		// ─── Pre-name block offsets (relative to nlp) ───────────────────────────────

		/// <summary>Offset from nlp to the precursor_global_index int32 field.</summary>
		public const int OffsetGlobalIndex = -48;

		/// <summary>Offset from nlp to the charge int32 field.</summary>
		public const int OffsetCharge = -44;

		/// <summary>Offset from nlp to the stripped_sequence_length int32 field.</summary>
		public const int OffsetStrippedSeqLen = -40;

		/// <summary>Offset from nlp to the precursor_mz float field.</summary>
		public const int OffsetPrecursorMz = -36;

		/// <summary>Offset from nlp to the iRT float field.</summary>
		public const int OffsetIRT = -32;

		/// <summary>Offset from nlp to the ion_mobility float field.</summary>
		public const int OffsetIonMobility = -28;

		/// <summary>Offset from nlp to the marker float field (4 bytes, one of the six patterns).</summary>
		public const int OffsetMarker = -8;

		/// <summary>Offset from nlp to the protein_group_index int32 field.</summary>
		public const int OffsetProteinGroupIndex = -4;

		// ─── Post-name block offsets (relative to name_end = nlp + 4 + name_length) ─

		/// <summary>Offset from name_end to the n_fragments int32 field.</summary>
		public const int OffsetNFragments = 16;

		/// <summary>Offset from name_end to the top_fragment_mz float field.</summary>
		public const int OffsetTopFragmentMz = 20;

		/// <summary>Offset from name_end to the top_fragment_intensity float field.</summary>
		public const int OffsetTopFragmentIntensity = 24;

		/// <summary>Total post-name block size in bytes (before fragment records begin).</summary>
		public const int PostNameBlockSize = 28;

		// ─── Fragment record layout (12 bytes each) ───────────────────────────────────

		/// <summary>Size in bytes of one fragment record.</summary>
		public const int FragmentRecordSize = 12;

		/// <summary>Total number of fragment slots per precursor (always 12, including 1 null).</summary>
		public const int FragmentSlotsPerPrecursor = 12;

		/// <summary>Byte offset within a fragment record for fragment_mz (float).</summary>
		public const int FragmentOffsetMz = 4;

		/// <summary>Byte offset within a fragment record for fragment_intensity (float).</summary>
		public const int FragmentOffsetIntensity = 8;

		// ─── Name validation bounds ───────────────────────────────────────────────────

		/// <summary>Minimum valid name_length value (shortest possible peptide name).</summary>
		public const int NameLengthMin = 5;

		/// <summary>
		/// Maximum valid name_length value. Long missed-cleavage peptides with multiple
		/// modifications can exceed 60 characters; 300 is a safe upper bound.
		/// </summary>
		public const int NameLengthMax = 300;

		// ─── Precursor validation bounds ──────────────────────────────────────────────

		/// <summary>Maximum plausible global_index (upper validation bound).</summary>
		public const int GlobalIndexMax = 100_000;

		/// <summary>Minimum plausible precursor m/z (Da).</summary>
		public const float PrecursorMzMin = 100.0f;

		/// <summary>Maximum plausible precursor m/z (Da).</summary>
		public const float PrecursorMzMax = 5000.0f;

		/// <summary>Minimum valid n_fragments value.</summary>
		public const int NFragmentsMin = 1;

		/// <summary>Maximum valid n_fragments value (generous upper bound).</summary>
		public const int NFragmentsMax = 50;

		// ─── Writer constants ─────────────────────────────────────────────────────────

		/// <summary>
		/// Opaque bytes written for a real fragment record (bytes 0–3).
		/// These are DIA-NN internal indices; safe placeholder values are used on write.
		/// </summary>
		public static readonly byte[] RealFragmentOpaqueBytes = { 1, 2, 1, 0 };

		/// <summary>
		/// Opaque bytes written for a null (placeholder) fragment record (bytes 0–3).
		/// </summary>
		public static readonly byte[] NullFragmentOpaqueBytes = { 1, 1, 1, 0 };
	}
}