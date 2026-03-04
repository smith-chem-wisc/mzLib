using System.Runtime.InteropServices;

namespace Readers.SpectralLibrary.DiaNNSpectralLibrary
{
    /// <summary>
    /// Binary struct definitions for DIA-NN .speclib format.
    /// 
    /// IMPORTANT: These layouts are ESTIMATES based on:
    /// - DIA-NN TSV export schema
    /// - ProteoWizard BiblioSpec reverse-engineered reader (PR #1097)
    /// - DIA-NN GitHub issues (#77, #1184, #1807)
    /// 
    /// The actual byte layouts MUST be validated against real .speclib binary files.
    /// The structs use Pack=1 to match C++ struct packing with #pragma pack(1).
    /// All multi-byte fields are little-endian (x86/x64 native).
    /// 
    /// Version history (from GitHub issues and BiblioSpec reader):
    ///   Version -2 : Very early legacy format
    ///   Version -1 : First versioned format (int32 = 0xFFFFFFFF at offset 0)
    ///   Version  1 : Added structured header
    ///   Version  2 : Added ion mobility fields
    ///   Version  3 : Added scores, additional metadata
    ///   Version  8 : DIA-NN 2.0+, major restructuring
    /// </summary>
    public static class DiaNNBinaryStructs
    {
        /// <summary>
        /// Known .speclib format versions. The version int32 is the first 4 bytes of the file.
        /// Negative values were used in early versions; positive values in later ones.
        /// </summary>
        public enum SpecLibVersion : int
        {
            /// <summary>Very early legacy format, no explicit version field</summary>
            LegacyV2 = -2,
            /// <summary>First versioned format (pre-1.8)</summary>
            V1Legacy = -1,
            /// <summary>~DIA-NN 1.7, added structured header</summary>
            V1 = 1,
            /// <summary>~DIA-NN 1.8, added ion mobility</summary>
            V2 = 2,
            /// <summary>~DIA-NN 1.8.1, added scores and metadata</summary>
            V3 = 3,
            /// <summary>DIA-NN 2.0+, current format with major restructuring</summary>
            V8 = 8,
        }

        /// <summary>
        /// File header for .speclib files. The first field identifies the format version.
        /// Additional header fields are version-dependent and will be populated during
        /// validation against real binary files (Prompt 8).
        /// </summary>
        [StructLayout(LayoutKind.Sequential, Pack = 1)]
        public struct SpecLibHeader
        {
            /// <summary>
            /// Format version. First 4 bytes of the file.
            /// Values: -2, -1, 1, 2, 3, 8 (see SpecLibVersion enum).
            /// NOT a SQLite magic number — this was confirmed by Skyline developers.
            /// </summary>
            public int Version;

            /// <summary>
            /// Whether the library contains generated decoys.
            /// 0 = no decoys, 1 = decoys included.
            /// </summary>
            public int GenDecoys;

            /// <summary>
            /// Total number of precursor entries in the library.
            /// Used to pre-allocate arrays during reading.
            /// </summary>
            public int PrecursorCount;

            // NOTE: Additional version-dependent header fields will be added here
            // after validation with real .speclib files. Possible fields include:
            //   - FragmentCount (total fragments across all precursors)
            //   - StringTableOffset (offset to the sequence/protein string table)
            //   - StringTableSize
            //   - Flags (e.g., whether ion mobility is stored)
        }

        /// <summary>
        /// Precursor-level record in the .speclib binary format.
        /// Each precursor represents one peptide species at a specific charge state.
        /// 
        /// The fragment data for this precursor starts at FragmentOffset in the file
        /// and contains FragmentCount consecutive DiaNNFragmentRecord structs.
        /// 
        /// The modified sequence is stored in a separate string table at SequenceOffset
        /// with SequenceLength bytes.
        /// 
        /// ESTIMATED LAYOUT — must be validated against real binary files.
        /// </summary>
        [StructLayout(LayoutKind.Sequential, Pack = 1)]
        public struct PrecursorRecord
        {
            /// <summary>Precursor m/z value (single-precision for compactness)</summary>
            public float PrecursorMz;

            /// <summary>Precursor charge state (1-6 typical)</summary>
            public short Charge;

            /// <summary>Retention time: iRT value or calibrated RT in minutes</summary>
            public float RetentionTime;

            /// <summary>Ion mobility value (1/K0 for timsTOF). 0 if not available.</summary>
            public float IonMobility;

            /// <summary>Byte offset into the string table for the modified sequence</summary>
            public int SequenceOffset;

            /// <summary>Length of the modified sequence string in bytes</summary>
            public short SequenceLength;

            /// <summary>Byte offset into the protein name string table</summary>
            public int ProteinOffset;

            /// <summary>Byte offset into the gene name string table</summary>
            public int GeneOffset;

            /// <summary>Byte offset in the file where this precursor's fragment data begins</summary>
            public long FragmentOffset;

            /// <summary>Number of fragment ions for this precursor</summary>
            public short FragmentCount;

            /// <summary>Decoy flag: 0 = target, 1 = decoy</summary>
            public byte IsDecoy;

            /// <summary>Proteotypic flag: 1 = unique to one protein, 0 = shared</summary>
            public byte IsProteotypic;

            /// <summary>Library-level q-value, if available. NaN or 0 if not set.</summary>
            public float QValue;
        }

        /// <summary>
        /// Fragment ion record in the .speclib binary format.
        /// Stored contiguously for each precursor, sorted by m/z.
        /// 
        /// The compact 12-byte layout enables efficient cache utilization:
        /// a 64-byte cache line holds 5 complete fragment records.
        /// 
        /// ESTIMATED LAYOUT — must be validated against real binary files.
        /// </summary>
        [StructLayout(LayoutKind.Sequential, Pack = 1)]
        public struct FragmentRecord
        {
            /// <summary>Fragment m/z value (single-precision)</summary>
            public float Mz;

            /// <summary>Relative intensity, pre-normalized to [0, 1] where 1 = most abundant</summary>
            public float Intensity;

            /// <summary>
            /// Ion series type encoded as byte:
            ///   0 = b-ion (N-terminal)
            ///   1 = y-ion (C-terminal)
            /// Additional types may be added in future versions.
            /// </summary>
            public byte IonType;

            /// <summary>Fragment series number (1-based, e.g., b3 = 3, y7 = 7)</summary>
            public byte SeriesNumber;

            /// <summary>Fragment charge state (typically 1 or 2)</summary>
            public byte Charge;

            /// <summary>
            /// Neutral loss type encoded as byte:
            ///   0 = no loss
            ///   1 = H2O loss (-18.0106 Da)
            ///   2 = NH3 loss (-17.0265 Da)
            ///   3 = H3PO4 loss (-97.9769 Da)
            /// </summary>
            public byte LossType;
        }

        /// <summary>
        /// Ion type byte values used in FragmentRecord.IonType.
        /// Maps to mzLib ProductType enum values.
        /// </summary>
        public static class IonTypeEncoding
        {
            public const byte B = 0;
            public const byte Y = 1;

            /// <summary>Convert a character ion type ('b', 'y') to the byte encoding</summary>
            public static byte FromChar(char ionType) => ionType switch
            {
                'b' or 'B' => B,
                'y' or 'Y' => Y,
                _ => throw new ArgumentException($"Unknown ion type: '{ionType}'. Expected 'b' or 'y'.")
            };

            /// <summary>Convert the byte encoding back to a character</summary>
            public static char ToChar(byte encoded) => encoded switch
            {
                B => 'b',
                Y => 'y',
                _ => throw new ArgumentException($"Unknown ion type encoding: {encoded}. Expected 0 (b) or 1 (y).")
            };

            /// <summary>Convert the byte encoding to the mzLib ProductType string name</summary>
            public static string ToProductTypeName(byte encoded) => encoded switch
            {
                B => "b",
                Y => "y",
                _ => throw new ArgumentException($"Unknown ion type encoding: {encoded}.")
            };
        }

        /// <summary>
        /// Neutral loss type byte values used in FragmentRecord.LossType.
        /// Maps to DIA-NN's string-based neutral loss names in TSV format.
        /// </summary>
        public static class LossTypeEncoding
        {
            public const byte NoLoss = 0;
            public const byte H2O = 1;
            public const byte NH3 = 2;
            public const byte H3PO4 = 3;
            public const byte HPO3 = 4;

            /// <summary>Convert a DIA-NN loss type string to the byte encoding</summary>
            public static byte FromString(string lossType)
            {
                if (string.IsNullOrEmpty(lossType) || lossType.Equals("noloss", StringComparison.OrdinalIgnoreCase))
                    return NoLoss;

                return lossType.ToUpperInvariant() switch
                {
                    "H2O" => H2O,
                    "NH3" => NH3,
                    "H3PO4" => H3PO4,
                    "HPO3" => HPO3,
                    _ => NoLoss // Default to no loss for unknown types
                };
            }

            /// <summary>Convert the byte encoding back to DIA-NN loss type string</summary>
            public static string ToString(byte encoded) => encoded switch
            {
                NoLoss => "noloss",
                H2O => "H2O",
                NH3 => "NH3",
                H3PO4 => "H3PO4",
                HPO3 => "HPO3",
                _ => "noloss"
            };

            /// <summary>Convert the byte encoding to a neutral loss mass in Daltons</summary>
            public static double ToMass(byte encoded) => encoded switch
            {
                NoLoss => 0.0,
                H2O => 18.010565,
                NH3 => 17.026549,
                H3PO4 => 97.976896,
                HPO3 => 79.966331,
                _ => 0.0
            };
        }

        /// <summary>
        /// Size constants for binary I/O. These reflect the ESTIMATED struct sizes
        /// and will be validated against real files.
        /// </summary>
        public static class SizeOf
        {
            /// <summary>Size of SpecLibHeader in bytes</summary>
            public static readonly int Header = Marshal.SizeOf<SpecLibHeader>();

            /// <summary>Size of PrecursorRecord in bytes</summary>
            public static readonly int Precursor = Marshal.SizeOf<PrecursorRecord>();

            /// <summary>Size of FragmentRecord in bytes</summary>
            public static readonly int Fragment = Marshal.SizeOf<FragmentRecord>();
        }
    }
}
