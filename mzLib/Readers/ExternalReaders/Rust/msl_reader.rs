//! msl_reader.rs — Reference reader for the .msl binary spectral library format.
//!
//! Implements the format specified in FinalProjectSummary.md Section 4.
//! No dependencies outside the Rust standard library.
//!
//! # Format overview
//!
//! ```text
//! File Header       64 bytes  (magic, version, counts, section offsets)
//! Protein Table     NProteins × 24 bytes
//! String Table      variable  (int32 count, int32 TotalBytes, length-prefixed UTF-8 entries)
//! Precursor Section NPrecursors × 56 bytes
//! Fragment Section  variable  (FragmentCount × 20 bytes per precursor)
//! Offset Table      NPrecursors × 8 bytes (int64[])
//! Footer            20 bytes  (last 20 bytes of file)
//! ```
//!
//! # Validation order (Section 4.10)
//!
//! 1. Header magic == `[0x4D, 0x5A, 0x4C, 0x42]`
//! 2. `FormatVersion` in `[1, 3]`
//! 3. Trailing footer magic (read as LE uint32) == `0x4D5A4C42`
//! 4. `footer.NPrecursors == header.NPrecursors`
//! 5. CRC-32/ISO-HDLC over bytes `[0, footer.OffsetTableOffset)` == `footer.DataCrc32`
//!
//! # Notes
//!
//! - Header magic: raw bytes compared as `[0x4D, 0x5A, 0x4C, 0x42]`.
//! - Footer trailing magic: read as LE uint32 and compared to `0x4D5A4C42u32`.
//!   On disk the bytes are `[0x42, 0x4C, 0x5A, 0x4D]` (LE representation).
//! - String table: after `NStrings` (int32) comes a `TotalBytes` (int32) field
//!   that is skipped (spec omission documented in Prompt 19 Note 7.3).
//! - NCE field: stored as `actual_nce × 10` on disk; divided by 10 on read.
//! - QValue: IEEE 754 NaN means unavailable; use `q_value.is_nan()` to detect.

use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::string::FromUtf8Error;

// ─────────────────────────────────────────────────────────────────────────────
// Constants
// ─────────────────────────────────────────────────────────────────────────────

const HEADER_MAGIC: [u8; 4] = [0x4D, 0x5A, 0x4C, 0x42];
const FOOTER_MAGIC_U32: u32 = 0x4D5A4C42;
const MIN_SUPPORTED_VERSION: i32 = 1;
const MAX_SUPPORTED_VERSION: i32 = 3;

const HEADER_SIZE: usize    = 64;
const PROTEIN_SIZE: usize   = 24;
const PRECURSOR_SIZE: usize = 56;
const FRAGMENT_SIZE: usize  = 20;
const FOOTER_SIZE: usize    = 20;

/// Neutral loss masses in daltons, indexed by 3-bit code (Section 4.11).
const NEUTRAL_LOSS_MASSES: [f64; 7] = [
    0.0,          // 0 — None
    -18.010565,   // 1 — H2O
    -17.026549,   // 2 — NH3
    -97.976895,   // 3 — H3PO4
    -79.966331,   // 4 — HPO3
    -115.987460,  // 5 — H3PO4 + H2O (PlusH2O)
    0.0,          // 6 — Custom (mass stored externally)
];

// ─────────────────────────────────────────────────────────────────────────────
// Error type
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug)]
pub enum MslError {
    Io(io::Error),
    InvalidMagic,
    UnsupportedVersion(i32),
    TrailingMagicMismatch,
    PrecursorCountMismatch { header: i32, footer: i32 },
    CrcMismatch { expected: u32, computed: u32 },
    InvalidUtf8(FromUtf8Error),
    InvalidData(String),
}

impl std::fmt::Display for MslError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MslError::Io(e)                    => write!(f, "I/O error: {}", e),
            MslError::InvalidMagic             => write!(f, "Invalid .msl file magic"),
            MslError::UnsupportedVersion(v)    => write!(f, "Unsupported .msl format version: {}", v),
            MslError::TrailingMagicMismatch    => write!(f, "Trailing footer magic mismatch (file may be truncated)"),
            MslError::PrecursorCountMismatch { header, footer } =>
                write!(f, "Precursor count mismatch: header={}, footer={}", header, footer),
            MslError::CrcMismatch { expected, computed } =>
                write!(f, "CRC-32 mismatch: stored=0x{:08X}, computed=0x{:08X}", expected, computed),
            MslError::InvalidUtf8(e)           => write!(f, "Invalid UTF-8 in string table: {}", e),
            MslError::InvalidData(msg)         => write!(f, "Invalid data: {}", msg),
        }
    }
}

impl From<io::Error> for MslError {
    fn from(e: io::Error) -> Self {
        MslError::Io(e)
    }
}

impl From<FromUtf8Error> for MslError {
    fn from(e: FromUtf8Error) -> Self {
        MslError::InvalidUtf8(e)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Data model
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct MslFragment {
    pub mz: f32,
    pub intensity: f32,
    pub product_type: i16,
    /// `None` when the on-disk value is `-1` (terminal ion).
    pub secondary_product_type: Option<i16>,
    pub fragment_number: i16,
    pub secondary_fragment_number: i16,
    pub residue_position: i16,
    pub charge: u8,
    pub is_internal: bool,
    pub is_diagnostic: bool,
    /// Raw 3-bit neutral loss code (0–6).
    pub neutral_loss_code: u8,
    /// Corresponding mass in Da from the lookup table.
    pub neutral_loss_mass: f64,
    pub exclude_from_quant: bool,
}

#[derive(Debug, Clone)]
pub struct MslPrecursor {
    pub precursor_mz: f32,
    pub irt: f32,
    pub ion_mobility: f32,
    pub charge: i16,
    pub fragment_count: i16,
    pub elution_group_id: i32,
    pub protein_idx: i32,
    pub modified_sequence: String,
    pub stripped_sequence: String,
    pub fragment_block_offset: i64,
    /// IEEE 754 NaN when unavailable.
    pub q_value: f32,
    pub stripped_seq_length: i32,
    pub molecule_type: i16,
    pub dissociation_type: i16,
    /// Actual NCE value (already ÷10 from on-disk × 10 encoding).
    pub nce: i16,
    pub is_decoy: bool,
    pub is_proteotypic: bool,
    pub rt_is_calibrated: bool,
    /// 0=Predicted, 1=Empirical, 2=EmpiricalRefined.
    pub source_type: u8,
    /// Resolved from protein table; empty string if no protein assigned.
    pub protein_accession: String,
    pub protein_name: String,
    pub gene_name: String,
    /// Empty in index-only mode.
    pub fragments: Vec<MslFragment>,
}

#[derive(Debug)]
pub struct MslLibrary {
    pub format_version: i32,
    pub n_precursors: i32,
    pub n_proteins: i32,
    pub has_ion_mobility: bool,
    pub has_protein_data: bool,
    pub has_gene_data: bool,
    pub is_predicted: bool,
    pub precursors: Vec<MslPrecursor>,
}

// ─────────────────────────────────────────────────────────────────────────────
// CRC-32/ISO-HDLC (no external crate)
// ─────────────────────────────────────────────────────────────────────────────

const CRC32_POLY: u32 = 0xEDB8_8320;

fn build_crc32_table() -> [u32; 256] {
    let mut table = [0u32; 256];
    for i in 0u32..256 {
        let mut entry = i;
        for _ in 0..8 {
            entry = if entry & 1 != 0 {
                (entry >> 1) ^ CRC32_POLY
            } else {
                entry >> 1
            };
        }
        table[i as usize] = entry;
    }
    table
}

fn compute_crc32(data: &[u8]) -> u32 {
    let table = build_crc32_table();
    let mut crc: u32 = 0xFFFF_FFFF;
    for &byte in data {
        crc = (crc >> 8) ^ table[((crc ^ byte as u32) & 0xFF) as usize];
    }
    crc ^ 0xFFFF_FFFF
}

/// Computes CRC-32/ISO-HDLC by streaming `len` bytes from the current file
/// position, without buffering the entire range in memory.
fn compute_crc32_streaming(file: &mut File, len: usize) -> io::Result<u32> {
    let table = build_crc32_table();
    let mut crc: u32 = 0xFFFF_FFFF;
    let mut remaining = len;
    let mut buf = [0u8; 65536];
    while remaining > 0 {
        let to_read = remaining.min(buf.len());
        file.read_exact(&mut buf[..to_read])?;
        for &byte in &buf[..to_read] {
            crc = (crc >> 8) ^ table[((crc ^ byte as u32) & 0xFF) as usize];
        }
        remaining -= to_read;
    }
    Ok(crc ^ 0xFFFF_FFFF)
}

// ─────────────────────────────────────────────────────────────────────────────
// Little-endian reading helpers
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn read_i32_le(buf: &[u8], offset: usize) -> i32 {
    i32::from_le_bytes(buf[offset..offset + 4].try_into().unwrap())
}

#[inline]
fn read_i64_le(buf: &[u8], offset: usize) -> i64 {
    i64::from_le_bytes(buf[offset..offset + 8].try_into().unwrap())
}

#[inline]
fn read_f32_le(buf: &[u8], offset: usize) -> f32 {
    f32::from_le_bytes(buf[offset..offset + 4].try_into().unwrap())
}

#[inline]
fn read_i16_le(buf: &[u8], offset: usize) -> i16 {
    i16::from_le_bytes(buf[offset..offset + 2].try_into().unwrap())
}

#[inline]
fn read_u32_le(buf: &[u8], offset: usize) -> u32 {
    u32::from_le_bytes(buf[offset..offset + 4].try_into().unwrap())
}

// ─────────────────────────────────────────────────────────────────────────────
// Flag decoding
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn decode_file_flags(flags: i32) -> (bool, bool, bool, bool) {
    (
        flags & 0x01 != 0, // has_ion_mobility
        flags & 0x02 != 0, // has_protein_data
        flags & 0x04 != 0, // has_gene_data
        flags & 0x08 != 0, // is_predicted
    )
}

#[inline]
fn decode_precursor_flags(flags: u8) -> (bool, bool, bool) {
    (
        flags & 0x01 != 0, // is_decoy
        flags & 0x02 != 0, // is_proteotypic
        flags & 0x04 != 0, // rt_is_calibrated
    )
}

/// Returns (is_internal, is_diagnostic, neutral_loss_code, exclude_from_quant).
#[inline]
fn decode_fragment_flags(flags: u8) -> (bool, bool, u8, bool) {
    (
        flags & 0x01 != 0,        // is_internal
        flags & 0x02 != 0,        // is_diagnostic
        (flags >> 2) & 0x07,      // neutral_loss_code (3 bits)
        flags & 0x20 != 0,        // exclude_from_quant
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// Parsing helpers
// ─────────────────────────────────────────────────────────────────────────────

fn parse_fragment(buf: &[u8]) -> MslFragment {
    // Fragment record layout (20 bytes, Section 4.6):
    //   0:  f32 Mz
    //   4:  f32 Intensity
    //   8:  i16 ProductType
    //  10:  i16 SecondaryProductType  (-1 = absent)
    //  12:  i16 FragmentNumber
    //  14:  i16 SecondaryFragmentNumber
    //  16:  i16 ResiduePosition
    //  18:  u8  Charge
    //  19:  u8  Flags
    let mz                          = read_f32_le(buf, 0);
    let intensity                   = read_f32_le(buf, 4);
    let product_type                = read_i16_le(buf, 8);
    let secondary_product_type_raw  = read_i16_le(buf, 10);
    let fragment_number             = read_i16_le(buf, 12);
    let secondary_fragment_number   = read_i16_le(buf, 14);
    let residue_position            = read_i16_le(buf, 16);
    let charge                      = buf[18];
    let flags_byte                  = buf[19];

    let secondary_product_type = if secondary_product_type_raw == -1 {
        None
    } else {
        Some(secondary_product_type_raw)
    };

    let (is_internal, is_diagnostic, nl_code, exclude_from_quant) =
        decode_fragment_flags(flags_byte);

    let nl_mass = if (nl_code as usize) < NEUTRAL_LOSS_MASSES.len() {
        NEUTRAL_LOSS_MASSES[nl_code as usize]
    } else {
        0.0
    };

    MslFragment {
        mz,
        intensity,
        product_type,
        secondary_product_type,
        fragment_number,
        secondary_fragment_number,
        residue_position,
        charge,
        is_internal,
        is_diagnostic,
        neutral_loss_code: nl_code,
        neutral_loss_mass: nl_mass,
        exclude_from_quant,
    }
}

fn read_fragments(data: &[u8], block_offset: i64, count: i16)
    -> Result<Vec<MslFragment>, MslError>
{
    let start  = block_offset as usize;
    let end    = start + (count as usize) * FRAGMENT_SIZE;
    if end > data.len() {
        return Err(MslError::InvalidData(format!(
            "Fragment block at offset {} size {} exceeds file length {}",
            block_offset, (count as usize) * FRAGMENT_SIZE, data.len()
        )));
    }
    let mut frags = Vec::with_capacity(count as usize);
    for i in 0..count as usize {
        let off = start + i * FRAGMENT_SIZE;
        frags.push(parse_fragment(&data[off..off + FRAGMENT_SIZE]));
    }
    Ok(frags)
}

fn read_string_table(data: &[u8], offset: i64, n_strings: i32)
    -> Result<Vec<String>, MslError>
{
    let off = offset as usize;
    // On-disk: int32 NStrings, int32 TotalBytes (skip), then entries.
    let n_read = read_i32_le(data, off);
    if n_read != n_strings {
        return Err(MslError::InvalidData(format!(
            "String table NStrings mismatch: header={}, table={}",
            n_strings, n_read
        )));
    }
    // Skip TotalBytes field at off+4
    let mut pos = off + 8; // past NStrings (4) and TotalBytes (4)

    let mut strings = Vec::with_capacity(n_strings as usize);
    for _ in 0..n_strings {
        let length = read_i32_le(data, pos) as usize;
        pos += 4;
        if pos + length > data.len() {
            return Err(MslError::InvalidData(
                "String table truncated while reading string data".to_string()
            ));
        }
        let s = String::from_utf8(data[pos..pos + length].to_vec())?;
        strings.push(s);
        pos += length;
    }
    Ok(strings)
}

fn read_protein_table(data: &[u8], offset: i64, n_proteins: i32, strings: &[String])
    -> Result<Vec<(String, String, String)>, MslError>
{
    // Protein record layout (24 bytes, Section 4.3):
    //   0:  i32 AccessionStringIdx
    //   4:  i32 NameStringIdx
    //   8:  i32 GeneStringIdx
    //  12:  i32 ProteinGroupId
    //  16:  i32 NPrecursors
    //  20:  i32 ProteinFlags
    let off = offset as usize;
    let total = n_proteins as usize * PROTEIN_SIZE;
    if off + total > data.len() {
        return Err(MslError::InvalidData("Protein table truncated".to_string()));
    }

    let resolve = |idx: i32| -> String {
        if idx >= 0 && (idx as usize) < strings.len() {
            strings[idx as usize].clone()
        } else {
            String::new()
        }
    };

    let mut proteins = Vec::with_capacity(n_proteins as usize);
    for i in 0..n_proteins as usize {
        let base = off + i * PROTEIN_SIZE;
        let acc_idx  = read_i32_le(data, base);
        let name_idx = read_i32_le(data, base + 4);
        let gene_idx = read_i32_le(data, base + 8);
        proteins.push((resolve(acc_idx), resolve(name_idx), resolve(gene_idx)));
    }
    Ok(proteins)
}

fn parse_precursor(buf: &[u8], strings: &[String],
                   proteins: &[(String, String, String)])
    -> MslPrecursor
{
    // Precursor record layout (56 bytes, Section 4.5):
    //   0:  f32 PrecursorMz
    //   4:  f32 Irt
    //   8:  f32 IonMobility
    //  12:  i16 Charge
    //  14:  i16 FragmentCount
    //  16:  i32 ElutionGroupId
    //  20:  i32 ProteinIdx
    //  24:  i32 ModifiedSeqStringIdx
    //  28:  i32 StrippedSeqStringIdx
    //  32:  i64 FragmentBlockOffset
    //  40:  f32 QValue
    //  44:  i32 StrippedSeqLength
    //  48:  i16 MoleculeType
    //  50:  i16 DissociationType
    //  52:  i16 Nce
    //  54:  u8  PrecursorFlags
    //  55:  u8  SourceType
    let precursor_mz        = read_f32_le(buf, 0);
    let irt                 = read_f32_le(buf, 4);
    let ion_mobility        = read_f32_le(buf, 8);
    let charge              = read_i16_le(buf, 12);
    let fragment_count      = read_i16_le(buf, 14);
    let elution_group_id    = read_i32_le(buf, 16);
    let protein_idx         = read_i32_le(buf, 20);
    let mod_seq_idx         = read_i32_le(buf, 24);
    let stripped_seq_idx    = read_i32_le(buf, 28);
    let fragment_block_offset = read_i64_le(buf, 32);
    let q_value             = read_f32_le(buf, 40);
    let stripped_seq_length = read_i32_le(buf, 44);
    let molecule_type       = read_i16_le(buf, 48);
    let dissociation_type   = read_i16_le(buf, 50);
    let nce_raw             = read_i16_le(buf, 52);
    let precursor_flags     = buf[54];
    let source_type         = buf[55];

    let (is_decoy, is_proteotypic, rt_is_calibrated) =
        decode_precursor_flags(precursor_flags);

    let resolve_str = |idx: i32| -> String {
        if idx >= 0 && (idx as usize) < strings.len() {
            strings[idx as usize].clone()
        } else {
            String::new()
        }
    };

    let modified_sequence = resolve_str(mod_seq_idx);
    let stripped_sequence = resolve_str(stripped_seq_idx);

    let (protein_accession, protein_name, gene_name) =
        if protein_idx >= 0 && (protein_idx as usize) < proteins.len() {
            proteins[protein_idx as usize].clone()
        } else {
            (String::new(), String::new(), String::new())
        };

    MslPrecursor {
        precursor_mz,
        irt,
        ion_mobility,
        charge,
        fragment_count,
        elution_group_id,
        protein_idx,
        modified_sequence,
        stripped_sequence,
        fragment_block_offset,
        q_value,
        stripped_seq_length,
        molecule_type,
        dissociation_type,
        nce: nce_raw / 10, // stored as actual_nce × 10
        is_decoy,
        is_proteotypic,
        rt_is_calibrated,
        source_type,
        protein_accession,
        protein_name,
        gene_name,
        fragments: Vec::new(),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Core loader (used by both load and load_index_only)
// ─────────────────────────────────────────────────────────────────────────────

fn load_impl(path: &str, load_fragment_data: bool)
    -> Result<MslLibrary, MslError>
{
    let mut file = File::open(path)?;
    let file_size = file.seek(SeekFrom::End(0))? as usize;

    if file_size < HEADER_SIZE + FOOTER_SIZE {
        return Err(MslError::InvalidData("File too short to be a valid .msl".to_string()));
    }

    // ── 1. Read header (64 bytes) ─────────────────────────────────────────────
    file.seek(SeekFrom::Start(0))?;
    let mut header_buf = [0u8; HEADER_SIZE];
    file.read_exact(&mut header_buf)?;

    if &header_buf[0..4] != &HEADER_MAGIC {
        return Err(MslError::InvalidMagic);
    }

    // Header layout (64 bytes, Section 4.2):
    //   0:  u8[4]  Magic
    //   4:  i32    FormatVersion
    //   8:  i32    FileFlags
    //  12:  i32    NPrecursors
    //  16:  i32    NProteins
    //  20:  i32    NElutionGroups
    //  24:  i32    NStrings
    //  28:  i32    Reserved (v1) / ExtAnnotationTableOffset (v2+)
    //  32:  i64    ProteinTableOffset
    //  40:  i64    StringTableOffset
    //  48:  i64    PrecursorSectionOffset
    //  56:  i64    FragmentSectionOffset
    let format_version          = read_i32_le(&header_buf, 4);
    if format_version < MIN_SUPPORTED_VERSION || format_version > MAX_SUPPORTED_VERSION {
        return Err(MslError::UnsupportedVersion(format_version));
    }

    let file_flags              = read_i32_le(&header_buf, 8);
    let n_precursors            = read_i32_le(&header_buf, 12);
    let n_proteins              = read_i32_le(&header_buf, 16);
    let n_strings               = read_i32_le(&header_buf, 24);
    // v2+ stores ExtAnnotationTableOffset at byte 28; v1 has Reserved here.
    // Either way, this i32 is not needed for basic reading — skip it.
    let _ext_annotation_table_offset = read_i32_le(&header_buf, 28);
    let protein_table_offset    = read_i64_le(&header_buf, 32);
    let string_table_offset     = read_i64_le(&header_buf, 40);
    let precursor_section_offset = read_i64_le(&header_buf, 48);

    let (has_ion_mobility, has_protein_data, has_gene_data, is_predicted) =
        decode_file_flags(file_flags);

    // ── 2. Read footer (last 20 bytes) ────────────────────────────────────────
    let footer_start = file_size - FOOTER_SIZE;
    file.seek(SeekFrom::Start(footer_start as u64))?;
    let mut footer_buf = [0u8; FOOTER_SIZE];
    file.read_exact(&mut footer_buf)?;

    let offset_table_offset = read_i64_le(&footer_buf, 0);
    let footer_n_precursors  = read_i32_le(&footer_buf, 8);
    let data_crc32           = read_u32_le(&footer_buf, 12);
    let trailing_magic       = read_u32_le(&footer_buf, 16);

    if trailing_magic != FOOTER_MAGIC_U32 {
        return Err(MslError::TrailingMagicMismatch);
    }

    // ── 3. Precursor count consistency ───────────────────────────────────────
    if footer_n_precursors != n_precursors {
        return Err(MslError::PrecursorCountMismatch {
            header: n_precursors,
            footer: footer_n_precursors,
        });
    }

    // ── 4. CRC-32 over [0, offset_table_offset) — streamed in chunks ─────────
    let crc_end = offset_table_offset as usize;
    if crc_end > file_size {
        return Err(MslError::InvalidData(format!(
            "offset_table_offset {} exceeds file size {}",
            crc_end, file_size
        )));
    }
    file.seek(SeekFrom::Start(0))?;
    let computed = compute_crc32_streaming(&mut file, crc_end)
        .map_err(MslError::Io)?;
    if computed != data_crc32 {
        return Err(MslError::CrcMismatch { expected: data_crc32, computed });
    }

    // ── 5. Read string table (targeted seek) ─────────────────────────────────
    file.seek(SeekFrom::Start(string_table_offset as u64))?;
    let mut st_header = [0u8; 8];
    file.read_exact(&mut st_header)?;
    let n_read = read_i32_le(&st_header, 0);
    if n_read != n_strings {
        return Err(MslError::InvalidData(format!(
            "String table NStrings mismatch: header={}, table={}",
            n_strings, n_read
        )));
    }
    let total_bytes = read_i32_le(&st_header, 4);
    let string_entry_size = (n_strings as usize) * 4 + total_bytes as usize;
    let mut string_entry_buf = vec![0u8; string_entry_size];
    file.read_exact(&mut string_entry_buf)?;

    let mut strings = Vec::with_capacity(n_strings as usize);
    let mut pos = 0usize;
    for _ in 0..n_strings {
        let length = read_i32_le(&string_entry_buf, pos) as usize;
        pos += 4;
        if pos + length > string_entry_buf.len() {
            return Err(MslError::InvalidData(
                "String table truncated while reading string data".to_string()
            ));
        }
        let s = String::from_utf8(string_entry_buf[pos..pos + length].to_vec())?;
        strings.push(s);
        pos += length;
    }

    // ── 6. Read protein table (targeted seek) ────────────────────────────────
    let proteins: Vec<(String, String, String)> =
        if n_proteins > 0 && protein_table_offset > 0 {
            let prot_total = n_proteins as usize * PROTEIN_SIZE;
            let mut prot_buf = vec![0u8; prot_total];
            file.seek(SeekFrom::Start(protein_table_offset as u64))?;
            file.read_exact(&mut prot_buf)?;

            let resolve = |idx: i32| -> String {
                if idx >= 0 && (idx as usize) < strings.len() {
                    strings[idx as usize].clone()
                } else {
                    String::new()
                }
            };

            let mut proteins = Vec::with_capacity(n_proteins as usize);
            for i in 0..n_proteins as usize {
                let base = i * PROTEIN_SIZE;
                let acc_idx  = read_i32_le(&prot_buf, base);
                let name_idx = read_i32_le(&prot_buf, base + 4);
                let gene_idx = read_i32_le(&prot_buf, base + 8);
                proteins.push((resolve(acc_idx), resolve(name_idx), resolve(gene_idx)));
            }
            proteins
        } else {
            Vec::new()
        };

    // ── 7. Read precursor section (targeted seek) ────────────────────────────
    let prec_total = n_precursors as usize * PRECURSOR_SIZE;
    let mut prec_buf = vec![0u8; prec_total];
    file.seek(SeekFrom::Start(precursor_section_offset as u64))?;
    file.read_exact(&mut prec_buf)?;

    let mut precursors: Vec<MslPrecursor> = Vec::with_capacity(n_precursors as usize);
    for i in 0..n_precursors as usize {
        let off = i * PRECURSOR_SIZE;
        let mut prec = parse_precursor(
            &prec_buf[off..off + PRECURSOR_SIZE],
            &strings,
            &proteins,
        );

        // ── 8. Fragment data: only read when requested ───────────────────────
        if load_fragment_data && prec.fragment_count > 0 {
            let frag_byte_count = prec.fragment_count as usize * FRAGMENT_SIZE;
            let mut frag_buf = vec![0u8; frag_byte_count];
            file.seek(SeekFrom::Start(prec.fragment_block_offset as u64))?;
            file.read_exact(&mut frag_buf)?;

            let mut frags = Vec::with_capacity(prec.fragment_count as usize);
            for j in 0..prec.fragment_count as usize {
                let foff = j * FRAGMENT_SIZE;
                frags.push(parse_fragment(&frag_buf[foff..foff + FRAGMENT_SIZE]));
            }
            prec.fragments = frags;
        }

        precursors.push(prec);
    }

    Ok(MslLibrary {
        format_version,
        n_precursors,
        n_proteins,
        has_ion_mobility,
        has_protein_data,
        has_gene_data,
        is_predicted,
        precursors,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// Full load: reads all precursors and all fragment blocks into memory.
///
/// Validates magic, version, trailing footer magic, precursor count, and CRC-32.
/// Returns `Err(MslError)` on any validation failure or I/O error.
pub fn load(path: &str) -> Result<MslLibrary, MslError> {
    load_impl(path, true)
}

/// Index-only load: reads precursor metadata only.
///
/// All `MslPrecursor.fragments` vecs are empty. Validation is identical to `load()`.
/// Use `load_fragments()` to fetch fragments for individual precursors on demand.
pub fn load_index_only(path: &str) -> Result<MslLibrary, MslError> {
    load_impl(path, false)
}

/// Loads and returns the fragment vec for a single precursor.
///
/// Opens the file, seeks to `precursor.fragment_block_offset`, and reads
/// `precursor.fragment_count` records. Does not re-validate the whole file.
pub fn load_fragments(path: &str, precursor: &MslPrecursor)
    -> Result<Vec<MslFragment>, MslError>
{
    if precursor.fragment_count <= 0 {
        return Err(MslError::InvalidData(format!(
            "Precursor has fragment_count={}; nothing to load.",
            precursor.fragment_count
        )));
    }

    let mut file = File::open(path)?;
    let file_size = file.seek(SeekFrom::End(0))?;
    let start = precursor.fragment_block_offset as u64;
    let end   = start + (precursor.fragment_count as u64) * FRAGMENT_SIZE as u64;
    if start >= file_size || end > file_size {
        return Err(MslError::InvalidData(format!(
            "fragment_block_offset {} out of range for file size {}",
            precursor.fragment_block_offset, file_size
        )));
    }

    file.seek(SeekFrom::Start(start))?;
    let byte_count = (precursor.fragment_count as usize) * FRAGMENT_SIZE;
    let mut buf = vec![0u8; byte_count];
    file.read_exact(&mut buf)?;

    let mut frags = Vec::with_capacity(precursor.fragment_count as usize);
    for i in 0..precursor.fragment_count as usize {
        let off = i * FRAGMENT_SIZE;
        frags.push(parse_fragment(&buf[off..off + FRAGMENT_SIZE]));
    }
    Ok(frags)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    // ── Synthetic .msl builder ────────────────────────────────────────────────

    struct MslBuilder {
        strings:    Vec<String>,
        proteins:   Vec<(i32, i32, i32)>,     // (acc_idx, name_idx, gene_idx)
        precursors: Vec<PrecursorSpec>,
    }

    struct PrecursorSpec {
        mz: f32, irt: f32, ion_mobility: f32,
        charge: i16, fragment_count: i16,
        protein_idx: i32,
        mod_seq_idx: i32, stripped_seq_idx: i32,
        fragment_block_offset: i64,
        q_value: f32,
        stripped_seq_length: i32,
        molecule_type: i16, dissociation_type: i16,
        nce_raw: i16, prec_flags: u8, source_type: u8,
        fragments: Vec<FragSpec>,
    }

    struct FragSpec {
        mz: f32, intensity: f32,
        product_type: i16, secondary_product_type: i16,
        fragment_number: i16, secondary_fragment_number: i16,
        residue_position: i16, charge: u8, flags: u8,
    }

    impl MslBuilder {
        fn new() -> Self {
            MslBuilder { strings: vec!["".to_string()], proteins: vec![], precursors: vec![] }
        }

        fn intern(&mut self, s: &str) -> i32 {
            if let Some(i) = self.strings.iter().position(|x| x == s) {
                return i as i32;
            }
            self.strings.push(s.to_string());
            (self.strings.len() - 1) as i32
        }

        fn add_protein(&mut self, acc: &str, name: &str, gene: &str) -> i32 {
            let (ai, ni, gi) = (self.intern(acc), self.intern(name), self.intern(gene));
            self.proteins.push((ai, ni, gi));
            (self.proteins.len() - 1) as i32
        }

        fn add_precursor(&mut self, mod_seq: &str, stripped: &str, mz: f32,
                         charge: i16, is_decoy: bool, q_value: f32,
                         protein_idx: i32, nce_actual: i16,
                         molecule_type: i16)
        {
            let mod_seq_idx     = self.intern(mod_seq);
            let stripped_seq_idx = self.intern(stripped);
            let mut flags = 0u8;
            if is_decoy { flags |= 0x01; }
            flags |= 0x02; // is_proteotypic
            self.precursors.push(PrecursorSpec {
                mz, irt: 10.0, ion_mobility: 0.0, charge,
                fragment_count: 0, protein_idx,
                mod_seq_idx, stripped_seq_idx,
                fragment_block_offset: 0,
                q_value, stripped_seq_length: stripped.len() as i32,
                molecule_type, dissociation_type: 0,
                nce_raw: nce_actual * 10, prec_flags: flags, source_type: 0,
                fragments: vec![],
            });
        }

        fn add_fragment(&mut self, prec_idx: usize, mz: f32, intensity: f32,
                        product_type: i16, fragment_number: i16,
                        charge: u8, neutral_loss_code: u8)
        {
            let mut flags = 0u8;
            flags |= (neutral_loss_code & 0x07) << 2;
            let prec = &mut self.precursors[prec_idx];
            prec.fragments.push(FragSpec {
                mz, intensity, product_type,
                secondary_product_type: -1,
                fragment_number, secondary_fragment_number: 0,
                residue_position: 0, charge, flags,
            });
            prec.fragment_count = prec.fragments.len() as i16;
        }

        fn build(&mut self) -> Vec<u8> {
            let n_prec    = self.precursors.len() as i32;
            let n_prot    = self.proteins.len() as i32;
            let n_strings = self.strings.len() as i32;

            let protein_size  = PROTEIN_SIZE;
            let precursor_size = PRECURSOR_SIZE;
            let fragment_size  = FRAGMENT_SIZE;

            let protein_offset: i64 = HEADER_SIZE as i64;
            let after_proteins = HEADER_SIZE + n_prot as usize * protein_size;

            // String table bytes
            let mut string_body = Vec::new();
            let mut total_body: i32 = 0;
            for s in &self.strings {
                let enc = s.as_bytes();
                total_body += enc.len() as i32;
                string_body.extend_from_slice(&(enc.len() as i32).to_le_bytes());
                string_body.extend_from_slice(enc);
            }
            let mut string_table_bytes = Vec::new();
            string_table_bytes.extend_from_slice(&n_strings.to_le_bytes());
            string_table_bytes.extend_from_slice(&total_body.to_le_bytes());
            string_table_bytes.extend_from_slice(&string_body);

            let string_table_offset = after_proteins as i64;
            let after_strings = after_proteins + string_table_bytes.len();

            let precursor_section_offset = after_strings as i64;
            let after_precursors = after_strings + n_prec as usize * precursor_size;

            let fragment_section_offset = after_precursors as i64;
            let mut frag_bytes = Vec::new();
            let mut cur_frag_off = fragment_section_offset;
            for prec in &mut self.precursors {
                prec.fragment_block_offset = cur_frag_off;
                for fr in &prec.fragments {
                    frag_bytes.extend_from_slice(&fr.mz.to_le_bytes());
                    frag_bytes.extend_from_slice(&fr.intensity.to_le_bytes());
                    frag_bytes.extend_from_slice(&fr.product_type.to_le_bytes());
                    frag_bytes.extend_from_slice(&fr.secondary_product_type.to_le_bytes());
                    frag_bytes.extend_from_slice(&fr.fragment_number.to_le_bytes());
                    frag_bytes.extend_from_slice(&fr.secondary_fragment_number.to_le_bytes());
                    frag_bytes.extend_from_slice(&fr.residue_position.to_le_bytes());
                    frag_bytes.push(fr.charge);
                    frag_bytes.push(fr.flags);
                }
                cur_frag_off += prec.fragments.len() as i64 * fragment_size as i64;
            }

            let after_fragments = fragment_section_offset as usize + frag_bytes.len();
            let offset_table_offset = after_fragments as i64;

            // Offset table
            let mut offset_table_bytes = Vec::new();
            for i in 0..n_prec as usize {
                let off = precursor_section_offset + (i * precursor_size) as i64;
                offset_table_bytes.extend_from_slice(&off.to_le_bytes());
            }

            // File flags
            let file_flags: i32 = if n_prot > 0 { 0x02 } else { 0x00 };

            // Header
            let mut header = Vec::with_capacity(HEADER_SIZE);
            header.extend_from_slice(&HEADER_MAGIC);
            header.extend_from_slice(&MIN_SUPPORTED_VERSION.to_le_bytes());
            header.extend_from_slice(&file_flags.to_le_bytes());
            header.extend_from_slice(&n_prec.to_le_bytes());
            header.extend_from_slice(&n_prot.to_le_bytes());
            header.extend_from_slice(&0i32.to_le_bytes()); // NElutionGroups
            header.extend_from_slice(&n_strings.to_le_bytes());
            header.extend_from_slice(&0i32.to_le_bytes()); // Reserved
            header.extend_from_slice(&protein_offset.to_le_bytes());
            header.extend_from_slice(&string_table_offset.to_le_bytes());
            header.extend_from_slice(&precursor_section_offset.to_le_bytes());
            header.extend_from_slice(&fragment_section_offset.to_le_bytes());
            assert_eq!(header.len(), HEADER_SIZE);

            // Protein table
            let mut protein_bytes = Vec::new();
            for (ai, ni, gi) in &self.proteins {
                protein_bytes.extend_from_slice(&ai.to_le_bytes());
                protein_bytes.extend_from_slice(&ni.to_le_bytes());
                protein_bytes.extend_from_slice(&gi.to_le_bytes());
                protein_bytes.extend_from_slice(&0i32.to_le_bytes()); // group_id
                protein_bytes.extend_from_slice(&0i32.to_le_bytes()); // n_prec
                protein_bytes.extend_from_slice(&0i32.to_le_bytes()); // flags
            }

            // Precursor section
            let mut precursor_bytes = Vec::new();
            for prec in &self.precursors {
                precursor_bytes.extend_from_slice(&prec.mz.to_le_bytes());
                precursor_bytes.extend_from_slice(&prec.irt.to_le_bytes());
                precursor_bytes.extend_from_slice(&prec.ion_mobility.to_le_bytes());
                precursor_bytes.extend_from_slice(&prec.charge.to_le_bytes());
                precursor_bytes.extend_from_slice(&prec.fragment_count.to_le_bytes());
                precursor_bytes.extend_from_slice(&0i32.to_le_bytes()); // elution_group_id
                precursor_bytes.extend_from_slice(&prec.protein_idx.to_le_bytes());
                precursor_bytes.extend_from_slice(&prec.mod_seq_idx.to_le_bytes());
                precursor_bytes.extend_from_slice(&prec.stripped_seq_idx.to_le_bytes());
                precursor_bytes.extend_from_slice(&prec.fragment_block_offset.to_le_bytes());
                precursor_bytes.extend_from_slice(&prec.q_value.to_le_bytes());
                precursor_bytes.extend_from_slice(&prec.stripped_seq_length.to_le_bytes());
                precursor_bytes.extend_from_slice(&prec.molecule_type.to_le_bytes());
                precursor_bytes.extend_from_slice(&prec.dissociation_type.to_le_bytes());
                precursor_bytes.extend_from_slice(&prec.nce_raw.to_le_bytes());
                precursor_bytes.push(prec.prec_flags);
                precursor_bytes.push(prec.source_type);
            }

            // Assemble body up to offset table
            let mut body = Vec::new();
            body.extend_from_slice(&header);
            body.extend_from_slice(&protein_bytes);
            body.extend_from_slice(&string_table_bytes);
            body.extend_from_slice(&precursor_bytes);
            body.extend_from_slice(&frag_bytes);
            body.extend_from_slice(&offset_table_bytes);

            // CRC over [0, offset_table_offset)
            let crc_data = &body[..offset_table_offset as usize];
            let data_crc32 = compute_crc32(crc_data);

            // Footer
            let n_prec_footer = n_prec;
            body.extend_from_slice(&offset_table_offset.to_le_bytes());
            body.extend_from_slice(&n_prec_footer.to_le_bytes());
            body.extend_from_slice(&data_crc32.to_le_bytes());
            body.extend_from_slice(&FOOTER_MAGIC_U32.to_le_bytes());

            body
        }
    }

    fn make_test_file() -> tempfile::NamedTempFile {
        let mut b = MslBuilder::new();
        let p0 = b.add_protein("P12345", "SOME_HUMAN", "GENE1");

        // Precursor 0 — PEPTIDER, non-decoy, q=0.01
        b.add_precursor("PEPTIDER", "PEPTIDER", 500.0, 2, false, 0.01f32, p0, 28, 0);
        b.add_fragment(0, 175.119, 1.0, 1, 1, 1, 0);
        b.add_fragment(0, 288.203, 0.8, 1, 2, 1, 0);
        b.add_fragment(0, 401.287, 0.6, 1, 3, 1, 0);

        // Precursor 1 — DECOY, is_decoy=true
        b.add_precursor("DECOY_PEPTIDER", "DECOY_PEPTIDER", 510.0, 2, true, 0.05f32, p0, 28, 0);
        b.add_fragment(1, 200.0, 1.0, 2, 1, 1, 0);

        // Precursor 2 — has H2O neutral loss fragment (code=1)
        b.add_precursor("ACDEFGHIK", "ACDEFGHIK", 700.0, 3, false, f32::NAN, -1, 35, 0);
        b.add_fragment(2, 175.119, 1.0, 1, 1, 1, 0);
        b.add_fragment(2, 175.119 - 18.011, 0.5, 1, 1, 1, 1); // H2O loss

        // Precursor 3 — Proteoform (molecule_type=1)
        b.add_precursor("TEIM[Ox]FDLK", "TEIMFDLK", 620.0, 2, false, f32::NAN, -1, 28, 1);
        b.add_fragment(3, 175.0, 1.0, 1, 1, 1, 0);

        let data = b.build();
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(&data).unwrap();
        f
    }

    // ── Test: valid file loads without error ──────────────────────────────────

    #[test]
    fn test_load_valid_file() {
        let f = make_test_file();
        let result = load(f.path().to_str().unwrap());
        assert!(result.is_ok(), "load() failed: {:?}", result.err());
    }

    // ── Test: wrong magic → InvalidMagic ─────────────────────────────────────

    #[test]
    fn test_magic_mismatch_returns_err() {
        let f = make_test_file();
        let mut data = std::fs::read(f.path()).unwrap();
        data[0] = 0xFF;
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.write_all(&data).unwrap();
        let result = load(tmp.path().to_str().unwrap());
        assert!(matches!(result, Err(MslError::InvalidMagic)), "Expected InvalidMagic, got {:?}", result.err());
    }

    // ── Test: wrong version → UnsupportedVersion ─────────────────────────────

    #[test]
    fn test_version_mismatch_returns_err() {
        let f = make_test_file();
        let mut data = std::fs::read(f.path()).unwrap();
        // FormatVersion at bytes 4–7
        let bad_ver: i32 = 99;
        data[4..8].copy_from_slice(&bad_ver.to_le_bytes());
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.write_all(&data).unwrap();
        let result = load(tmp.path().to_str().unwrap());
        assert!(matches!(result, Err(MslError::UnsupportedVersion(99))),
                "Expected UnsupportedVersion(99), got {:?}", result.err());
    }

    // ── Test: corrupt byte → CrcMismatch ─────────────────────────────────────

    #[test]
    fn test_crc_mismatch_returns_err() {
        let f = make_test_file();
        let mut data = std::fs::read(f.path()).unwrap();
        // Corrupt a byte in the fragment section (offset from header's FragmentSectionOffset)
        let frag_offset = i64::from_le_bytes(data[56..64].try_into().unwrap()) as usize;
        data[frag_offset] ^= 0xFF;
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        tmp.write_all(&data).unwrap();
        let result = load(tmp.path().to_str().unwrap());
        assert!(matches!(result, Err(MslError::CrcMismatch { .. })),
                "Expected CrcMismatch, got {:?}", result.err());
    }

    // ── Test: precursor count ─────────────────────────────────────────────────

    #[test]
    fn test_precursor_count() {
        let f = make_test_file();
        let lib = load(f.path().to_str().unwrap()).unwrap();
        assert_eq!(lib.n_precursors, 4);
        assert_eq!(lib.precursors.len(), 4);
    }

    // ── Test: precursor fields ────────────────────────────────────────────────

    #[test]
    fn test_precursor_fields() {
        let f = make_test_file();
        let lib = load(f.path().to_str().unwrap()).unwrap();
        let p = &lib.precursors[0];
        assert!((p.precursor_mz - 500.0).abs() < 0.001, "mz mismatch: {}", p.precursor_mz);
        assert_eq!(p.charge, 2);
        assert_eq!(p.modified_sequence, "PEPTIDER");
        assert!(!p.is_decoy);
        assert!((p.q_value - 0.01).abs() < 1e-5, "q_value mismatch: {}", p.q_value);
    }

    // ── Test: fragment fields ─────────────────────────────────────────────────

    #[test]
    fn test_fragment_fields() {
        let f = make_test_file();
        let lib = load(f.path().to_str().unwrap()).unwrap();
        let frags = &lib.precursors[0].fragments;
        assert_eq!(frags.len(), 3);
        assert!((frags[0].mz - 175.119).abs() < 0.001, "frag mz mismatch: {}", frags[0].mz);
        assert!((frags[0].intensity - 1.0).abs() < 1e-6);
        assert_eq!(frags[0].charge, 1);
    }

    // ── Test: neutral loss code 0 → mass 0.0 ─────────────────────────────────

    #[test]
    fn test_neutral_loss_none() {
        let f = make_test_file();
        let lib = load(f.path().to_str().unwrap()).unwrap();
        let frag = &lib.precursors[0].fragments[0];
        assert_eq!(frag.neutral_loss_code, 0);
        assert_eq!(frag.neutral_loss_mass, 0.0);
    }

    // ── Test: neutral loss code 1 → -18.0106 ─────────────────────────────────

    #[test]
    fn test_neutral_loss_h2o() {
        let f = make_test_file();
        let lib = load(f.path().to_str().unwrap()).unwrap();
        let frags = &lib.precursors[2].fragments;
        let h2o_frag = frags.iter().find(|fr| fr.neutral_loss_code == 1)
            .expect("No H2O-loss fragment found");
        assert!((h2o_frag.neutral_loss_mass - (-18.010565)).abs() < 1e-4,
                "H2O mass: {}", h2o_frag.neutral_loss_mass);
    }

    // ── Test: index-only load → fragments empty ───────────────────────────────

    #[test]
    fn test_index_only_fragments_empty() {
        let f = make_test_file();
        let lib = load_index_only(f.path().to_str().unwrap()).unwrap();
        for prec in &lib.precursors {
            assert!(prec.fragments.is_empty(),
                    "Expected empty fragments for {}", prec.modified_sequence);
        }
    }

    // ── Test: load_fragments on demand ───────────────────────────────────────

    #[test]
    fn test_load_fragments_on_demand() {
        let f = make_test_file();
        let lib = load_index_only(f.path().to_str().unwrap()).unwrap();
        let prec = &lib.precursors[0];
        assert!(prec.fragments.is_empty());

        let frags = load_fragments(f.path().to_str().unwrap(), prec).unwrap();
        assert_eq!(frags.len(), 3);
        assert!((frags[0].mz - 175.119).abs() < 0.001);
    }

    // ── Test: decoy flag ──────────────────────────────────────────────────────

    #[test]
    fn test_decoy_flag() {
        let f = make_test_file();
        let lib = load(f.path().to_str().unwrap()).unwrap();
        assert!(!lib.precursors[0].is_decoy);
        assert!(lib.precursors[1].is_decoy);
    }

    // ── Test: proteoform molecule_type ────────────────────────────────────────

    #[test]
    fn test_proteoform_molecule_type() {
        let f = make_test_file();
        let lib = load(f.path().to_str().unwrap()).unwrap();
        assert_eq!(lib.precursors[3].molecule_type, 1); // Proteoform
    }

    // ── Test: MslError variants implement Debug without panic ─────────────────

    #[test]
    fn test_error_display() {
        let errors: Vec<MslError> = vec![
            MslError::InvalidMagic,
            MslError::UnsupportedVersion(99),
            MslError::TrailingMagicMismatch,
            MslError::PrecursorCountMismatch { header: 10, footer: 9 },
            MslError::CrcMismatch { expected: 0xDEAD, computed: 0xBEEF },
            MslError::InvalidData("test message".to_string()),
        ];
        for e in &errors {
            // Debug formatting must not panic
            let s = format!("{:?}", e);
            assert!(!s.is_empty());
        }
    }

    // ── Test: protein table resolved correctly ────────────────────────────────

    #[test]
    fn test_protein_table_resolution() {
        let f = make_test_file();
        let lib = load(f.path().to_str().unwrap()).unwrap();
        assert_eq!(lib.precursors[0].protein_accession, "P12345");
        assert_eq!(lib.precursors[0].protein_name, "SOME_HUMAN");
        assert_eq!(lib.precursors[0].gene_name, "GENE1");
    }

    // ── Test: NCE decoding (×10 on disk → actual value) ──────────────────────

    #[test]
    fn test_nce_decoding() {
        let f = make_test_file();
        let lib = load(f.path().to_str().unwrap()).unwrap();
        assert_eq!(lib.precursors[0].nce, 28);  // stored as 280
        assert_eq!(lib.precursors[2].nce, 35);  // stored as 350
    }

    // ── Test: string table resolution ────────────────────────────────────────

    #[test]
    fn test_string_table_resolution() {
        let f = make_test_file();
        let lib = load(f.path().to_str().unwrap()).unwrap();
        assert_eq!(lib.precursors[1].modified_sequence, "DECOY_PEPTIDER");
        assert_eq!(lib.precursors[3].modified_sequence, "TEIM[Ox]FDLK");
    }
}
