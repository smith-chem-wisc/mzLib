"""
msl_reader.py — Reference reader for the .msl binary spectral library format.

Implements the format specified in FinalProjectSummary.md Section 4.
No dependencies outside the Python standard library.

Format overview:
  File Header       64 bytes  (magic, version, counts, section offsets)
  Protein Table     NProteins × 24 bytes
  String Table      variable  (int32 count, int32 TotalBytes, then length-prefixed UTF-8)
  Precursor Section NPrecursors × 56 bytes
  Fragment Section  variable  (FragmentCount × 20 bytes per precursor)
  Offset Table      NPrecursors × 8 bytes (int64[])
  Footer            20 bytes  (last 20 bytes of file)

Validation order (Section 4.10):
  1. Header magic == 0x4D 0x5A 0x4C 0x42
  2. FormatVersion in [1, 3]  (v2 adds ExtAnnotationTableOffset; v3 adds optional
     zstd compression of the fragment section)
  3. Trailing footer magic (read as little-endian uint32) == 0x4D5A4C42
  4. footer.NPrecursors == header.NPrecursors
  5. CRC-32/ISO-HDLC over bytes [0, footer.OffsetTableOffset) == footer.DataCrc32
"""

import struct
import zlib
import math
import io
from dataclasses import dataclass, field
from typing import Optional, List


# ---------------------------------------------------------------------------
# Struct format strings (all little-endian, prefix '<')
# ---------------------------------------------------------------------------

# File Header: 64 bytes
# 4s  = Magic (4 raw bytes)
# i   = FormatVersion
# i   = FileFlags
# i   = NPrecursors
# i   = NProteins
# i   = NElutionGroups
# i   = NStrings
# i   = Reserved
# i   = (padding to reach offset 32 — last i is Reserved; no extra pad needed)
# q   = ProteinTableOffset  (int64)
# q   = StringTableOffset   (int64)
# q   = PrecursorSectionOffset (int64)
# q   = FragmentSectionOffset  (int64)
HEADER_FMT = '<4s i i i i i i i q q q q'

# Protein Record: 24 bytes
# i   = AccessionStringIdx
# i   = NameStringIdx
# i   = GeneStringIdx
# i   = ProteinGroupId
# i   = NPrecursors
# i   = ProteinFlags
PROTEIN_FMT = '<i i i i i i'

# Precursor Record: 56 bytes
# f   = PrecursorMz
# f   = Irt
# f   = IonMobility
# h   = Charge  (int16)
# h   = FragmentCount (int16)
# i   = ElutionGroupId
# i   = ProteinIdx
# i   = ModifiedSeqStringIdx
# i   = StrippedSeqStringIdx
# q   = FragmentBlockOffset (int64)
# f   = QValue
# i   = StrippedSeqLength
# h   = MoleculeType (int16)
# h   = DissociationType (int16)
# h   = Nce (int16)
# B   = PrecursorFlags (uint8)
# B   = SourceType (uint8)
PRECURSOR_FMT = '<f f f h h i i i i q f i h h h B B'

# Fragment Record: 20 bytes
# f   = Mz
# f   = Intensity
# h   = ProductType (int16)
# h   = SecondaryProductType (int16)
# h   = FragmentNumber (int16)
# h   = SecondaryFragmentNumber (int16)
# h   = ResiduePosition (int16)
# B   = Charge (uint8)
# B   = Flags (uint8)
FRAGMENT_FMT = '<f f h h h h h B B'

# Footer: 20 bytes
# q   = OffsetTableOffset (int64)
# i   = NPrecursors (int32)
# I   = DataCrc32 (uint32)
# I   = TrailingMagic (uint32)
FOOTER_FMT = '<q i I I'

# Expected sizes — verify at import time
_EXPECTED_SIZES = {
    'HEADER':    (HEADER_FMT, 64),
    'PROTEIN':   (PROTEIN_FMT, 24),
    'PRECURSOR': (PRECURSOR_FMT, 56),
    'FRAGMENT':  (FRAGMENT_FMT, 20),
    'FOOTER':    (FOOTER_FMT, 20),
}
for _name, (_fmt, _expected) in _EXPECTED_SIZES.items():
    _actual = struct.calcsize(_fmt)
    if _actual != _expected:
        raise ImportError(
            f'msl_reader: {_name}_FMT calcsize={_actual}, expected={_expected}. '
            f'Format string is wrong.'
        )

# ---------------------------------------------------------------------------
# Neutral loss mass table (Section 4.11)
# ---------------------------------------------------------------------------

NEUTRAL_LOSS_MASSES = {
    0: 0.0,
    1: -18.010565,   # H2O
    2: -17.026549,   # NH3
    3: -97.976895,   # H3PO4
    4: -79.966331,   # HPO3
    5: -115.987460,  # H3PO4 + H2O  (H3PO4AndH2O — formerly misnamed PlusH2O)
    6: 0.0,          # Custom — exact mass stored externally (Prompt 11)
}

# ---------------------------------------------------------------------------
# Magic constants
# ---------------------------------------------------------------------------

HEADER_MAGIC = bytes([0x4D, 0x5A, 0x4C, 0x42])   # raw bytes: MZLB
FOOTER_MAGIC_UINT32 = 0x4D5A4C42                   # read as LE uint32 (bytes on disk: 42 4C 5A 4D)

# Version 1: original format. Header byte 28-31 is Reserved (always 0).
# Version 2: header byte 28-31 repurposed as ExtAnnotationTableOffset for
#             custom neutral-loss masses (FileFlagHasExtAnnotations).
#             Version-1 readers safely ignore this field (it was 0 in all v1 files).

MIN_SUPPORTED_VERSION = 1
MAX_SUPPORTED_VERSION = 3  # must equal MslFormat.CurrentVersion in MslFormat.cs
                            # Version 3 adds optional zstd compression of the
                            # fragment section (FileFlagIsCompressed, bit 5).
                            # Uncompressed v3 files are identical to v2 files
                            # except for the FormatVersion field.


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class MslFragment:
    mz: float
    intensity: float
    product_type: int               # raw int16
    secondary_product_type: Optional[int]  # None when -1 on disk
    fragment_number: int
    secondary_fragment_number: int
    residue_position: int
    charge: int
    is_internal: bool
    is_diagnostic: bool
    neutral_loss_code: int
    neutral_loss_mass: float        # Da; 0.0 for code 0
    exclude_from_quant: bool


@dataclass
class MslPrecursor:
    precursor_mz: float
    irt: float
    ion_mobility: float
    charge: int
    fragment_count: int
    elution_group_id: int
    protein_idx: int
    modified_sequence: str
    stripped_sequence: str
    fragment_block_offset: int
    q_value: float
    stripped_seq_length: int
    molecule_type: int              # 0=Peptide, 1=Proteoform, 2=Oligonucleotide, 3=Glycopeptide
    dissociation_type: int
    nce: int                        # actual NCE (already ÷10 from on-disk × 10 encoding)
    is_decoy: bool
    is_proteotypic: bool
    rt_is_calibrated: bool
    source_type: int                # 0=Predicted, 1=Empirical, 2=EmpiricalRefined
    protein_accession: str          # resolved from protein table; '' if no protein
    protein_name: str
    gene_name: str
    fragments: List[MslFragment] = field(default_factory=list)


@dataclass
class MslLibrary:
    format_version: int
    n_precursors: int
    n_proteins: int
    has_ion_mobility: bool
    has_protein_data: bool
    has_gene_data: bool
    is_predicted: bool
    is_compressed: bool = False   # True when FileFlagIsCompressed (bit 5) is set
    precursors: List[MslPrecursor] = field(default_factory=list)


# ---------------------------------------------------------------------------
# CRC-32 helper
# ---------------------------------------------------------------------------

def _compute_crc32(data: bytes) -> int:
    """CRC-32/ISO-HDLC — same polynomial as mzLib's implementation."""
    return zlib.crc32(data) & 0xFFFFFFFF


# ---------------------------------------------------------------------------
# Internal parsing helpers
# ---------------------------------------------------------------------------

def _decode_file_flags(file_flags: int):
    has_ion_mobility  = bool(file_flags & 0x01)
    has_protein_data  = bool(file_flags & 0x02)
    has_gene_data     = bool(file_flags & 0x04)
    is_predicted      = bool(file_flags & 0x08)
    # bit 4: has_ext_annotations (handled implicitly via ExtAnnotationTableOffset)
    is_compressed     = bool(file_flags & 0x20)  # bit 5: FileFlagIsCompressed (v3+)
    return has_ion_mobility, has_protein_data, has_gene_data, is_predicted, is_compressed


def _decode_precursor_flags(flags_byte: int):
    """
    Precursor flags byte layout (Section 4.5):
      bit 0 = is_decoy
      bit 1 = is_proteotypic
      bit 2 = rt_is_calibrated
    """
    is_decoy         = bool(flags_byte & 0x01)
    is_proteotypic   = bool(flags_byte & 0x02)
    rt_is_calibrated = bool(flags_byte & 0x04)
    return is_decoy, is_proteotypic, rt_is_calibrated


def _decode_file_flags(file_flags: int):
    has_ion_mobility  = bool(file_flags & 0x01)
    has_protein_data  = bool(file_flags & 0x02)
    has_gene_data     = bool(file_flags & 0x04)
    is_predicted      = bool(file_flags & 0x08)
    has_ion_mobility, has_protein_data, has_gene_data, is_predicted, is_compressed = \
        _decode_file_flags(file_flags)


def _parse_fragment(raw: bytes) -> MslFragment:
    (mz, intensity, product_type, secondary_product_type_raw,
     fragment_number, secondary_fragment_number, residue_position,
     charge, flags_byte) = struct.unpack(FRAGMENT_FMT, raw)

    secondary_product_type = (None if secondary_product_type_raw == -1
                              else secondary_product_type_raw)

    is_internal, is_diagnostic, nl_code, exclude_from_quant = \
        _decode_fragment_flags(flags_byte)

    nl_mass = NEUTRAL_LOSS_MASSES.get(nl_code, 0.0)

    return MslFragment(
        mz=mz,
        intensity=intensity,
        product_type=product_type,
        secondary_product_type=secondary_product_type,
        fragment_number=fragment_number,
        secondary_fragment_number=secondary_fragment_number,
        residue_position=residue_position,
        charge=charge,
        is_internal=is_internal,
        is_diagnostic=is_diagnostic,
        neutral_loss_code=nl_code,
        neutral_loss_mass=nl_mass,
        exclude_from_quant=exclude_from_quant,
    )


def _read_fragments(f: io.RawIOBase, block_offset: int, count: int) -> List[MslFragment]:
    """Seek to block_offset and read count fragment records."""
    frag_size = struct.calcsize(FRAGMENT_FMT)  # 20
    f.seek(block_offset)
    raw = f.read(frag_size * count)
    if len(raw) < frag_size * count:
        raise ValueError(
            f'Fragment block truncated: expected {frag_size * count} bytes '
            f'at offset {block_offset}, got {len(raw)}'
        )
    fragments = []
    for i in range(count):
        chunk = raw[i * frag_size:(i + 1) * frag_size]
        fragments.append(_parse_fragment(chunk))
    return fragments


def _read_string_table(f: io.RawIOBase, offset: int, n_strings: int) -> List[str]:
    """
    Deserialize the interned string table.

    On-disk layout (Note 7.3 — spec omission):
        int32   NStrings       (= header.NStrings)
        int32   TotalBytes     (total UTF-8 body bytes — skip this field)
        for each string:
            int32   length
            bytes   data       (UTF-8, no null terminator)
    """
    f.seek(offset)
    n_read = struct.unpack('<i', f.read(4))[0]
    if n_read != n_strings:
        raise ValueError(
            f'String table NStrings mismatch: header says {n_strings}, '
            f'table header says {n_read}'
        )
    _total_body_bytes = struct.unpack('<i', f.read(4))[0]  # skip TotalBytes field

    strings: List[str] = []
    for _ in range(n_strings):
        length = struct.unpack('<i', f.read(4))[0]
        data = f.read(length)
        if len(data) < length:
            raise ValueError('String table truncated while reading string data')
        strings.append(data.decode('utf-8'))
    return strings


def _read_protein_table(f: io.RawIOBase, offset: int, n_proteins: int,
                         strings: List[str]):
    """
    Read NProteins × 24-byte protein records and resolve string indices.
    Returns list of (accession, name, gene) tuples.
    """
    prot_size = struct.calcsize(PROTEIN_FMT)  # 24
    f.seek(offset)
    raw = f.read(prot_size * n_proteins)
    if len(raw) < prot_size * n_proteins:
        raise ValueError('Protein table truncated')

    proteins = []
    for i in range(n_proteins):
        chunk = raw[i * prot_size:(i + 1) * prot_size]
        (acc_idx, name_idx, gene_idx,
         _group_id, _n_prec, _flags) = struct.unpack(PROTEIN_FMT, chunk)

        def resolve(idx):
            if idx < 0 or idx >= len(strings):
                return ''
            return strings[idx]

        proteins.append((resolve(acc_idx), resolve(name_idx), resolve(gene_idx)))
    return proteins


def _validate_and_load(path: str, load_fragments_data: bool) -> MslLibrary:
    """
    Core loading function used by both load() and load_index_only().
    Performs all five validation checks from Section 4.10.
    """
    with open(path, 'rb') as f:
        # ── 1. Read and validate header magic ────────────────────────────────
        magic_bytes = f.read(4)
        if magic_bytes != HEADER_MAGIC:
            raise ValueError(
                f'Invalid .msl magic: expected {list(HEADER_MAGIC)}, '
                f'got {list(magic_bytes)}'
            )

        # ── 2. Read full header ───────────────────────────────────────────────
        f.seek(0)
        header_size = struct.calcsize(HEADER_FMT)  # 64
        header_raw = f.read(header_size)
        if len(header_raw) < header_size:
            raise ValueError('File too short to contain a valid .msl header')

        (magic, format_version, file_flags, n_precursors, n_proteins,
         n_elution_groups, n_strings, reserved,
         protein_table_offset, string_table_offset,
         precursor_section_offset, fragment_section_offset) = \
            struct.unpack(HEADER_FMT, header_raw)

        # ── 2. Validate format version ────────────────────────────────────────
        if not (MIN_SUPPORTED_VERSION <= format_version <= MAX_SUPPORTED_VERSION):
            raise ValueError(
                f'Unsupported .msl format version: {format_version} '
                f'(supported: {MIN_SUPPORTED_VERSION}–{MAX_SUPPORTED_VERSION})'
            )

        # ── 3. Read footer (last 20 bytes) and validate trailing magic ────────
        footer_size = struct.calcsize(FOOTER_FMT)  # 20
        f.seek(-footer_size, 2)
        footer_raw = f.read(footer_size)
        if len(footer_raw) < footer_size:
            raise ValueError('File too short to contain a valid .msl footer')

        (offset_table_offset, footer_n_precursors,
         data_crc32, trailing_magic) = struct.unpack(FOOTER_FMT, footer_raw)

        if trailing_magic != FOOTER_MAGIC_UINT32:
            raise ValueError(
                f'Trailing magic mismatch: expected 0x{FOOTER_MAGIC_UINT32:08X}, '
                f'got 0x{trailing_magic:08X}. File may be truncated or corrupt.'
            )

        # ── 4. Validate precursor count ───────────────────────────────────────
        if footer_n_precursors != n_precursors:
            raise ValueError(
                f'Precursor count mismatch: header={n_precursors}, '
                f'footer={footer_n_precursors}. File may be truncated.'
            )

        # ── 5. CRC-32 validation ──────────────────────────────────────────────
        f.seek(0)
        crc_data = f.read(offset_table_offset)
        if len(crc_data) < offset_table_offset:
            raise ValueError(
                f'File truncated: could not read {offset_table_offset} bytes '
                f'for CRC computation'
            )
        computed_crc = _compute_crc32(crc_data)
        if computed_crc != data_crc32:
            raise ValueError(
                f'CRC-32 mismatch: file stores 0x{data_crc32:08X}, '
                f'computed 0x{computed_crc:08X}. Data is corrupt.'
            )

        # ── Decode file flags ─────────────────────────────────────────────────
        has_ion_mobility, has_protein_data, has_gene_data, is_predicted = \
            _decode_file_flags(file_flags)

        # ── String table ──────────────────────────────────────────────────────
        strings = _read_string_table(f, string_table_offset, n_strings)

        # ── Protein table ─────────────────────────────────────────────────────
        proteins = []
        if n_proteins > 0 and protein_table_offset > 0:
            proteins = _read_protein_table(f, protein_table_offset, n_proteins, strings)

        # ── Precursor section ─────────────────────────────────────────────────
        prec_size = struct.calcsize(PRECURSOR_FMT)  # 56
        f.seek(precursor_section_offset)
        prec_raw = f.read(prec_size * n_precursors)
        if len(prec_raw) < prec_size * n_precursors:
            raise ValueError('Precursor section truncated')

        precursors: List[MslPrecursor] = []
        for i in range(n_precursors):
            chunk = prec_raw[i * prec_size:(i + 1) * prec_size]
            (prec_mz, irt, ion_mobility, charge, fragment_count,
             elution_group_id, protein_idx,
             modified_seq_idx, stripped_seq_idx,
             fragment_block_offset, q_value, stripped_seq_length,
             molecule_type, dissociation_type, nce_raw,
             precursor_flags, source_type) = struct.unpack(PRECURSOR_FMT, chunk)

            is_decoy, is_proteotypic, rt_is_calibrated = \
                _decode_precursor_flags(precursor_flags)

            modified_sequence = (strings[modified_seq_idx]
                                 if 0 <= modified_seq_idx < len(strings) else '')
            stripped_sequence = (strings[stripped_seq_idx]
                                 if 0 <= stripped_seq_idx < len(strings) else '')

            acc, pname, gname = '', '', ''
            if 0 <= protein_idx < len(proteins):
                acc, pname, gname = proteins[protein_idx]

            nce = nce_raw // 10  # stored as actual_nce × 10 on disk

            prec = MslPrecursor(
                precursor_mz=prec_mz,
                irt=irt,
                ion_mobility=ion_mobility,
                charge=charge,
                fragment_count=fragment_count,
                elution_group_id=elution_group_id,
                protein_idx=protein_idx,
                modified_sequence=modified_sequence,
                stripped_sequence=stripped_sequence,
                fragment_block_offset=fragment_block_offset,
                q_value=q_value,
                stripped_seq_length=stripped_seq_length,
                molecule_type=molecule_type,
                dissociation_type=dissociation_type,
                nce=nce,
                is_decoy=is_decoy,
                is_proteotypic=is_proteotypic,
                rt_is_calibrated=rt_is_calibrated,
                source_type=source_type,
                protein_accession=acc,
                protein_name=pname,
                gene_name=gname,
                fragments=[],
            )
            precursors.append(prec)

                # ── Compression descriptor and fragment decompression (version 3+) ────────
        # When FileFlagIsCompressed is set (bit 5), the fragment section on disk is
        # a single zstd frame. A 16-byte compression descriptor at file offset 64
        # (immediately after the header) records the compressed and uncompressed sizes.
        #
        # FragmentBlockOffset values in the precursor records are byte offsets into
        # the DECOMPRESSED buffer, not absolute file positions.  We decompress into
        # a BytesIO buffer and pass that to _read_fragments instead of the file handle.
        #
        # For uncompressed files (is_compressed == False), fragment_source remains
        # the real file handle and FragmentBlockOffset values are absolute positions.
        fragment_source = f  # default: read directly from the file
        if is_compressed:
            # Read compression descriptor at offset 64 (immediately after 64-byte header)
            COMPRESSION_DESCRIPTOR_OFFSET = 64
            COMPRESSION_DESCRIPTOR_SIZE   = 16  # int64 compressed + int64 uncompressed
            f.seek(COMPRESSION_DESCRIPTOR_OFFSET)
            desc_raw = f.read(COMPRESSION_DESCRIPTOR_SIZE)
            if len(desc_raw) < COMPRESSION_DESCRIPTOR_SIZE:
                raise ValueError(
                    'Compressed .msl file is missing the 16-byte compression descriptor '
                    f'at offset {COMPRESSION_DESCRIPTOR_OFFSET}.'
                )
            compressed_size, uncompressed_size = struct.unpack('<qq', desc_raw)
 
            # The compressed fragment frame begins at fragment_section_offset
            f.seek(fragment_section_offset)
            compressed_data = f.read(compressed_size)
            if len(compressed_data) < compressed_size:
                raise ValueError(
                    f'Compressed fragment section truncated: expected {compressed_size} '
                    f'bytes at offset {fragment_section_offset}, got {len(compressed_data)}.'
                )
 
            # zstd decompression — requires the zstandard package
            try:
                import zstandard as zstd
                dctx = zstd.ZstdDecompressor()
                decompressed_data = dctx.decompress(compressed_data,
                                                     max_output_size=uncompressed_size)
            except ImportError:
                raise ImportError(
                    "Reading compressed .msl files requires the 'zstandard' package. "
                    "Install it with: pip install zstandard"
                )
 
            if len(decompressed_data) != uncompressed_size:
                raise ValueError(
                    f'Decompressed fragment data length {len(decompressed_data)} does not '
                    f'match expected uncompressed size {uncompressed_size}.'
                )
 
            fragment_source = io.BytesIO(decompressed_data)

        # ── Fragment data (full load only) ────────────────────────────────────
        if load_fragments_data:
            for prec in precursors:
                if prec.fragment_count > 0:
                    prec.fragments = _read_fragments(
                        fragment_source,
                        prec.fragment_block_offset,
                        prec.fragment_count
                    )

        return MslLibrary(
            format_version=format_version,
            n_precursors=n_precursors,
            n_proteins=n_proteins,
            has_ion_mobility=has_ion_mobility,
            has_protein_data=has_protein_data,
            has_gene_data=has_gene_data,
            is_predicted=is_predicted,
            is_compressed=is_compressed,
            precursors=precursors,
        )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def load(path: str) -> MslLibrary:
    """
    Full load: reads all precursors and all fragment blocks into memory.

    Validates magic, version, trailing footer magic, precursor count, and CRC-32.
    Raises ValueError on any validation failure.
    Raises FileNotFoundError if path does not exist.

    Returns an MslLibrary with each MslPrecursor.fragments fully populated.
    """
    return _validate_and_load(path, load_fragments_data=True)


def load_index_only(path: str) -> MslLibrary:
    """
    Index-only load: reads precursor metadata only.

    Fragment blocks are not loaded; each MslPrecursor.fragments list is empty.
    Validation is identical to load() — all five checks are performed.

    Returns an MslLibrary with empty fragment lists. Use load_fragments() to
    fetch fragments for individual precursors on demand.
    """
    return _validate_and_load(path, load_fragments_data=False)


def load_fragments(path: str, precursor: MslPrecursor) -> List[MslFragment]:
    \"\"\"
    Loads and returns the fragment list for a single precursor.
 
    NOTE: This function is not supported for compressed .msl files (FileFlagIsCompressed).
    For compressed files, use load() which decompresses the fragment section automatically.
    An ImportError will be raised if zstandard is not installed and the file is compressed.
 
    Opens the file, seeks to precursor.fragment_block_offset, and reads
    precursor.fragment_count records. Does not re-validate the whole file.
 
    Raises ValueError if fragment_count == 0 or block_offset is out of range.
    Raises FileNotFoundError if path does not exist.
    \"\"\"
    if precursor.fragment_count <= 0:
        raise ValueError(
            f'Precursor has fragment_count={precursor.fragment_count}; nothing to load.'
        )
 
    # Check whether this file is compressed — if so, we cannot seek directly.
    # Read only the file flags byte (offset 8, 4 bytes).
    with open(path, 'rb') as f:
        f.seek(8)
        file_flags = struct.unpack('<i', f.read(4))[0]
    is_compressed = bool(file_flags & 0x20)
    if is_compressed:
        raise ValueError(
            'load_fragments() does not support compressed .msl files. '
            'Use load() instead, which decompresses the fragment section automatically.'
        )
 
    import os
    file_size = os.path.getsize(path)
    frag_size = struct.calcsize(FRAGMENT_FMT)
    end_offset = precursor.fragment_block_offset + frag_size * precursor.fragment_count
    if precursor.fragment_block_offset < 0 or end_offset > file_size:
        raise ValueError(
            f'fragment_block_offset {precursor.fragment_block_offset} is out of range '
            f'for file size {file_size}'
        )
    with open(path, 'rb') as f:
        return _read_fragments(f, precursor.fragment_block_offset, precursor.fragment_count)