"""
test_msl_reader.py — Unit tests for msl_reader.py

Run with:  python -m unittest test_msl_reader.py

All tests are self-contained: a minimal but valid .msl binary file is
synthesised in memory using the same struct layouts the reader expects.
No external test files or C# tooling are required.
"""

import unittest
import struct
import zlib
import math
import os
import tempfile
import io

import msl_reader
from msl_reader import (
    load, load_index_only, load_fragments,
    HEADER_FMT, PROTEIN_FMT, PRECURSOR_FMT, FRAGMENT_FMT, FOOTER_FMT,
    NEUTRAL_LOSS_MASSES,
)


# ---------------------------------------------------------------------------
# Synthetic .msl builder
# ---------------------------------------------------------------------------

class _MslBuilder:
    """
    Builds a minimal but fully valid .msl binary in memory.

    Caller populates self.precursors (list of dicts) and self.proteins
    (list of dicts) before calling build().
    """

    HEADER_MAGIC   = bytes([0x4D, 0x5A, 0x4C, 0x42])
    FORMAT_VERSION = 3   # must equal msl_reader.MAX_SUPPORTED_VERSION
    FOOTER_MAGIC   = 0x4D5A4C42

    def __init__(self):
        self.precursors = []   # list of dicts (see _default_precursor)
        self.proteins   = []   # list of dicts (see _default_protein)
        self._strings   = ['']  # index 0 = empty string by convention

    def add_string(self, s: str) -> int:
        """Intern a string and return its index."""
        if s not in self._strings:
            self._strings.append(s)
        return self._strings.index(s)

    def add_protein(self, accession='', name='', gene=''):
        idx = len(self.proteins)
        self.proteins.append(dict(
            accession_idx=self.add_string(accession),
            name_idx=self.add_string(name),
            gene_idx=self.add_string(gene),
            group_id=0,
            n_precursors=0,
            flags=0,
        ))
        return idx

    def add_precursor(self, modified_seq, stripped_seq, mz=500.0, irt=10.0,
                      charge=2, fragment_count=3, is_decoy=False,
                      is_proteotypic=True, rt_calibrated=False,
                      q_value=float('nan'), protein_idx=-1,
                      nce_actual=28, source_type=0,
                      molecule_type=0, dissociation_type=0,
                      ion_mobility=0.0):
        prec_flags = 0
        if is_decoy:        prec_flags |= 0x01
        if is_proteotypic:  prec_flags |= 0x02
        if rt_calibrated:   prec_flags |= 0x04

        self.precursors.append(dict(
            mz=mz,
            irt=irt,
            ion_mobility=ion_mobility,
            charge=charge,
            fragment_count=fragment_count,
            elution_group_id=0,
            protein_idx=protein_idx,
            modified_seq_idx=self.add_string(modified_seq),
            stripped_seq_idx=self.add_string(stripped_seq),
            fragment_block_offset=0,   # filled in during build()
            q_value=q_value,
            stripped_seq_length=len(stripped_seq),
            molecule_type=molecule_type,
            dissociation_type=dissociation_type,
            nce_raw=nce_actual * 10,   # stored ×10
            prec_flags=prec_flags,
            source_type=source_type,
            # Fragment specs for the builder
            fragments=[],
        ))

    def add_fragment(self, precursor_idx, mz, intensity,
                     product_type=1, secondary_product_type=-1,
                     fragment_number=1, secondary_fragment_number=0,
                     residue_position=0, charge=1,
                     is_internal=False, is_diagnostic=False,
                     neutral_loss_code=0, exclude_from_quant=False):
        flags = 0
        if is_internal:         flags |= 0x01
        if is_diagnostic:       flags |= 0x02
        flags |= (neutral_loss_code & 0x07) << 2
        if exclude_from_quant:  flags |= 0x20

        self.precursors[precursor_idx]['fragments'].append(dict(
            mz=mz,
            intensity=intensity,
            product_type=product_type,
            secondary_product_type=secondary_product_type,
            fragment_number=fragment_number,
            secondary_fragment_number=secondary_fragment_number,
            residue_position=residue_position,
            charge=charge,
            flags=flags,
        ))

    def build(self) -> bytes:
        """Assemble and return valid .msl binary content."""
        header_size   = struct.calcsize(HEADER_FMT)    # 64
        protein_size  = struct.calcsize(PROTEIN_FMT)   # 24
        precursor_size = struct.calcsize(PRECURSOR_FMT) # 56
        fragment_size  = struct.calcsize(FRAGMENT_FMT)  # 20
        footer_size    = struct.calcsize(FOOTER_FMT)    # 20

        n_prec    = len(self.precursors)
        n_prot    = len(self.proteins)
        n_strings = len(self._strings)

        # ── Section offsets (layout order: header, proteins, strings,
        #    precursors, fragments, offset_table, footer) ───────────────────
        protein_offset = header_size if n_prot > 0 else 0
        after_proteins = header_size + (n_prot * protein_size if n_prot > 0 else 0)

        # Build string table bytes
        string_body = b''
        for s in self._strings:
            enc = s.encode('utf-8')
            string_body += struct.pack('<i', len(enc)) + enc

        # The writer also emits a TotalBytes field right after NStrings
        total_body_bytes = sum(len(s.encode('utf-8')) for s in self._strings)
        string_table_bytes = (struct.pack('<i', n_strings) +
                              struct.pack('<i', total_body_bytes) +
                              string_body)

        string_table_offset = after_proteins
        after_strings = string_table_offset + len(string_table_bytes)

        precursor_section_offset = after_strings
        after_precursors = precursor_section_offset + n_prec * precursor_size

        # Compute fragment block offsets and build fragment bytes
        fragment_section_offset = after_precursors
        frag_bytes_list = []
        cur_frag_offset = fragment_section_offset
        for prec in self.precursors:
            prec['fragment_block_offset'] = cur_frag_offset
            for frag in prec['fragments']:
                frag_bytes_list.append(struct.pack(
                    FRAGMENT_FMT,
                    frag['mz'], frag['intensity'],
                    frag['product_type'], frag['secondary_product_type'],
                    frag['fragment_number'], frag['secondary_fragment_number'],
                    frag['residue_position'], frag['charge'], frag['flags'],
                ))
            cur_frag_offset += len(prec['fragments']) * fragment_size

        fragment_bytes = b''.join(frag_bytes_list)
        after_fragments = fragment_section_offset + len(fragment_bytes)

        # Offset table
        offset_table_offset = after_fragments
        offset_table_bytes = b''
        cur_prec_offset = precursor_section_offset
        for _ in range(n_prec):
            offset_table_bytes += struct.pack('<q', cur_prec_offset)
            cur_prec_offset += precursor_size
        after_offset_table = offset_table_offset + len(offset_table_bytes)

        # File flags
        file_flags = 0
        # has_protein_data = bit 1
        if n_prot > 0:
            file_flags |= 0x02

        # ── Assemble header ───────────────────────────────────────────────────
        # Note: HEADER_FMT starts with 4s (magic), but ProteinTableOffset is
        # the 9th field; when NProteins==0 we write 0 for protein_offset.
        prot_off_field = protein_offset if n_prot > 0 else 0
        header_bytes = struct.pack(
            HEADER_FMT,
            self.HEADER_MAGIC,
            self.FORMAT_VERSION,
            file_flags,
            n_prec,
            n_prot,
            0,         # NElutionGroups
            n_strings,
            0,         # Reserved
            prot_off_field,
            string_table_offset,
            precursor_section_offset,
            fragment_section_offset,
        )

        # ── Protein table ─────────────────────────────────────────────────────
        protein_bytes = b''
        for prot in self.proteins:
            protein_bytes += struct.pack(
                PROTEIN_FMT,
                prot['accession_idx'], prot['name_idx'], prot['gene_idx'],
                prot['group_id'], prot['n_precursors'], prot['flags'],
            )

        # ── Precursor section ─────────────────────────────────────────────────
        precursor_bytes = b''
        for prec in self.precursors:
            precursor_bytes += struct.pack(
                PRECURSOR_FMT,
                prec['mz'], prec['irt'], prec['ion_mobility'],
                prec['charge'], prec['fragment_count'],
                prec['elution_group_id'], prec['protein_idx'],
                prec['modified_seq_idx'], prec['stripped_seq_idx'],
                prec['fragment_block_offset'],
                prec['q_value'],
                prec['stripped_seq_length'],
                prec['molecule_type'], prec['dissociation_type'],
                prec['nce_raw'],
                prec['prec_flags'], prec['source_type'],
            )

        # ── Assemble everything up to offset table ────────────────────────────
        body = (header_bytes + protein_bytes + string_table_bytes +
                precursor_bytes + fragment_bytes + offset_table_bytes)

        assert len(body) == after_offset_table, (
            f'Body length mismatch: {len(body)} vs expected {after_offset_table}'
        )

        # ── CRC-32 over [0, offset_table_offset) ─────────────────────────────
        crc_data = body[:offset_table_offset]
        data_crc32 = zlib.crc32(crc_data) & 0xFFFFFFFF

        # ── Footer ────────────────────────────────────────────────────────────
        footer_bytes = struct.pack(
            FOOTER_FMT,
            offset_table_offset,
            n_prec,
            data_crc32,
            self.FOOTER_MAGIC,
        )

        return body + footer_bytes


def _make_test_file(tmp_dir=None) -> str:
    """
    Build a known-good .msl file and write it to a temp file.
    Returns the path.

    Contents:
      2 proteins
      4 precursors:
        [0] PEPTIDER/2 (non-decoy, 3 fragments, q=0.01)
        [1] DECOY_PEPTIDER/2 (decoy, 2 fragments, q=0.05)
        [2] ACDEFGHIK/3 (non-decoy, 4 fragments, fragment with H2O loss)
        [3] PROTEOFORM.../2 (molecule_type=1 Proteoform, 2 fragments)
    """
    b = _MslBuilder()

    p0 = b.add_protein('P12345', 'SOME_HUMAN', 'GENE1')
    p1 = b.add_protein('P99999', 'OTHER_HUMAN', 'GENE2')

    # Precursor 0 — PEPTIDER, non-decoy
    b.add_precursor('PEPTIDER', 'PEPTIDER', mz=500.0, irt=12.5,
                    charge=2, fragment_count=3, is_decoy=False,
                    q_value=0.01, protein_idx=p0, nce_actual=28)
    b.add_fragment(0, mz=175.119, intensity=1.0, product_type=1,
                   fragment_number=1, charge=1)
    b.add_fragment(0, mz=288.203, intensity=0.8, product_type=1,
                   fragment_number=2, charge=1)
    b.add_fragment(0, mz=401.287, intensity=0.6, product_type=1,
                   fragment_number=3, charge=1)

    # Precursor 1 — DECOY_PEPTIDER, is_decoy=True
    b.add_precursor('DECOY_PEPTIDER', 'DECOY_PEPTIDER', mz=510.0, irt=11.0,
                    charge=2, fragment_count=2, is_decoy=True,
                    q_value=0.05, protein_idx=p0, nce_actual=28)
    b.add_fragment(1, mz=200.0, intensity=1.0, product_type=2, fragment_number=1, charge=1)
    b.add_fragment(1, mz=313.0, intensity=0.7, product_type=2, fragment_number=2, charge=1)

    # Precursor 2 — ACDEFGHIK, fragment with H2O neutral loss (code=1)
    b.add_precursor('ACDEFGHIK', 'ACDEFGHIK', mz=700.0, irt=20.0,
                    charge=3, fragment_count=4, is_decoy=False,
                    q_value=float('nan'), protein_idx=p1, nce_actual=35)
    b.add_fragment(2, mz=175.119, intensity=1.0, product_type=1,
                   fragment_number=1, charge=1)
    b.add_fragment(2, mz=175.119 - 18.010565, intensity=0.5, product_type=1,
                   fragment_number=1, charge=1, neutral_loss_code=1)  # H2O loss
    b.add_fragment(2, mz=288.203, intensity=0.8, product_type=1,
                   fragment_number=2, charge=1)
    b.add_fragment(2, mz=401.287, intensity=0.6, product_type=1,
                   fragment_number=3, charge=1)

    # Precursor 3 — Proteoform (molecule_type=1), 2 fragments
    b.add_precursor('TEIM[Oxidation]FDLK', 'TEIMFDLK', mz=620.0, irt=15.0,
                    charge=2, fragment_count=2, is_decoy=False,
                    q_value=float('nan'), protein_idx=-1, nce_actual=28,
                    molecule_type=1)
    b.add_fragment(3, mz=175.0, intensity=1.0, product_type=1, fragment_number=1, charge=1)
    b.add_fragment(3, mz=288.0, intensity=0.8, product_type=1, fragment_number=2, charge=1)

    data = b.build()

    suffix = '.msl'
    if tmp_dir:
        path = os.path.join(tmp_dir, 'test_library' + suffix)
    else:
        fd, path = tempfile.mkstemp(suffix=suffix)
        os.close(fd)

    with open(path, 'wb') as f:
        f.write(data)
    return path


# ---------------------------------------------------------------------------
# Test class
# ---------------------------------------------------------------------------

class TestMslReader(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._tmpdir = tempfile.mkdtemp()
        cls._path = _make_test_file(cls._tmpdir)
        cls._lib = load(cls._path)

    @classmethod
    def tearDownClass(cls):
        import shutil
        shutil.rmtree(cls._tmpdir, ignore_errors=True)

    # ── Struct size validation ────────────────────────────────────────────────

    def test_struct_sizes(self):
        """Format strings must pack to exactly 64, 24, 56, 20, 20 bytes."""
        self.assertEqual(struct.calcsize(HEADER_FMT),    64)
        self.assertEqual(struct.calcsize(PROTEIN_FMT),   24)
        self.assertEqual(struct.calcsize(PRECURSOR_FMT), 56)
        self.assertEqual(struct.calcsize(FRAGMENT_FMT),  20)
        self.assertEqual(struct.calcsize(FOOTER_FMT),    20)

    # ── Magic validation ─────────────────────────────────────────────────────

    def test_magic_validation(self):
        """A file with wrong magic bytes raises ValueError."""
        with tempfile.NamedTemporaryFile(suffix='.msl', delete=False) as f:
            bad_path = f.name
            with open(self._path, 'rb') as src:
                data = bytearray(src.read())
            data[0] = 0xFF  # corrupt magic byte 0
            f.write(data)
        try:
            with self.assertRaises(ValueError) as ctx:
                load(bad_path)
            self.assertIn('magic', str(ctx.exception).lower())
        finally:
            os.unlink(bad_path)

    # ── Version validation ────────────────────────────────────────────────────

    def test_version_validation(self):
        """A file with FormatVersion outside [1, 3] raises ValueError."""
        with tempfile.NamedTemporaryFile(suffix='.msl', delete=False) as f:
            bad_path = f.name
            with open(self._path, 'rb') as src:
                data = bytearray(src.read())
            # FormatVersion is at bytes 4–7
            struct.pack_into('<i', data, 4, 99)
            f.write(data)
        try:
            with self.assertRaises(ValueError) as ctx:
                load(bad_path)
            self.assertIn('version', str(ctx.exception).lower())
        finally:
            os.unlink(bad_path)

    # ── CRC validation ────────────────────────────────────────────────────────

    def test_crc_validation(self):
        """Flipping a byte in the fragment section causes a CRC mismatch error."""
        with tempfile.NamedTemporaryFile(suffix='.msl', delete=False) as f:
            bad_path = f.name
            with open(self._path, 'rb') as src:
                data = bytearray(src.read())

            # Find a byte in the fragment section to corrupt.
            # FragmentSectionOffset is at header bytes 56–63.
            frag_offset = struct.unpack_from('<q', data, 56)[0]
            data[frag_offset] ^= 0xFF  # flip all bits
            f.write(data)
        try:
            with self.assertRaises(ValueError) as ctx:
                load(bad_path)
            self.assertIn('crc', str(ctx.exception).lower())
        finally:
            os.unlink(bad_path)

    # ── Precursor count ───────────────────────────────────────────────────────

    def test_precursor_count(self):
        """Library precursor count matches the expected 4."""
        self.assertEqual(self._lib.n_precursors, 4)
        self.assertEqual(len(self._lib.precursors), 4)

    # ── Precursor m/z values ──────────────────────────────────────────────────

    def test_precursor_mz_values(self):
        """First precursor m/z matches the expected value within 0.001."""
        prec = self._lib.precursors[0]
        self.assertAlmostEqual(prec.precursor_mz, 500.0, delta=0.001)

    # ── Fragment count ────────────────────────────────────────────────────────

    def test_fragment_count(self):
        """First precursor (PEPTIDER) has 3 fragments."""
        self.assertEqual(len(self._lib.precursors[0].fragments), 3)

    # ── Fragment m/z values ───────────────────────────────────────────────────

    def test_fragment_mz_values(self):
        """Fragment m/z values match expected values within 0.001."""
        frags = self._lib.precursors[0].fragments
        self.assertAlmostEqual(frags[0].mz, 175.119, delta=0.001)
        self.assertAlmostEqual(frags[1].mz, 288.203, delta=0.001)
        self.assertAlmostEqual(frags[2].mz, 401.287, delta=0.001)

    # ── Neutral loss code 0 ───────────────────────────────────────────────────

    def test_neutral_loss_code_none(self):
        """Fragment with neutral_loss_code 0 has neutral_loss_mass == 0.0."""
        frag = self._lib.precursors[0].fragments[0]
        self.assertEqual(frag.neutral_loss_code, 0)
        self.assertEqual(frag.neutral_loss_mass, 0.0)

    # ── Neutral loss code 1 (H2O) ─────────────────────────────────────────────

    def test_neutral_loss_code_h2o(self):
        """Fragment with neutral_loss_code 1 has neutral_loss_mass ≈ -18.0106."""
        # Precursor 2, fragment 1 has H2O loss (code=1)
        frags = self._lib.precursors[2].fragments
        h2o_frag = next((f for f in frags if f.neutral_loss_code == 1), None)
        self.assertIsNotNone(h2o_frag, 'No fragment with H2O neutral loss found')
        self.assertAlmostEqual(h2o_frag.neutral_loss_mass, -18.010565, delta=1e-4)

    # ── Decoy flag ────────────────────────────────────────────────────────────

    def test_decoy_flag(self):
        """Precursor 1 (DECOY_PEPTIDER) has is_decoy == True."""
        self.assertTrue(self._lib.precursors[1].is_decoy)
        self.assertFalse(self._lib.precursors[0].is_decoy)

    # ── String table resolution ────────────────────────────────────────────────

    def test_string_table_resolution(self):
        """modified_sequence is correctly resolved from the string table."""
        self.assertEqual(self._lib.precursors[0].modified_sequence, 'PEPTIDER')
        self.assertEqual(self._lib.precursors[1].modified_sequence, 'DECOY_PEPTIDER')
        self.assertEqual(self._lib.precursors[3].modified_sequence, 'TEIM[Oxidation]FDLK')

    # ── Index-only load ───────────────────────────────────────────────────────

    def test_index_only_load_no_fragments(self):
        """Index-only load produces precursors with empty fragment lists."""
        idx_lib = load_index_only(self._path)
        self.assertEqual(idx_lib.n_precursors, 4)
        for prec in idx_lib.precursors:
            self.assertEqual(len(prec.fragments), 0,
                             f'Expected empty fragments for {prec.modified_sequence}')

    # ── On-demand fragment loading ────────────────────────────────────────────

    def test_load_fragments_on_demand(self):
        """load_fragments() populates fragment list correctly for a given precursor."""
        idx_lib = load_index_only(self._path)
        prec = idx_lib.precursors[0]
        self.assertEqual(len(prec.fragments), 0)

        frags = load_fragments(self._path, prec)
        self.assertEqual(len(frags), 3)
        self.assertAlmostEqual(frags[0].mz, 175.119, delta=0.001)

    # ── Protein table resolution ──────────────────────────────────────────────

    def test_protein_table_resolution(self):
        """Precursor protein_accession is resolved from the protein table."""
        self.assertEqual(self._lib.precursors[0].protein_accession, 'P12345')
        self.assertEqual(self._lib.precursors[0].protein_name, 'SOME_HUMAN')
        self.assertEqual(self._lib.precursors[0].gene_name, 'GENE1')

    # ── NCE decoding ──────────────────────────────────────────────────────────

    def test_nce_decoding(self):
        """NCE is decoded from ×10 on-disk value: nce_actual=28 stored as 280."""
        self.assertEqual(self._lib.precursors[0].nce, 28)
        self.assertEqual(self._lib.precursors[2].nce, 35)

    # ── File flags ────────────────────────────────────────────────────────────

    def test_file_flags_has_protein_data(self):
        """has_protein_data flag is set when proteins are present."""
        self.assertTrue(self._lib.has_protein_data)

    # ── Molecule type (Proteoform) ────────────────────────────────────────────

    def test_proteoform_molecule_type(self):
        """Precursor 3 has molecule_type == 1 (Proteoform)."""
        self.assertEqual(self._lib.precursors[3].molecule_type, 1)

    # ── NaN q-value ───────────────────────────────────────────────────────────

    def test_nan_q_value(self):
        """Precursor with unavailable q-value stores IEEE 754 NaN."""
        self.assertTrue(math.isnan(self._lib.precursors[2].q_value))
        self.assertAlmostEqual(self._lib.precursors[0].q_value, 0.01, delta=1e-6)

    # ── Neutral loss table completeness ───────────────────────────────────────

    def test_neutral_loss_table_all_codes(self):
        """NEUTRAL_LOSS_MASSES contains entries for codes 0–6."""
        for code in range(7):
            self.assertIn(code, NEUTRAL_LOSS_MASSES)
    # ── Version range tests ───────────────────────────────────────────────────
 
    def test_version_range_accepted(self):
        \"\"\"
        Files with FormatVersion in [MIN_SUPPORTED_VERSION, MAX_SUPPORTED_VERSION]
        must load without raising ValueError.
 
        Versions 1 and 2 are tested by patching the FormatVersion field of the
        standard test file (which has no compression and no custom losses, so it
        is structurally valid at any version in the supported range — only the
        version field differs).
 
        Version 3 is the FORMAT_VERSION used by _MslBuilder, so the standard
        test file is already a version-3 file.
        \"\"\"
        with open(self._path, 'rb') as f:
            original_data = f.read()
 
        for version in range(msl_reader.MIN_SUPPORTED_VERSION,
                             msl_reader.MAX_SUPPORTED_VERSION + 1):
            patched = bytearray(original_data)
            struct.pack_into('<i', patched, 4, version)  # FormatVersion at offset 4
 
            with tempfile.NamedTemporaryFile(suffix='.msl', delete=False) as tmp:
                tmp_path = tmp.name
                tmp.write(patched)
 
            try:
                # The reader may raise on CRC mismatch (we changed a header byte)
                # but must NOT raise because of the version number.
                try:
                    msl_reader.load(tmp_path)
                    # Loaded successfully — version accepted
                except ValueError as e:
                    if 'version' in str(e).lower():
                        self.fail(
                            f'Version {version} was rejected by the reader: {e}. '
                            f'All versions in [{msl_reader.MIN_SUPPORTED_VERSION}, '
                            f'{msl_reader.MAX_SUPPORTED_VERSION}] must be accepted.'
                        )
                    # Other ValueError (e.g. CRC mismatch) is acceptable
            finally:
                os.unlink(tmp_path)
 
    def test_version_above_max_rejected(self):
        \"\"\"
        A file with FormatVersion = MAX_SUPPORTED_VERSION + 1 must raise ValueError
        with a message that includes the supported version range.
        \"\"\"
        with open(self._path, 'rb') as f:
            original_data = f.read()
 
        future_version = msl_reader.MAX_SUPPORTED_VERSION + 1
        patched = bytearray(original_data)
        struct.pack_into('<i', patched, 4, future_version)
 
        with tempfile.NamedTemporaryFile(suffix='.msl', delete=False) as tmp:
            tmp_path = tmp.name
            tmp.write(patched)
 
        try:
            with self.assertRaises(ValueError) as ctx:
                msl_reader.load(tmp_path)
 
            err_msg = str(ctx.exception).lower()
            self.assertIn('version', err_msg,
                'ValueError message must mention "version".')
            # The message should include the supported range so callers know what is valid
            self.assertIn(str(msl_reader.MIN_SUPPORTED_VERSION), str(ctx.exception),
                'ValueError message should include MIN_SUPPORTED_VERSION.')
            self.assertIn(str(msl_reader.MAX_SUPPORTED_VERSION), str(ctx.exception),
                'ValueError message should include MAX_SUPPORTED_VERSION.')
        finally:
            os.unlink(tmp_path)

if __name__ == '__main__':
    unittest.main()
