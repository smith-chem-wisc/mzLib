using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests.VariantCorpus
{
    /// <summary>
    /// Data-driven sequence-variant corpus. Each <see cref="CorpusCase"/> is a node in the possibility tree:
    /// a tiny protein spec (base sequence + mods + variants) paired with the proteoforms it MUST digest to,
    /// hand-derived from the domain "bible" (DatabaseTests/VariantCorpus/README.md) — never captured from the
    /// code. The cases live inline as records (not an external table) so the columns are compiler-checked and
    /// the input sits next to its expected answer, matching the house style in Transcriptomics/TestDigestion.cs.
    ///
    /// Two families run over the same oracle:
    ///   Corpus_Node           — build XML -> LoadProteinXML (variants expand at load) -> digest -> compare.
    ///   Corpus_Node_RoundTrip — same, but the consensus protein is serialized through the real ProteinDbWriter
    ///                           and re-read first, so encode+decode is exercised too (invariant 6 / #1083).
    /// L0 (substitutions + PTMs) + L1 (deletions: single/multi-base, before/on/after a PTM, out-of-range, and the
    /// VCF depth cutoff = minAlleleDepth). The S-series deepens the substitution axis (before/on/after a PTM,
    /// protein start/end, VCF depth cutoff, multi-residue MNV, and double substitutions) and pulls the bottom-up
    /// DIGESTION axis forward for the "a substitution moves the knife" cases (installment 5): trypsin cut-site
    /// create/destroy and the not-before-proline rule (trypsin|P). The I-series adds anchored insertions (L1).
    /// Processing is a later layer.
    /// </summary>
    [TestFixture]
    internal class VariantCorpusTests
    {
        /// <summary>
        /// One corpus node. <see cref="ExpectedForms"/> is the oracle (the exact proteoforms, in canonical
        /// least-modified-first order, in mzLib full-sequence notation). <see cref="ExpectedCount"/> is kept
        /// alongside as an independent over/under-generation guard (mirrors RnaDigestionTestCase, which keeps
        /// both DigestionProductCount and Sequences). Verdict/Reason document WHY, tied to a bible installment.
        /// </summary>
        internal record CorpusCase(
            string Id, string Layer, string Tests,
            string Base, string Mods, string Variants, string Protease,
            int MaxIsoforms, int MaxMods,
            int ExpectedCount, string Verdict, string Reason,
            string[] ExpectedForms,
            int MinAlleleDepth = 1, int MaxHeterozygous = 4,
            string Opts = "-");

        // Canonical test-mod registry (see README "Canonical test mods"). name -> (ModificationType, monoisotopicMass).
        private static readonly Dictionary<string, (string Type, double Mass)> ModRegistry = new()
        {
            ["Phosphorylation"] = ("Biological", 79.966331),
        };

        /// <summary>
        /// The corpus. L0 foundation: substitutions + PTMs only (encoding-unambiguous). Grow in layers per the
        /// README; simple nodes must keep passing as complexity is added.
        /// </summary>
        internal static IEnumerable<CorpusCase> GetCases()
        {
            // F00 — the floor: no variant, no mod -> one consensus proteoform. Caps non-binding.
            yield return new CorpusCase(
                Id: "F00", Layer: "L0a", Tests: "baseline-identity",
                Base: "PEPTIDE", Mods: "-", Variants: "-", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 1, Verdict: "applied",
                Reason: "No variant, no mod -> one consensus proteoform. Caps non-binding. The floor.",
                ExpectedForms: new[] { "PEPTIDE" });

            // F01 — variable PTM -> unmodified AND modified (0..MaxMods). 1 mod avail, cap 2 -> 2 forms.
            yield return new CorpusCase(
                Id: "F01", Layer: "L0b", Tests: "mod-only-variable",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "-", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 2, Verdict: "applied",
                Reason: "Variable PTM -> unmodified AND modified (0..MaxMods). 1 mod avail, cap 2 -> 2 forms.",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE" });

            // F02 — substitution -> consensus + variant sequence (installment 2).
            yield return new CorpusCase(
                Id: "F02", Layer: "L0a", Tests: "sub-only",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=T VAR=V POS=4 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 2, Verdict: "applied",
                Reason: "Substitution -> consensus + variant sequence (installment 2).",
                ExpectedForms: new[] { "PEPTIDE", "PEPVIDE" });

            // F03 — sub x mod interaction: T->V drops the phospho at the edited position; no phosphorylated
            // variant produced. 3 forms NOT 4 (installment 2; invariants 1,2).
            yield return new CorpusCase(
                Id: "F03", Layer: "L0-int", Tests: "sub-x-mod-legal-consensus-only",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=T VAR=V POS=4 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 3, Verdict: "pruned:illegal-mod",
                Reason: "T->V drops the phospho at the edited position; no phosphorylated variant produced. 3 forms NOT 4 (installment 2; invariants 1,2).",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE", "PEPVIDE" });

            // F04 — boundary: position 9 > length 7 -> variant pruned; only consensus (reader ~line 776;
            // outOfRangeVariant.xml).
            yield return new CorpusCase(
                Id: "F04", Layer: "L0a", Tests: "sub-out-of-range",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=- VAR=A POS=9 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 1, Verdict: "pruned:out-of-range",
                Reason: "Position 9 > length 7 -> variant pruned; only consensus (reader ~line 776; outOfRangeVariant.xml).",
                ExpectedForms: new[] { "PEPTIDE" });

            // ---- L1: deletions (installments 3, 10, 12) ---------------------------------------------------
            // UniProt-native deletions are encoded as an EMPTY <variation> (VAR omitted here). Until the reader
            // gap at ProteinXmlEntry.cs:588 is fixed these D0x nodes are RED (the deletion is silently dropped,
            // so only the consensus forms are produced); the oracle below is the biologically-correct answer per
            // the bible, so the fix turns them green. D05/D06 (out-of-range) and D07/D08 (VCF depth) are green
            // regardless of that fix.

            // D00 — single-base deletion, unmodified. Delete P3 -> PETIDE.
            yield return new CorpusCase(
                Id: "D00", Layer: "L1a", Tests: "del-single-unmod",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=P POS=3 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 2, Verdict: "applied",
                Reason: "Single-base deletion of P3 -> consensus + PETIDE. UniProt-native (empty <variation>); RED until ProteinXmlEntry.cs:588 fix.",
                ExpectedForms: new[] { "PEPTIDE", "PETIDE" });

            // D01 — multi-base deletion, unmodified. Delete E2-P3 -> PTIDE.
            yield return new CorpusCase(
                Id: "D01", Layer: "L1a", Tests: "del-multi-unmod",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=EP BEGIN=2 END=3 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 2, Verdict: "applied",
                Reason: "Multi-base deletion of E2-P3 -> consensus + PTIDE (begin/end range). RED until reader fix.",
                ExpectedForms: new[] { "PEPTIDE", "PTIDE" });

            // D02 — deletion BEFORE the PTM residue. Phospho@T4; delete P3 -> T slides 4->3, phospho follows
            // (installment 3). 4 forms.
            yield return new CorpusCase(
                Id: "D02", Layer: "L1-int", Tests: "del-before-ptm",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=P POS=3 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 4, Verdict: "applied",
                Reason: "Deletion upstream of the mod re-anchors the phospho with its residue (T 4->3). Consensus{0,1 phospho} + variant{0,1 phospho} = 4. Installment 3; AdjustModificationIndices k>end shift.",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE", "PETIDE", "PET[Biological:Phosphorylation on T]IDE" });

            // D03 — deletion ON the PTM residue. Phospho@T4; delete T4 -> anchor gone, no phospho on the variant
            // (installment 3; D7 positional drop). 3 forms NOT 4.
            yield return new CorpusCase(
                Id: "D03", Layer: "L1-int", Tests: "del-on-ptm",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=T POS=4 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 3, Verdict: "pruned:illegal-mod",
                Reason: "Deleting the modified residue drops the phospho on the variant (anchor gone), kept on consensus. 3 forms NOT 4. Installment 3; AdjustModificationIndices begin<=k<=end -> dropped (D7).",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE", "PEPIDE" });

            // D04 — deletion AFTER the PTM residue. Phospho@T4; delete I5 -> T4 unmoved, phospho stays. 4 forms.
            yield return new CorpusCase(
                Id: "D04", Layer: "L1-int", Tests: "del-after-ptm",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=I POS=5 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 4, Verdict: "applied",
                Reason: "Deletion downstream of the mod leaves the phospho anchored at T4. Consensus{0,1} + variant PEPTDE{0,1} = 4. AdjustModificationIndices k<begin -> unaffected.",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE", "PEPTDE", "PEPT[Biological:Phosphorylation on T]DE" });

            // D05 — deletion OUT OF BOUNDS (single). Position 9 > length 7 -> pruned; consensus only. Green regardless.
            yield return new CorpusCase(
                Id: "D05", Layer: "L1a", Tests: "del-out-of-range",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=A POS=9 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 1, Verdict: "pruned:out-of-range",
                Reason: "Deletion position 9 > length 7 -> PruneOutOfRangeSequenceVariants removes it; consensus only.",
                ExpectedForms: new[] { "PEPTIDE" });

            // D06 — deletion straddling the C-terminus. Original "E" = residue 7 (the last real residue) and end 8
            // is one past length 7, so the SOLE reason for pruning is end > len (the original matches residues in
            // range; it is not malformed on any other axis). This isolates PruneOutOfRangeSequenceVariants' end-
            // beyond-length branch cleanly. Consensus only.
            yield return new CorpusCase(
                Id: "D06", Layer: "L1a", Tests: "del-straddle-end",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=E BEGIN=7 END=8 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 1, Verdict: "pruned:out-of-range",
                Reason: "Deletion end 8 > length 7 (begin 7 = last residue E; original matches residues[begin,end] in range) -> pruned by end>len alone; consensus only.",
                ExpectedForms: new[] { "PEPTIDE" });

            // D07 — VCF deletion PASSING the depth cutoff. Heterozygous 0/1, alt AD=6 >= minAlleleDepth 5 ->
            // keep consensus AND apply the deletion (installment 12; ApplyVariants het path). 2 forms.
            yield return new CorpusCase(
                Id: "D07", Layer: "L1-vcf", Tests: "del-depth-pass",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=EP VAR=E BEGIN=2 END=3 SRC=vcf GT=0/1 AD=8,6 DP=14", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 2, Verdict: "applied",
                Reason: "VCF-anchored deletion of P3; het alt AD=6 >= minAlleleDepth=5 -> consensus + PETIDE. Depth filter = minAlleleDepth (VCF AD).",
                ExpectedForms: new[] { "PEPTIDE", "PETIDE" },
                MinAlleleDepth: 5);

            // D08 — same VCF deletion FAILING the depth cutoff. alt AD=6 < minAlleleDepth 7 -> allele filtered;
            // consensus only. The depth cutoff, encoded from the sequencing data (AD), drops the deletion.
            yield return new CorpusCase(
                Id: "D08", Layer: "L1-vcf", Tests: "del-depth-fail",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=EP VAR=E BEGIN=2 END=3 SRC=vcf GT=0/1 AD=8,6 DP=14", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 1, Verdict: "filtered:below-depth-cutoff",
                Reason: "Same VCF deletion; alt AD=6 < minAlleleDepth=7 -> isDeepAlternateAllele false, deletion not applied; consensus only.",
                ExpectedForms: new[] { "PEPTIDE" },
                MinAlleleDepth: 7);

            // D09 — anchored deletion whose KEPT anchor carries the PTM. Phospho@T4; TI->T over [4,5] deletes I5 and
            // keeps T4 unchanged, so the phospho survives on the variant (invariant 2; D7 retained-anchor half,
            // settled 2026-09-16). Consensus{0,1} + PEPTDE{0,1} = 4. Same proteoforms as D04, reached through the
            // anchored encoding, which puts the modified residue inside [begin,end].
            yield return new CorpusCase(
                Id: "D09", Layer: "L1-int", Tests: "del-anchored-on-ptm-kept",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=TI VAR=T BEGIN=4 END=5 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 4, Verdict: "applied",
                Reason: "Anchored deletion TI->T keeps T4 and removes I5; the mod on the kept anchor survives even though position 4 lies inside [begin,end]. Consensus{0,1} + PEPTDE{0,1} = 4, the same forms as D04. Deletion half of the I01 retained-anchor rule.",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE", "PEPTDE", "PEPT[Biological:Phosphorylation on T]DE" });

            // ---- L1 multi: two SEPARATE deletions in one protein (combinatorial path; invariants 3-4) -------
            // Two genotype-less deletions expand to consensus + each single + the pair (least-modified-first,
            // applied descending-position). Kept to PAIRS so base+2+1 = 4 stays under the default cap (no
            // throttle/determinism-under-cap yet — that's a later layer). Each combo re-anchors mods independently.

            // MD00 — two single-base deletions, unmodified. Delete P3 and D6 -> {consensus, PETIDE, PEPTIE, PETIE}.
            yield return new CorpusCase(
                Id: "MD00", Layer: "L1-multi", Tests: "del-x2-unmod",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=P POS=3 SRC=uniprot; OP=D POS=6 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 4, Verdict: "applied",
                Reason: "Two independent single-base deletions -> consensus + each single + the pair. base+2+1=4 (under cap).",
                ExpectedForms: new[] { "PEPTIDE", "PETIDE", "PEPTIE", "PETIE" });

            // MD01 — two deletions FLANKING a PTM. Phospho@T4; delete P3 (before) and I5 (after). Phospho
            // re-anchors to T per combo (4->3 whenever P3 is deleted). 4 isoforms x variable phospho = 8 forms.
            yield return new CorpusCase(
                Id: "MD01", Layer: "L1-multi", Tests: "del-x2-flanking-ptm",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=P POS=3 SRC=uniprot; OP=I POS=5 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 8, Verdict: "applied",
                Reason: "Deletions before (P3) and after (I5) the mod; phospho survives all 4 combos, re-anchored to T (pos 4 or 3). 4 isoforms x {0,1 phospho} = 8.",
                ExpectedForms: new[]
                {
                    "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE",
                    "PETIDE", "PET[Biological:Phosphorylation on T]IDE",
                    "PEPTDE", "PEPT[Biological:Phosphorylation on T]DE",
                    "PETDE", "PET[Biological:Phosphorylation on T]DE"
                });

            // MD02 — two deletions where ONE removes the PTM residue. Phospho@T4; delete E2 (before) and T4 (on).
            // Phospho survives combos WITHOUT the T4 deletion (consensus, delE2), dropped in combos WITH it
            // (delT4, delE2+delT4). 2+2+1+1 = 6 forms.
            yield return new CorpusCase(
                Id: "MD02", Layer: "L1-multi", Tests: "del-x2-one-on-ptm",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=E POS=2 SRC=uniprot; OP=T POS=4 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 6, Verdict: "applied",
                Reason: "Phospho kept where T4 survives (consensus{0,1}, delE2{0,1}), dropped where T4 is deleted (delT4, delE2+delT4). 2+2+1+1 = 6.",
                ExpectedForms: new[]
                {
                    "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE",
                    "PPTIDE", "PPT[Biological:Phosphorylation on T]IDE",
                    "PEPIDE", "PPIDE"
                });

            // ---- L1: insertions (installment 3; invariants 2-3) --------------------------------------------
            // UniProt writes an insertion ANCHORED on a retained residue: <original>P</original><variation>PAG
            // </variation> = "insert AG after P". The variant is longer than the original, so everything after the
            // anchor shifts RIGHT by the inserted length (the mirror of a deletion) and surviving mods ride along.
            // I00-I02 walk a two-residue insertion past a PTM at both boundaries: anchored immediately before the
            // modified residue, ON it (the modified residue is the retained anchor), and immediately after it.
            // I03-I04 put the insertion at the protein's start and end; I05-I06 are the VCF depth-cutoff twins;
            // I07 moves the trypsin knife; I08 carries a variant-borne mod; MI00-MI01 are indel pairs.

            // I00 — insertion immediately BEFORE the PTM residue. Phospho@T4; P3->PAG inserts AG between P3 and T4.
            // T slides 4->6 and the phospho follows it (installment 3). Consensus{0,1} + PEPAGTIDE{0,1} = 4 forms.
            yield return new CorpusCase(
                Id: "I00", Layer: "L1-int", Tests: "ins-before-ptm",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=P VAR=PAG POS=3 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 4, Verdict: "applied",
                Reason: "Insertion upstream of the mod re-anchors the phospho with its residue (T 4->6, +2). Mod sits at end+1 of the edit, the shift boundary. Consensus{0,1} + PEPAGTIDE{0,1} = 4. Installment 3; mirror of D02.",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE", "PEPAGTIDE", "PEPAGT[Biological:Phosphorylation on T]IDE" });

            // I01 — insertion ANCHORED ON the PTM residue. Phospho@T4; T4->TAG keeps T4 where it is and inserts AG
            // after it. The anchor residue survives unchanged (same residue, same position), so the phospho survives
            // on the variant (invariant 2). 4 forms. AdjustModificationIndices used to drop every mod inside the edited
            // [begin,end], including a retained anchor; it now keeps mods on the unchanged prefix/suffix.
            yield return new CorpusCase(
                Id: "I01", Layer: "L1-int", Tests: "ins-anchored-on-ptm",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=T VAR=TAG POS=4 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 4, Verdict: "applied",
                Reason: "T4->TAG retains T4 (same residue, same position) and inserts AG after it, so the phospho's anchor survives and the mod stays at 4 (invariant 2). Consensus{0,1} + PEPTAGIDE{0,1} = 4. Pins the retained-anchor half of D7: only residues the variant actually changes lose their mods.",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE", "PEPTAGIDE", "PEPT[Biological:Phosphorylation on T]AGIDE" });

            // I02 — insertion immediately AFTER the PTM residue. Phospho@T4; I5->IAG inserts AG after I5. T4 is
            // before the edit, so it is unmoved and keeps its phospho. Consensus{0,1} + PEPTIAGDE{0,1} = 4 forms.
            yield return new CorpusCase(
                Id: "I02", Layer: "L1-int", Tests: "ins-after-ptm",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=I VAR=IAG POS=5 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 4, Verdict: "applied",
                Reason: "Insertion downstream of the mod leaves the phospho anchored at T4 (mod at begin-1, the unaffected boundary). Consensus{0,1} + PEPTIAGDE{0,1} = 4. AdjustModificationIndices k<begin -> unaffected; mirror of D04.",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE", "PEPTIAGDE", "PEPT[Biological:Phosphorylation on T]IAGDE" });

            // I03 — insertion at the protein START, ahead of the first residue. Phospho@T4; P1->AGP keeps P1 as the
            // right-hand anchor and puts AG in front of it, so the new protein begins AG. Every residue shifts +2 and
            // T4->6 carries its phospho. M-free base + Retain, so no initiator-Met interplay (that node is parked).
            // Consensus{0,1} + AGPEPTIDE{0,1} = 4 forms.
            yield return new CorpusCase(
                Id: "I03", Layer: "L1-int", Tests: "ins-at-start",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=P VAR=AGP POS=1 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 4, Verdict: "applied",
                Reason: "N-terminal insertion (P1->AGP, AG before the first residue) shifts the whole protein +2; the phospho follows T 4->6. Consensus{0,1} + AGPEPTIDE{0,1} = 4. Start boundary of the insertion axis (mirror of S03); M-free base + Retain.",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE", "AGPEPTIDE", "AGPEPT[Biological:Phosphorylation on T]IDE" });

            // I04 — insertion at the protein END, after the last residue. Phospho@T4; E7->EAG extends the C-terminus
            // by AG. Nothing sits downstream of the edit, so T4 and its phospho are unmoved and the variant is two
            // residues longer than the consensus. Consensus{0,1} + PEPTIDEAG{0,1} = 4 forms.
            yield return new CorpusCase(
                Id: "I04", Layer: "L1-int", Tests: "ins-at-end",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=E VAR=EAG POS=7 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 4, Verdict: "applied",
                Reason: "C-terminal insertion (E7->EAG) extends the protein past its consensus length; the anchor is the last residue so nothing shifts and the phospho stays at T4. Consensus{0,1} + PEPTIDEAG{0,1} = 4. End boundary of the insertion axis (mirror of S04; contrast D06, where end > length is pruned).",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE", "PEPTIDEAG", "PEPT[Biological:Phosphorylation on T]IDEAG" });

            // I05 — VCF insertion PASSING the depth cutoff. VCF anchors on the reference base: REF=E ALT=EAG at 2
            // inserts AG after E2 -> PEAGPTIDE. Het 0/1, alt AD=6 >= minAlleleDepth 5 -> keep consensus AND apply
            // (installment 12; ApplyVariants het path). Phospho@T4 sits downstream, so it shifts 4->6 through the
            // VCF path, which re-anchors mods separately from the UniProt combinatorial path I00 covers.
            // Consensus{0,1} + PEAGPTIDE{0,1} = 4 forms.
            yield return new CorpusCase(
                Id: "I05", Layer: "L1-vcf", Tests: "ins-depth-pass",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=E VAR=EAG POS=2 SRC=vcf GT=0/1 AD=8,6 DP=14", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 4, Verdict: "applied",
                Reason: "VCF-anchored insertion (REF=E ALT=EAG at 2); het alt AD=6 >= minAlleleDepth=5 -> consensus + PEAGPTIDE, with the phospho re-anchored T 4->6 on the variant. Consensus{0,1} + variant{0,1} = 4. Depth filter = minAlleleDepth (VCF AD); mirror of D07.",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE", "PEAGPTIDE", "PEAGPT[Biological:Phosphorylation on T]IDE" },
                MinAlleleDepth: 5);

            // I06 — same VCF insertion FAILING the depth cutoff. alt AD=6 < minAlleleDepth 7 -> allele filtered; only
            // the consensus and its phospho form survive. The variant is still READ, so the round-trip family checks
            // the writer preserves an unapplied insertion (as D08/S06 do for their variant types).
            yield return new CorpusCase(
                Id: "I06", Layer: "L1-vcf", Tests: "ins-depth-fail",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=E VAR=EAG POS=2 SRC=vcf GT=0/1 AD=8,6 DP=14", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 2, Verdict: "filtered:below-depth-cutoff",
                Reason: "Same VCF insertion; alt AD=6 < minAlleleDepth=7 -> isDeepAlternateAllele false, insertion not applied; consensus{0,1} only. Mirror of D08.",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE" },
                MinAlleleDepth: 7);

            // I07 — insertion CREATES a trypsin cut site (installment 5; invariant 7). A4->AK inserts a K after A4.
            // Consensus PEPAIDE (no K/R) -> {PEPAIDE}; variant PEPAKIDE cuts after K5 -> {PEPAK, IDE}. Union = 3.
            yield return new CorpusCase(
                Id: "I07", Layer: "L3-digest", Tests: "ins-creates-cut-site",
                Base: "PEPAIDE", Mods: "-", Variants: "OP=A VAR=AK POS=4 SRC=uniprot", Protease: "trypsin",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 3, Verdict: "applied",
                Reason: "Inserting K after A4 gives PEPAKIDE, which trypsin cuts after K5 -> PEPAK + IDE; consensus PEPAIDE has no K/R -> one peptide. An inserted K/R adds a knife (installment 5). Insertion analog of S08.",
                ExpectedForms: new[] { "PEPAIDE", "PEPAK", "IDE" });

            // I08 — insertion that CARRIES A MOD (variant-borne; installment 4, #1083). P3->PT inserts a T after P3,
            // turning PEPIDE into PEPTIDE, and a phospho is stored ON THE VARIANT at the inserted T (subposition 4,
            // in the applied frame). The consensus PEPIDE has no T and stays bare; only the applied variant carries
            // the variable phospho. 3 forms. Unlike S11 the mod sits past the variant's own position (4 vs 3), so it
            // also checks the stored subposition is read in the applied frame. The round-trip family checks the
            // writer keeps a mod that lives on an insertion.
            yield return new CorpusCase(
                Id: "I08", Layer: "L1-int", Tests: "ins-carries-variant-borne-mod",
                Base: "PEPIDE", Mods: "-", Variants: "OP=P VAR=PT POS=3 SRC=uniprot VMOD=Phosphorylation@4", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 3, Verdict: "applied",
                Reason: "P3->PT inserts a phosphorylatable T at applied position 4, and the phospho is stored on the variant (SequenceVariation.OneBasedModifications). Consensus PEPIDE is bare; applied PEPTIDE carries a variable phospho. 3 forms. Insertion analog of S11 (installment 4; AdjustModificationIndices merges variant.OneBasedModifications).",
                ExpectedForms: new[] { "PEPIDE", "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE" });

            // ---- L1 multi: two indels in one protein (combinatorial path; invariants 3-4) ----------------------
            // Genotype-less pairs expand to consensus + each single + the pair, applied descending-position (D6),
            // so each step re-frames the coordinates the next one works in. Kept to pairs (base+2+1 = 4 isoforms).

            // MI00 — two insertions FLANKING a PTM. Phospho@T4; E2->EAG (before) and I5->IAG (after). T sits at 6
            // whenever the E2 insertion is applied, else at 4. 4 isoforms x {0,1 phospho} = 8 forms. Mirror of MD01.
            yield return new CorpusCase(
                Id: "MI00", Layer: "L1-multi", Tests: "ins-x2-flanking-ptm",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=E VAR=EAG POS=2 SRC=uniprot; OP=I VAR=IAG POS=5 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 8, Verdict: "applied",
                Reason: "Insertions before (E2) and after (I5) the mod; the phospho survives all 4 combos, at T6 when E2->EAG is applied and T4 otherwise. Pair applied I5 first, then E2 -> PEAGPTIAGDE. 4 isoforms x {0,1 phospho} = 8.",
                ExpectedForms: new[]
                {
                    "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE",
                    "PEAGPTIDE", "PEAGPT[Biological:Phosphorylation on T]IDE",
                    "PEPTIAGDE", "PEPT[Biological:Phosphorylation on T]IAGDE",
                    "PEAGPTIAGDE", "PEAGPT[Biological:Phosphorylation on T]IAGDE"
                });

            // MI01 — an insertion AND a deletion, both upstream of the PTM, with opposite shifts. Phospho@T4;
            // P1->PA (+1) and delete P3 (-1). Singly they move T to 5 and 3; together they cancel and T ends at 4
            // again, but only if each step re-anchors in the frame the previous step left (applied P3 first -> PETIDE,
            // T3; then P1->PA -> PAETIDE, T4). 4 isoforms x {0,1 phospho} = 8 forms.
            yield return new CorpusCase(
                Id: "MI01", Layer: "L1-multi", Tests: "ins-plus-del-net-zero",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=P VAR=PA POS=1 SRC=uniprot; OP=P POS=3 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 8, Verdict: "applied",
                Reason: "Insertion P1->PA (+1) and deletion of P3 (-1), both before T4. Singles: PAEPTIDE (T5), PETIDE (T3). Pair: PAETIDE with T back at 4, after two opposite re-anchorings (invariants 3-4). 4 isoforms x {0,1 phospho} = 8.",
                ExpectedForms: new[]
                {
                    "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE",
                    "PAEPTIDE", "PAEPT[Biological:Phosphorylation on T]IDE",
                    "PETIDE", "PET[Biological:Phosphorylation on T]IDE",
                    "PAETIDE", "PAET[Biological:Phosphorylation on T]IDE"
                });

            // ---- S-series: the substitution axis, deepened (installments 2, 5, 11; invariants 1-2, 5) ---------
            // Substitutions HOLD position (no downstream renumbering, unlike indels), so a mod survives iff the
            // substitution does not land ON its anchor residue. S00-S02 walk a sub past a PTM (before/on/after);
            // S03-S04 hit the protein termini; S05-S06 are the VCF depth-cutoff twins; S07-S09 move the trypsin
            // knife (bottom-up); S10 is a multi-residue MNV; SS00-SS01 are double substitutions.

            // S00 — substitution BEFORE the PTM. Phospho@T4; sub E2->A. T4 unmoved (subs hold position), phospho
            // survives on the variant. Consensus{0,1} + variant PAPTIDE{0,1} = 4 forms.
            yield return new CorpusCase(
                Id: "S00", Layer: "L0-sub-int", Tests: "sub-before-ptm",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=E VAR=A POS=2 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 4, Verdict: "applied",
                Reason: "Substitution upstream of the mod holds position (no renumbering); phospho stays anchored at T4. Consensus{0,1} + PAPTIDE{0,1} = 4. Installment 2; AdjustModificationIndices leaves k outside [begin,end] untouched.",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE", "PAPTIDE", "PAPT[Biological:Phosphorylation on T]IDE" });

            // S01 — substitution ON the PTM residue, motif-preserving-LOOKING (T->S). The phospho is dropped on the
            // variant anyway: positionally (D7 — position 4 lies in the edited [4,4]) AND motif-wise (our annotated
            // mod is "Phosphorylation on T"; serine does not satisfy a T-specific motif). Consensus keeps it. 3
            // forms NOT 4. This pins the D7 case for a RESIDUE-SPECIFIC motif (both rules agree -> drop); the truly
            // ambiguous D7 (a mod whose motif would still fit post-sub, e.g. a phospho annotated "on S or T") stays
            // parked for a later node/decision (installment 11).
            yield return new CorpusCase(
                Id: "S01", Layer: "L0-sub-int", Tests: "sub-on-ptm-motif-preserving",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=T VAR=S POS=4 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 3, Verdict: "pruned:illegal-mod",
                Reason: "T->S drops the phospho on the variant: positional drop (pos 4 in edited [4,4]) and our motif is T-specific so S does not fit either. Consensus keeps it. 3 forms NOT 4. Pins D7 for a residue-specific motif; the motif-still-fits D7 case is parked (installment 11).",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE", "PEPSIDE" });

            // S02 — substitution AFTER the PTM. Phospho@T4; sub D6->A. T4 unmoved, phospho survives. 4 forms.
            yield return new CorpusCase(
                Id: "S02", Layer: "L0-sub-int", Tests: "sub-after-ptm",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=D VAR=A POS=6 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 4, Verdict: "applied",
                Reason: "Substitution downstream of the mod holds position; phospho stays anchored at T4. Consensus{0,1} + PEPTIAE{0,1} = 4.",
                ExpectedForms: new[] { "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE", "PEPTIAE", "PEPT[Biological:Phosphorylation on T]IAE" });

            // S03 — substitution at the protein START (position 1). P1->A. No mod. Consensus + AEPTIDE = 2 forms.
            yield return new CorpusCase(
                Id: "S03", Layer: "L0-sub", Tests: "sub-at-start",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=P VAR=A POS=1 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 2, Verdict: "applied",
                Reason: "Substitution at the first residue (P1->A). No initiator-Met interplay (M-free base, Retain behavior). Consensus + AEPTIDE = 2.",
                ExpectedForms: new[] { "PEPTIDE", "AEPTIDE" });

            // S04 — substitution at the protein END (position 7 = last residue). E7->A. Consensus + PEPTIDA = 2 forms.
            yield return new CorpusCase(
                Id: "S04", Layer: "L0-sub", Tests: "sub-at-end",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=E VAR=A POS=7 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 2, Verdict: "applied",
                Reason: "Substitution at the last residue (E7->A). In-range boundary (contrast F04, pos 9 > len). Consensus + PEPTIDA = 2.",
                ExpectedForms: new[] { "PEPTIDE", "PEPTIDA" });

            // S05 — VCF substitution PASSING the depth cutoff. T4->V, het 0/1, alt AD=6 >= minAlleleDepth 5 ->
            // consensus AND variant (installment 12; ApplyVariants het path). 2 forms.
            yield return new CorpusCase(
                Id: "S05", Layer: "L0-sub-vcf", Tests: "sub-depth-pass",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=T VAR=V POS=4 SRC=vcf GT=0/1 AD=8,6 DP=14", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 2, Verdict: "applied",
                Reason: "VCF substitution T4->V; het alt AD=6 >= minAlleleDepth=5 -> consensus + PEPVIDE. Depth filter = minAlleleDepth (VCF AD).",
                ExpectedForms: new[] { "PEPTIDE", "PEPVIDE" },
                MinAlleleDepth: 5);

            // S06 — same VCF substitution FAILING the depth cutoff. alt AD=6 < minAlleleDepth 7 -> allele filtered;
            // consensus only. The variant is still READ (SequenceVariation created), so the round-trip family asserts
            // the writer preserves it (teeth, as with D08).
            yield return new CorpusCase(
                Id: "S06", Layer: "L0-sub-vcf", Tests: "sub-depth-fail",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=T VAR=V POS=4 SRC=vcf GT=0/1 AD=8,6 DP=14", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 1, Verdict: "filtered:below-depth-cutoff",
                Reason: "Same VCF substitution; alt AD=6 < minAlleleDepth=7 -> isDeepAlternateAllele false, not applied; consensus only.",
                ExpectedForms: new[] { "PEPTIDE" },
                MinAlleleDepth: 7);

            // ---- Bottom-up: a substitution moves the trypsin knife (installment 5; invariant 5) ---------------
            // These are the first DIGESTION-axis nodes: protease = trypsin, so consensus and variant digest
            // independently and the peptide SET is a function of which variant is applied. FullSequence of an
            // unmodified peptide == its base sequence.

            // S07 — substitution DESTROYS a cut site (K4->A). trypsin cuts after K/R: consensus PEPKIDE -> {PEPK, IDE};
            // variant PEPAIDE has no K/R -> {PEPAIDE}. Union = 3 peptides.
            yield return new CorpusCase(
                Id: "S07", Layer: "L3-digest", Tests: "sub-destroys-cut-site",
                Base: "PEPKIDE", Mods: "-", Variants: "OP=K VAR=A POS=4 SRC=uniprot", Protease: "trypsin",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 3, Verdict: "applied",
                Reason: "K4->A removes trypsin's cut after K4: consensus PEPKIDE -> PEPK + IDE; variant PEPAIDE (no K/R) -> one peptide. K/R->X merges two peptides into one (installment 5).",
                ExpectedForms: new[] { "PEPK", "IDE", "PEPAIDE" });

            // S08 — substitution CREATES a cut site (A4->K), the inverse of S07. Consensus PEPAIDE (no K/R) ->
            // {PEPAIDE}; variant PEPKIDE -> {PEPK, IDE}. Union = 3 peptides.
            yield return new CorpusCase(
                Id: "S08", Layer: "L3-digest", Tests: "sub-creates-cut-site",
                Base: "PEPAIDE", Mods: "-", Variants: "OP=A VAR=K POS=4 SRC=uniprot", Protease: "trypsin",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 3, Verdict: "applied",
                Reason: "A4->K introduces a trypsin cut after position 4: consensus PEPAIDE -> one peptide; variant PEPKIDE -> PEPK + IDE. X->K/R splits one peptide into two (installment 5). Inverse of S07.",
                ExpectedForms: new[] { "PEPAIDE", "PEPK", "IDE" });

            // S09 — substitution BLOCKS a cut site via the not-before-proline rule (I5->P), under trypsin|P
            // (K[P]|,R[P]|). Consensus PEPKIDE -> {PEPK, IDE} (K4 before I, cuts); variant PEPKPDE -> {PEPKPDE}
            // (K4 now before P, cut suppressed). Under plain trypsin the variant would still cut (PEPK + PDE), so
            // this node specifically pins trypsin|P's proline restriction, not just any K-adjacent edit.
            yield return new CorpusCase(
                Id: "S09", Layer: "L3-digest", Tests: "sub-blocks-cut-before-proline",
                Base: "PEPKIDE", Mods: "-", Variants: "OP=I VAR=P POS=5 SRC=uniprot", Protease: "trypsin|P",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 3, Verdict: "applied",
                Reason: "I5->P puts a proline immediately after K4; trypsin|P (K[P]|) does not cut before P, so variant PEPKPDE stays whole while consensus PEPKIDE -> PEPK + IDE. A variant destroys a cut site WITHOUT touching K/R (installment 5). Requires trypsin|P, not trypsin.",
                ExpectedForms: new[] { "PEPK", "IDE", "PEPKPDE" });

            // S10 — MULTI-RESIDUE substitution (MNV), same length. TI->VL over [4,5] -> PEPVLDE. Exercises the
            // begin/end range path for substitutions (equal length, no coordinate shift). 2 forms.
            yield return new CorpusCase(
                Id: "S10", Layer: "L0-sub", Tests: "sub-multi-residue-mnv",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=TI VAR=VL BEGIN=4 END=5 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 2, Verdict: "applied",
                Reason: "Two-residue substitution TI->VL over [4,5] (equal length -> no renumbering). Consensus + PEPVLDE = 2. Tests the range (begin/end) apply path for a same-length variant.",
                ExpectedForms: new[] { "PEPTIDE", "PEPVLDE" });

            // S11 — substitution that CARRIES A MOD into the proteoform (variant-borne mod; installment 4, #1083).
            // The reciprocal of F03: there the consensus bears the phospho and the T->V variant loses it; here the
            // CONSENSUS residue (V) cannot bear phospho, and a phospho is stored ON THE VARIANT's new T (V4->T). So
            // the consensus yields one bare form, and only the applied variant carries the (variable) phospho.
            // 3 forms: PEPVIDE, PEPTIDE, PEPT[phospho]IDE. Encoded via the nested <subfeature>/<subposition>; the
            // round-trip family gives this real teeth on variant-borne-mod SERIALIZATION (a writer that drops the
            // sub-feature would produce only 2 forms on reload).
            yield return new CorpusCase(
                Id: "S11", Layer: "L0-sub-int", Tests: "sub-carries-variant-borne-mod",
                Base: "PEPVIDE", Mods: "-", Variants: "OP=V VAR=T POS=4 SRC=uniprot VMOD=Phosphorylation@4", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 3, Verdict: "applied",
                Reason: "V4->T introduces a phosphorylatable T; the phospho is stored on the VARIANT (SequenceVariation.OneBasedModifications), so the consensus PEPVIDE is bare while the applied PEPTIDE carries a variable phospho. Reciprocal of F03 (installment 4; AdjustModificationIndices merges variant.OneBasedModifications). 3 forms.",
                ExpectedForms: new[] { "PEPVIDE", "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE" });

            // ---- L4 multi: double substitutions (combinatorial path; invariants 3-4) --------------------------

            // SS00 — two INDEPENDENT substitutions. E2->A and D6->A -> consensus + each single + the pair. Subs hold
            // position so no re-anchoring; base+2+1 = 4 forms (under cap). Mirror of the deletion MD00.
            yield return new CorpusCase(
                Id: "SS00", Layer: "L4-multi", Tests: "sub-x2-independent",
                Base: "PEPTIDE", Mods: "-", Variants: "OP=E VAR=A POS=2 SRC=uniprot; OP=D VAR=A POS=6 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 4, Verdict: "applied",
                Reason: "Two independent substitutions -> consensus + each single + the pair. base+2+1 = 4 (under cap). No re-anchoring (subs hold position).",
                ExpectedForms: new[] { "PEPTIDE", "PAPTIDE", "PEPTIAE", "PAPTIAE" });

            // SS01 — two substitutions where ONE lands on the PTM residue. Phospho@T4; sub E2->A (before) and
            // T4->V (on). Phospho survives the combos WITHOUT the T4 sub (consensus, subE2), dropped in the combos
            // WITH it (subT4, both). 2+2+1+1 = 6 forms. Mirror of the deletion MD02.
            yield return new CorpusCase(
                Id: "SS01", Layer: "L4-multi", Tests: "sub-x2-one-on-ptm",
                Base: "PEPTIDE", Mods: "Phosphorylation@4", Variants: "OP=E VAR=A POS=2 SRC=uniprot; OP=T VAR=V POS=4 SRC=uniprot", Protease: "top-down",
                MaxIsoforms: 1024, MaxMods: 2,
                ExpectedCount: 6, Verdict: "applied",
                Reason: "Phospho kept where T4 survives (consensus{0,1}, subE2{0,1}), dropped where T4->V (subT4, both). 2+2+1+1 = 6.",
                ExpectedForms: new[]
                {
                    "PEPTIDE", "PEPT[Biological:Phosphorylation on T]IDE",
                    "PAPTIDE", "PAPT[Biological:Phosphorylation on T]IDE",
                    "PEPVIDE", "PAPVIDE"
                });
        }

        private static IEnumerable<TestCaseData> Cases()
        {
            foreach (var c in GetCases())
                yield return new TestCaseData(c).SetName($"Corpus_{c.Id}_{c.Tests}");
        }

        [TestCaseSource(nameof(Cases))]
        public void Corpus_Node(CorpusCase c)
        {
            AssertCaseSelfConsistent(c);

            List<string> produced = RunNode(c);
            List<string> producedAgain = RunNode(c);

            // Determinism (invariant 11): identical output, same order, across runs.
            Assert.That(producedAgain, Is.EqualTo(produced),
                $"{c.Id}: nondeterministic output across two runs.");

            // Count (fast over/under-generation check).
            Assert.That(produced.Count, Is.EqualTo(c.ExpectedCount),
                $"{c.Id} ({c.Tests}): expected {c.ExpectedCount} forms, got {produced.Count}.\n  got: {string.Join(" | ", produced)}");

            // Membership (order-independent for now; canonical-order assertion is a later tightening).
            Assert.That(produced.OrderBy(x => x, StringComparer.Ordinal),
                Is.EqualTo(c.ExpectedForms.OrderBy(x => x, StringComparer.Ordinal)),
                $"{c.Id} ({c.Tests}): form set mismatch.\n  expected: {string.Join(" | ", c.ExpectedForms)}\n  got:      {string.Join(" | ", produced)}");
        }

        /// <summary>
        /// Round-trip variant: same oracle, but the CONSENSUS protein is serialized through the REAL
        /// ProteinDbWriter and re-read before digestion. This exercises encode+decode (invariant 6 / #1083)
        /// WITHOUT the writer ever sitting upstream of the ground truth — the input XML is still the
        /// hand-authored, bible-grounded one; the writer only has to preserve what we already loaded.
        /// A writer that drops a variant or a variant-adjacent mod fails here while Corpus_Node stays green.
        /// </summary>
        [TestCaseSource(nameof(Cases))]
        public void Corpus_Node_RoundTrip(CorpusCase c)
        {
            AssertCaseSelfConsistent(c);

            List<string> produced = RunNodeRoundTrip(c);

            Assert.That(produced.Count, Is.EqualTo(c.ExpectedCount),
                $"{c.Id} ({c.Tests}) [round-trip]: expected {c.ExpectedCount} forms, got {produced.Count}.\n  got: {string.Join(" | ", produced)}");

            Assert.That(produced.OrderBy(x => x, StringComparer.Ordinal),
                Is.EqualTo(c.ExpectedForms.OrderBy(x => x, StringComparer.Ordinal)),
                $"{c.Id} ({c.Tests}) [round-trip]: form set mismatch after write/read cycle.\n  expected: {string.Join(" | ", c.ExpectedForms)}\n  got:      {string.Join(" | ", produced)}");
        }

        // The case's own two statements of size must agree — catches a typo where ExpectedForms and
        // ExpectedCount drift (mirrors the redundancy the TSV used to enforce across cases.tsv + expected/).
        private static void AssertCaseSelfConsistent(CorpusCase c) =>
            Assert.That(c.ExpectedForms.Length, Is.EqualTo(c.ExpectedCount),
                $"{c.Id}: case defines ExpectedCount={c.ExpectedCount} but lists {c.ExpectedForms.Length} ExpectedForms.");

        private static List<string> RunNodeRoundTrip(CorpusCase c)
        {
            List<Modification> knownMods = BuildMods(c);
            string xml = BuildXml(c);
            string tmp1 = Path.Combine(Path.GetTempPath(), $"corpus_{c.Id}_{Guid.NewGuid():N}.xml");
            string tmp2 = Path.Combine(Path.GetTempPath(), $"corpus_{c.Id}_rt_{Guid.NewGuid():N}.xml");
            try
            {
                // First load: the hand-authored (bible-grounded) XML -> consensus + expanded variants.
                File.WriteAllText(tmp1, xml);
                List<Protein> loaded = ProteinDbLoader.LoadProteinXML(tmp1,
                    generateTargets: true, decoyType: DecoyType.None,
                    allKnownModifications: knownMods, isContaminant: false,
                    modTypesToExclude: new List<string>(), out _,
                    maxHeterozygousVariants: c.MaxHeterozygous, minAlleleDepth: c.MinAlleleDepth);

                // Write ONLY the consensus protein(s): they retain SequenceVariations, so re-expansion
                // happens on reload exactly as it did on the first load. Writing the applied isoforms
                // instead would double-count and pin the expansion to write time.
                List<Protein> consensus = loaded.Where(p => !p.IsDecoy && !p.AppliedSequenceVariations.Any()).ToList();
                ProteinDbWriter.WriteXmlDatabase(
                    new Dictionary<string, HashSet<Tuple<int, Modification>>>(), consensus, tmp2);

                List<Protein> reread = ProteinDbLoader.LoadProteinXML(tmp2,
                    generateTargets: true, decoyType: DecoyType.None,
                    allKnownModifications: knownMods, isContaminant: false,
                    modTypesToExclude: new List<string>(), out _,
                    maxHeterozygousVariants: c.MaxHeterozygous, minAlleleDepth: c.MinAlleleDepth);

                // Writer fidelity on the variant metadata itself, independent of the digested forms. When a node's
                // oracle is consensus-only because the variant is read-but-not-applied (e.g. D08, depth-filtered),
                // the digested output equals the base sequence whether or not the writer preserved the
                // SequenceVariation, so the form-set assertion in Corpus_Node_RoundTrip cannot catch a writer that
                // dropped it. Assert the re-read consensus carries exactly the variations the loader produced.
                // (Out-of-range nodes like D05/D06 are pruned from the consensus at load, so both sides are
                // legitimately empty — there is nothing there for the writer to drop.)
                AssertConsensusVariationsRoundTrip(c, loaded, reread);

                return Digest(reread, c);
            }
            finally { File.Delete(tmp1); File.Delete(tmp2); }
        }

        // Compares the SequenceVariations carried by the consensus protein(s) before write vs after re-read, so a
        // ProteinDbWriter that silently drops or mangles a read-but-unapplied variant fails the round-trip even
        // when the digested proteoform set is identical (invariant 6 / #1083).
        private static void AssertConsensusVariationsRoundTrip(CorpusCase c, List<Protein> loaded, List<Protein> reread)
        {
            static List<(int Begin, int End, string Original, string Variant)> ConsensusVariations(List<Protein> proteins) => proteins
                .Where(p => !p.IsDecoy && !p.AppliedSequenceVariations.Any())
                .SelectMany(p => p.SequenceVariations)
                .Select(v => (v.OneBasedBeginPosition, v.OneBasedEndPosition, v.OriginalSequence, v.VariantSequence))
                .OrderBy(v => v.Item1).ThenBy(v => v.Item2).ToList();

            Assert.That(ConsensusVariations(reread), Is.EqualTo(ConsensusVariations(loaded)),
                $"{c.Id} ({c.Tests}) [round-trip]: writer did not preserve the consensus SequenceVariations across the write/read cycle.");
        }

        private static List<string> RunNode(CorpusCase c)
        {
            List<Modification> knownMods = BuildMods(c);
            string xml = BuildXml(c);
            string tmp = Path.Combine(Path.GetTempPath(), $"corpus_{c.Id}_{Guid.NewGuid():N}.xml");
            try
            {
                File.WriteAllText(tmp, xml);
                List<Protein> proteins = ProteinDbLoader.LoadProteinXML(tmp,
                    generateTargets: true, decoyType: DecoyType.None,
                    allKnownModifications: knownMods, isContaminant: false,
                    modTypesToExclude: new List<string>(), out _,
                    maxHeterozygousVariants: c.MaxHeterozygous, minAlleleDepth: c.MinAlleleDepth);
                return Digest(proteins, c);
            }
            finally { File.Delete(tmp); }
        }

        private static List<string> Digest(List<Protein> proteins, CorpusCase c)
        {
            var dp = new DigestionParams(protease: c.Protease, maxMissedCleavages: 0, minPeptideLength: 1,
                maxModificationIsoforms: c.MaxIsoforms, maxModsForPeptides: c.MaxMods,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            return proteins
                .SelectMany(p => p.Digest(dp, new List<Modification>(), new List<Modification>()))
                .Select(pw => pw.FullSequence)
                .ToList();
        }

        // --- helpers: parse the case's mods/variants DSL and emit the matching mzLib XML ------------------

        // Every mod the reader might need in allKnownModifications: the consensus mods (Mods DSL, motif = the
        // consensus residue) AND any variant-borne mods (a VMOD token on a variant, motif = the VARIANT residue,
        // since the mod sits on the post-substitution proteoform). De-duplicated by IdWithMotif so a mod shared by
        // both sides is registered once.
        private static List<Modification> BuildMods(CorpusCase c)
        {
            var byId = new Dictionary<string, Modification>();
            void Register(string name, char motifChar)
            {
                ModificationMotif.TryGetMotif(motifChar.ToString(), out ModificationMotif motif);
                var (type, mass) = ModRegistry[name];
                var mod = new Modification(_originalId: name, _modificationType: type, _target: motif,
                    _locationRestriction: "Anywhere.", _monoisotopicMass: mass);
                byId[mod.IdWithMotif] = mod;
            }

            if (c.Mods != "-")
                foreach (var token in c.Mods.Split(';', StringSplitOptions.RemoveEmptyEntries))
                {
                    var (name, pos) = ParseNameAtPos(token.Trim());
                    Register(name, c.Base[pos - 1]);
                }

            if (c.Variants != "-")
                foreach (var v in c.Variants.Split(';', StringSplitOptions.RemoveEmptyEntries))
                {
                    var kv = ParseKeyVals(v.Trim());
                    if (!kv.TryGetValue("VMOD", out var vmodSpec)) continue;
                    string variation = kv.TryGetValue("VAR", out var vv) && vv != "-" ? vv : "";
                    int begin = kv.ContainsKey("BEGIN") ? int.Parse(kv["BEGIN"]) : int.Parse(kv["POS"]);
                    foreach (var m in vmodSpec.Split(',', StringSplitOptions.RemoveEmptyEntries))
                    {
                        var (name, subpos) = ParseNameAtPos(m);
                        Register(name, variation[subpos - begin]);  // motif = variant residue at the mod site
                    }
                }

            return byId.Values.ToList();
        }

        private static string BuildXml(CorpusCase c)
        {
            var features = new System.Text.StringBuilder();

            if (c.Mods != "-")
            {
                foreach (var token in c.Mods.Split(';', StringSplitOptions.RemoveEmptyEntries))
                {
                    var (name, pos) = ParseNameAtPos(token.Trim());
                    char motifChar = c.Base[pos - 1];
                    // The reader binds a "modified residue" feature by description == Modification.IdWithMotif.
                    string idWithMotif = $"{name} on {motifChar}";
                    features.Append(
                        $"    <feature type=\"modified residue\" description=\"{idWithMotif}\">\n" +
                        $"      <location><position position=\"{pos}\" /></location>\n" +
                        $"    </feature>\n");
                }
            }

            if (c.Variants != "-")
            {
                foreach (var v in c.Variants.Split(';', StringSplitOptions.RemoveEmptyEntries))
                {
                    var kv = ParseKeyVals(v.Trim());
                    // OP = original residue(s); OP=- -> placeholder "X". VAR = variation residue(s);
                    // VAR omitted or VAR=- -> EMPTY variation = a UniProt-native deletion (installment 10).
                    string op = kv.TryGetValue("OP", out var o) && o != "-" ? o : "X";
                    string variation = kv.TryGetValue("VAR", out var vv) && vv != "-" ? vv : "";
                    // Location is either a single POS or a BEGIN/END range (multi-residue indel).
                    bool hasRange = kv.ContainsKey("BEGIN") && kv.ContainsKey("END");
                    int begin = hasRange ? int.Parse(kv["BEGIN"]) : int.Parse(kv["POS"]);
                    int end = hasRange ? int.Parse(kv["END"]) : begin;

                    // SRC=vcf -> the description IS a VCF record (so VariantCallFormat parses genotype + AD,
                    // routing through the minAlleleDepth-aware ApplyVariants path). Otherwise a prose label
                    // (UniProt-native: no genotype, goes through the combinatorial path).
                    string description;
                    if (kv.TryGetValue("SRC", out var src) && src == "vcf")
                    {
                        string gt = kv["GT"], ad = kv["AD"], dp = kv.TryGetValue("DP", out var d) ? d : "0";
                        // Effect label is cosmetic for the reader (only ANN's Allele field, == ALT, drives
                        // AlleleIndex); label it faithfully anyway: empty variation = deletion, equal-length = sub,
                        // longer = in-frame insertion (a protein-level variant cannot express a frameshift).
                        string effect = variation == "" ? "inframe_deletion"
                            : variation.Length == op.Length ? "missense_variant"
                            : variation.Length > op.Length ? "inframe_insertion"
                            : "frameshift_variant";
                        // Tabs as literal "\t" (VariantCallFormat normalizes them); ANN carries the 16
                        // pipe-delimited fields SnpEffAnnotation indexes, with Allele == ALT so AlleleIndex resolves.
                        description =
                            $"1\\t{begin}\\t.\\t{op}\\t{variation}\\t.\\tPASS\\t" +
                            $"ANN={variation}|{effect}|MODERATE|TEST|TESTID|transcript|tx1|protein_coding|1/1|c.1del|p.1del||1/7|1/7||\\t" +
                            $"GT:AD:DP\\t{gt}:{ad}:{dp}";
                    }
                    else
                    {
                        description = variation == "" ? $"deletion {op}{begin}" : $"{op}{begin}{variation} variant";
                    }

                    string locInner = hasRange
                        ? $"<begin position=\"{begin}\" /><end position=\"{end}\" />"
                        : $"<position position=\"{begin}\" />";

                    // VMOD=Name@subpos[,Name@subpos...] -> variant-BORNE mod(s): a mod that lives on the variant,
                    // not the consensus. Encoded (per the reader/writer) as a <subfeature type="modified residue">
                    // NESTED inside the variant's <location>, keyed by <subposition subposition="N"/> (distinct from
                    // the variant's own <position>). Only when the variant is applied does the mod reach the
                    // proteoform (installment 4 / #1083 reciprocal storage).
                    var subfeatures = new System.Text.StringBuilder();
                    if (kv.TryGetValue("VMOD", out var vmodSpec))
                    {
                        foreach (var m in vmodSpec.Split(',', StringSplitOptions.RemoveEmptyEntries))
                        {
                            var (mname, subpos) = ParseNameAtPos(m);
                            char motifChar = variation[subpos - begin];
                            string idWithMotif = $"{mname} on {motifChar}";
                            subfeatures.Append(
                                $"        <subfeature type=\"modified residue\" description=\"{idWithMotif}\">\n" +
                                $"          <location><subposition subposition=\"{subpos}\" /></location>\n" +
                                $"        </subfeature>\n");
                        }
                    }

                    string location = subfeatures.Length == 0
                        ? $"      <location>{locInner}</location>\n"
                        : $"      <location>\n        {locInner}\n{subfeatures}      </location>\n";

                    features.Append(
                        $"    <feature type=\"sequence variant\" description=\"{description}\">\n" +
                        $"      <original>{op}</original>\n      <variation>{variation}</variation>\n" +
                        location +
                        $"    </feature>\n");
                }
            }

            return
                "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n" +
                "<mzLibProteinDb>\n" +
                "  <entry>\n" +
                $"    <accession>P_{c.Id}</accession>\n" +
                $"    <name>TEST_{c.Id}</name>\n" +
                "    <protein><recommendedName><fullName>Test protein</fullName></recommendedName></protein>\n" +
                features +
                $"    <sequence length=\"{c.Base.Length}\">{c.Base}</sequence>\n" +
                "  </entry>\n" +
                "</mzLibProteinDb>\n";
        }

        private static (string name, int pos) ParseNameAtPos(string token)
        {
            var parts = token.Split('@');
            return (parts[0], int.Parse(parts[1]));
        }

        private static Dictionary<string, string> ParseKeyVals(string s)
        {
            var d = new Dictionary<string, string>();
            foreach (var kv in s.Split(' ', StringSplitOptions.RemoveEmptyEntries))
            {
                var i = kv.IndexOf('=');
                if (i > 0) d[kv.Substring(0, i)] = kv.Substring(i + 1);
            }
            return d;
        }
    }
}
