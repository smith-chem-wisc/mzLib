using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// Locks down what "semi-specific digestion" produces, for every way mzLib can be asked for it.
    /// </summary>
    /// <remarks>
    /// <para><b>What semi-specific means.</b> A peptide is semi-specific when AT LEAST ONE of its two ends
    /// was made by the protease. For trypsin that is either orientation:</para>
    /// <code>
    ///   ...K | ABCDER | F...    fully tryptic      both ends tryptic
    ///   ...K | ABCDE    ...     semi-tryptic       tryptic N-terminus, ragged C-terminus
    ///   ...     BCDER | F...    semi-tryptic       ragged N-terminus, tryptic C-terminus
    ///   ...     BCDE    ...     NOT semi-tryptic   both ends ragged: that is non-specific
    /// </code>
    /// <para>The ragged end is made by something other than the protease: endogenous processing (urine and
    /// plasma proteins arrive already cut), a second enzyme such as the O-glycoprotease StcE, or sample
    /// degradation. A semi-specific search is therefore the normal choice for biofluids and for
    /// glycoprotease + trypsin glycoproteomics (Lu et al., Nat. Methods 2020, the O-Pair Search paper,
    /// searched StcE + trypsin mucin digests with "semi-trypsin"). Fully specific peptides are part of the
    /// semi-specific set.</para>
    ///
    /// <para><b>Three code paths can be asked for it, and they must agree.</b></para>
    /// <list type="number">
    /// <item><description>A protease whose own specificity is <see cref="CleavageSpecificity.Semi"/>
    /// (the "classic" path, <c>Protease.SemiProteolyticDigestion</c>). mzLib shipped such a protease as
    /// <c>semi-trypsin</c> until #1005 removed it in February 2026; users can still define one in a custom
    /// proteases.tsv.</description></item>
    /// <item><description>A fully specific protease with <see cref="DigestionParams.SearchModeType"/> =
    /// <see cref="CleavageSpecificity.Semi"/> and <see cref="DigestionParams.FragmentationTerminus"/> =
    /// <see cref="FragmentationTerminus.Both"/>, the default. This must return the real semi-specific
    /// set. Before the fix it silently returned only the C-terminal "seed" peptides described next (38
    /// peptides instead of 1,726 for AMBP), and MetaMorpheus Glyco/Classic/Modern searches configured this
    /// way lost most of their identifications without any error.</description></item>
    /// <item><description>The same, but with FragmentationTerminus <see cref="FragmentationTerminus.N"/>
    /// or <see cref="FragmentationTerminus.C"/>. This deliberately does NOT enumerate semi-specific
    /// peptides. It returns "seeds" (<c>ProteinDigestion.SpeedySemiSpecificDigestion</c>): the longest
    /// fully specific stretch from each fixed terminus, which MetaMorpheus's non-specific search engine
    /// scores with ions from that terminus only and then trims to the length the precursor mass supports.
    /// These tests do not demand that seeds equal the semi set; they demand that every semi peptide is
    /// reachable by trimming a seed (coverage) and that no seed could trim to a peptide that cannot exist
    /// (tightness).</description></item>
    /// </list>
    ///
    /// <para><b>How the tests know the right answer.</b> <see cref="EnumerateReferencePeptides"/> is a
    /// deliberately simple brute-force enumeration written from the definition above and sharing NO code
    /// with mzLib's digestion. It is itself checked against mzLib's fully specific digestion (group A),
    /// which has been correct for years, and against a hand-worked example. Every other group compares a
    /// digestion path to that reference over a matrix of proteins chosen for their edge cases and every
    /// combination of missed cleavages, length limits and initiator-methionine behavior.</para>
    ///
    /// <para><b>Test proteases.</b> The tests register their own trypsin, "cleave after K or R, with no
    /// proline rule" (motifs <c>K|,R|</c>), under names unique to this fixture, in both a Full and a Semi
    /// version. That keeps the tests independent of the embedded protease table, whose trypsin naming is
    /// being changed by an open PR (#1186). The one test that uses the embedded <c>trypsin</c> checks first
    /// that its protein has no K/R followed by P, so the result is the same under either naming.</para>
    ///
    /// <para><b>Not covered by the reference:</b> truncation products (signal peptides, processed chains).
    /// Digestion adds extra peptides at their boundaries; group E pins only the invariant that those
    /// additions never remove a semi-specific peptide and always sit on a product boundary.</para>
    /// </remarks>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class SemiSpecificDigestionTests
    {
        #region Test proteins and settings

        /// <summary>
        /// Proteins chosen so that, between them, every branch of the three digestion paths is exercised.
        /// The comment on each says which edge it exists for; do not remove one without replacing its edge.
        /// </summary>
        private static readonly Dictionary<string, string> Proteins = new()
        {
            // Short enough to work through by hand; see Reference_SmallProtein_MatchesHandWorkedPeptides.
            ["HandWorkedExample"] = "PEKTIDR",

            // The protein the existing TestNonAndSemiSpecificDigests uses: initiator Met, evenly spaced K.
            ["ExistingToyProtein"] = "MAAKCCKDDKEEKFFKGG",

            // Initiator Met followed directly by a cleavage residue, and a K at the protein C-terminus.
            ["InitiatorMetThenLysine"] = "MKPEPTIDEKAAAARGGGGGK",

            // No initiator Met. The classic path used to treat residue 1 as a removable Met anyway and
            // dropped semi peptides that start at residue 2.
            ["NoInitiatorMethionine"] = "APEPTIDEKGGGG",

            // Runs of consecutive cleavage sites (every one-residue peptide K, KK, KKR... is an edge).
            ["AdjacentCleavageSites"] = "AKKRRPEPTIDEKRAAK",

            // One internal site but up to three missed cleavages allowed: the classic path's main loop
            // never runs and everything comes from its "finish the protein ends" loops.
            ["FewerSitesThanMissedCleavages"] = "MAAAKPEPTIDEGG",

            // Two internal sites: the boundary case for MaxMissedCleavages = 2 versus 3.
            ["TwoInternalSites"] = "MPEPTIDEKPEPTIDERGG",

            // No internal cleavage site at all: only protein termini are specific.
            ["NoCleavageSites"] = "MPEPTIDEPEPTIDE",

            // Residue 1 is itself a cleavage residue, so residue 2 is a specific N-terminus without any Met.
            ["StartsWithCleavageResidue"] = "KPEPTIDEAAAAGGRAA",

            // Smallest meaningful protein: Met plus one cleavage residue.
            ["TinyProtein"] = "MK",

            // A real protein (human AMBP, P02760, 380 residues) so the matrix is not only toy sequences.
            ["RealProteinAmbp"] = AmbpSequence,

            // A real mucin-domain stretch (human PSGL-1, Q14242, residues 41-140): dense Ser/Thr, few K/R.
            // This is the kind of sequence a glycoproteomics semi-tryptic search exists for.
            ["MucinFragmentPsgl1"] = Psgl1MucinFragment,
        };

        private const string AmbpSequence =
            "MRSLGALLLLLSACLAVSAGPVPTPPDNIQVQENFNISRVGTPAAHLLPGFAVMPETPLLGIWPPLKIYGKWYNLAIGSTCPWLKKIMDRMTVSTLVLGEG" +
            "ATEAEISMTSTRWRKGVCEETSGAYEKTDTDGKFLYHKSKWNITMESYVVHTNYDEYAIFLTKKFSRHHGPTITAKLYGRAPQLRETLLQDFRVVAQGVGIPE" +
            "DSIFTMADRGECVPGEQEPEPILIPRVRRAVLPQEEEGSGGGQLVTEVTKKEDSCQLGYSAGPCMGMTSRYFYNGTSMACETFQYGGCMGNGNNFVTEKECL" +
            "QTCRTVAACNLPIVRGPCRAFIQLWAFDAVKGKCVLFPYGGCQGNGNKFYSEKECREYCGVPGDGDEELLRFSN";

        private const string Psgl1MucinFragment =
            "RQATEYEYLDYDFLPETEPPEMLRNSTDTTPLTGPGTPESTTVEPAARRSTGLDAGGAVTELTTELANMGNLSTDSAAMEIQTTQPAATEAQTTQPVPTE";

        private static IEnumerable<string> ProteinNames => Proteins.Keys;

        /// <summary>One combination of the digestion settings that change which peptides exist.</summary>
        public readonly record struct DigestionSettings(int MaxMissedCleavages, int MinLength, int MaxLength, InitiatorMethionineBehavior InitiatorMethionine)
        {
            public override string ToString() =>
                $"missedCleavages={MaxMissedCleavages} minLength={MinLength} maxLength={(MaxLength == int.MaxValue ? "unlimited" : MaxLength.ToString())} initiatorMet={InitiatorMethionine}";
        }

        /// <summary>
        /// Every combination the tests sweep: 0-3 missed cleavages, a permissive and a typical minimum
        /// length, a restrictive and an unlimited maximum length, and all three initiator-Met behaviors.
        /// </summary>
        private static IEnumerable<DigestionSettings> AllSettings()
        {
            foreach (int missed in new[] { 0, 1, 2, 3 })
            foreach (int min in new[] { 1, 7 })
            foreach (int max in new[] { 25, int.MaxValue })
            foreach (InitiatorMethionineBehavior met in new[] { InitiatorMethionineBehavior.Variable, InitiatorMethionineBehavior.Retain, InitiatorMethionineBehavior.Cleave })
                yield return new DigestionSettings(missed, min, max, met);
        }

        #endregion

        #region The reference: semi-specific tryptic peptides from the definition alone

        /// <summary>A peptide as the reference defines it: 1-based inclusive residues, plus why it exists.</summary>
        public readonly record struct ReferencePeptide(int Start, int End, bool NTerminusSpecific, bool CTerminusSpecific, int MissedCleavages)
        {
            public bool FullySpecific => NTerminusSpecific && CTerminusSpecific;
            public override string ToString() => $"[{Start}-{End}]";
        }

        /// <summary>
        /// Brute-force enumeration of every peptide with at least one specific terminus, for a protease that
        /// cleaves after every K and R (no proline rule). Written from the definition, sharing no code with
        /// mzLib's digestion, so that it can judge it.
        /// </summary>
        /// <remarks>
        /// The rules, each matching mzLib's long-standing fully specific digestion (group A proves it):
        /// <list type="bullet">
        /// <item><description>A cleavage site sits after every K or R. The protein's own start and end are
        /// always specific termini.</description></item>
        /// <item><description>The N-terminus at residue <c>s</c> is specific when <c>s</c> is 1, when residue
        /// <c>s-1</c> is K/R, or when <c>s</c> is 2, residue 1 is Met, and initiator Met may be removed
        /// (Variable or Cleave).</description></item>
        /// <item><description>With initiator Met behavior Cleave and residue 1 Met, no peptide starts at
        /// residue 1: that Met does not exist in the sample.</description></item>
        /// <item><description>The C-terminus at residue <c>e</c> is specific when <c>e</c> is the last residue
        /// or residue <c>e</c> is K/R.</description></item>
        /// <item><description>Missed cleavages are the K/R residues INSIDE the peptide, i.e. at positions
        /// <c>s..e-1</c>; a K/R at the peptide's last residue is its cleavage site, not a missed one. The
        /// same count applies whether the peptide is fully or semi-specific.</description></item>
        /// <item><description>Length is <c>e - s + 1</c> and must lie within [MinLength, MaxLength].</description></item>
        /// </list>
        /// </remarks>
        public static List<ReferencePeptide> EnumerateReferencePeptides(string sequence, DigestionSettings settings)
        {
            int length = sequence.Length;
            bool IsCleavageResidue(int oneBasedPosition) => sequence[oneBasedPosition - 1] is 'K' or 'R';
            bool startsWithMet = sequence[0] == 'M';
            bool metMayBeRemoved = startsWithMet && settings.InitiatorMethionine != InitiatorMethionineBehavior.Retain;
            bool metMustBeRemoved = startsWithMet && settings.InitiatorMethionine == InitiatorMethionineBehavior.Cleave;

            var peptides = new List<ReferencePeptide>();
            for (int start = 1; start <= length; start++)
            {
                if (start == 1 && metMustBeRemoved)
                    continue;

                bool nSpecific = start == 1 || IsCleavageResidue(start - 1) || (start == 2 && metMayBeRemoved);
                int missed = 0;
                for (int end = start; end <= length; end++)
                {
                    // Extending the peptide by one residue pulls residue end-1 inside it; if that residue is
                    // K/R it becomes a missed cleavage. Both limits only get worse as the peptide grows.
                    if (end > start && IsCleavageResidue(end - 1))
                        missed++;
                    if (missed > settings.MaxMissedCleavages || end - start + 1 > settings.MaxLength)
                        break;
                    if (end - start + 1 < settings.MinLength)
                        continue;

                    bool cSpecific = end == length || IsCleavageResidue(end);
                    if (nSpecific || cSpecific)
                        peptides.Add(new ReferencePeptide(start, end, nSpecific, cSpecific, missed));
                }
            }
            return peptides;
        }

        #endregion

        #region Helpers

        private const string TestFullTrypsinName = "SemiSpecificDigestionTests-trypsin-full-no-proline-rule";
        private const string TestSemiTrypsinName = "SemiSpecificDigestionTests-trypsin-semi-no-proline-rule";

        /// <summary>
        /// Registers this fixture's own Full and Semi trypsin (cleave after K or R, no proline rule) for the
        /// lifetime of a <c>using</c> block, and removes them afterwards. <see cref="DigestionParams"/> looks
        /// proteases up by name in the static <see cref="ProteaseDictionary.Dictionary"/>, so a test
        /// protease has to be registered there; removing it in Dispose keeps other fixtures unaffected.
        /// </summary>
        private sealed class TestProteaseRegistration : IDisposable
        {
            public TestProteaseRegistration()
            {
                List<DigestionMotif> afterLysineOrArginine = DigestionMotif.ParseDigestionMotifsFromString("K|,R|");
                ProteaseDictionary.Dictionary[TestFullTrypsinName] = new Protease(TestFullTrypsinName, CleavageSpecificity.Full, null, null, afterLysineOrArginine);
                ProteaseDictionary.Dictionary[TestSemiTrypsinName] = new Protease(TestSemiTrypsinName, CleavageSpecificity.Semi, null, null, afterLysineOrArginine);
            }

            public void Dispose()
            {
                ProteaseDictionary.Dictionary.Remove(TestFullTrypsinName);
                ProteaseDictionary.Dictionary.Remove(TestSemiTrypsinName);
            }
        }

        /// <summary>What a digestion path actually returned, reduced to the facts the tests judge.</summary>
        private readonly record struct DigestedPeptide(int Start, int End, CleavageSpecificity Label, int MissedCleavages, string Description)
        {
            public override string ToString() => $"[{Start}-{End}]";
        }

        private static List<DigestedPeptide> Digest(Protein protein, DigestionParams digestionParams, List<Modification> variableModifications = null) =>
            protein.Digest(digestionParams, new List<Modification>(), variableModifications ?? new List<Modification>())
                .Select(p => new DigestedPeptide(p.OneBasedStartResidue, p.OneBasedEndResidue, p.CleavageSpecificityForFdrCategory, p.MissedCleavages, p.PeptideDescription))
                .ToList();

        private static DigestionParams Params(string protease, DigestionSettings s,
            CleavageSpecificity searchModeType = CleavageSpecificity.Full,
            FragmentationTerminus terminus = FragmentationTerminus.Both,
            bool keepOGlycopeptide = false) =>
            new(protease, s.MaxMissedCleavages, s.MinLength, s.MaxLength,
                initiatorMethionineBehavior: s.InitiatorMethionine,
                searchModeType: searchModeType, fragmentationTerminus: terminus,
                keepOGlycopeptide: keepOGlycopeptide);

        /// <summary>
        /// Compares a digestion result with the reference as sets of (start, end), reports duplicates, and
        /// returns human-readable failure lines (empty when they agree).
        /// </summary>
        private static List<string> CompareSets(string what, string proteinName, DigestionSettings settings,
            IReadOnlyCollection<(int Start, int End)> expected, IReadOnlyCollection<DigestedPeptide> actual)
        {
            var failures = new List<string>();
            var actualSet = actual.Select(p => (p.Start, p.End)).ToHashSet();
            var expectedSet = expected.ToHashSet();
            var missing = expectedSet.Except(actualSet).OrderBy(x => x).ToList();
            var extra = actualSet.Except(expectedSet).OrderBy(x => x).ToList();
            var duplicated = actual.GroupBy(p => (p.Start, p.End)).Where(g => g.Count() > 1).Select(g => g.Key).OrderBy(x => x).ToList();

            if (missing.Count > 0 || extra.Count > 0 || duplicated.Count > 0)
            {
                failures.Add($"{what} | {proteinName} | {settings}: expected {expectedSet.Count}, got {actualSet.Count} distinct"
                    + (missing.Count > 0 ? $"; MISSING {missing.Count} e.g. {Show(missing)}" : "")
                    + (extra.Count > 0 ? $"; EXTRA {extra.Count} e.g. {Show(extra)}" : "")
                    + (duplicated.Count > 0 ? $"; DUPLICATED {duplicated.Count} e.g. {Show(duplicated)}" : ""));
            }
            return failures;
        }

        private static string Show(IEnumerable<(int Start, int End)> peptides) =>
            string.Join(" ", peptides.Take(8).Select(p => $"[{p.Start}-{p.End}]"));

        private static void AssertNoFailures(List<string> failures, string explanation)
        {
            Assert.That(failures, Is.Empty,
                explanation + Environment.NewLine + $"{failures.Count} failing combination(s); first 25:" + Environment.NewLine +
                string.Join(Environment.NewLine, failures.Take(25)));
        }

        #endregion

        #region A. The reference itself is trustworthy

        /// <summary>
        /// The reference's fully specific peptides must equal mzLib's fully specific digestion for every
        /// protein and setting. Fully specific digestion has been correct for years, so agreement here shows
        /// the reference encodes mzLib's conventions for missed cleavages, length limits and initiator Met
        /// correctly, which is what lets the other groups use it as the judge for semi-specific digestion.
        /// </summary>
        [Test]
        [TestCaseSource(nameof(ProteinNames))]
        public static void Reference_FullySpecificSubset_MatchesFullTrypsinDigestion(string proteinName)
        {
            using var proteases = new TestProteaseRegistration();
            var protein = new Protein(Proteins[proteinName], proteinName);
            var failures = new List<string>();

            foreach (DigestionSettings settings in AllSettings())
            {
                var expected = EnumerateReferencePeptides(protein.BaseSequence, settings).Where(p => p.FullySpecific).Select(p => (p.Start, p.End)).ToList();
                var actual = Digest(protein, Params(TestFullTrypsinName, settings));
                failures.AddRange(CompareSets("full digestion vs reference", proteinName, settings, expected, actual));
            }

            AssertNoFailures(failures, "The reference disagrees with fully specific digestion, so it cannot be used to judge semi-specific digestion. Fix the reference first.");
        }

        /// <summary>
        /// The reference against a list worked out by hand, so it rests on something a person can check and
        /// not only on agreement with mzLib.
        /// </summary>
        /// <remarks>
        /// PEKTIDR has one internal cleavage site, after K3.
        /// <para>0 missed cleavages: peptides may not contain K3 unless it is their last residue.
        /// Tryptic N-terminus (start 1 or 4): P, PE, PEK, T, TI, TID, TIDR.
        /// Tryptic C-terminus (end 3 or 7) with a ragged start: EK, K, IDR, DR, R.
        /// Fully tryptic: PEK and TIDR. Total 12.</para>
        /// <para>1 missed cleavage adds the peptides that span K3:
        /// start 1 ending at 4-7 (PEKT, PEKTI, PEKTID, PEKTIDR) and end 7 starting at 2-3 (EKTIDR, KTIDR).
        /// Total 18.</para>
        /// </remarks>
        [Test]
        public static void Reference_SmallProtein_MatchesHandWorkedPeptides()
        {
            const string sequence = "PEKTIDR";
            string[] Sequences(int missed) => EnumerateReferencePeptides(sequence, new DigestionSettings(missed, 1, int.MaxValue, InitiatorMethionineBehavior.Variable))
                .Select(p => sequence.Substring(p.Start - 1, p.End - p.Start + 1)).ToArray();

            Assert.That(Sequences(0), Is.EquivalentTo(new[] { "P", "PE", "PEK", "T", "TI", "TID", "TIDR", "EK", "K", "IDR", "DR", "R" }));
            Assert.That(Sequences(1), Is.EquivalentTo(new[] { "P", "PE", "PEK", "T", "TI", "TID", "TIDR", "EK", "K", "IDR", "DR", "R",
                                                              "PEKT", "PEKTI", "PEKTID", "PEKTIDR", "EKTIDR", "KTIDR" }));
        }

        /// <summary>
        /// The reference hard-codes "cleave after K or R". This pins that the test protease really has that
        /// rule, so a change to motif parsing cannot make the reference and the protease quietly disagree.
        /// </summary>
        [Test]
        [TestCaseSource(nameof(ProteinNames))]
        public static void Reference_CleavageRule_MatchesTestProteaseCleavageSites(string proteinName)
        {
            using var proteases = new TestProteaseRegistration();
            string sequence = Proteins[proteinName];
            var expectedSites = new List<int> { 0 };
            for (int position = 1; position < sequence.Length; position++)
                if (sequence[position - 1] is 'K' or 'R')
                    expectedSites.Add(position);
            expectedSites.Add(sequence.Length);

            Assert.That(ProteaseDictionary.Dictionary[TestFullTrypsinName].GetDigestionSiteIndices(sequence), Is.EqualTo(expectedSites));
        }

        #endregion

        #region B. The classic path: a protease whose own specificity is Semi

        /// <summary>
        /// A Semi-specificity protease must return exactly the reference's semi-specific set, once each.
        /// </summary>
        /// <remarks>
        /// Two defects in <c>Protease.SemiProteolyticDigestion</c> used to break this, both in the loops that
        /// finish the protein ends:
        /// (1) residue 1 was treated as a removable initiator Met even when it is not M, so for proteins that
        /// do not start with M, semi peptides starting at residue 2 and ending at one of the first cleavage
        /// sites were dropped (NoInitiatorMethionine, AdjacentCleavageSites);
        /// (2) when a protein has fewer internal cleavage sites than MaxMissedCleavages the main loop never
        /// runs, so every peptide starting after a removed Met was lost, and with initiator Met Cleave the
        /// Met-retaining peptides were emitted anyway (FewerSitesThanMissedCleavages, NoCleavageSites).
        /// </remarks>
        [Test]
        [TestCaseSource(nameof(ProteinNames))]
        public static void ClassicSemiProtease_ReturnsExactlyTheReferencePeptides(string proteinName)
        {
            using var proteases = new TestProteaseRegistration();
            var protein = new Protein(Proteins[proteinName], proteinName);
            var failures = new List<string>();

            foreach (DigestionSettings settings in AllSettings())
            {
                var expected = EnumerateReferencePeptides(protein.BaseSequence, settings).Select(p => (p.Start, p.End)).ToList();
                var actual = Digest(protein, Params(TestSemiTrypsinName, settings));
                failures.AddRange(CompareSets("classic semi protease vs reference", proteinName, settings, expected, actual));
            }

            AssertNoFailures(failures, "A Semi-specificity protease must produce every peptide with at least one specific terminus, and nothing else.");
        }

        /// <summary>
        /// Each semi-specific peptide must be labelled Full when both ends are specific and Semi otherwise.
        /// MetaMorpheus uses this label to split PSMs into separate FDR categories, so a wrong label moves a
        /// peptide into the wrong target-decoy competition.
        /// </summary>
        [Test]
        [TestCaseSource(nameof(ProteinNames))]
        public static void ClassicSemiProtease_LabelsFullAndSemiPeptidesCorrectly(string proteinName)
        {
            using var proteases = new TestProteaseRegistration();
            var protein = new Protein(Proteins[proteinName], proteinName);
            AssertNoFailures(CompareLabelsAndMissedCleavages(protein, proteinName, settings => Params(TestSemiTrypsinName, settings), checkMissedCleavages: false),
                "Peptide labels (CleavageSpecificityForFdrCategory) disagree with which termini are specific.");
        }

        /// <summary>
        /// Each semi-specific peptide must report its true number of missed cleavages (internal K/R).
        /// </summary>
        /// <remarks>
        /// The classic path used to pass the peptide LENGTH where the missed-cleavage count belongs, for every
        /// peptide it made. MetaMorpheus feeds <c>MissedCleavages</c> to its PEP model as the
        /// MissedCleavagesCount feature and writes it to the glyco PSM output, so every semi-trypsin search
        /// trained PEP on a meaningless number.
        /// </remarks>
        [Test]
        [TestCaseSource(nameof(ProteinNames))]
        public static void ClassicSemiProtease_ReportsTrueMissedCleavages(string proteinName)
        {
            using var proteases = new TestProteaseRegistration();
            var protein = new Protein(Proteins[proteinName], proteinName);
            AssertNoFailures(CompareLabelsAndMissedCleavages(protein, proteinName, settings => Params(TestSemiTrypsinName, settings), checkLabels: false),
                "MissedCleavages must count the cleavage residues inside the peptide.");
        }

        /// <summary>
        /// Shared check for labels and missed cleavages against the reference, matched by (start, end).
        /// Peptides missing from the result are reported by the set tests, not here.
        /// </summary>
        private static List<string> CompareLabelsAndMissedCleavages(Protein protein, string proteinName, Func<DigestionSettings, DigestionParams> makeParams,
            bool checkLabels = true, bool checkMissedCleavages = true)
        {
            var failures = new List<string>();
            foreach (DigestionSettings settings in AllSettings())
            {
                var reference = EnumerateReferencePeptides(protein.BaseSequence, settings).ToDictionary(p => (p.Start, p.End));
                var wrong = new List<string>();
                foreach (DigestedPeptide peptide in Digest(protein, makeParams(settings)))
                {
                    if (!reference.TryGetValue((peptide.Start, peptide.End), out ReferencePeptide expected))
                        continue;
                    CleavageSpecificity expectedLabel = expected.FullySpecific ? CleavageSpecificity.Full : CleavageSpecificity.Semi;
                    if (checkLabels && peptide.Label != expectedLabel)
                        wrong.Add($"{peptide} labelled {peptide.Label}, expected {expectedLabel}");
                    if (checkMissedCleavages && peptide.MissedCleavages != expected.MissedCleavages)
                        wrong.Add($"{peptide} reports {peptide.MissedCleavages} missed cleavages, expected {expected.MissedCleavages}");
                }
                if (wrong.Count > 0)
                    failures.Add($"{proteinName} | {settings}: {wrong.Count} wrong, e.g. {string.Join("; ", wrong.Take(5))}");
            }
            return failures;
        }

        #endregion

        #region E. Truncation products only ever add peptides at their own boundaries

        /// <summary>
        /// A protein with a truncation product (for example a signal peptide removed at residue 20) gains
        /// peptides that start or end at the product boundary. Whatever else digestion does with them, it must
        /// never lose a semi-specific peptide the protein would have without the product, and every extra
        /// peptide must start at a product's first residue or end at its last. Checked for both ways of asking
        /// for a semi-specific digest.
        /// </summary>
        [Test]
        [TestCase(TestSemiTrypsinName, CleavageSpecificity.Full)]
        public static void SemiDigestion_WithTruncationProduct_OnlyAddsPeptidesOnTheProductBoundary(string protease, CleavageSpecificity searchModeType)
        {
            using var proteases = new TestProteaseRegistration();
            var products = new List<TruncationProduct> { new TruncationProduct(20, 150, "chain") };
            var protein = new Protein(AmbpSequence, "P02760", proteolysisProducts: products);
            var failures = new List<string>();

            foreach (DigestionSettings settings in AllSettings())
            {
                var reference = EnumerateReferencePeptides(AmbpSequence, settings).Select(p => (p.Start, p.End)).ToHashSet();
                List<DigestedPeptide> digested = Digest(protein, Params(protease, settings, searchModeType, FragmentationTerminus.Both));
                var actual = digested.Select(p => (p.Start, p.End)).ToHashSet();

                var lost = reference.Except(actual).OrderBy(x => x).ToList();
                var offBoundary = actual.Except(reference).Where(p => !products.Any(t => p.Start == t.OneBasedBeginPosition || p.End == t.OneBasedEndPosition)).OrderBy(x => x).ToList();

                // Boundary peptides are still trypsin peptides inside: their missed cleavages are the K/R
                // residues inside them, exactly as for every other peptide (not their length, and not
                // length - 1, which is the count only a non-specific protease would give).
                var wrongMissedCleavages = digested
                    .Where(p => p.MissedCleavages != AmbpSequence.Substring(p.Start - 1, p.End - p.Start).Count(c => c is 'K' or 'R'))
                    .Select(p => $"{p} reports {p.MissedCleavages}").Take(5).ToList();

                if (lost.Count > 0 || offBoundary.Count > 0 || wrongMissedCleavages.Count > 0)
                    failures.Add($"{settings}: lost {lost.Count} e.g. {Show(lost)}; extra peptides off the product boundary {offBoundary.Count} e.g. {Show(offBoundary)}; wrong missed cleavages e.g. {string.Join("; ", wrongMissedCleavages)}");
            }

            AssertNoFailures(failures, "Truncation products may only add peptides that start or end on the product boundary, and every peptide must report its true missed cleavages.");
        }

        #endregion
    }
}
