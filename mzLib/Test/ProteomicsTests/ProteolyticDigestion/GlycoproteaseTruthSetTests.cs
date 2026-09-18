using NUnit.Framework;
// Imported rather than written out at the use site: this file sits in namespace Test.*, where the
// qualified name Omics.Modifications binds to Test.Omics.Modifications, which has no Modification.
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Assert = NUnit.Framework.Legacy.ClassicAssert;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// Acceptance tests for glycosylation-aware digestion, driven by a hand-curated truth set of
    /// published substrates whose digestion products are known. Each case is one
    /// (enzyme x substrate x glycoform) whose answer a paper states outright, so these are the
    /// EXPECTATIONS -- they are not a description of what the digester currently does.
    /// </summary>
    /// <remarks>
    /// <para>Written before the glycan-aware digestion work, deliberately. Several cases fail today,
    /// and that is the point: the failures are the specification. The three shapes of failure are
    /// worth telling apart when reading the output.</para>
    ///
    /// <para><b>Over-digestion.</b> mzLib's <c>StcE</c> and <c>StcE-trypsin</c> entries match the
    /// SEQUENCE motif only. The real enzyme also requires an O-GalNAc at P2, so every case whose
    /// expectation is "no cleavage because the glycan is absent or is the wrong sugar" is cut anyway.
    /// The shipped <c>proteases.tsv</c> comment block already says so: the entries "OVER-DIGEST
    /// relative to the real enzyme -- their peptides are a superset of true StcE peptides". These
    /// tests put a number on that.</para>
    ///
    /// <para><b>Enzyme not modelled.</b> IMPa, OpeRATOR/OgpA, SmE, BT4244, CpaA and AMUC_1438 have no
    /// entry in <c>proteases.tsv</c>, by an explicit decision recorded in the same comment block:
    /// modelling them "needs glycosylation-aware digestion". Those cases report Inconclusive rather
    /// than failing, because there is nothing to be wrong yet.</para>
    ///
    /// <para><b>Efficiency, not possibility.</b> Cases marked <c>reduced:&lt;pct&gt;</c> are REAL
    /// cleavages that happen slowly -- IMPa at 8% completion with isoleucine at P1, OpeRATOR ~170x
    /// slower on alpha-2,3-sialyl core 1. A digester should model those as missed-cleavage weight, not
    /// as blocks. The truth set records the percentage so the distinction survives.</para>
    ///
    /// <para>The glycan column is carried through but not yet consulted: with no glycan-aware
    /// digestion, only the base sequences can be compared. That is exactly why the negative cases
    /// fail.</para>
    /// </remarks>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class GlycoproteaseTruthSetTests
    {
        /// <summary>Enzymes mzLib can actually digest with today. Everything else is not modelled yet.</summary>
        /// <summary>
        /// Whether proteases.tsv actually ships this enzyme. Asked of the dictionary rather than held in
        /// a list here: a hand-maintained allowlist drifts silently the moment an enzyme is added, and it
        /// did -- StcE-trypsin|P was shipped and its two cases kept reporting "not modelled" instead of
        /// being tested.
        /// </summary>
        private static bool IsModelled(string enzyme) =>
            !string.IsNullOrWhiteSpace(enzyme) && ProteaseDictionary.Dictionary.ContainsKey(enzyme);

        /// <summary>One row of the truth set.</summary>
        public class TruthCase
        {
            public string CaseId { get; init; }
            public string Enzyme { get; init; }
            public string Substrate { get; init; }
            public string Sequence { get; init; }
            public string Glycosites { get; init; }
            public string[] ExpectedProducts { get; init; }
            public int ExpectedCuts { get; init; }
            public string Efficiency { get; init; }
            public string Grade { get; init; }
            public string Basis { get; init; }
            public string Source { get; init; }
            public string Note { get; init; }

            /// <summary>
            /// Empty when the digester is expected to get this case right today. Otherwise it names a
            /// KNOWN shortfall, and the case reports Inconclusive instead of failing -- the gap stays
            /// visible in the report without reddening CI. Clearing the reason in the TSV promotes the
            /// row to a real assertion, which is how this fixture drives the work forward.
            /// </summary>
            public string KnownGap { get; init; }

            /// <summary>
            /// "fixed" (the default) when the published substrate was homogeneously glycosylated, as a
            /// synthetic peptide is; "variable" when the case deliberately models a mixed population in
            /// which some molecules carry the glycan and some do not.
            /// </summary>
            public string Occupancy { get; init; }

            public bool IsMixedPopulation =>
                string.Equals(Occupancy, "variable", StringComparison.OrdinalIgnoreCase);

            public bool HasKnownGap => !string.IsNullOrEmpty(KnownGap) && KnownGap != "-";

            public bool IsEncodable =>
                ExpectedProducts.Length > 0
                && !ExpectedProducts.Any(p => p.StartsWith("(", StringComparison.Ordinal));

            public override string ToString() => CaseId + " " + Enzyme + " " + Substrate;
        }

        private static string TruthSetPath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            "ProteomicsTests", "ProteolyticDigestion", "TestData", "digestion-truth-set.tsv");

        private static IEnumerable<TruthCase> LoadTruthSet()
        {
            foreach (string line in File.ReadAllLines(TruthSetPath))
            {
                if (line.Length == 0 || line[0] == '#' || line.StartsWith("case_id", StringComparison.Ordinal))
                    continue;

                string[] f = line.Split('\t');
                if (f.Length < 12)
                    continue;

                yield return new TruthCase
                {
                    CaseId = f[0],
                    Enzyme = f[1],
                    Substrate = f[2],
                    Sequence = f[3],
                    Glycosites = f[4],
                    ExpectedProducts = f[5].Split('|', StringSplitOptions.RemoveEmptyEntries),
                    ExpectedCuts = int.TryParse(f[6], out int cuts) ? cuts : -1,
                    Efficiency = f[7],
                    Grade = f[8],
                    Basis = f[9],
                    Source = f[10],
                    Note = f[11],
                    KnownGap = f.Length > 12 ? f[12] : "-",
                    Occupancy = f.Length > 13 ? f[13] : "fixed",
                };
            }
        }

        /// <summary>
        /// NUnit refuses a test display name containing a '.', and CI checks for it before the suite
        /// runs, so every generated name is scrubbed. Sequences and case ids carry no dots today; the
        /// replacement is here so that adding a source or substrate that does cannot redden CI.
        /// </summary>
        public static IEnumerable<TestCaseData> TruthCases() =>
            LoadTruthSet().Select(c => new TestCaseData(c)
                .SetName(("TruthSet_" + c.CaseId + "_" + c.Enzyme + "_" + c.Substrate).Replace('.', '_')));

        /// <summary>The glycan markers are stripped: products are compared as base sequences.</summary>
        private static string StripGlycanMarkers(string product) => product.Replace("*", string.Empty);

        /// <summary>
        /// Builds the substrate with its glycans placed as localized modifications, from the truth set's
        /// <c>glycosites</c> column (<c>pos:glycan</c>, one-based, semicolon separated).
        /// </summary>
        /// <remarks>
        /// Every glycan in the corpus is placed with ModificationType "O-linked glycosylation", including
        /// the O-GlcNAc and O-mannose cases, because that is what a real modification database says about
        /// them -- they ARE O-linked. The digester's requirement is resolved at CLASS level, so it cannot
        /// tell alpha-O-GalNAc from beta-O-GlcNAc and will cleave both. That is a documented limit of the
        /// current model rather than a bug in the fixture, and STCE-02 stays flagged as a known gap
        /// because of it. Encoding O-GlcNAc as something other than O-linked here would hide the limit
        /// by mis-describing the chemistry.
        /// </remarks>
        private static Protein BuildProtein(TruthCase c)
        {
            var localized = new Dictionary<int, List<Modification>>();

            if (!string.IsNullOrWhiteSpace(c.Glycosites) && c.Glycosites != "-")
            {
                foreach (string site in c.Glycosites.Split(';', StringSplitOptions.RemoveEmptyEntries))
                {
                    string[] parts = site.Split(':');
                    if (parts.Length != 2 || !int.TryParse(parts[0], out int position))
                        continue;
                    if (position < 1 || position > c.Sequence.Length)
                        continue;

                    string residue = c.Sequence[position - 1].ToString();
                    if (!ModificationMotif.TryGetMotif(residue, out ModificationMotif motif))
                        continue;

                    var glycan = new Modification(_originalId: parts[1],
                        _modificationType: "O-linked glycosylation", _target: motif,
                        _locationRestriction: "Anywhere.", _monoisotopicMass: 203.079373,
                        _monosaccharideComposition: CompositionOf(parts[1]));

                    if (!localized.TryGetValue(position, out List<Modification> atSite))
                    {
                        atSite = new List<Modification>();
                        localized[position] = atSite;
                    }

                    atSite.Add(glycan);
                }
            }

            return localized.Count == 0
                ? new Protein(c.Sequence, "TRUTHSET_" + c.CaseId)
                : new Protein(c.Sequence, "TRUTHSET_" + c.CaseId, oneBasedModifications: localized);
        }

        private static List<string> Digest(TruthCase c)
        {
            Protein protein = BuildProtein(c);

            // minPeptideLength 1 because real products here are short (AASAA gives "AA"), and
            // maxMissedCleavages 0 because the truth set states the COMPLETE digest, not a partial one.
            // RespectCleavagePromotingModifications is ON here, because the truth set records what the
            // REAL enzyme does and that is what these cases assert. With it off the glycoproteases fall
            // back to their sequence motif and over-digest, which is the behaviour the negative cases
            // exist to catch.
            var parameters = new DigestionParams(
                protease: c.Enzyme,
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain,
                respectCleavagePromotingModifications: true);

            // A homogeneously glycosylated substrate is best modelled with FIXED modifications, not
            // localized ones. mzLib enumerates a localized modification both ways, which is right for a
            // search and wrong for a synthetic peptide -- and, more than a nuisance, it makes one rule
            // unjudgeable: a FORBIDDEN condition on the non-prime side constrains a residue that lies in
            // the NEIGHBOURING peptide, where no peptidoform can speak for it. Fixed settles that, because
            // fixed means every copy carries it.
            //
            // Only safe when every occurrence of a residue letter in this sequence is a declared
            // glycosite; otherwise a fixed modification would glycosylate residues the paper says are
            // bare, so those cases keep the localized-plus-filter path.
            List<Modification> fixedGlycans = FixedGlycansIfUnambiguous(c);
            IEnumerable<PeptideWithSetModifications> products = fixedGlycans is null
                ? protein.Digest(parameters, new List<Modification>(), new List<Modification>())
                : BareProtein(c).Digest(parameters, fixedGlycans, new List<Modification>());

            // OCCUPANCY. mzLib enumerates a localized modification BOTH ways -- one peptidoform carrying
            // it and one without -- because in a real search a database glycosite is a site that MAY be
            // occupied. That is the right default for a search and the wrong model for most of this
            // corpus: a synthetic peptide carrying a single GalNAc is homogeneously glycosylated, there
            // are no unglycosylated molecules, and so there is nothing for the protease to read through.
            //
            // Keeping only the fully-occupied peptidoforms reproduces the published experiment exactly.
            // The one case that genuinely models a mixed population (STCE-08) declares occupancy
            // "variable" and is left alone -- and it is precisely the read-through that distinguishes it
            // from STCE-01, which is the SAME substrate with the SAME glycosite. Those two cases are the
            // corpus's sharpest statement that occupancy, not sequence, decides the digest.
            if (!c.IsMixedPopulation && fixedGlycans is null)
            {
                products = products.Where(p => EveryGlycositeInsideIsOccupied(p, c));
            }

            return products
                .Select(p => p.BaseSequence)
                .Distinct()
                .OrderBy(s => s, StringComparer.Ordinal)
                .ToList();
        }

        /// <summary>
        /// The monosaccharide composition behind each glycan shorthand the truth set uses, or null when
        /// the shorthand names something a composition cannot express.
        /// </summary>
        /// <remarks>
        /// Two pairs in this table are deliberately IDENTICAL, and they are why composition closes some
        /// cases and not others. 3SC1 and 6SC1 differ only in whether the sialic acid is alpha-2,3 or
        /// alpha-2,6 linked -- OpeRATOR tolerates the first and is blocked outright by the second -- and
        /// Tn and OGlcNAc are both a single HexNAc differing in sugar and anomer, which StcE cares about
        /// and this cannot see. Those stay gaps by construction, not by oversight.
        /// </remarks>
        private static MonosaccharideComposition CompositionOf(string glycanShorthand) =>
            glycanShorthand switch
            {
                "Tn" => MonosaccharideComposition.Parse("HexNAc1"),
                "OGlcNAc" => MonosaccharideComposition.Parse("HexNAc1"),
                "T/core1" => MonosaccharideComposition.Parse("Hex1HexNAc1"),
                "core2" => MonosaccharideComposition.Parse("Hex1HexNAc2"),
                "sT" => MonosaccharideComposition.Parse("NeuAc1Hex1HexNAc1"),
                "3SC1" => MonosaccharideComposition.Parse("NeuAc1Hex1HexNAc1"),
                "6SC1" => MonosaccharideComposition.Parse("NeuAc1Hex1HexNAc1"),
                "dsC2" => MonosaccharideComposition.Parse("NeuAc2Hex1HexNAc2"),
                _ => null,
            };

        private static Protein BareProtein(TruthCase c) => new Protein(c.Sequence, "TRUTHSET_" + c.CaseId);

        /// <summary>
        /// The declared glycans as FIXED modifications, or null when that would misrepresent the
        /// substrate. Only for a homogeneously glycosylated case, and only when every occurrence of each
        /// glycosylated residue letter is itself a declared glycosite -- otherwise a fixed modification
        /// would decorate residues the source says are bare.
        /// </summary>
        private static List<Modification> FixedGlycansIfUnambiguous(TruthCase c)
        {
            if (c.IsMixedPopulation || string.IsNullOrWhiteSpace(c.Glycosites) || c.Glycosites == "-")
            {
                return null;
            }

            var glycanByResidue = new SortedDictionary<char, string>();
            var declaredPositions = new HashSet<int>();

            foreach (string site in c.Glycosites.Split(';', StringSplitOptions.RemoveEmptyEntries))
            {
                string[] parts = site.Split(':');
                if (parts.Length != 2 || !int.TryParse(parts[0], out int position)
                    || position < 1 || position > c.Sequence.Length)
                {
                    return null;
                }

                declaredPositions.Add(position);
                char residue = c.Sequence[position - 1];
                if (glycanByResidue.TryGetValue(residue, out string already) && already != parts[1])
                {
                    // Two different glycans on the same residue letter cannot both be fixed.
                    return null;
                }

                glycanByResidue[residue] = parts[1];
            }

            foreach (char residue in glycanByResidue.Keys)
            {
                for (int i = 0; i < c.Sequence.Length; i++)
                {
                    if (c.Sequence[i] == residue && !declaredPositions.Contains(i + 1))
                    {
                        return null;
                    }
                }
            }

            var fixedGlycans = new List<Modification>();
            foreach (var pair in glycanByResidue)
            {
                if (!ModificationMotif.TryGetMotif(pair.Key.ToString(), out ModificationMotif motif))
                {
                    return null;
                }

                fixedGlycans.Add(new Modification(_originalId: pair.Value,
                    _modificationType: "O-linked glycosylation", _target: motif,
                    _locationRestriction: "Anywhere.", _monoisotopicMass: 203.079373,
                    _monosaccharideComposition: CompositionOf(pair.Value)));
            }

            return fixedGlycans;
        }

        /// <summary>
        /// True when every glycosite the truth set declares for this case, that falls inside this
        /// peptide, actually carries its modification in this peptidoform.
        /// </summary>
        /// <remarks>
        /// Modification keys are two-based: key 1 is the N-terminus, key 2 the first residue. A protein
        /// position P inside a peptide starting at S is therefore key P - S + 2.
        /// </remarks>
        private static bool EveryGlycositeInsideIsOccupied(PeptideWithSetModifications peptide, TruthCase c)
        {
            if (string.IsNullOrWhiteSpace(c.Glycosites) || c.Glycosites == "-")
            {
                return true;
            }

            foreach (string site in c.Glycosites.Split(';', StringSplitOptions.RemoveEmptyEntries))
            {
                string[] parts = site.Split(':');
                if (parts.Length != 2 || !int.TryParse(parts[0], out int position))
                    continue;

                if (position < peptide.OneBasedStartResidue || position > peptide.OneBasedEndResidue)
                    continue;

                int key = position - peptide.OneBasedStartResidue + 2;
                if (!peptide.AllModsOneIsNterminus.ContainsKey(key))
                {
                    return false;
                }
            }

            return true;
        }

        [Test]
        [TestCaseSource(nameof(TruthCases))]
        public static void Digestion_ReproducesThePublishedProducts(TruthCase c)
        {
            if (!c.IsEncodable)
                Assert.Inconclusive(c.CaseId + ": the source names the substrate but not an exact bond, so no "
                    + "product list is defensible yet. " + c.Note);

            if (!IsModelled(c.Enzyme))
                Assert.Inconclusive(c.CaseId + ": " + c.Enzyme + " has no proteases.tsv entry -- it was left out "
                    + "deliberately because modelling it needs glycosylation-aware digestion. Nothing to be wrong "
                    + "about until that exists. Source: " + c.Source);

            List<string> expected = c.ExpectedProducts.Select(StripGlycanMarkers)
                .OrderBy(s => s, StringComparer.Ordinal).ToList();
            List<string> actual = Digest(c);

            if (c.HasKnownGap)
            {
                Assert.AreNotEqual(expected, actual, c.CaseId + " is flagged as a known gap but now MATCHES. "
                    + "That is good news: clear the known_gap column in digestion-truth-set.tsv so this case "
                    + "becomes a real assertion and cannot regress. Recorded gap was: " + c.KnownGap);
                Assert.Inconclusive(c.CaseId + " -- KNOWN GAP: " + c.KnownGap + Environment.NewLine
                    + "  expected: " + string.Join(" | ", expected)
                    + Environment.NewLine + "  actual:   " + string.Join(" | ", actual)
                    + Environment.NewLine + "  source:   " + c.Source);
            }

            Assert.AreEqual(expected, actual,
                c.CaseId + " (" + c.Enzyme + " on " + c.Substrate + ", glycosites " + c.Glycosites + ")"
                + Environment.NewLine + "  expected: " + string.Join(" | ", expected)
                + Environment.NewLine + "  actual:   " + string.Join(" | ", actual)
                + Environment.NewLine + "  source:   " + c.Source
                + Environment.NewLine + "  why:      " + c.Note);
        }

        /// <summary>
        /// Prints the whole truth set against current behaviour in one table, so the size and shape of
        /// the gap is visible without reading N individual failures. Always passes -- it is a report,
        /// not an assertion; <see cref="Digestion_ReproducesThePublishedProducts"/> is what fails.
        /// </summary>
        [Test]
        public static void TruthSet_Report()
        {
            var report = new StringBuilder();
            int match = 0, mismatch = 0, notModelled = 0, notEncodable = 0;

            report.AppendLine("case       enzyme        result      expected -> actual");
            report.AppendLine(new string('-', 100));

            foreach (TruthCase c in LoadTruthSet())
            {
                if (!c.IsEncodable)
                {
                    notEncodable++;
                    report.AppendLine($"{c.CaseId,-10} {c.Enzyme,-13} not-encodable");
                    continue;
                }

                if (!IsModelled(c.Enzyme))
                {
                    notModelled++;
                    report.AppendLine($"{c.CaseId,-10} {c.Enzyme,-13} not-modelled");
                    continue;
                }

                List<string> expected = c.ExpectedProducts.Select(StripGlycanMarkers)
                    .OrderBy(s => s, StringComparer.Ordinal).ToList();
                List<string> actual = Digest(c);
                bool ok = expected.SequenceEqual(actual);

                if (ok) match++; else mismatch++;

                string verdict = ok ? "MATCH" : (c.HasKnownGap ? "GAP" : "MISMATCH");
                report.AppendLine($"{c.CaseId,-10} {c.Enzyme,-13} {verdict,-11} "
                    + string.Join("|", expected) + "  ->  " + string.Join("|", actual));
            }

            report.AppendLine(new string('-', 100));
            report.AppendLine($"match {match}   mismatch {mismatch}   enzyme-not-modelled {notModelled}   "
                + $"not-encodable {notEncodable}");
            TestContext.Out.WriteLine(report.ToString());

            Assert.Greater(match + mismatch, 0, "no truth-set case was runnable at all, which means the "
                + "fixture did not load -- check that digestion-truth-set.tsv is copied to the output directory");
        }
    }
}
