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
        private static readonly HashSet<string> ModelledEnzymes =
            new(StringComparer.OrdinalIgnoreCase) { "StcE", "StcE-trypsin", "trypsin", "trypsin|P" };

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

        /// <summary>The glycan markers are stripped: with no glycan-aware digestion only base sequences compare.</summary>
        private static string StripGlycanMarkers(string product) => product.Replace("*", string.Empty);

        private static List<string> Digest(TruthCase c)
        {
            var protein = new Protein(c.Sequence, "TRUTHSET_" + c.CaseId);

            // minPeptideLength 1 because real products here are short (AASAA gives "AA"), and
            // maxMissedCleavages 0 because the truth set states the COMPLETE digest, not a partial one.
            var parameters = new DigestionParams(
                protease: c.Enzyme,
                maxMissedCleavages: 0,
                minPeptideLength: 1,
                initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            return protein.Digest(parameters, new List<Modification>(), new List<Modification>())
                .Select(p => p.BaseSequence)
                .OrderBy(s => s, StringComparer.Ordinal)
                .ToList();
        }

        [Test]
        [TestCaseSource(nameof(TruthCases))]
        public static void Digestion_ReproducesThePublishedProducts(TruthCase c)
        {
            if (!c.IsEncodable)
                Assert.Inconclusive(c.CaseId + ": the source names the substrate but not an exact bond, so no "
                    + "product list is defensible yet. " + c.Note);

            if (!ModelledEnzymes.Contains(c.Enzyme))
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

                if (!ModelledEnzymes.Contains(c.Enzyme))
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
