using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class PeptideWithSetModsTruncationTests
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        /// <summary>
        /// CRITICAL: Tests top-down proteomics digestion with proteolysis products.
        /// Validates that full-length proteins and their truncation products are correctly
        /// generated. Essential for intact protein and proteoform identification.
        /// </summary>
        [Test]
        public static void TestTopDownDigestion()
        {
            List<TruncationProduct> proteolysisProducts = new List<TruncationProduct>
            {
                new TruncationProduct(5, 20, "asdf")
            };
            Protein protein = new Protein("MACDEFGHIKLMNOPQRSTVWYMACDEFGHIKLMNOPQRSTVWYMACDEFGHIKLMNOPQRSTVWY", "testProtein", "Mus", proteolysisProducts: proteolysisProducts);
            DigestionParams topdownParams = new DigestionParams("top-down");
            List<PeptideWithSetModifications> peptides = protein.Digest(topdownParams, null, null).ToList();
            Assert.IsTrue(peptides.Count == 3);
        }

        /// <summary>
        /// CRITICAL: Tests that top-down digestion returns expected truncation products.
        /// For insulin, expects 68 truncation products. Essential for top-down and
        /// middle-down proteomics workflows that search for proteoforms.
        /// </summary>
        [Test]
        public static void TestPeptideWithSetModsReturnsTruncationsInTopDown()
        {
            string xmlDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.xml");

            Protein insulin = ProteinDbLoader.LoadProteinXML(xmlDatabase, true,
                DecoyType.None, null, false, null, out var unknownModifications, addTruncations: true)[0];

            Protease protease = new Protease("top-down", CleavageSpecificity.None, "", "", new List<DigestionMotif>(), null);
            List<PeptideWithSetModifications> insulinTruncations = insulin.Digest(new DigestionParams(protease: protease.Name), new List<Modification>(), new List<Modification>(), topDownTruncationSearch: true).ToList();
            Assert.AreEqual(68, insulinTruncations.Count);
        }
        /// <summary>
        /// CRITICAL: Regression test for top-down truncation product generation.
        /// Generates a comprehensive table of all truncation peptides for target/decoy
        /// insulin and compares against expected baseline. Ensures truncation logic
        /// remains stable across code changes. Essential for top-down proteomics.
        /// </summary>
        [Test]
        public static void TestTopDownTruncationTableMatchesExpected()
        {
            // PURPOSE
            // Generate a comprehensive TSV of all top-down truncation peptides for target and decoy insulin entries,
            // then compare the result against the stored expected file for regression checking.
            //
            // OUTPUT COLUMNS
            // - Sequence: peptide/proteoform base sequence
            // - Type: "Target" or "Decoy" based on parent protein
            // - Begin: one-based start residue within the parent protein
            // - End: one-based end residue within the parent protein
            // - RetainedMethionine: TRUE when the peptide includes the protein's N-terminal Met (Begin == 1 and Parent.BaseSequence[0] == 'M'), else FALSE

            // Arrange: load insulin with reverse decoys and truncations enabled
            string xmlDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "humanInsulin.xml");
            List<Protein> insulinProteins = ProteinDbLoader.LoadProteinXML(
                xmlDatabase,
                generateTargets: true,
                decoyType: DecoyType.Reverse,
                allKnownModifications: null,
                isContaminant: false,
                modTypesToExclude: null,
                unknownModifications: out var unknownModifications,
                addTruncations: true);

            Assert.That(insulinProteins.Any(p => !p.IsDecoy), "Expected at least one target protein");
            Assert.That(insulinProteins.Any(p => p.IsDecoy), "Expected at least one decoy protein");
            Assert.That(unknownModifications == null || unknownModifications.Count == 0, "No unknown modifications expected from insulin XML");

            // Digest: enumerate truncations for a representative target/decoy pair (parity sanity)
            static string Reverse(string s) => new string(s.Reverse().ToArray());
            var target = insulinProteins.First(p => !p.IsDecoy);
            string expectedDecoySeq = target.BaseSequence.Length > 0 && target.BaseSequence[0] == 'M'
                ? "M" + Reverse(target.BaseSequence.Substring(1))
                : Reverse(target.BaseSequence);
            var decoy = insulinProteins.FirstOrDefault(p => p.IsDecoy && p.BaseSequence == expectedDecoySeq)
                        ?? insulinProteins.First(p => p.IsDecoy && p.Length == target.Length);

            Assert.IsFalse(string.IsNullOrWhiteSpace(target.Accession));
            Assert.IsTrue(decoy.Accession.StartsWith("DECOY_"));
            Assert.AreEqual(target.Length, decoy.Length);
            Assert.AreNotEqual(target.BaseSequence, decoy.BaseSequence);
            Assert.AreEqual(expectedDecoySeq, decoy.BaseSequence, "Decoy must follow 'retain M, reverse remainder' rule");

            var dp = new DigestionParams(protease: "top-down");
            List<PeptideWithSetModifications> targetTruncs = target
                .Digest(dp, new List<Modification>(), new List<Modification>(), topDownTruncationSearch: true)
                .Cast<PeptideWithSetModifications>().ToList();
            List<PeptideWithSetModifications> decoyTruncs = decoy
                .Digest(dp, new List<Modification>(), new List<Modification>(), topDownTruncationSearch: true)
                .Cast<PeptideWithSetModifications>().ToList();

            // Parity and sanity checks for the selected pair
            Assert.AreEqual(68, targetTruncs.Count, "Target should yield 68 truncation products in top-down mode");
            Assert.AreEqual(68, decoyTruncs.Count, "Decoy should yield 68 truncation products in top-down mode");
            Assert.That(targetTruncs.All(p => p.DigestionParams?.DigestionAgent?.Name == "top-down"));
            Assert.That(decoyTruncs.All(p => p.DigestionParams?.DigestionAgent?.Name == "top-down"));
            Assert.AreEqual(targetTruncs.Count, targetTruncs.Select(p => p.BaseSequence).Distinct().Count());
            Assert.AreEqual(decoyTruncs.Count, decoyTruncs.Select(p => p.BaseSequence).Distinct().Count());
            Assert.IsTrue(targetTruncs.Any(p => p.BaseSequence == target.BaseSequence));
            Assert.IsTrue(decoyTruncs.Any(p => p.BaseSequence == decoy.BaseSequence));

            // Build the table rows
            static bool HasRetainedMet(PeptideWithSetModifications p) =>
                p.OneBasedStartResidueInProtein == 1 &&
                p.Parent?.BaseSequence?.Length > 0 &&
                p.Parent.BaseSequence[0] == 'M';

            // We only compare the combined truncations for the chosen target/decoy in this test (68 + 68 rows)
            var rows = targetTruncs.Concat(decoyTruncs)
                .Select(pep =>
                {
                    bool isDecoy = (pep.Parent as Protein)?.IsDecoy == true;
                    string type = isDecoy ? "Decoy" : "Target";
                    string retained = HasRetainedMet(pep) ? "TRUE" : "FALSE"; // normalize to match expected
                    return string.Join("\t", new[]
                    {
                        pep.BaseSequence,
                        type,
                        pep.OneBasedStartResidueInProtein.ToString(),
                        pep.OneBasedEndResidueInProtein.ToString(),
                        retained
                    });
                })
                // Sort deterministically to avoid platform/iteration-order differences
                .OrderBy(r => r, StringComparer.Ordinal)
                .ToList();

            var header = "Sequence\tType\tBegin\tEnd\tRetainedMethionine";
            var outputLines = new List<string>(capacity: rows.Count + 1) { header };
            outputLines.AddRange(rows);

            // Persist the generated table for inspection
            string workPath = Path.Combine(TestContext.CurrentContext.WorkDirectory, "topdown_truncations_table.tsv");
            File.WriteAllLines(workPath, outputLines);
            Console.WriteLine($"Generated truncation table: {workPath}");

            // Load expected and compare
            string expectedPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "truncationsExpected.tsv");
            Assert.That(File.Exists(expectedPath), $"Expected file not found: {expectedPath}");
            var expectedAll = File.ReadAllLines(expectedPath)
                                  .Where(l => l is not null)
                                  .Select(l => l.TrimEnd('\r', '\n'))
                                  .ToList();

            Assert.That(expectedAll.Count > 0, "Expected file is empty");
            string expectedHeader = expectedAll[0];
            var expectedRows = expectedAll.Skip(1)
                                          .Where(l => !string.IsNullOrWhiteSpace(l))
                                          .OrderBy(l => l, StringComparer.Ordinal)
                                          .ToList();

            // Header check
            if (!string.Equals(header, expectedHeader, StringComparison.Ordinal))
            {
                TestContext.Out.WriteLine($"Header mismatch:");
                TestContext.Out.WriteLine($"  Expected: {expectedHeader}");
                TestContext.Out.WriteLine($"  Actual  : {header}");
            }

            // Multiset comparison for rows (counts of duplicates matter)
            static Dictionary<string, int> ToCounts(IEnumerable<string> lines)
                => lines.GroupBy(x => x, StringComparer.Ordinal)
                        .ToDictionary(g => g.Key, g => g.Count(), StringComparer.Ordinal);

            var expCounts = ToCounts(expectedRows);
            var gotCounts = ToCounts(rows);

            var missing = new List<string>();   // in expected more times than in actual
            var extra = new List<string>();     // in actual more times than in expected

            foreach (var kv in expCounts)
            {
                gotCounts.TryGetValue(kv.Key, out int got);
                if (got < kv.Value)
                {
                    int deficit = kv.Value - got;
                    for (int i = 0; i < deficit; i++) missing.Add(kv.Key);
                }
            }
            foreach (var kv in gotCounts)
            {
                expCounts.TryGetValue(kv.Key, out int exp);
                if (kv.Value > exp)
                {
                    int surplus = kv.Value - exp;
                    for (int i = 0; i < surplus; i++) extra.Add(kv.Key);
                }
            }

            if (missing.Count == 0 && extra.Count == 0)
            {
                TestContext.Out.WriteLine("Top-down truncation table matches expected.");
                TestContext.Out.WriteLine("Sample (first 5 rows):");
                foreach (var l in outputLines.Take(6)) TestContext.Out.WriteLine(l);
            }
            else
            {
                TestContext.Out.WriteLine("Top-down truncation table differs from expected.");
                TestContext.Out.WriteLine($"Missing rows (expected but not found or under-counted): {missing.Count}");
                foreach (var l in missing.Take(20)) TestContext.Out.WriteLine($"  MISSING: {l}");
                if (missing.Count > 20) TestContext.Out.WriteLine($"  ...and {missing.Count - 20} more");

                TestContext.Out.WriteLine($"Extra rows (found but not expected or over-counted): {extra.Count}");
                foreach (var l in extra.Take(20)) TestContext.Out.WriteLine($"  EXTRA:   {l}");
                if (extra.Count > 20) TestContext.Out.WriteLine($"  ...and {extra.Count - 20} more");

                Assert.Fail($"Generated top-down truncation table does not match expected.\nExpected file: {expectedPath}\nActual file: {workPath}\nMissing: {missing.Count}, Extra: {extra.Count}");
            }
        }
    }
}
