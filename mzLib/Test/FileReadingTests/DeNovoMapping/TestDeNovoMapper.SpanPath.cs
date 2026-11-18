using NUnit.Framework;
using Omics;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers.DeNovoMapping;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Test.FileReadingTests.DeNovoMapping
{
    /// <summary>
    /// The mapper looks digestion windows up as spans instead of digesting each protein into peptide objects. These tests
    /// hold it to what <see cref="IBioPolymer.Digest"/> gives.
    /// </summary>
    internal partial class TestDeNovoMapper
    {
        private static readonly string[] SpanPathProteins =
        [
            ProteinA,
            ProteinB,
            "MIKLLRPAAAKIILLRGGGGGGGK",      // initiator Met, I/L, K before P
            "MKAAAAAAAAR",                    // Met then an immediate cleavage site
            "M",
            "MK",
            "K",
            "KKKKRRRRPK",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAA",  // no cleavage site
            "DAAAAAADKAAAAAAAMAAAAAAWPAAAAAAEK",
        ];

        private static IEnumerable<TestCaseData> SpanPathDigestions()
        {
            yield return new TestCaseData(new DigestionParams("trypsin", 2, 7), true).SetName("SpanPath(trypsin, mc2, len7+)");
            yield return new TestCaseData(new DigestionParams("trypsin", 0, 1), true).SetName("SpanPath(trypsin, mc0, len1+)");
            yield return new TestCaseData(new DigestionParams("trypsin", 3, 4, 12), false).SetName("SpanPath(trypsin, mc3, len4-12, no I/L)");
            yield return new TestCaseData(new DigestionParams("trypsin|P", 2, 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Cleave), true)
                .SetName("SpanPath(trypsin|P, Met cleaved)");
            yield return new TestCaseData(new DigestionParams("trypsin", 1, 2, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain), true)
                .SetName("SpanPath(trypsin, Met retained)");
            yield return new TestCaseData(new DigestionParams("Asp-N", 2, 3), true).SetName("SpanPath(Asp-N)");
            yield return new TestCaseData(new DigestionParams("Lys-N", 1, 3), true).SetName("SpanPath(Lys-N)");
            yield return new TestCaseData(new DigestionParams("CNBr", 2, 3), true).SetName("SpanPath(CNBr, cleavage modification)");
            yield return new TestCaseData(new DigestionParams("non-specific", 3, 2, 6), true).SetName("SpanPath(non-specific protease)");
            yield return new TestCaseData(new DigestionParams("StcE-trypsin", 2, 3), true).SetName("SpanPath(StcE-trypsin)");
        }

        private static IEnumerable<TestCaseData> FallbackDigestions()
        {
            yield return new TestCaseData(new DigestionParams("trypsin", 2, 7, searchModeType: CleavageSpecificity.Semi)).SetName("Fallback(semi)");
            yield return new TestCaseData(new DigestionParams("trypsin", 2, 7, searchModeType: CleavageSpecificity.Semi,
                fragmentationTerminus: FragmentationTerminus.N)).SetName("Fallback(semi seeds)");
            yield return new TestCaseData(new DigestionParams("trypsin", 2, 7, 20, searchModeType: CleavageSpecificity.None)).SetName("Fallback(non-specific search)");
            yield return new TestCaseData(new DigestionParams("top-down", 0, 1)).SetName("Fallback(top-down)");
            yield return new TestCaseData(new DigestionParams("trypsin", 2, 5, keepNGlycopeptide: true)).SetName("Fallback(N-glycopeptides)");
        }

        [TestCaseSource(nameof(SpanPathDigestions))]
        public void WindowsAreTheDigestPeptides(DigestionParams digestionParams, bool replaceIWithL)
        {
            var digester = new MappingSequenceDigester(digestionParams, replaceIWithL);
            foreach (var protein in SpanPathProteinsWithDecoysAndRandom())
            {
                Assert.That(digester.EnumeratesWindows(protein), Is.True, protein.BaseSequence);
                Assert.That(Windows(digester, protein), Is.EquivalentTo(DigestPeptides(protein, digestionParams, replaceIWithL)), protein.BaseSequence);
            }
        }

        [TestCaseSource(nameof(FallbackDigestions))]
        public void OtherDigestionsGoThroughDigest(DigestionParams digestionParams)
        {
            var digester = new MappingSequenceDigester(digestionParams, true);
            foreach (var protein in SpanPathProteinsWithDecoysAndRandom().Take(20))
            {
                Assert.That(digester.EnumeratesWindows(protein), Is.False);
                Assert.That(Windows(digester, protein), Is.EquivalentTo(DigestPeptides(protein, digestionParams, true)), protein.BaseSequence);
            }
        }

        [TestCase(0, 7, InitiatorMethionineBehavior.Variable)]
        [TestCase(2, 4, InitiatorMethionineBehavior.Variable)]
        [TestCase(3, 1, InitiatorMethionineBehavior.Cleave)]
        [TestCase(1, 5, InitiatorMethionineBehavior.Retain)]
        public void TruncationProductWindowsAreTheDigestPeptides(int maxMissedCleavages, int minLength, InitiatorMethionineBehavior initiatorMethionine)
        {
            // UniProt XML gives most proteins chains, signal peptides and propeptides, so the host's are enumerated too.
            var digestionParams = new DigestionParams("trypsin", maxMissedCleavages, minLength, 40, initiatorMethionineBehavior: initiatorMethionine);
            var digester = new MappingSequenceDigester(digestionParams, true);
            var random = new Random(maxMissedCleavages * 100 + minLength);

            for (int p = 0; p < 300; p++)
            {
                string sequence = RandomProteinSequence(random);
                var truncations = new List<TruncationProduct>();
                for (int t = random.Next(1, 4); t > 0; t--)
                {
                    int begin = random.Next(1, sequence.Length + 1);
                    int end = random.Next(begin, sequence.Length + 1);
                    // Ends on a cleavage site, at the termini, and the whole protein all occur in UniProt.
                    if (random.Next(5) == 0) begin = 1;
                    if (random.Next(5) == 0) end = sequence.Length;
                    if (random.Next(10) == 0 && end >= 2) begin = 2;
                    truncations.Add(new TruncationProduct(begin, end, "chain"));
                }
                if (random.Next(10) == 0)
                    truncations.Add(new TruncationProduct(null, sequence.Length, "unknown start"));

                var protein = new Protein(sequence, $"TRUNC{p}", proteolysisProducts: truncations);
                Assert.That(digester.EnumeratesWindows(protein), Is.True);
                Assert.That(Windows(digester, protein), Is.EquivalentTo(DigestPeptides(protein, digestionParams, true)),
                    $"{sequence} {string.Join(";", truncations.Select(t => $"{t.OneBasedBeginPosition}-{t.OneBasedEndPosition}"))}");
            }
        }

        [Test]
        public void ProteinsDigestCannotEnumerateGoThroughDigest()
        {
            var digestionParams = new DigestionParams("trypsin", 2, 4, respectCleavageBlockingModifications: true);
            ModificationMotif.TryGetMotif("K", out var motif);
            var succinyl = new Modification(_originalId: "N6-succinyllysine", _modificationType: "Test", _target: motif,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 100.01604);

            var truncated = new Protein("MAAAAKAAAAAAAKAAAAAAAAAAAKAAAAAAR", "TRUNC",
                proteolysisProducts: [new TruncationProduct(3, 20, "chain"), new TruncationProduct(12, 33, "peptide")]);
            var reversedTruncation = new Protein("MAAAAKAAAAAAAKAAAAAAAAAAAKAAAAAAR", "TRUNC_REVERSED",
                proteolysisProducts: [new TruncationProduct(20, 3, "chain")]);
            var blocked = new Protein("AAAAKAAAAKAAAAKAAAAR", "BLOCKED",
                oneBasedModifications: new Dictionary<int, List<Modification>> { [10] = [succinyl] });
            var plain = new Protein("AAAAKAAAAKAAAAKAAAAR", "PLAIN");
            var digester = new MappingSequenceDigester(digestionParams, true);

            Assert.Multiple(() =>
            {
                Assert.That(digester.EnumeratesWindows(truncated), Is.True);
                Assert.That(digester.EnumeratesWindows(reversedTruncation), Is.False);
                Assert.That(digester.EnumeratesWindows(blocked), Is.False);
                Assert.That(digester.EnumeratesWindows(plain), Is.True);
                foreach (var protein in new[] { truncated, reversedTruncation, blocked, plain })
                    Assert.That(Windows(digester, protein), Is.EquivalentTo(DigestPeptides(protein, digestionParams, true)), protein.Accession);

                // The policy's read-through peptides are why the blocked protein cannot be enumerated from its sites.
                Assert.That(Windows(digester, blocked), Is.Not.EquivalentTo(Windows(digester, plain)));
            });
        }

        [Test]
        public void TruncationProductOutsideTheProteinGoesThroughDigest()
        {
            var digester = new MappingSequenceDigester(new DigestionParams("trypsin", 2, 4), true);
            var outside = new Protein("MAAAAKAAAAAAAKAAAAAAAAAAAKAAAAAAR", "OUTSIDE", proteolysisProducts: [new TruncationProduct(3, 50, "chain")]);

            Assert.That(digester.EnumeratesWindows(outside), Is.False);
        }

        [TestCase(1, false)]
        [TestCase(1, true)]
        [TestCase(4, false)]
        [TestCase(4, true)]
        public void MapDatabaseMatchesTheDigestPath(int maxThreads, bool withHost)
        {
            var random = new Random(1506);
            var digestionParams = new DigestionParams("trypsin", 2, 5, initiatorMethionineBehavior: InitiatorMethionineBehavior.Variable);
            var targets = SpanPathProteins.Select((s, i) => new Protein(s, $"P{i}"))
                .Concat(Enumerable.Range(0, 60).Select(i => new Protein(RandomProteinSequence(random), $"R{i}")))
                .Append(new Protein("MKVLSAADKGGGGGGGR", "SAME_ACCESSION"))
                .Append(new Protein("TPEPTIDERGGGGGGGK", "SAME_ACCESSION"))
                .Append(new Protein("MAAAAKAAAAAAAKAAAAAAAAAAAKAAAAAAR", "TRUNC", proteolysisProducts: [new TruncationProduct(3, 20, "chain")]))
                .ToList();
            var decoys = targets.Select(p => new Protein(Reverse(p.BaseSequence), "DECOY_" + p.Accession, isDecoy: true))
                // A decoy producing a target peptide, which must not count as a decoy.
                .Append(new Protein("GGGGGGGRVLSAADKWWWWWWWK", "DECOY_CLASH", isDecoy: true))
                .ToList();
            var host = Enumerable.Range(0, 20).Select(i => new Protein(RandomProteinSequence(random), $"H{i}"))
                .Concat(targets.Skip(12).Take(10).Select(p => new Protein(p.BaseSequence, "HOST_" + p.Accession)))
                .Append(new Protein(ProteinA, "HOST_A"))
                .Append(new Protein("GGGGGGGRTPEPTLDER", "DECOY_HOST", isDecoy: true))
                .ToList();
            var database = targets.Concat(decoys).Cast<IBioPolymer>().ToList();

            // Predictions: peptides of the targets, decoys and host (some written with I), duplicates, and misses.
            var peptides = database.Concat(host).SelectMany(p => DigestPeptides(p, digestionParams, false)).Distinct().OrderBy(s => s, StringComparer.Ordinal).ToList();
            var psms = new List<(string, double, double?)>();
            for (int i = 0; i < 400; i++)
            {
                string sequence = random.Next(10) == 0 ? RandomProteinSequence(random)[..7] : peptides[random.Next(peptides.Count)];
                if (random.Next(4) == 0)
                    sequence = sequence.Replace('L', 'I');
                psms.Add((sequence, Math.Round(random.NextDouble(), 4), 60.0 * i));
            }
            string mzTab = WriteMzTab("RunSpan.mztab", "file:///Z:/RunSpan.mgf", psms.ToArray());

            var parameters = new DeNovoMappingParameters
            {
                DigestionParams = digestionParams,
                MapRetentionTimes = false,
                MaxThreads = maxThreads,
                HostProteins = withHost ? host : null,
            };
            var mapper = new DeNovoMapper([mzTab], parameters);
            var actual = mapper.MapDatabase("Span", database);
            var expected = MapByDigestion(mapper, parameters, "Span", database, withHost ? host : null);

            Assert.That(actual.TargetPredictions, Is.GreaterThan(50), "the test should exercise many hits");
            Assert.That(actual.DecoyPredictions, Is.GreaterThan(5));
            if (withHost)
                Assert.That(actual.AmbiguousTargetPredictions, Is.GreaterThan(5));
            AssertSameResult(actual, expected);
        }

        #region Span path helpers

        private static IEnumerable<IBioPolymer> SpanPathProteinsWithDecoysAndRandom()
        {
            var random = new Random(42);
            var targets = SpanPathProteins.Concat(Enumerable.Range(0, 40).Select(_ => RandomProteinSequence(random))).ToList();
            foreach (var (sequence, i) in targets.Select((s, i) => (s, i)))
            {
                yield return new Protein(sequence, $"T{i}");
                yield return new Protein(Reverse(sequence), $"DECOY_T{i}", isDecoy: true);
            }
        }

        private static string Reverse(string sequence) => new(sequence.Reverse().ToArray());

        /// <summary>
        /// A protein rich in cleavage sites, prolines, I and L, starting with M half the time.
        /// </summary>
        private static string RandomProteinSequence(Random random)
        {
            const string residues = "ACDEFGHIKLMNPQRSTVWYKRKRPILIL";
            var builder = new StringBuilder(random.Next(2) == 0 ? "M" : "");
            int length = random.Next(8, 90);
            while (builder.Length < length)
                builder.Append(residues[random.Next(residues.Length)]);
            return builder.ToString();
        }

        private static List<string> Windows(MappingSequenceDigester digester, IBioPolymer protein)
        {
            var sink = new WindowList();
            digester.Digest(protein, sink);
            return sink.Windows.Distinct().ToList();
        }

        private static List<string> DigestPeptides(IBioPolymer protein, DigestionParams digestionParams, bool replaceIWithL) =>
            protein.Digest(digestionParams, new List<Modification>(), new List<Modification>())
                .Select(p => DeNovoPrediction.NormalizeSequence(p.BaseSequence, replaceIWithL))
                .Distinct()
                .ToList();

        private sealed class WindowList : IMappingWindowSink
        {
            public List<string> Windows { get; } = new();
            public void OnWindow(string normalizedSource, int start, int length) => Windows.Add(normalizedSource.Substring(start, length));
        }

        /// <summary>
        /// The mapping as it was done before the span path: every protein digested into peptides, every prediction visited.
        /// </summary>
        private static DeNovoMappingResult MapByDigestion(DeNovoMapper mapper, DeNovoMappingParameters parameters, string identifier,
            List<IBioPolymer> proteins, List<Protein>? host)
        {
            HashSet<string> Digest(IBioPolymer protein) =>
                DigestPeptides(protein, (DigestionParams)parameters.DigestionParams, parameters.ReplaceIWithL).ToHashSet();

            var predictions = mapper.Predictions;
            var bySequence = predictions.Select((p, i) => (p.MappingSequence, i)).GroupBy(t => t.MappingSequence)
                .ToDictionary(g => g.Key, g => g.Select(t => t.i).ToArray());

            var hostShared = host?.Where(p => !p.IsDecoy).SelectMany(Digest).Where(bySequence.ContainsKey).ToHashSet() ?? [];
            Assert.That(predictions.Where(p => p.IsHostShared).Select(p => p.MappingSequence).ToHashSet(), Is.EquivalentTo(hostShared),
                "host marking");
            bool hasHost = host is not null;

            var targetHit = new bool[predictions.Count];
            var decoyHit = new bool[predictions.Count];
            var databasePeptides = new HashSet<string>();
            var mappedPeptides = new HashSet<string>();
            var mappedProteins = new HashSet<string>();
            var unambiguousPeptides = new HashSet<string>();
            var unambiguousProteins = new HashSet<string>();
            var targets = proteins.Where(p => !p.IsDecoy).ToList();

            foreach (var protein in targets)
            {
                bool matched = false, matchedUnambiguously = false;
                foreach (var peptide in Digest(protein))
                {
                    databasePeptides.Add(peptide);
                    if (!bySequence.TryGetValue(peptide, out var indices))
                        continue;
                    foreach (int p in indices)
                        targetHit[p] = true;
                    mappedPeptides.Add(peptide);
                    matched = true;
                    if (hasHost && !hostShared.Contains(peptide))
                    {
                        unambiguousPeptides.Add(peptide);
                        matchedUnambiguously = true;
                    }
                }
                if (matched)
                    mappedProteins.Add(protein.Accession);
                if (matchedUnambiguously)
                    unambiguousProteins.Add(protein.Accession);
            }

            foreach (var protein in proteins.Where(p => p.IsDecoy))
            foreach (var peptide in Digest(protein))
            {
                if (databasePeptides.Contains(peptide) || !bySequence.TryGetValue(peptide, out var indices))
                    continue;
                foreach (int p in indices)
                    decoyHit[p] = true;
            }

            var result = new DeNovoMappingResult
            {
                DatabaseIdentifier = identifier,
                DatabaseProteinCount = targets.Count,
                DatabasePeptideCount = databasePeptides.Count,
                UniquePeptidesMapped = mappedPeptides.Count,
                UniqueProteinsMapped = mappedProteins.Count,
                HostProteinCount = hasHost ? host!.Count(p => !p.IsDecoy) : -1,
                UniqueUnambiguousPeptidesMapped = unambiguousPeptides.Count,
                UniqueUnambiguousProteinsMapped = unambiguousProteins.Count,
            };

            for (int i = 0; i < predictions.Count; i++)
            {
                var prediction = predictions[i];
                bool unambiguous = hasHost && !hostShared.Contains(prediction.MappingSequence);
                if (targetHit[i])
                {
                    result.TargetPredictions++;
                    result.TargetPredictionScores.Add(prediction.Score);
                    if (prediction is { ObservedRetentionTime: { } observed, CalibratedRetentionTime: { } calibrated })
                    {
                        result.RetentionTimeErrors.Add(observed - calibrated);
                        if (unambiguous)
                            result.UnambiguousRetentionTimeErrors.Add(observed - calibrated);
                    }
                    if (unambiguous)
                    {
                        result.UnambiguousTargetPredictions++;
                        result.UnambiguousTargetPredictionScores.Add(prediction.Score);
                    }
                    else if (hasHost)
                    {
                        result.AmbiguousTargetPredictions++;
                    }
                }
                if (decoyHit[i])
                {
                    result.DecoyPredictions++;
                    result.DecoyPredictionScores.Add(prediction.Score);
                    if (unambiguous)
                    {
                        result.UnambiguousDecoyPredictions++;
                        result.UnambiguousDecoyPredictionScores.Add(prediction.Score);
                    }
                }
                if (targetHit[i] || decoyHit[i])
                    result.PredictionCountsByFile[prediction.FileKey] = result.PredictionCountsByFile.GetValueOrDefault(prediction.FileKey) + 1;
            }

            result.FinalizeValues();
            return result;
        }

        private static void AssertSameResult(DeNovoMappingResult actual, DeNovoMappingResult expected)
        {
            Assert.Multiple(() =>
            {
                // Every column of the summary, with the lists in the same order.
                foreach (var property in typeof(DeNovoMappingResult).GetProperties())
                {
                    object? a = property.GetValue(actual), e = property.GetValue(expected);
                    if (a is IEnumerable<double> actualList)
                        Assert.That(actualList, Is.EqualTo((IEnumerable<double>)e!), property.Name);
                    else if (a is double actualValue)
                        Assert.That(actualValue, double.IsNaN((double)e!) ? Is.NaN : Is.EqualTo((double)e!), property.Name);
                    else
                        Assert.That(a, Is.EqualTo(e), property.Name);
                }
            });
        }

        #endregion
    }
}
