using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.FlashLFQ
{
    /// <summary>
    /// Real-data match-between-runs tests that share a single FlashLFQ run.
    ///
    /// Both <see cref="RealDataMbrTest"/> and <see cref="MbrTargetDecoyScoreDistributionTest"/> analyze the
    /// same MBR run over the sliced f1r1/f1r2 data. Running the engine (which crosses the PEP path) is by far
    /// the most expensive thing here, so it is executed exactly once in <see cref="RunSharedMbrEngineOnce"/>
    /// and both tests consume the stashed <see cref="FlashLfqResults"/>. RealDataMbrTest additionally runs the
    /// engine with two other configurations to exercise the requireMsmsIdInCondition and
    /// peptideSequencesToQuantify features; those are unrelated to the shared analysis and remain in the test.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal class MbrRealDataTests
    {
        // Inputs + results of the shared MBR run, built once in OneTimeSetUp and consumed by both tests.
        private static List<Identification> _ids;
        private static SpectraFileInfo _f1r1;
        private static SpectraFileInfo _f1r2;
        private static FlashLfqResults _sharedResults;

        // Config shared by both tests. This crosses the thresholds (>100 MBR peaks, >20 random-RT decoys)
        // that trigger the PEP path, so the scorer is fully exercised.
        [OneTimeSetUp]
        public static void RunSharedMbrEngineOnce()
        {
            string psmFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "FlashLFQ", "TestData", @"PSMsForMbrTest.psmtsv");

            _f1r1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "FlashLFQ", "TestData", @"f1r1_sliced_mbr.raw"), "a", 0, 0, 0);
            _f1r2 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "FlashLFQ", "TestData", @"f1r2_sliced_mbr.raw"), "a", 1, 0, 0);

            _ids = LoadMbrIdentifications(psmFile, _f1r1, _f1r2);

            // Remove any stale PEP model.zip left in the input folder by an older build, so the regression
            // assertion in RealDataMbrTest reflects only what this run produced.
            DeletePepModelArtifact();

            var engine = new FlashLfqEngine(_ids, matchBetweenRuns: true, requireMsmsIdInCondition: false, maxThreads: 1, matchBetweenRunsFdrThreshold: 0.15, maxMbrWindow: 1);
            _sharedResults = engine.Run();
        }

        [OneTimeTearDown]
        public static void CleanUpPepModelArtifact()
        {
            // FlashLFQ's PEP engine used to drop a scratch "model.zip" into the input data folder
            // (mzLib#1124). The fix routes it to a temp dir; clean up here as well so the test suite
            // never accumulates that artifact. Done once, after both tests, so it can't mask the
            // RealDataMbrTest regression assertion regardless of test execution order.
            DeletePepModelArtifact();
        }

        [Test]
        public static void RealDataMbrTest()
        {
            var ids = _ids;
            var f1r1 = _f1r1;
            var f1r2 = _f1r2;
            var results = _sharedResults;

            // REGRESSION (mzLib#1124): this data triggers the PEP path, which must write its scratch
            // model.zip to a temp dir, NOT the input data folder. Assert before TearDown cleans up.
            Assert.That(File.Exists(Path.Combine(TestContext.CurrentContext.TestDirectory, "FlashLFQ", "TestData", "model.zip")), Is.False,
                "FlashLFQ left its PEP model.zip in the input data folder.");

            // Count the number of MBR results in each file
            var f1r1MbrResults = results
                .PeptideModifiedSequences
                .Where(p => p.Value.GetDetectionType(f1r1) == DetectionType.MBR && p.Value.GetDetectionType(f1r2) == DetectionType.MSMS)
                .ToList();
            var f1r2MbrResults = results
                .PeptideModifiedSequences
                .Where(p => p.Value.GetDetectionType(f1r1) == DetectionType.MSMS && p.Value.GetDetectionType(f1r2) == DetectionType.MBR)
                .ToList();

            // Due to the small number of results in the test data, the counts are machine-dependent: the PEP
            // path trains an ML.NET model whose output varies across platforms, thread counts, and runtime
            // versions, so the exact count shifts by a few from run to run (e.g. 142/77 on CI, 141/76 elsewhere).
            // Assert a tolerant range around the expected values rather than an exact equality that is red on
            // some machines. See smith-chem-wisc/mzLib#1283 for the planned top-N-transfers rewrite of this test.
            Console.WriteLine("r1 PIP event count: " + f1r1MbrResults.Count);
            Console.WriteLine("r2 PIP event count: " + f1r2MbrResults.Count);
            Assert.That(f1r1MbrResults.Count, Is.EqualTo(142).Within(10));
            Assert.That(f1r2MbrResults.Count, Is.EqualTo(77).Within(10));

            // Check that MS/MS identified peaks and MBR identified peaks have similar intensities
            List<(double, double)> peptideIntensities = f1r1MbrResults.Select(pep => (Math.Log(pep.Value.GetIntensity(f1r1)), Math.Log(pep.Value.GetIntensity(f1r2)))).ToList();
            double corrRun1 = Correlation.Pearson(peptideIntensities.Select(p => p.Item1), peptideIntensities.Select(p => p.Item2));

            peptideIntensities = f1r2MbrResults.Select(pep => (Math.Log(pep.Value.GetIntensity(f1r1)), Math.Log(pep.Value.GetIntensity(f1r2)))).ToList();
            double corrRun2 = Correlation.Pearson(peptideIntensities.Select(p => p.Item1), peptideIntensities.Select(p => p.Item2));

            // These values are also sensitive, changes can cause them to dip as low as 0.6 (specifically the corrRun2 value)
            Console.WriteLine("r1 correlation: " + corrRun1);
            Console.WriteLine("r2 correlation: " + corrRun2);
            Assert.Greater(corrRun1, 0.75);
            Assert.Greater(corrRun2, 0.65);

            // the "requireMsmsIdInCondition" field requires that at least one MS/MS identification from a protein
            // has to be observed in a condition for match-between-runs

            f1r1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "FlashLFQ", "TestData", @"f1r1_sliced_mbr.raw"), "b", 0, 0, 0);
            var engine = new FlashLfqEngine(ids, matchBetweenRuns: true, requireMsmsIdInCondition: true, maxThreads: 5);
            results = engine.Run();
            var proteinsObservedInF1 = ids.Where(id => !id.IsDecoy).Where(p => p.FileInfo == f1r1).SelectMany(p => p.ProteinGroups).Distinct().ToList();
            var proteinsObservedInF2 = ids.Where(id => !id.IsDecoy).Where(p => p.FileInfo == f1r2).SelectMany(p => p.ProteinGroups).Distinct().ToList();
            var proteinsObservedInF1ButNotF2 = proteinsObservedInF1.Except(proteinsObservedInF2).ToList();
            foreach (ProteinGroup protein in proteinsObservedInF1ButNotF2)
            {
                Assert.That(results.ProteinGroups[protein.ProteinGroupName].GetIntensity(f1r2) == 0);
            }

            // Test that no decoys are reported in the final resultsw
            Assert.AreEqual(0, ids.Where(id => id.IsDecoy).Count(id => results.ProteinGroups.ContainsKey(id.ProteinGroups.First().ProteinGroupName)));

            List<string> peptidesToUse = ids.Where(id => id.QValue <= 0.007 & !id.IsDecoy).Select(id => id.ModifiedSequence).Distinct().ToList();
            engine = new FlashLfqEngine(ids, matchBetweenRuns: true, requireMsmsIdInCondition: true, maxThreads: 1, matchBetweenRunsFdrThreshold: 0.5, maxMbrWindow: 1, peptideSequencesToQuantify: peptidesToUse);
            results = engine.Run();

            CollectionAssert.AreEquivalent(results.PeptideModifiedSequences.Select(kvp => kvp.Key), peptidesToUse);
        }

        /// <summary>
        /// Uses the same real MBR run as <see cref="RealDataMbrTest"/>, but instead of checking
        /// intensity correlations it inspects the score distributions of target vs. decoy MBR transfers.
        /// After a run, results.Peaks retains every scored MBR peak — real targets plus two flavors of
        /// decoy: random-retention-time decoys (RandomRt) and decoy-peptide transfers (DecoyPeptide).
        /// For each score component (the overall MbrScore and the five sub-scores that compose it) we
        /// report how well targets separate from decoys using a rank-based AUC (the probability that a
        /// randomly chosen target outscores a randomly chosen decoy; 0.5 = no separation, 1 = perfect).
        /// This makes it easy to see which individual score components actually discriminate.
        /// </summary>
        [Test]
        public static void MbrTargetDecoyScoreDistributionTest()
        {
            var results = _sharedResults;

            // Gather every scored MBR peak across both acceptor files. The MbrQValue threshold only
            // affects which peaks are quantified into PeptideModifiedSequences; results.Peaks still holds
            // the full target+decoy population that the scorer produced.
            List<MbrChromatographicPeak> mbrPeaks = results.Peaks
                .SelectMany(kvp => kvp.Value)
                .OfType<MbrChromatographicPeak>()
                .ToList();

            // Four categories, per the FDR model in CalculateFdrForMbrPeaks:
            //   (DecoyPeptide, RandomRt) => target / decoy-peptide / random-RT decoy / double decoy
            var targets = mbrPeaks.Where(p => !p.DecoyPeptide && !p.RandomRt).ToList();
            var randomRtDecoys = mbrPeaks.Where(p => !p.DecoyPeptide && p.RandomRt).ToList();
            var decoyPeptides = mbrPeaks.Where(p => p.DecoyPeptide && !p.RandomRt).ToList();
            var doubleDecoys = mbrPeaks.Where(p => p.DecoyPeptide && p.RandomRt).ToList();
            var allDecoys = mbrPeaks.Where(p => p.DecoyPeptide || p.RandomRt).ToList();

            Console.WriteLine($"Total MBR peaks: {mbrPeaks.Count}");
            Console.WriteLine($"  Targets:            {targets.Count}");
            Console.WriteLine($"  Random-RT decoys:   {randomRtDecoys.Count}");
            Console.WriteLine($"  Decoy peptides:     {decoyPeptides.Count}");
            Console.WriteLine($"  Double decoys:      {doubleDecoys.Count}");
            Console.WriteLine($"  All decoys:         {allDecoys.Count}");

            // Need a meaningful population of each class for the comparison to mean anything. The decoy
            // count is modest because results.Peaks is heavily deduplicated after scoring, so keep the
            // floor low enough to survive run-to-run ML.NET variability in the PEP path.
            Assert.That(targets.Count, Is.GreaterThan(50), "Not enough target MBR transfers to compare distributions.");
            Assert.That(allDecoys.Count, Is.GreaterThan(10), "Not enough decoy MBR transfers to compare distributions.");

            var scoreComponents = new (string Name, Func<MbrChromatographicPeak, double> Selector)[]
            {
                ("MbrScore (overall)",        p => p.MbrScore),
                ("IntensityScore",            p => p.IntensityScore),
                ("RtScore",                   p => p.RtScore),
                ("PpmScore",                  p => p.PpmScore),
                ("ScanCountScore",            p => p.ScanCountScore),
                ("IsotopicDistributionScore", p => p.IsotopicDistributionScore),
            };

            Console.WriteLine();
            Console.WriteLine($"{"Score component",-28}{"targetMed",12}{"decoyMed",12}{"rtDecoyMed",12}{"pepDecoyMed",12}{"AUC(t>d)",10}");

            double overallScoreAuc = double.NaN;
            foreach (var (name, selector) in scoreComponents)
            {
                List<double> targetScores = targets.Select(selector).ToList();
                List<double> decoyScores = allDecoys.Select(selector).ToList();

                double targetMedian = targetScores.Median();
                double decoyMedian = decoyScores.Median();
                double rtDecoyMedian = randomRtDecoys.Count > 0 ? randomRtDecoys.Select(selector).Median() : double.NaN;
                double pepDecoyMedian = decoyPeptides.Count > 0 ? decoyPeptides.Select(selector).Median() : double.NaN;
                double auc = RankSumAuc(targetScores, decoyScores);

                Console.WriteLine($"{name,-28}{targetMedian,12:F4}{decoyMedian,12:F4}{rtDecoyMedian,12:F4}{pepDecoyMedian,12:F4}{auc,10:F3}");

                if (name.StartsWith("MbrScore"))
                    overallScoreAuc = auc;
            }

            // The composite MBR score is the discriminant the FDR model relies on, so at minimum it must
            // rank targets above decoys. (AUC > 0.5 means a random target beats a random decoy more often
            // than not.) Individual sub-scores are reported above for inspection but not asserted, since
            // the point of this test is to observe which components separate and which don't.
            Assert.That(overallScoreAuc, Is.GreaterThan(0.5),
                "Targets did not rank above decoys on the overall MBR score.");
        }

        /// <summary>
        /// Loads the Identifications used by the real-data MBR tests from a MetaMorpheus psmtsv file,
        /// assigning each PSM to its acceptor file by name and preserving the target/decoy flag.
        /// </summary>
        private static List<Identification> LoadMbrIdentifications(string psmFile, SpectraFileInfo f1r1, SpectraFileInfo f1r2)
        {
            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });

                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains("f1r1"))
                {
                    file = f1r1;
                }
                else if (split[0].Contains("f1r2"))
                {
                    file = f1r2;
                }

                string baseSequence = split[12];
                string fullSequence = split[13];
                double monoMass = double.Parse(split[21]);
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[24].Split(new char[] { '|' });
                bool decoyPeptide = split[39].Equals("D");
                List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups, decoy: decoyPeptide);
                ids.Add(id);
            }

            return ids;
        }

        /// <summary>
        /// Rank-based AUC: the probability that a randomly chosen target score exceeds a randomly chosen
        /// decoy score (ties count as half). 0.5 means the two distributions are indistinguishable, 1.0
        /// means perfect separation with targets on top. Equivalent to the normalized Mann-Whitney U.
        /// </summary>
        private static double RankSumAuc(IReadOnlyList<double> targetScores, IReadOnlyList<double> decoyScores)
        {
            if (targetScores.Count == 0 || decoyScores.Count == 0)
                return double.NaN;

            // Pool the scores, sort ascending, and assign fractional (average) ranks to ties.
            var pooled = targetScores.Select(s => (score: s, isTarget: true))
                .Concat(decoyScores.Select(s => (score: s, isTarget: false)))
                .OrderBy(x => x.score)
                .ToList();

            double[] ranks = new double[pooled.Count];
            int i = 0;
            while (i < pooled.Count)
            {
                int j = i;
                while (j + 1 < pooled.Count && pooled[j + 1].score == pooled[i].score)
                    j++;
                double averageRank = ((i + 1) + (j + 1)) / 2.0; // ranks are 1-based
                for (int k = i; k <= j; k++)
                    ranks[k] = averageRank;
                i = j + 1;
            }

            double targetRankSum = 0;
            for (int k = 0; k < pooled.Count; k++)
                if (pooled[k].isTarget)
                    targetRankSum += ranks[k];

            int nTarget = targetScores.Count;
            int nDecoy = decoyScores.Count;
            double u = targetRankSum - (double)nTarget * (nTarget + 1) / 2.0;
            return u / ((double)nTarget * nDecoy);
        }

        private static void DeletePepModelArtifact()
        {
            string pepModelArtifact = Path.Combine(TestContext.CurrentContext.TestDirectory, "FlashLFQ", "TestData", "model.zip");
            if (File.Exists(pepModelArtifact)) File.Delete(pepModelArtifact);
        }
    }
}
