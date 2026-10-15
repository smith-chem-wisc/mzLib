using FlashLFQ;
using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using ChromatographicPeak = FlashLFQ.ChromatographicPeak;
using IsotopicEnvelope = FlashLFQ.IsotopicEnvelope;
using Peptide = FlashLFQ.Peptide;

namespace Test.FlashLFQ
{
    /// <summary>
    /// A design may skip a replicate number, when a sample was lost, and the numbers are the user's: they are never
    /// closed up, since replicate 4 of one condition can be the same subject as replicate 4 of another. Normalization
    /// and the Bayesian protein step must therefore give, on a design with a gap, exactly what they give on the same
    /// data numbered without one. Every test here runs both and compares them.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class TestReplicateGaps
    {
        private const int NumPeptides = 12;

        /// <summary>A peptide's intensity in a file is its own abundance times the file's loading factor.</summary>
        private static double Abundance(int peptide) => 1000.0 * (peptide + 1) * (1 + 0.1 * (peptide % 3));

        /// <summary>
        /// One synthetic peak per peptide per file, with intensity <see cref="Abundance"/> times the file's loading,
        /// so normalization has a known answer: every file's loading divided out. A peptide listed in
        /// <paramref name="absent"/> for a file has no peak there.
        /// </summary>
        private static FlashLfqResults Build(List<(SpectraFileInfo File, double Loading)> files,
            Func<SpectraFileInfo, int, bool> absent = null)
        {
            var ids = new List<Identification>();
            var peaks = new List<ChromatographicPeak>();
            foreach (var (file, loading) in files)
            {
                for (int p = 0; p < NumPeptides; p++)
                {
                    if (absent != null && absent(file, p))
                    {
                        continue;
                    }
                    string sequence = "PEPTIDE" + (char)('A' + p);
                    var id = new Identification(file, sequence, sequence, 1000 + p, 10, 2,
                        new List<ProteinGroup> { new ProteinGroup("Protein", "Gene", "Organism") });
                    ids.Add(id);

                    var peak = new ChromatographicPeak(id, file);
                    peak.ResolveIdentifications();
                    peak.IsotopicEnvelopes.Add(new IsotopicEnvelope(new IndexedMassSpectralPeak(0, 0, 0, 0), 2, Abundance(p) * loading, 1));
                    peak.CalculateIntensityForThisFeature(false);
                    peaks.Add(peak);
                }
            }

            var results = new FlashLfqResults(files.Select(f => f.File).ToList(), ids);
            foreach (var peak in peaks)
            {
                results.Peaks[peak.SpectraFileInfo].Add(peak);
            }
            results.CalculatePeptideResults(quantifyAmbiguousPeptides: false);
            return results;
        }

        private static FlashLfqResults Normalize(FlashLfqResults results, bool silent = true)
        {
            new IntensityNormalizationEngine(results, integrate: false, silent: silent, maxThreads: 1).NormalizeResults();
            return results;
        }

        /// <summary>Every peptide's intensity in every file, files in the order given.</summary>
        private static double[][] Intensities(FlashLfqResults results) =>
            results.SpectraFiles.Select(file =>
                Enumerable.Range(0, NumPeptides)
                    .Select(p => results.PeptideModifiedSequences.TryGetValue("PEPTIDE" + (char)('A' + p), out var peptide) ? peptide.GetIntensity(file) : 0)
                    .ToArray())
                .ToArray();

        /// <summary>Every file's intensities equal the reference file's: its loading has been divided out.</summary>
        private static void AssertEveryFileMatches(double[][] intensities, double[] reference, string because)
        {
            foreach (var file in intensities)
            {
                Assert.That(file, Is.EqualTo(reference).Within(1e-9).Percent, because);
            }
        }

        private static SpectraFileInfo File(string name, string condition, int biorep, int fraction = 0, int techrep = 0) =>
            new SpectraFileInfo(name, condition, biorep, techrep, fraction);

        [Test]
        public static void BiorepNormalizationOnAGappedDesignMatchesTheSameDesignNumberedWithoutTheGap()
        {
            // condition a lost biorep 2 (index 1); condition b lost biorep 1 (index 0)
            var gapped = Normalize(Build(new()
            {
                (File("a1", "a", 0), 1.0), (File("a3", "a", 2), 2.0),
                (File("b2", "b", 1), 0.5), (File("b3", "b", 2), 3.0),
            }));
            var dense = Normalize(Build(new()
            {
                (File("a1", "a", 0), 1.0), (File("a3", "a", 1), 2.0),
                (File("b2", "b", 0), 0.5), (File("b3", "b", 1), 3.0),
            }));

            Assert.That(Intensities(gapped), Is.EqualTo(Intensities(dense)));
            AssertEveryFileMatches(Intensities(gapped), Intensities(gapped)[0], "every file's loading is divided out");
        }

        /// <summary>
        /// With the first condition's biorep 1 lost there was no reference at all: biorep normalization quietly
        /// normalized nothing. The reference is now the lowest biorep that exists.
        /// </summary>
        [Test]
        public static void LosingTheReferenceBiorepStillNormalizesAgainstTheLowestBiorepLeft()
        {
            var results = Normalize(Build(new()
            {
                (File("a2", "a", 1), 1.0), (File("a3", "a", 2), 4.0),
                (File("b1", "b", 0), 0.25), (File("b2", "b", 1), 2.0),
            }));

            AssertEveryFileMatches(Intensities(results), Intensities(results)[0], "normalized against a2, the lowest biorep left");
        }

        /// <summary>
        /// A biorep that shares no peptide with the reference cannot be normalized, but it used to stop biorep
        /// normalization for every biorep after it too. Now only that biorep is left as it is.
        /// </summary>
        [Test]
        public static void ABiorepWithNoPeptideInCommonWithTheReferenceDoesNotStopTheRest()
        {
            var lonely = File("a2", "a", 1);
            List<(SpectraFileInfo, double)> design = new()
            {
                (File("a1", "a", 0), 1.0), (lonely, 5.0), (File("a3", "a", 2), 2.0), (File("b1", "b", 0), 3.0),
            };
            Func<SpectraFileInfo, int, bool> absent = (file, p) => file == lonely ? p < NumPeptides / 2 : p >= NumPeptides / 2;
            var before = Intensities(Build(design, absent));
            var results = Normalize(Build(design, absent));

            var intensities = Intensities(results);
            AssertEveryFileMatches(new[] { intensities[2], intensities[3] }, intensities[0],
                "a3, after the lonely biorep, and b1, in the next condition, are normalized");
            Assert.That(intensities[1], Is.EqualTo(before[1]), "the lonely biorep is left as it was");
        }

        [Test]
        public static void FractionNormalizationOnAGappedDesignMatchesTheSameDesignNumberedWithoutTheGap()
        {
            // two fractions per biorep; condition a lost biorep 2 (index 1). This used to throw.
            List<(SpectraFileInfo, double)> Design(int secondBiorep) => new()
            {
                (File("a1f1", "a", 0, 0), 1.0), (File("a1f2", "a", 0, 1), 1.0),
                (File("aXf1", "a", secondBiorep, 0), 2.0), (File("aXf2", "a", secondBiorep, 1), 0.5),
                (File("b1f1", "b", 0, 0), 3.0), (File("b1f2", "b", 0, 1), 1.5),
            };

            var gapped = Normalize(Build(Design(2)));
            var dense = Normalize(Build(Design(1)));

            Assert.That(Intensities(gapped), Is.EqualTo(Intensities(dense)));
        }

        /// <summary>A fraction whose first technical replicate was lost used to throw in biorep normalization.</summary>
        [Test]
        public static void ALostFirstTechrepMatchesTheSameDesignNumberedWithoutTheGap()
        {
            List<(SpectraFileInfo, double)> Design(int techrep) => new()
            {
                (File("a1", "a", 0), 1.0),
                (File("a2t", "a", 1, 0, techrep), 2.0),
                (File("b1", "b", 0), 0.5),
            };

            Assert.That(Intensities(Normalize(Build(Design(1)))), Is.EqualTo(Intensities(Normalize(Build(Design(0))))));
        }

        /// <summary>
        /// The first technical replicate is the lowest techrep number, not the first file listed: the same files listed
        /// in another order normalize the same way, and the warning that names "the first technical replicate" is true.
        /// </summary>
        [Test]
        public static void TechrepNormalizationDoesNotDependOnTheOrderFilesAreListed()
        {
            var inOrder = new List<(SpectraFileInfo, double)>
            {
                (File("a1t1", "a", 0, 0, 0), 1.0), (File("a1t2", "a", 0, 0, 1), 2.0), (File("b1", "b", 0), 0.5),
            };
            var reversed = new List<(SpectraFileInfo, double)> { inOrder[1], inOrder[0], inOrder[2] };

            Dictionary<string, double[]> ByName(FlashLfqResults results) => results.SpectraFiles
                .Zip(Intensities(results), (file, intensities) => (file.FullFilePathWithExtension, intensities))
                .ToDictionary(x => x.FullFilePathWithExtension, x => x.intensities);

            var expected = ByName(Normalize(Build(inOrder)));
            var actual = ByName(Normalize(Build(reversed)));

            foreach (var name in expected.Keys)
            {
                Assert.That(actual[name], Is.EqualTo(expected[name]).Within(1e-9).Percent, name);
            }
        }

        /// <summary>
        /// A technical replicate that shares no peptide with the first one cannot be normalized. It is left as it is,
        /// with a warning naming it, and the technical replicates after it are still normalized.
        /// </summary>
        [Test]
        public static void ATechrepWithNoPeptideInCommonWithTheFirstIsWarnedAboutAndDoesNotStopTheRest()
        {
            var lonely = File("a1t2", "a", 0, 0, 1);
            List<(SpectraFileInfo, double)> design = new()
            {
                (File("a1t1", "a", 0, 0, 0), 1.0), (lonely, 5.0), (File("a1t3", "a", 0, 0, 2), 2.0), (File("b1", "b", 0), 0.5),
            };
            Func<SpectraFileInfo, int, bool> absent = (file, p) => file == lonely ? p < NumPeptides / 2 : file.Condition == "a" && p >= NumPeptides / 2;
            var before = Intensities(Build(design, absent));

            var console = Console.Out;
            var output = new System.IO.StringWriter();
            FlashLfqResults results;
            try
            {
                Console.SetOut(output);
                results = Normalize(Build(design, absent), silent: false);
            }
            finally
            {
                Console.SetOut(console);
            }

            var intensities = Intensities(results);
            Assert.That(intensities[1], Is.EqualTo(before[1]), "the lonely techrep is left as it was");
            Assert.That(intensities[2].Take(NumPeptides / 2), Is.EqualTo(intensities[0].Take(NumPeptides / 2)).Within(1e-9).Percent,
                "a1t3, after the lonely techrep, is normalized to the first");
            Assert.That(output.ToString(), Does.Contain("Warning: technical replicate 2 of condition \"a\" biorep 1 shares no peptides with the first technical replicate, so it was not normalized"));
        }

        /// <summary>
        /// Two peptides of one protein, three bioreps per condition, condition b about twice condition a. Returns the
        /// Bayesian result for b against a, and its text as written to the fold-change table.
        /// </summary>
        private static (UnpairedProteinQuantResult Result, string Text) Bayesian(int[] aBioreps, int[] bBioreps)
        {
            var files = aBioreps.Select((b, i) => File("a" + i, "a", b))
                .Concat(bBioreps.Select((b, i) => File("b" + i, "b", b))).ToList();
            var ids = new List<Identification>
            {
                new(null, "PEPTIDEA", "PEPTIDEA", 0, 0, 0, new List<ProteinGroup> { new("Protein", "Gene", "Organism") }),
                new(null, "PEPTIDEB", "PEPTIDEB", 0, 0, 0, new List<ProteinGroup> { new("Protein", "Gene", "Organism") }),
            };
            var results = new FlashLfqResults(files, ids);

            double[] aIntensities = { 900, 1000, 1100 };
            double[] bIntensities = { 1950, 2000, 2050 };
            for (int i = 0; i < aBioreps.Length; i++)
            {
                results.PeptideModifiedSequences["PEPTIDEA"].SetIntensity(files[i], aIntensities[i]);
                results.PeptideModifiedSequences["PEPTIDEB"].SetIntensity(files[i], 3 * aIntensities[i]);
            }
            for (int i = 0; i < bBioreps.Length; i++)
            {
                results.PeptideModifiedSequences["PEPTIDEA"].SetIntensity(files[aBioreps.Length + i], bIntensities[i]);
                results.PeptideModifiedSequences["PEPTIDEB"].SetIntensity(files[aBioreps.Length + i], 3 * bIntensities[i]);
            }

            new ProteinQuantificationEngine(results, maxThreads: 1, controlCondition: "a", randomSeed: 0, foldChangeCutoff: 0.1).Run();
            var result = (UnpairedProteinQuantResult)results.ProteinGroups["Protein"].ConditionToQuantificationResults["b"];
            return (result, result.ToString());
        }

        [Test]
        public static void TheBayesianStepOnAGappedDesignMatchesTheSameDesignNumberedWithoutTheGap()
        {
            var gapped = Bayesian(new[] { 0, 1, 3 }, new[] { 0, 2, 3 });
            var dense = Bayesian(new[] { 0, 1, 2 }, new[] { 0, 1, 2 });

            Assert.That(gapped.Text, Is.EqualTo(dense.Text));
            Assert.That(gapped.Result.FoldChangePointEstimate, Is.EqualTo(dense.Result.FoldChangePointEstimate));
            Assert.That(gapped.Result.PosteriorErrorProbability, Is.EqualTo(dense.Result.PosteriorErrorProbability));
            Assert.That(Math.Round(gapped.Result.FoldChangePointEstimate, 1), Is.EqualTo(1.0), "b is about twice a");
        }

        /// <summary>A control condition with one biorep takes its own path, which looped 0..max over the bioreps too.</summary>
        [Test]
        public static void AControlWithOneBiorepThatIsNotBiorep1MatchesOneThatIs()
        {
            var gapped = Bayesian(new[] { 2 }, new[] { 0, 1, 2 });
            var dense = Bayesian(new[] { 0 }, new[] { 0, 1, 2 });

            Assert.That(gapped.Text, Is.EqualTo(dense.Text));
        }
    }
}
