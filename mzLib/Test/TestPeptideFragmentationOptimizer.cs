using NUnit.Framework;
using PredictionClients.AI;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace Test
{
    [TestFixture]
    [Category("Koina")]
    public class TestPeptideFragmentationOptimizer
    {
        private static string OutputDir => Path.Combine(@"C:\Users\trish\Downloads", "FragmentationOptimizer");
        private static string ClaudeAnalysisPath => Path.Combine(OutputDir, "claude_analysis.json");

        [OneTimeSetUp]
        public void Setup()
        {
            if (!Directory.Exists(OutputDir)) Directory.CreateDirectory(OutputDir);
        }

        /// <summary>
        /// Linear peptide optimization - finds peptides that fragment evenly in a single orientation.
        /// </summary>
        [Test]
        [Explicit("Requires Koina API access")]
        public async Task OptimizeLinearFragmentation()
        {
            var optimizer = new PeptideFragmentationOptimizer(
                peptideLength: 10,
                populationSize: 30,
                collisionEnergy: 35,
                precursorCharge: 2,
                randomSeed: 42);

            // --- Load Claude's analysis if available ---
            List<string>? claudeSeeds = null;
            if (File.Exists(ClaudeAnalysisPath))
            {
                Console.WriteLine($"Claude analysis found at: {ClaudeAnalysisPath}");
                try
                {
                    var analysis = ClaudeAnalysisReader.LoadFromFile(ClaudeAnalysisPath);
                    claudeSeeds = optimizer.InjectClaudeAnalysis(analysis);

                    Console.WriteLine($"  Injected {claudeSeeds.Count} seed sequences from Claude.");
                    Console.WriteLine($"  Dipeptide adjustments: {analysis.DipeptideAdjustments.Count}");
                    Console.WriteLine($"  Terminal adjustments: {analysis.TerminalAdjustments.Count}");

                    if (analysis.ProtectedCores.Count > 0)
                        Console.WriteLine($"  Protected cores: {string.Join(", ", analysis.ProtectedCores)}");
                    if (analysis.ChampionSequences.Count > 0)
                        Console.WriteLine($"  Champions: {string.Join(", ", analysis.ChampionSequences)}");
                    if (analysis.PrematureConvergenceWarning)
                        Console.WriteLine("  WARNING: Premature convergence flagged.");
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"  Could not load analysis: {ex.Message}");
                }
            }
            else
            {
                Console.WriteLine("No Claude analysis — cold start.");
            }

            var best = await optimizer.OptimizeAsync(
                generations: 20,
                eliteRatio: 0.05,
                mutationRate: 0.2,
                claudeSeeds: claudeSeeds,
                onGenerationComplete: (gen, fitness, diversity) =>
                    Console.WriteLine($"Gen {gen}: {fitness.Sequence}  Fitness={fitness.OverallFitness:F4}  Diversity={diversity:P0}"));

            File.WriteAllText(Path.Combine(OutputDir, "linear_progress.txt"), optimizer.GetProgressReport());
            File.WriteAllText(Path.Combine(OutputDir, "linear_patterns.txt"), optimizer.GetLearnedPatternsSummary());
            File.WriteAllText(Path.Combine(OutputDir, "linear_claude_prompt.txt"), optimizer.GenerateClaudePrompt());

            Console.WriteLine($"\nBEST LINEAR: {best.Sequence} (Fitness: {best.OverallFitness:F4})");
            Console.WriteLine($"Results saved to: {OutputDir}");

            Assert.That(best.Sequence.Length, Is.EqualTo(10));
            Assert.That(best.Sequence.Distinct().Count(), Is.EqualTo(10));
        }

        /// <summary>
        /// Circular peptide optimization - finds peptides that fragment evenly 
        /// regardless of where the ring is opened (all rotations tested).
        /// </summary>
        [Test]
        [Explicit("Requires Koina API access - intensive (N rotations per peptide)")]
        public async Task OptimizeCircularFragmentation()
        {
            // Smaller population because each peptide = 10 API calls (one per rotation)
            var optimizer = new PeptideFragmentationOptimizer(
                peptideLength: 10,
                populationSize: 20,  // 20 peptides × 10 rotations = 200 API calls/generation
                collisionEnergy: 35,
                precursorCharge: 2,
                randomSeed: 42);

            // --- Load Claude's analysis if available ---
            List<string>? claudeSeeds = null;
            string circularAnalysisPath = Path.Combine(OutputDir, "circular_claude_analysis.json");
            if (File.Exists(circularAnalysisPath))
            {
                Console.WriteLine($"Claude circular analysis found at: {circularAnalysisPath}");
                try
                {
                    var analysis = ClaudeAnalysisReader.LoadFromFile(circularAnalysisPath);
                    claudeSeeds = optimizer.InjectClaudeAnalysis(analysis);
                    Console.WriteLine($"  Injected {claudeSeeds.Count} seed sequences.");
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"  Could not load: {ex.Message}");
                }
            }
            else
            {
                Console.WriteLine("No circular analysis — cold start.");
            }

            Console.WriteLine($"\nStarting CIRCULAR optimization...");
            Console.WriteLine($"Each peptide will be tested at all {optimizer.GetPeptideLength()} rotation points.\n");

            var best = await optimizer.OptimizeCircularAsync(
                generations: 15,
                eliteRatio: 0.1,
                mutationRate: 0.2,
                claudeSeeds: claudeSeeds,
                onGenerationComplete: (gen, fitness, diversity) =>
                {
                    Console.WriteLine($"Gen {gen}: {fitness.CanonicalSequence}  " +
                        $"CircFit={fitness.OverallCircularFitness:F4}  " +
                        $"Min={fitness.MinRotationFitness:F4}  " +
                        $"Worst@{fitness.WorstRotation.Sequence}  " +
                        $"Diversity={diversity:P0}");
                });

            // Save outputs
            File.WriteAllText(Path.Combine(OutputDir, "circular_progress.txt"),
                optimizer.GetCircularProgressReport());
            File.WriteAllText(Path.Combine(OutputDir, "circular_patterns.txt"),
                optimizer.GetLearnedPatternsSummary());
            File.WriteAllText(Path.Combine(OutputDir, "circular_claude_prompt.txt"),
                optimizer.GenerateCircularClaudePrompt());

            Console.WriteLine($"\n=== BEST CIRCULAR PEPTIDE ===");
            Console.WriteLine($"Sequence: {best.CanonicalSequence}");
            Console.WriteLine($"Circular Fitness: {best.OverallCircularFitness:F4}");
            Console.WriteLine($"Min Rotation Fitness: {best.MinRotationFitness:F4}");
            Console.WriteLine($"Avg Rotation Fitness: {best.AvgRotationFitness:F4}");
            Console.WriteLine($"Fitness StdDev: {best.FitnessStdDev:F4}");
            Console.WriteLine($"Worst Rotation: {best.WorstRotation.Sequence} (Fit={best.WorstRotation.OverallFitness:F4})");
            Console.WriteLine($"Best Rotation: {best.BestRotation.Sequence} (Fit={best.BestRotation.OverallFitness:F4})");

            Console.WriteLine($"\nAll rotations of best peptide:");
            foreach (var rot in best.RotationFitnesses.OrderByDescending(r => r.OverallFitness))
            {
                string marker = rot.Sequence == best.WorstRotation.Sequence ? " ? WORST" :
                               rot.Sequence == best.BestRotation.Sequence ? " ? BEST" : "";
                Console.WriteLine($"  {rot.Sequence}  Fit={rot.OverallFitness:F4}  B-CV={rot.BIonCV:F3}  Y-CV={rot.YIonCV:F3}{marker}");
            }

            Console.WriteLine($"\nResults saved to: {OutputDir}");
            Console.WriteLine($"Send circular_claude_prompt.txt to Claude for analysis.");

            Assert.That(best.CanonicalSequence.Length, Is.EqualTo(10));
            Assert.That(best.CanonicalSequence.Distinct().Count(), Is.EqualTo(10));
            Assert.That(best.RotationFitnesses.Count, Is.EqualTo(10)); // All rotations tested
        }

        /// <summary>
        /// Quick test of the rotation helper methods.
        /// </summary>
        [Test]
        public void TestRotationHelpers()
        {
            string peptide = "ABCDEFGHIJ";

            var rotations = PeptideFragmentationOptimizer.GetAllRotations(peptide);

            Assert.That(rotations.Count, Is.EqualTo(10));
            Assert.That(rotations[0], Is.EqualTo("ABCDEFGHIJ")); // Original
            Assert.That(rotations[1], Is.EqualTo("BCDEFGHIJA")); // Rotate by 1
            Assert.That(rotations[5], Is.EqualTo("FGHIJABCDE")); // Rotate by 5

            // All rotations should have same length and same characters
            foreach (var rot in rotations)
            {
                Assert.That(rot.Length, Is.EqualTo(10));
                Assert.That(rot.OrderBy(c => c).SequenceEqual(peptide.OrderBy(c => c)), Is.True);
            }

            // Test specific rotation
            Assert.That(PeptideFragmentationOptimizer.GetRotation(peptide, 3), Is.EqualTo("DEFGHIJABC"));

            Console.WriteLine("All rotations of ABCDEFGHIJ:");
            for (int i = 0; i < rotations.Count; i++)
                Console.WriteLine($"  Rotation {i}: {rotations[i]}");
        }

        /// <summary>
        /// Compares a single peptide's linear vs circular fitness.
        /// Useful for understanding how much variance exists across rotations.
        /// </summary>
        [Test]
        [Explicit("Requires Koina API access")]
        public async Task AnalyzeSinglePeptideCircularFitness()
        {
            string testPeptide = "FSAHRYTGNE"; // Replace with your peptide of interest

            var optimizer = new PeptideFragmentationOptimizer(
                peptideLength: testPeptide.Length,
                collisionEnergy: 35,
                precursorCharge: 2);

            Console.WriteLine($"Analyzing circular fitness for: {testPeptide}\n");

            // Evaluate all rotations
            var rotations = PeptideFragmentationOptimizer.GetAllRotations(testPeptide);
            var fitnessResults = await optimizer.EvaluatePeptidesAsync(rotations, generation: 0);

            Console.WriteLine("| Rotation   | Fitness | B-CV  | Y-CV  | MinInt | Evenness |");
            Console.WriteLine("|------------|---------|-------|-------|--------|----------|");

            foreach (var (rotation, fitness) in rotations.Zip(fitnessResults))
            {
                Console.WriteLine($"| {rotation} | {fitness.OverallFitness:F4}  | {fitness.BIonCV:F3} | {fitness.YIonCV:F3} | {fitness.MinIntensity:F4} | {fitness.EvennessScore:F4}   |");
            }

            double minFit = fitnessResults.Min(f => f.OverallFitness);
            double avgFit = fitnessResults.Average(f => f.OverallFitness);
            double maxFit = fitnessResults.Max(f => f.OverallFitness);
            double stdDev = Math.Sqrt(fitnessResults.Select(f => Math.Pow(f.OverallFitness - avgFit, 2)).Average());

            Console.WriteLine($"\nSummary:");
            Console.WriteLine($"  Min fitness: {minFit:F4} @ {fitnessResults.OrderBy(f => f.OverallFitness).First().Sequence}");
            Console.WriteLine($"  Avg fitness: {avgFit:F4}");
            Console.WriteLine($"  Max fitness: {maxFit:F4} @ {fitnessResults.OrderByDescending(f => f.OverallFitness).First().Sequence}");
            Console.WriteLine($"  Std dev: {stdDev:F4}");
            Console.WriteLine($"  Circular fitness (0.6×min + 0.3×avg + 0.1×max): {0.6 * minFit + 0.3 * avgFit + 0.1 * maxFit:F4}");

            Assert.That(fitnessResults.Count, Is.EqualTo(testPeptide.Length));
        }
    }
}