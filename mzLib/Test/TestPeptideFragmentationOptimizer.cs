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
        private static string CircularAnalysisPath => Path.Combine(OutputDir, "circular_claude_analysis.json");

        [OneTimeSetUp]
        public void Setup()
        {
            if (!Directory.Exists(OutputDir)) Directory.CreateDirectory(OutputDir);
        }

        [Test]
        [Explicit("Requires Koina API access - intensive (N rotations per peptide)")]
        public async Task OptimizeCircularFragmentation()
        {
            var optimizer = new PeptideFragmentationOptimizer(
                peptideLength: 10,
                populationSize: 20,
                collisionEnergy: 35,
                precursorCharge: 2,
                randomSeed: 42);

            List<string>? claudeSeeds = null;
            if (File.Exists(CircularAnalysisPath))
            {
                Console.WriteLine($"Loading Claude analysis from: {CircularAnalysisPath}");
                try
                {
                    var analysis = ClaudeAnalysisReader.LoadFromFile(CircularAnalysisPath);
                    claudeSeeds = optimizer.InjectClaudeAnalysis(analysis);
                    Console.WriteLine($"  Injected {claudeSeeds.Count} seeds");
                }
                catch (Exception ex)
                {
                    Console.WriteLine($"  Failed to load: {ex.Message}");
                }
            }
            else
            {
                Console.WriteLine("No analysis file - cold start");
            }

            Console.WriteLine($"\nOptimizing circular peptides ({optimizer.GetPeptideLength()} rotations each)...\n");

            var best = await optimizer.OptimizeCircularAsync(
                generations: 15,
                eliteRatio: 0.1,
                mutationRate: 0.2,
                claudeSeeds: claudeSeeds,
                onGenerationComplete: (gen, fitness, diversity) =>
                {
                    string canonical = PeptideFragmentationOptimizer.GetCanonicalRotation(fitness.CanonicalSequence);
                    Console.WriteLine($"Gen {gen}: {canonical}  CircFit={fitness.OverallCircularFitness:F4}  Min={fitness.MinRotationFitness:F4}  Diversity={diversity:P0}");
                });

            File.WriteAllText(Path.Combine(OutputDir, "circular_progress.txt"), optimizer.GetCircularProgressReport());
            File.WriteAllText(Path.Combine(OutputDir, "circular_patterns.txt"), optimizer.GetLearnedPatternsSummary());
            File.WriteAllText(Path.Combine(OutputDir, "circular_claude_prompt.txt"), optimizer.GenerateCircularClaudePrompt());

            string bestCanonical = PeptideFragmentationOptimizer.GetCanonicalRotation(best.CanonicalSequence);

            Console.WriteLine($"\n=== BEST CIRCULAR PEPTIDE ===");
            Console.WriteLine($"Canonical: {bestCanonical}");
            Console.WriteLine($"CircularFitness: {best.OverallCircularFitness:F4}");
            Console.WriteLine($"Min: {best.MinRotationFitness:F4} Avg: {best.AvgRotationFitness:F4} Max: {best.MaxRotationFitness:F4}");
            Console.WriteLine($"Worst rotation: {best.WorstRotation.Sequence} (Fit={best.WorstRotation.OverallFitness:F4})");
            Console.WriteLine($"Best rotation:  {best.BestRotation.Sequence} (Fit={best.BestRotation.OverallFitness:F4})");

            Console.WriteLine($"\nAll rotations:");
            foreach (var rot in best.RotationFitnesses.OrderByDescending(r => r.OverallFitness))
            {
                string marker = rot.Sequence == best.WorstRotation.Sequence ? " <- WORST" :
                               rot.Sequence == best.BestRotation.Sequence ? " <- BEST" : "";
                Console.WriteLine($"  {rot.Sequence} Fit={rot.OverallFitness:F4} B-CV={rot.BIonCV:F3} Y-CV={rot.YIonCV:F3}{marker}");
            }

            Console.WriteLine($"\nOutputs saved to: {OutputDir}");

            Assert.That(best.CanonicalSequence.Length, Is.EqualTo(10));
            Assert.That(best.CanonicalSequence.Distinct().Count(), Is.EqualTo(10));
            Assert.That(best.RotationFitnesses.Count, Is.EqualTo(10));
        }

        [Test]
        public void TestCanonicalRotation()
        {
            Assert.That(PeptideFragmentationOptimizer.GetCanonicalRotation("BCDA"), Is.EqualTo("ABCD"));
            Assert.That(PeptideFragmentationOptimizer.GetCanonicalRotation("CDAB"), Is.EqualTo("ABCD"));
            Assert.That(PeptideFragmentationOptimizer.GetCanonicalRotation("DABC"), Is.EqualTo("ABCD"));
            Assert.That(PeptideFragmentationOptimizer.GetCanonicalRotation("ABCD"), Is.EqualTo("ABCD"));

            var rotations = PeptideFragmentationOptimizer.GetAllRotations("ABCD");
            Assert.That(rotations.Count, Is.EqualTo(4));
            Assert.That(rotations, Does.Contain("ABCD"));
            Assert.That(rotations, Does.Contain("BCDA"));
            Assert.That(rotations, Does.Contain("CDAB"));
            Assert.That(rotations, Does.Contain("DABC"));
        }
    }
}