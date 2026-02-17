using System;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using NUnit.Framework;
using PredictionClients.AI;

namespace Test
{
    [TestFixture]
    [Category("Koina")]
    public class TestPeptideFragmentationOptimizer
    {
        private static string OutputDir => Path.Combine(@"C:\Users\trish\Downloads", "FragmentationOptimizer");

        [OneTimeSetUp]
        public void Setup()
        {
            if (!Directory.Exists(OutputDir)) Directory.CreateDirectory(OutputDir);
        }

        [Test]
        [Explicit("Requires Koina API access")]
        public async Task OptimizeEvenFragmentation()
        {
            var optimizer = new PeptideFragmentationOptimizer(
                peptideLength: 12,
                populationSize: 30,
                collisionEnergy: 25,
                precursorCharge: 2,
                randomSeed: 42);

            var best = await optimizer.OptimizeAsync(
                generations: 5,
                onGenerationComplete: (gen, fitness) =>
                    Console.WriteLine($"Gen {gen}: {fitness.Sequence} Fitness={fitness.OverallFitness:F4}"));

            File.WriteAllText(Path.Combine(OutputDir, "progress.txt"), optimizer.GetProgressReport());
            File.WriteAllText(Path.Combine(OutputDir, "patterns.txt"), optimizer.GetLearnedPatternsSummary());
            File.WriteAllText(Path.Combine(OutputDir, "claude_prompt.txt"), optimizer.GenerateClaudePrompt());

            Console.WriteLine($"\nBEST: {best.Sequence} (Fitness: {best.OverallFitness:F4})");
            Console.WriteLine($"\nResults saved to: {OutputDir}");

            Assert.That(best.Sequence.Length, Is.EqualTo(12));
            Assert.That(best.Sequence.Distinct().Count(), Is.EqualTo(12));
        }
    }
}