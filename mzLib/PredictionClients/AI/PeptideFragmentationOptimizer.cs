using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;

namespace PredictionClients.AI
{
    /// <summary>
    /// Fitness metrics for evaluating peptide fragmentation quality.
    /// </summary>
    public class FragmentationFitness
    {
        public string Sequence { get; set; } = "";
        public double[] BIonIntensities { get; set; } = Array.Empty<double>();
        public double[] YIonIntensities { get; set; } = Array.Empty<double>();
        public double BIonCV { get; set; }
        public double YIonCV { get; set; }
        public double EvennessScore { get; set; }
        public double MinIntensity { get; set; }
        public double OverallFitness { get; set; }
        public int Generation { get; set; }

        public override string ToString() =>
            $"{Sequence}: Fitness={OverallFitness:F4}, Evenness={EvennessScore:F4}, MinInt={MinIntensity:F4}";
    }

    /// <summary>
    /// Optimizes peptide sequences for even HCD fragmentation using Prosit predictions
    /// and an evolutionary algorithm that learns favorable dipeptide patterns.
    /// </summary>
    public class PeptideFragmentationOptimizer
    {
        //private static readonly char[] AvailableAminoAcids =
        //    { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' };
        private static readonly char[] AvailableAminoAcids =
            { 'A', 'D', 'E', 'F', 'G', 'H', 'K', 'L', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' };
        private readonly int _peptideLength;
        private readonly int _populationSize;
        private readonly int _collisionEnergy;
        private readonly int _precursorCharge;
        private readonly Random _random;

        public List<FragmentationFitness> AllEvaluatedPeptides { get; } = new();
        public List<FragmentationFitness> BestPerGeneration { get; } = new();
        public Dictionary<string, double> LearnedPatterns { get; } = new();

        public PeptideFragmentationOptimizer(
            int peptideLength = 10,
            int populationSize = 1000,
            int collisionEnergy = 35,
            int precursorCharge = 2,
            int? randomSeed = null)
        {
            _peptideLength = peptideLength;
            _populationSize = populationSize;
            _collisionEnergy = collisionEnergy;
            _precursorCharge = precursorCharge;
            _random = randomSeed.HasValue ? new Random(randomSeed.Value) : new Random();
        }

        public string GenerateRandomPeptide()
        {
            var available = AvailableAminoAcids.ToList();
            var sequence = new char[_peptideLength];
            for (int i = 0; i < _peptideLength; i++)
            {
                int index = _random.Next(available.Count);
                sequence[i] = available[index];
                available.RemoveAt(index);
            }
            return new string(sequence);
        }

        public List<string> GenerateInitialPopulation()
        {
            var population = new HashSet<string>();
            while (population.Count < _populationSize)
                population.Add(GenerateRandomPeptide());
            return population.ToList();
        }

        public async Task<List<FragmentationFitness>> EvaluatePeptidesAsync(List<string> peptides, int generation)
        {
            var charges = Enumerable.Repeat(_precursorCharge, peptides.Count).ToList();
            var energies = Enumerable.Repeat(_collisionEnergy, peptides.Count).ToList();
            var retentionTimes = Enumerable.Repeat<double?>(null, peptides.Count).ToList();

            var model = new Prosit2020IntensityHCD(peptides, charges, energies, retentionTimes, out var warnings);

            if (model.PeptideSequences.Count == 0)
                throw new InvalidOperationException($"No valid peptides. Warnings: {warnings?.Message}");

            await model.RunInferenceAsync();

            var fitnessResults = new List<FragmentationFitness>();
            for (int i = 0; i < model.Predictions.Count; i++)
            {
                var fitness = CalculateFitness(peptides[i], model.Predictions[i], generation);
                fitnessResults.Add(fitness);
                AllEvaluatedPeptides.Add(fitness);
            }
            return fitnessResults;
        }

        private FragmentationFitness CalculateFitness(string sequence, PeptideFragmentIntensityPrediction prediction, int generation)
        {
            var bIons = new List<double>();
            var yIons = new List<double>();

            for (int i = 0; i < prediction.FragmentAnnotations.Count; i++)
            {
                string annotation = prediction.FragmentAnnotations[i];
                double intensity = prediction.FragmentIntensities[i];
                if (intensity < 0) continue;

                if (annotation.StartsWith("b")) bIons.Add(intensity);
                else if (annotation.StartsWith("y")) yIons.Add(intensity);
            }

            double bCV = CalculateCV(bIons);
            double yCV = CalculateCV(yIons);
            double evennessScore = 2.0 / (1.0 + bCV + yCV);

            double minIntensity = 0;
            if (bIons.Count > 0 && yIons.Count > 0) minIntensity = Math.Min(bIons.Min(), yIons.Min());
            else if (bIons.Count > 0) minIntensity = bIons.Min();
            else if (yIons.Count > 0) minIntensity = yIons.Min();

            var allIons = bIons.Concat(yIons).ToList();
            double avgIntensity = allIons.Count > 0 ? allIons.Average() : 0;
            double overallFitness = 0.5 * evennessScore + 0.3 * (minIntensity * 10) + 0.2 * avgIntensity;

            return new FragmentationFitness
            {
                Sequence = sequence,
                BIonIntensities = bIons.ToArray(),
                YIonIntensities = yIons.ToArray(),
                BIonCV = bCV, YIonCV = yCV,
                EvennessScore = evennessScore,
                MinIntensity = minIntensity,
                OverallFitness = overallFitness,
                Generation = generation
            };
        }

        private double CalculateCV(List<double> values)
        {
            if (values.Count == 0) return 1.0;
            double mean = values.Average();
            if (mean == 0) return 1.0;
            double std = Math.Sqrt(values.Select(v => Math.Pow(v - mean, 2)).Average());
            return std / mean;
        }

        public string Crossover(string parent1, string parent2)
        {
            int crossPoint = _random.Next(1, _peptideLength - 1);
            var child = new char[_peptideLength];
            var usedAAs = new HashSet<char>();

            for (int i = 0; i < crossPoint; i++) { child[i] = parent1[i]; usedAAs.Add(parent1[i]); }

            int childPos = crossPoint;
            foreach (char aa in parent2)
                if (childPos < _peptideLength && !usedAAs.Contains(aa)) { child[childPos++] = aa; usedAAs.Add(aa); }

            foreach (char aa in AvailableAminoAcids)
                if (childPos < _peptideLength && !usedAAs.Contains(aa)) { child[childPos++] = aa; usedAAs.Add(aa); }

            return new string(child);
        }

        public string Mutate(string peptide, double mutationRate = 0.1)
        {
            var sequence = peptide.ToCharArray();
            var usedAAs = new HashSet<char>(sequence);
            var unusedAAs = AvailableAminoAcids.Where(aa => !usedAAs.Contains(aa)).ToList();

            for (int i = 0; i < _peptideLength; i++)
            {
                if (_random.NextDouble() < mutationRate)
                {
                    if (_random.NextDouble() < 0.5 || unusedAAs.Count == 0)
                    {
                        int j = _random.Next(_peptideLength);
                        (sequence[i], sequence[j]) = (sequence[j], sequence[i]);
                    }
                    else
                    {
                        int idx = _random.Next(unusedAAs.Count);
                        char old = sequence[i];
                        sequence[i] = unusedAAs[idx];
                        unusedAAs.RemoveAt(idx);
                        unusedAAs.Add(old);
                    }
                }
            }
            return new string(sequence);
        }

        public async Task<FragmentationFitness> OptimizeAsync(
            int generations = 50, double eliteRatio = 0.05, double mutationRate = 0.2,
            Action<int, FragmentationFitness>? onGenerationComplete = null)
        {
            var population = GenerateInitialPopulation();
            FragmentationFitness? bestEver = null;

            for (int gen = 0; gen < generations; gen++)
            {
                var fitnessResults = (await EvaluatePeptidesAsync(population, gen))
                    .OrderByDescending(f => f.OverallFitness).ToList();

                var bestThisGen = fitnessResults.First();
                BestPerGeneration.Add(bestThisGen);
                if (bestEver == null || bestThisGen.OverallFitness > bestEver.OverallFitness)
                    bestEver = bestThisGen;

                LearnPatterns(fitnessResults);
                onGenerationComplete?.Invoke(gen, bestThisGen);

                var nextPopulation = new List<string>();
                int eliteCount = (int)(_populationSize * eliteRatio);
                foreach (var elite in fitnessResults.Take(eliteCount))
                    nextPopulation.Add(elite.Sequence);

                while (nextPopulation.Count < _populationSize)
                {
                    var p1 = TournamentSelect(fitnessResults);
                    var p2 = TournamentSelect(fitnessResults);
                    var child = Mutate(Crossover(p1.Sequence, p2.Sequence), mutationRate);
                    if (_random.NextDouble() < 0.3) child = ApplyLearnedPatterns(child);
                    nextPopulation.Add(child);
                }
                population = nextPopulation;
            }
            return bestEver!;
        }

        private FragmentationFitness TournamentSelect(List<FragmentationFitness> pop, int size = 3)
        {
            FragmentationFitness? best = null;
            for (int i = 0; i < size; i++)
            {
                var c = pop[_random.Next(pop.Count)];
                if (best == null || c.OverallFitness > best.OverallFitness) best = c;
            }
            return best!;
        }

        private void LearnPatterns(List<FragmentationFitness> results)
        {
            var top = results.Take(10).Select(f => f.Sequence).ToList();
            var bottom = results.TakeLast(10).Select(f => f.Sequence).ToList();

            foreach (var pep in top)
            {
                for (int i = 0; i < pep.Length - 1; i++)
                {
                    string di = pep.Substring(i, 2);
                    LearnedPatterns.TryAdd(di, 0); LearnedPatterns[di] += 0.1;
                }
                LearnedPatterns.TryAdd($"N-{pep[0]}", 0); LearnedPatterns[$"N-{pep[0]}"] += 0.05;
                LearnedPatterns.TryAdd($"C-{pep[^1]}", 0); LearnedPatterns[$"C-{pep[^1]}"] += 0.05;
            }
            foreach (var pep in bottom)
                for (int i = 0; i < pep.Length - 1; i++)
                { string di = pep.Substring(i, 2); LearnedPatterns.TryAdd(di, 0); LearnedPatterns[di] -= 0.1; }
        }

        private string ApplyLearnedPatterns(string peptide)
        {
            var seq = peptide.ToCharArray();
            var good = LearnedPatterns.Where(kv => kv.Value > 0.3 && kv.Key.Length == 2)
                .OrderByDescending(kv => kv.Value).Take(3).Select(kv => kv.Key);

            foreach (var di in good)
            {
                if (_random.NextDouble() < 0.3)
                {
                    int p1 = Array.IndexOf(seq, di[0]), p2 = Array.IndexOf(seq, di[1]);
                    if (p1 >= 0 && p2 >= 0 && p1 != p2 - 1 && p1 + 1 < _peptideLength && p1 + 1 != p2)
                        (seq[p2], seq[p1 + 1]) = (seq[p1 + 1], seq[p2]);
                }
            }
            return new string(seq);
        }

        public string GetProgressReport()
        {
            var sb = new StringBuilder();
            sb.AppendLine($"=== OPTIMIZATION PROGRESS ===\nPeptides evaluated: {AllEvaluatedPeptides.Count}\nGenerations: {BestPerGeneration.Count}\n");
            sb.AppendLine("Best per generation:");
            for (int i = 0; i < BestPerGeneration.Count; i++)
                sb.AppendLine($"  Gen {i,3}: {BestPerGeneration[i].Sequence} - Fitness: {BestPerGeneration[i].OverallFitness:F4}");
            sb.AppendLine("\nTop 5 overall:");
            foreach (var p in AllEvaluatedPeptides.OrderByDescending(f => f.OverallFitness).Take(5))
                sb.AppendLine($"  {p}");
            return sb.ToString();
        }

        public string GetLearnedPatternsSummary()
        {
            var sb = new StringBuilder();
            sb.AppendLine("=== LEARNED PATTERNS ===\n\nFavorable dipeptides:");
            foreach (var p in LearnedPatterns.Where(kv => kv.Key.Length == 2).OrderByDescending(kv => kv.Value).Take(10))
                sb.AppendLine($"  {p.Key}: +{p.Value:F2}");
            sb.AppendLine("\nUnfavorable dipeptides:");
            foreach (var p in LearnedPatterns.Where(kv => kv.Key.Length == 2).OrderBy(kv => kv.Value).Take(10))
                sb.AppendLine($"  {p.Key}: {p.Value:F2}");
            return sb.ToString();
        }

        public string GenerateClaudePrompt()
        {
            var sb = new StringBuilder();
            sb.AppendLine("# PEPTIDE FRAGMENTATION OPTIMIZATION RESULTS");
            sb.AppendLine();
            sb.AppendLine("## GOAL");
            sb.AppendLine($"Design {_peptideLength}-AA peptides that fragment EVENLY across all backbone positions in HCD.");
            sb.AppendLine("Even fragmentation means: low coefficient of variation (CV) for both b and y ion series,");
            sb.AppendLine("and no 'dead spots' (positions with very low intensity).");
            sb.AppendLine();
            sb.AppendLine("## CONSTRAINTS");
            sb.AppendLine($"- Exactly {_peptideLength} amino acids");
            sb.AppendLine("- Each amino acid used only ONCE (no repeats)");
            sb.AppendLine($"- Available AAs: {string.Join(", ", AvailableAminoAcids)}");
            sb.AppendLine($"- Collision energy: {_collisionEnergy} eV, Charge: +{_precursorCharge}");
            sb.AppendLine();

            sb.AppendLine("## OPTIMIZATION STATISTICS");
            sb.AppendLine($"- Total peptides evaluated: {AllEvaluatedPeptides.Count:N0}");
            sb.AppendLine($"- Generations completed: {BestPerGeneration.Count}");
            sb.AppendLine($"- Fitness range: {AllEvaluatedPeptides.Min(p => p.OverallFitness):F4} to {AllEvaluatedPeptides.Max(p => p.OverallFitness):F4}");
            sb.AppendLine();

            sb.AppendLine("## TOP 10 BEST FRAGMENTING PEPTIDES");
            sb.AppendLine("(Lower CV = more even; Higher MinInt = no dead spots)");
            sb.AppendLine();
            sb.AppendLine("| Rank | Sequence     | Fitness | B-ion CV | Y-ion CV | Min Intensity | Evenness |");
            sb.AppendLine("|------|--------------|---------|----------|----------|---------------|----------|");
            int rank = 1;
            foreach (var p in AllEvaluatedPeptides.OrderByDescending(f => f.OverallFitness).Take(10))
            {
                sb.AppendLine($"| {rank,4} | {p.Sequence} | {p.OverallFitness:F4}  | {p.BIonCV:F4}   | {p.YIonCV:F4}   | {p.MinIntensity:F4}        | {p.EvennessScore:F4}   |");
                rank++;
            }
            sb.AppendLine();

            sb.AppendLine("## BOTTOM 10 WORST FRAGMENTING PEPTIDES");
            sb.AppendLine("(Analyze what makes these fragment poorly)");
            sb.AppendLine();
            sb.AppendLine("| Rank | Sequence     | Fitness | B-ion CV | Y-ion CV | Min Intensity |");
            sb.AppendLine("|------|--------------|---------|----------|----------|---------------|");
            rank = 1;
            foreach (var p in AllEvaluatedPeptides.OrderBy(f => f.OverallFitness).Take(10))
            {
                sb.AppendLine($"| {rank,4} | {p.Sequence} | {p.OverallFitness:F4}  | {p.BIonCV:F4}   | {p.YIonCV:F4}   | {p.MinIntensity:F4}        |");
                rank++;
            }
            sb.AppendLine();

            sb.AppendLine("## LEARNED DIPEPTIDE PATTERNS");
            sb.AppendLine("(Positive = appears more in good peptides; Negative = appears more in bad peptides)");
            sb.AppendLine();

            var dipeptides = LearnedPatterns.Where(kv => kv.Key.Length == 2).OrderByDescending(kv => kv.Value).ToList();

            sb.AppendLine("### Most Favorable (top 15):");
            foreach (var p in dipeptides.Take(15))
                sb.AppendLine($"  {p.Key}: {p.Value:+0.00;-0.00}");

            sb.AppendLine();
            sb.AppendLine("### Most Unfavorable (bottom 15):");
            foreach (var p in dipeptides.TakeLast(15).Reverse())
                sb.AppendLine($"  {p.Key}: {p.Value:+0.00;-0.00}");
            sb.AppendLine();

            sb.AppendLine("## TERMINAL PREFERENCES");
            sb.AppendLine("(Which amino acids work best at N-terminus and C-terminus)");
            sb.AppendLine();
            var terminals = LearnedPatterns.Where(kv => kv.Key.StartsWith("N-") || kv.Key.StartsWith("C-"))
                .OrderByDescending(kv => kv.Value).ToList();
            sb.AppendLine("### N-terminal preferences:");
            foreach (var p in terminals.Where(kv => kv.Key.StartsWith("N-")).Take(5))
                sb.AppendLine($"  {p.Key}: {p.Value:+0.00;-0.00}");
            sb.AppendLine("### C-terminal preferences:");
            foreach (var p in terminals.Where(kv => kv.Key.StartsWith("C-")).Take(5))
                sb.AppendLine($"  {p.Key}: {p.Value:+0.00;-0.00}");
            sb.AppendLine();

            sb.AppendLine("## CONVERGENCE HISTORY");
            sb.AppendLine("(Best fitness per generation - is it still improving?)");
            sb.AppendLine();
            for (int i = 0; i < BestPerGeneration.Count; i++)
                sb.AppendLine($"  Gen {i,3}: {BestPerGeneration[i].Sequence} = {BestPerGeneration[i].OverallFitness:F4}");
            sb.AppendLine();

            sb.AppendLine("## QUESTIONS FOR ANALYSIS");
            sb.AppendLine("1. What chemical/structural patterns distinguish good from bad fragmenting peptides?");
            sb.AppendLine("2. Are there specific amino acids that should be at N-terminus, C-terminus, or middle?");
            sb.AppendLine("3. What dipeptide combinations should be encouraged or avoided?");
            sb.AppendLine("4. Is Proline (P) causing dead spots? Where should it be placed if at all?");
            sb.AppendLine("5. Are basic residues (K, R, H) helping the y-ion series?");
            sb.AppendLine();

            sb.AppendLine($"## PLEASE SUGGEST 10 NEW {_peptideLength}-AA SEQUENCES");
            sb.AppendLine($"Each must use exactly {_peptideLength} UNIQUE amino acids from: {string.Join("", AvailableAminoAcids)}");
            sb.AppendLine("Optimize for even fragmentation based on the patterns above.");

            return sb.ToString();
        }
    }
}