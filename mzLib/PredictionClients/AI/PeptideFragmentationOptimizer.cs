using System.Text;
using System.Text.Json;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;

namespace PredictionClients.AI
{
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

    public class CircularFragmentationFitness
    {
        public string CanonicalSequence { get; set; } = "";
        public int Generation { get; set; }
        public List<FragmentationFitness> RotationFitnesses { get; set; } = new();
        public double MinRotationFitness { get; set; }
        public double AvgRotationFitness { get; set; }
        public double MaxRotationFitness { get; set; }
        public double FitnessStdDev { get; set; }
        public FragmentationFitness WorstRotation { get; set; } = null!;
        public FragmentationFitness BestRotation { get; set; } = null!;
        public double OverallCircularFitness { get; set; }

        public override string ToString() =>
            $"{CanonicalSequence}: CircFit={OverallCircularFitness:F4} (Min={MinRotationFitness:F4}, Avg={AvgRotationFitness:F4})";
    }

    public class ClaudeAnalysisReader
    {
        public List<string> SeedSequences { get; private set; } = new();
        public Dictionary<string, double> DipeptideAdjustments { get; private set; } = new();
        public Dictionary<string, double> TerminalAdjustments { get; private set; } = new();
        public bool PrematureConvergenceWarning { get; private set; }
        public double RecommendedMutationMultiplier { get; private set; } = 1.0;
        public List<string> ProtectedCores { get; private set; } = new();
        public List<string> ChampionSequences { get; private set; } = new();

        public static ClaudeAnalysisReader LoadFromFile(string jsonPath)
        {
            var reader = new ClaudeAnalysisReader();
            var text = File.ReadAllText(jsonPath);
            using var doc = JsonDocument.Parse(text);
            var root = doc.RootElement;

            if (root.TryGetProperty("flags", out var flags))
            {
                if (flags.TryGetProperty("premature_convergence", out var pc))
                    reader.PrematureConvergenceWarning = pc.GetBoolean();
                if (flags.TryGetProperty("recommended_mutation_rate_multiplier", out var mm))
                    reader.RecommendedMutationMultiplier = mm.GetDouble();
            }

            if (root.TryGetProperty("seed_sequences", out var seeds))
                foreach (var seed in seeds.EnumerateArray())
                    if (seed.TryGetProperty("sequence", out var seq))
                        reader.SeedSequences.Add(seq.GetString() ?? "");

            if (root.TryGetProperty("dipeptide_weight_adjustments", out var dwa))
            {
                foreach (var section in new[] { "reinforce", "penalise_more" })
                    if (dwa.TryGetProperty(section, out var block))
                        foreach (var kv in block.EnumerateObject())
                            if (kv.Name.Length == 2)
                                reader.DipeptideAdjustments[kv.Name] =
                                    reader.DipeptideAdjustments.GetValueOrDefault(kv.Name) + kv.Value.GetDouble();
            }

            if (root.TryGetProperty("terminal_preference_adjustments", out var tpa))
            {
                foreach (var termSection in new[] { "n_terminal", "c_terminal" })
                    if (tpa.TryGetProperty(termSection, out var block))
                    {
                        string prefix = termSection == "n_terminal" ? "N-" : "C-";
                        foreach (var kv in block.EnumerateObject())
                            if (kv.Name.Length == 1)
                                reader.TerminalAdjustments[$"{prefix}{kv.Name}"] =
                                    reader.TerminalAdjustments.GetValueOrDefault($"{prefix}{kv.Name}") + kv.Value.GetDouble();
                    }
            }

            if (root.TryGetProperty("protected_cores", out var cores))
                foreach (var core in cores.EnumerateArray())
                    if (core.GetString() is string s && !string.IsNullOrWhiteSpace(s))
                        reader.ProtectedCores.Add(s);

            if (root.TryGetProperty("champion_sequences", out var champs))
                foreach (var champ in champs.EnumerateArray())
                    if (champ.GetString() is string s && !string.IsNullOrWhiteSpace(s))
                        reader.ChampionSequences.Add(s);

            return reader;
        }
    }

    public class PeptideFragmentationOptimizer
    {
        private static readonly char[] AvailableAminoAcids =
            { 'A', 'E', 'F', 'H', 'Q', 'R', 'S', 'T', 'Y', 'G' };

        private readonly int _peptideLength;
        private readonly int _populationSize;
        private readonly int _collisionEnergy;
        private readonly int _precursorCharge;
        private readonly Random _random;
        private double _lastGenerationDiversity = 1.0;
        private readonly List<string> _protectedCores = new();
        private readonly List<string> _championSequences = new();

        public List<FragmentationFitness> AllEvaluatedPeptides { get; } = new();
        public List<FragmentationFitness> BestPerGeneration { get; } = new();
        public List<CircularFragmentationFitness> AllCircularEvaluatedPeptides { get; } = new();
        public List<CircularFragmentationFitness> BestCircularPerGeneration { get; } = new();
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

        public int GetPeptideLength() => _peptideLength;

        public List<string> InjectClaudeAnalysis(ClaudeAnalysisReader analysis)
        {
            foreach (var kv in analysis.DipeptideAdjustments)
            {
                LearnedPatterns.TryAdd(kv.Key, 0);
                LearnedPatterns[kv.Key] += kv.Value;
            }

            foreach (var kv in analysis.TerminalAdjustments)
            {
                LearnedPatterns.TryAdd(kv.Key, 0);
                LearnedPatterns[kv.Key] += kv.Value;
            }

            foreach (var core in analysis.ProtectedCores)
                if (!_protectedCores.Contains(core))
                    _protectedCores.Add(core);
            _protectedCores.Sort((a, b) => b.Length.CompareTo(a.Length));

            foreach (var champ in analysis.ChampionSequences)
                if (!_championSequences.Contains(champ) &&
                    champ.Length == _peptideLength &&
                    champ.Distinct().Count() == _peptideLength &&
                    champ.All(c => AvailableAminoAcids.Contains(c)))
                    _championSequences.Add(champ);

            return analysis.SeedSequences
                .Where(s => s.Length == _peptideLength &&
                            s.Distinct().Count() == _peptideLength &&
                            s.All(c => AvailableAminoAcids.Contains(c)))
                .ToList();
        }

        private int FindProtectedCore(char[] sequence, out string foundCore)
        {
            var s = new string(sequence);
            foreach (var core in _protectedCores)
            {
                int idx = s.IndexOf(core, StringComparison.Ordinal);
                if (idx >= 0) { foundCore = core; return idx; }
            }
            foundCore = "";
            return -1;
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

        public List<string> GenerateInitialPopulation(List<string>? seeds = null)
        {
            var population = new HashSet<string>();
            if (seeds != null)
                foreach (var s in seeds)
                    population.Add(s);

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

            double deadSpotScore = minIntensity < 1e-6 ? 0.0 : 1.0 - Math.Exp(-minIntensity / 0.01);

            var allIons = bIons.Concat(yIons).ToList();
            double avgIntensity = allIons.Count > 0 ? allIons.Average() : 0;

            double overallFitness = 0.45 * evennessScore + 0.35 * deadSpotScore + 0.20 * avgIntensity;

            return new FragmentationFitness
            {
                Sequence = sequence,
                BIonIntensities = bIons.ToArray(),
                YIonIntensities = yIons.ToArray(),
                BIonCV = bCV,
                YIonCV = yCV,
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
            int minCross = 1, maxCross = _peptideLength - 1;

            foreach (var core in _protectedCores)
            {
                int idx = parent1.IndexOf(core, StringComparison.Ordinal);
                if (idx < 0) continue;
                int coreEnd = idx + core.Length;
                if (idx < _peptideLength / 2) minCross = Math.Max(minCross, coreEnd);
                else maxCross = Math.Min(maxCross, idx);
            }

            if (minCross >= maxCross) minCross = 1;
            int crossPoint = _random.Next(minCross, Math.Max(minCross + 1, maxCross));

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

            var locked = new bool[_peptideLength];
            int coreStart = FindProtectedCore(sequence, out string foundCore);
            if (coreStart >= 0)
                for (int k = coreStart; k < coreStart + foundCore.Length; k++)
                    locked[k] = true;

            for (int i = 0; i < _peptideLength; i++)
            {
                if (locked[i] || _random.NextDouble() >= mutationRate) continue;

                if (_random.NextDouble() < 0.5 || unusedAAs.Count == 0)
                {
                    int attempts = 0, j;
                    do { j = _random.Next(_peptideLength); attempts++; } while (locked[j] && attempts < 10);
                    if (!locked[j]) (sequence[i], sequence[j]) = (sequence[j], sequence[i]);
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
            return new string(sequence);
        }

        private FragmentationFitness TournamentSelect(List<FragmentationFitness> pop, int size = 5)
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
            int sampleSize = Math.Max(10, results.Count / 5);
            var top = results.Take(sampleSize).Select(f => f.Sequence).ToList();
            var bottom = results.TakeLast(sampleSize).Select(f => f.Sequence).ToList();

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
            {
                for (int i = 0; i < pep.Length - 1; i++)
                {
                    string di = pep.Substring(i, 2);
                    LearnedPatterns.TryAdd(di, 0); LearnedPatterns[di] -= 0.1;
                }
                LearnedPatterns.TryAdd($"N-{pep[0]}", 0); LearnedPatterns[$"N-{pep[0]}"] -= 0.05;
                LearnedPatterns.TryAdd($"C-{pep[^1]}", 0); LearnedPatterns[$"C-{pep[^1]}"] -= 0.05;
            }
        }

        private string ApplyLearnedPatterns(string peptide)
        {
            var seq = peptide.ToCharArray();

            var good = LearnedPatterns.Where(kv => kv.Value > 0.3 && kv.Key.Length == 2)
                .OrderByDescending(kv => kv.Value).Take(5).Select(kv => kv.Key);

            foreach (var di in good)
            {
                if (_random.NextDouble() < 0.5)
                {
                    int p1 = Array.IndexOf(seq, di[0]), p2 = Array.IndexOf(seq, di[1]);
                    if (p1 >= 0 && p2 >= 0 && p1 != p2 - 1 && p1 + 1 < _peptideLength && p1 + 1 != p2)
                        (seq[p2], seq[p1 + 1]) = (seq[p1 + 1], seq[p2]);
                }
            }
            return new string(seq);
        }

        #region Circular Peptide Support

        public static string GetCanonicalRotation(string peptide)
        {
            string canonical = peptide;
            for (int i = 1; i < peptide.Length; i++)
            {
                string rotation = peptide.Substring(i) + peptide.Substring(0, i);
                if (string.Compare(rotation, canonical, StringComparison.Ordinal) < 0)
                    canonical = rotation;
            }
            return canonical;
        }

        public static List<string> GetAllRotations(string peptide)
        {
            var rotations = new List<string>();
            for (int i = 0; i < peptide.Length; i++)
                rotations.Add(peptide.Substring(i) + peptide.Substring(0, i));
            return rotations;
        }

        public List<CircularFragmentationFitness> GetUniqueTopCircularPeptides(int count = 5)
        {
            var unique = new List<CircularFragmentationFitness>();
            var seenCanonical = new HashSet<string>();

            foreach (var p in AllCircularEvaluatedPeptides.OrderByDescending(f => f.OverallCircularFitness))
            {
                string canonical = GetCanonicalRotation(p.CanonicalSequence);
                if (!seenCanonical.Contains(canonical))
                {
                    seenCanonical.Add(canonical);
                    unique.Add(p);
                    if (unique.Count >= count) break;
                }
            }
            return unique;
        }

        private async Task<CircularFragmentationFitness> CalculateCircularFitnessAsync(string peptide, int generation)
        {
            var rotations = GetAllRotations(peptide);
            var charges = Enumerable.Repeat(_precursorCharge, rotations.Count).ToList();
            var energies = Enumerable.Repeat(_collisionEnergy, rotations.Count).ToList();
            var retentionTimes = Enumerable.Repeat<double?>(null, rotations.Count).ToList();

            var model = new Prosit2020IntensityHCD(rotations, charges, energies, retentionTimes, out var warnings);

            if (model.PeptideSequences.Count == 0)
                throw new InvalidOperationException($"No valid rotations. Warnings: {warnings?.Message}");

            await model.RunInferenceAsync();

            var rotationFitnesses = new List<FragmentationFitness>();
            for (int i = 0; i < model.Predictions.Count; i++)
                rotationFitnesses.Add(CalculateFitness(rotations[i], model.Predictions[i], generation));

            double minFitness = rotationFitnesses.Min(f => f.OverallFitness);
            double avgFitness = rotationFitnesses.Average(f => f.OverallFitness);
            double maxFitness = rotationFitnesses.Max(f => f.OverallFitness);
            double stdDev = Math.Sqrt(rotationFitnesses.Select(f => Math.Pow(f.OverallFitness - avgFitness, 2)).Average());

            return new CircularFragmentationFitness
            {
                CanonicalSequence = peptide,
                Generation = generation,
                RotationFitnesses = rotationFitnesses,
                MinRotationFitness = minFitness,
                AvgRotationFitness = avgFitness,
                MaxRotationFitness = maxFitness,
                FitnessStdDev = stdDev,
                WorstRotation = rotationFitnesses.OrderBy(f => f.OverallFitness).First(),
                BestRotation = rotationFitnesses.OrderByDescending(f => f.OverallFitness).First(),
                OverallCircularFitness = 0.6 * minFitness + 0.3 * avgFitness + 0.1 * maxFitness
            };
        }

        public async Task<List<CircularFragmentationFitness>> EvaluateCircularPeptidesAsync(List<string> peptides, int generation)
        {
            var results = new List<CircularFragmentationFitness>();
            foreach (var peptide in peptides)
            {
                var result = await CalculateCircularFitnessAsync(peptide, generation);
                results.Add(result);
                AllCircularEvaluatedPeptides.Add(result);
            }
            return results;
        }

        public async Task<CircularFragmentationFitness> OptimizeCircularAsync(
            int generations = 50,
            double eliteRatio = 0.05,
            double mutationRate = 0.2,
            List<string>? claudeSeeds = null,
            Action<int, CircularFragmentationFitness, double>? onGenerationComplete = null)
        {
            var population = GenerateInitialPopulation(claudeSeeds);
            CircularFragmentationFitness? bestEver = null;
            int stallCount = 0;
            double prevBestFitness = 0;

            for (int gen = 0; gen < generations; gen++)
            {
                var fitnessResults = (await EvaluateCircularPeptidesAsync(population, gen))
                    .OrderByDescending(f => f.OverallCircularFitness).ToList();

                var bestThisGen = fitnessResults.First();
                BestCircularPerGeneration.Add(bestThisGen);

                if (bestEver == null || bestThisGen.OverallCircularFitness > bestEver.OverallCircularFitness)
                    bestEver = bestThisGen;

                double diversity = (double)fitnessResults.Select(f => GetCanonicalRotation(f.CanonicalSequence)).Distinct().Count() / fitnessResults.Count;
                _lastGenerationDiversity = diversity;

                if (Math.Abs(bestThisGen.OverallCircularFitness - prevBestFitness) < 1e-5) stallCount++;
                else { stallCount = 0; prevBestFitness = bestThisGen.OverallCircularFitness; }

                double activeMutationRate = mutationRate;
                if (stallCount >= 3 || diversity < 0.5) activeMutationRate = Math.Min(mutationRate * 2.0, 0.6);
                if (stallCount >= 6 || diversity < 0.2) activeMutationRate = Math.Min(mutationRate * 4.0, 0.8);

                LearnPatterns(fitnessResults.Select(f => f.WorstRotation).ToList());
                onGenerationComplete?.Invoke(gen, bestThisGen, diversity);

                var nextPopulation = new List<string>();

                int eliteCount = (int)(_populationSize * eliteRatio);
                foreach (var elite in fitnessResults.Take(eliteCount))
                    nextPopulation.Add(elite.CanonicalSequence);

                foreach (var champ in _championSequences)
                    if (!nextPopulation.Contains(champ))
                        nextPopulation.Add(champ);

                bool collapsed = fitnessResults.Take(10).Select(f => GetCanonicalRotation(f.CanonicalSequence)).Distinct().Count() == 1;
                if (collapsed)
                    for (int i = 0; i < (int)(_populationSize * 0.20); i++)
                        nextPopulation.Add(GenerateRandomPeptide());

                while (nextPopulation.Count < _populationSize)
                {
                    var p1 = fitnessResults[_random.Next(Math.Min(5, fitnessResults.Count))];
                    var p2 = fitnessResults[_random.Next(Math.Min(5, fitnessResults.Count))];
                    var child = Mutate(Crossover(p1.CanonicalSequence, p2.CanonicalSequence), activeMutationRate);
                    if (_random.NextDouble() < 0.6) child = ApplyLearnedPatterns(child);
                    nextPopulation.Add(child);
                }

                population = nextPopulation;
            }

            return bestEver!;
        }

        #endregion

        #region Reporting

        public string GetCircularProgressReport()
        {
            var sb = new StringBuilder();
            sb.AppendLine("=== CIRCULAR PEPTIDE OPTIMIZATION ===");
            sb.AppendLine($"Peptides evaluated: {AllCircularEvaluatedPeptides.Count}");
            sb.AppendLine($"Generations: {BestCircularPerGeneration.Count}");
            sb.AppendLine();
            sb.AppendLine("Top 5 UNIQUE circular peptides:");
            foreach (var p in GetUniqueTopCircularPeptides(5))
            {
                string canonical = GetCanonicalRotation(p.CanonicalSequence);
                sb.AppendLine($"  {canonical} CircFit={p.OverallCircularFitness:F4} Min={p.MinRotationFitness:F4}");
                sb.AppendLine($"    Worst: {p.WorstRotation.Sequence} Best: {p.BestRotation.Sequence}");
            }
            return sb.ToString();
        }

        public string GetLearnedPatternsSummary()
        {
            var sb = new StringBuilder();
            sb.AppendLine("=== LEARNED PATTERNS ===");
            sb.AppendLine("Favorable:");
            foreach (var p in LearnedPatterns.Where(kv => kv.Key.Length == 2).OrderByDescending(kv => kv.Value).Take(10))
                sb.AppendLine($"  {p.Key}: +{p.Value:F2}");
            sb.AppendLine("Unfavorable:");
            foreach (var p in LearnedPatterns.Where(kv => kv.Key.Length == 2).OrderBy(kv => kv.Value).Take(10))
                sb.AppendLine($"  {p.Key}: {p.Value:F2}");
            return sb.ToString();
        }

        public string GenerateCircularClaudePrompt()
        {
            var sb = new StringBuilder();
            sb.AppendLine("# CIRCULAR PEPTIDE FRAGMENTATION OPTIMIZATION");
            sb.AppendLine();
            sb.AppendLine($"## CONSTRAINTS: {_peptideLength} AA ring, each used once");
            sb.AppendLine($"Available: {string.Join("", AvailableAminoAcids)}");
            sb.AppendLine("Rotations of same ring = IDENTICAL");
            sb.AppendLine();
            sb.AppendLine($"## STATS: {AllCircularEvaluatedPeptides.Count} peptides, {BestCircularPerGeneration.Count} generations");
            sb.AppendLine();
            sb.AppendLine("## TOP 5 UNIQUE RINGS");
            foreach (var p in GetUniqueTopCircularPeptides(5))
            {
                string canonical = GetCanonicalRotation(p.CanonicalSequence);
                sb.AppendLine($"{canonical} CircFit={p.OverallCircularFitness:F4} Min={p.MinRotationFitness:F4} Worst@{p.WorstRotation.Sequence}");
            }
            sb.AppendLine();
            sb.AppendLine("## LEARNED PATTERNS");
            var dip = LearnedPatterns.Where(kv => kv.Key.Length == 2).OrderByDescending(kv => kv.Value).ToList();
            sb.AppendLine("Good: " + string.Join(", ", dip.Take(8).Select(p => $"{p.Key}({p.Value:+0.0;-0.0})")));
            sb.AppendLine("Bad: " + string.Join(", ", dip.TakeLast(8).Reverse().Select(p => $"{p.Key}({p.Value:+0.0;-0.0})")));
            sb.AppendLine();
            sb.AppendLine($"## SUGGEST 10 NEW {_peptideLength}-AA CIRCULAR SEQUENCES (canonical form)");
            sb.AppendLine("Return JSON: {\"seed_sequences\": [{\"sequence\": \"...\"}], \"champion_sequences\": [...]}");
            return sb.ToString();
        }

        #endregion
    }
}