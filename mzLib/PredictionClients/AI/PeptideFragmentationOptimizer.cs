using System.Text;
using System.Text.Json;
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
    /// Fitness metrics for circular peptides, considering all possible linearization points.
    /// </summary>
    public class CircularFragmentationFitness
    {
        public string CanonicalSequence { get; set; } = "";
        public int Generation { get; set; }

        /// <summary>Fitness results for each rotation of the circular peptide.</summary>
        public List<FragmentationFitness> RotationFitnesses { get; set; } = new();

        /// <summary>Minimum fitness across all rotations (the weak link).</summary>
        public double MinRotationFitness { get; set; }

        /// <summary>Average fitness across all rotations.</summary>
        public double AvgRotationFitness { get; set; }

        /// <summary>Maximum fitness across all rotations.</summary>
        public double MaxRotationFitness { get; set; }

        /// <summary>Standard deviation of fitness across rotations (lower = more consistent).</summary>
        public double FitnessStdDev { get; set; }

        /// <summary>The rotation with the worst fragmentation.</summary>
        public FragmentationFitness WorstRotation { get; set; } = null!;

        /// <summary>The rotation with the best fragmentation.</summary>
        public FragmentationFitness BestRotation { get; set; } = null!;

        /// <summary>
        /// Overall circular fitness score. 
        /// Heavily weighted toward minimum to penalize any bad rotation.
        /// </summary>
        public double OverallCircularFitness { get; set; }

        public override string ToString() =>
            $"{CanonicalSequence}: CircularFitness={OverallCircularFitness:F4} (Min={MinRotationFitness:F4}, Avg={AvgRotationFitness:F4}, StdDev={FitnessStdDev:F4})";
    }

    /// <summary>
    /// Deserialises the structured JSON file produced by Claude analysis
    /// and exposes seed sequences, weight adjustments, and flags in a form
    /// the optimizer can consume directly.
    /// </summary>
    public class ClaudeAnalysisReader
    {
        public List<string> SeedSequences { get; private set; } = new();
        public Dictionary<string, double> DipeptideAdjustments { get; private set; } = new();
        public Dictionary<string, double> TerminalAdjustments { get; private set; } = new();
        public bool PrematureConvergenceWarning { get; private set; }
        public double RecommendedMutationMultiplier { get; private set; } = 1.0;
        public int RecommendedPopulationIncrease { get; private set; } = 0;
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
                if (flags.TryGetProperty("recommended_population_increase", out var pi))
                    reader.RecommendedPopulationIncrease = pi.GetInt32();
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
                {
                    var s = core.GetString();
                    if (!string.IsNullOrWhiteSpace(s)) reader.ProtectedCores.Add(s);
                }

            if (root.TryGetProperty("champion_sequences", out var champs))
                foreach (var champ in champs.EnumerateArray())
                {
                    var s = champ.GetString();
                    if (!string.IsNullOrWhiteSpace(s)) reader.ChampionSequences.Add(s);
                }

            return reader;
        }
    }

    /// <summary>
    /// Optimizes peptide sequences for even HCD fragmentation using Prosit predictions
    /// and an evolutionary algorithm that learns favorable dipeptide patterns.
    /// </summary>
    public class PeptideFragmentationOptimizer
    {
        private static readonly char[] AvailableAminoAcids =
            { 'A', 'E', 'F', 'H', 'N', 'R', 'S', 'T', 'Y', 'K' };

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

        /// <summary>Gets the peptide length this optimizer is configured for.</summary>
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

            var validSeeds = analysis.SeedSequences
                .Where(s => s.Length == _peptideLength &&
                            s.Distinct().Count() == _peptideLength &&
                            s.All(c => AvailableAminoAcids.Contains(c)))
                .ToList();

            return validSeeds;
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

            double deadSpotScore = minIntensity < 1e-6
                ? 0.0
                : 1.0 - Math.Exp(-minIntensity / 0.01);

            var allIons = bIons.Concat(yIons).ToList();
            double avgIntensity = allIons.Count > 0 ? allIons.Average() : 0;

            double overallFitness = 0.45 * evennessScore
                                  + 0.35 * deadSpotScore
                                  + 0.20 * avgIntensity;

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
            int minCross = 1;
            int maxCross = _peptideLength - 1;

            foreach (var core in _protectedCores)
            {
                int idx = parent1.IndexOf(core, StringComparison.Ordinal);
                if (idx < 0) continue;
                int coreEnd = idx + core.Length;
                if (idx < _peptideLength / 2)
                    minCross = Math.Max(minCross, coreEnd);
                else
                    maxCross = Math.Min(maxCross, idx);
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
                if (locked[i]) continue;
                if (_random.NextDouble() >= mutationRate) continue;

                if (_random.NextDouble() < 0.5 || unusedAAs.Count == 0)
                {
                    int attempts = 0;
                    int j;
                    do { j = _random.Next(_peptideLength); attempts++; }
                    while (locked[j] && attempts < 10);
                    if (!locked[j])
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
            return new string(sequence);
        }

        public async Task<FragmentationFitness> OptimizeAsync(
            int generations = 50,
            double eliteRatio = 0.05,
            double mutationRate = 0.2,
            List<string>? claudeSeeds = null,
            Action<int, FragmentationFitness, double>? onGenerationComplete = null)
        {
            var population = GenerateInitialPopulation(claudeSeeds);
            FragmentationFitness? bestEver = null;
            int stallCount = 0;
            double prevBestFitness = 0;

            for (int gen = 0; gen < generations; gen++)
            {
                var fitnessResults = (await EvaluatePeptidesAsync(population, gen))
                    .OrderByDescending(f => f.OverallFitness).ToList();

                var bestThisGen = fitnessResults.First();
                BestPerGeneration.Add(bestThisGen);
                if (bestEver == null || bestThisGen.OverallFitness > bestEver.OverallFitness)
                    bestEver = bestThisGen;

                double diversity = (double)fitnessResults.Select(f => f.Sequence).Distinct().Count()
                                   / fitnessResults.Count;
                _lastGenerationDiversity = diversity;

                if (Math.Abs(bestThisGen.OverallFitness - prevBestFitness) < 1e-5) stallCount++;
                else { stallCount = 0; prevBestFitness = bestThisGen.OverallFitness; }

                double activeMutationRate = mutationRate;
                if (stallCount >= 3 || diversity < 0.5)
                    activeMutationRate = Math.Min(mutationRate * 2.0, 0.6);
                if (stallCount >= 6 || diversity < 0.2)
                    activeMutationRate = Math.Min(mutationRate * 4.0, 0.8);

                LearnPatterns(fitnessResults);
                onGenerationComplete?.Invoke(gen, bestThisGen, diversity);

                var nextPopulation = new List<string>();

                int eliteCount = (int)(_populationSize * eliteRatio);
                foreach (var elite in fitnessResults.Take(eliteCount))
                    nextPopulation.Add(elite.Sequence);

                foreach (var champ in _championSequences)
                    if (!nextPopulation.Contains(champ))
                        nextPopulation.Add(champ);

                bool collapsed = fitnessResults.Take(10).Select(f => f.Sequence).Distinct().Count() == 1;
                int injectionCount = collapsed ? (int)(_populationSize * 0.20) : 0;
                for (int i = 0; i < injectionCount; i++)
                    nextPopulation.Add(GenerateRandomPeptide());

                while (nextPopulation.Count < _populationSize)
                {
                    var p1 = TournamentSelect(fitnessResults, size: 5);
                    var p2 = TournamentSelect(fitnessResults, size: 5);
                    var child = Mutate(Crossover(p1.Sequence, p2.Sequence), activeMutationRate);

                    if (_random.NextDouble() < 0.6)
                        child = ApplyLearnedPatterns(child);

                    nextPopulation.Add(child);
                }

                population = nextPopulation;
            }

            return bestEver!;
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
                    LearnedPatterns.TryAdd(di, 0);
                    LearnedPatterns[di] += 0.1;
                }
                LearnedPatterns.TryAdd($"N-{pep[0]}", 0);
                LearnedPatterns[$"N-{pep[0]}"] += 0.05;
                LearnedPatterns.TryAdd($"C-{pep[^1]}", 0);
                LearnedPatterns[$"C-{pep[^1]}"] += 0.05;
            }

            foreach (var pep in bottom)
            {
                for (int i = 0; i < pep.Length - 1; i++)
                {
                    string di = pep.Substring(i, 2);
                    LearnedPatterns.TryAdd(di, 0);
                    LearnedPatterns[di] -= 0.1;
                }
                LearnedPatterns.TryAdd($"N-{pep[0]}", 0);
                LearnedPatterns[$"N-{pep[0]}"] -= 0.05;
                LearnedPatterns.TryAdd($"C-{pep[^1]}", 0);
                LearnedPatterns[$"C-{pep[^1]}"] -= 0.05;
            }
        }

        private string ApplyLearnedPatterns(string peptide)
        {
            var seq = peptide.ToCharArray();

            var good = LearnedPatterns
                .Where(kv => kv.Value > 0.3 && kv.Key.Length == 2)
                .OrderByDescending(kv => kv.Value)
                .Take(5)
                .Select(kv => kv.Key);

            foreach (var di in good)
            {
                if (_random.NextDouble() < 0.5)
                {
                    int p1 = Array.IndexOf(seq, di[0]);
                    int p2 = Array.IndexOf(seq, di[1]);
                    if (p1 >= 0 && p2 >= 0 && p1 != p2 - 1 && p1 + 1 < _peptideLength && p1 + 1 != p2)
                        (seq[p2], seq[p1 + 1]) = (seq[p1 + 1], seq[p2]);
                }
            }

            var nTermScore = LearnedPatterns.GetValueOrDefault($"N-{seq[0]}", 0);
            if (nTermScore < -0.1 && _random.NextDouble() < 0.4)
            {
                var bestNTerm = LearnedPatterns
                    .Where(kv => kv.Key.StartsWith("N-") && kv.Value > 0)
                    .OrderByDescending(kv => kv.Value)
                    .Select(kv => kv.Key[2])
                    .FirstOrDefault(c => seq.Contains(c) && c != seq[0]);

                if (bestNTerm != default)
                {
                    int swapIdx = Array.IndexOf(seq, bestNTerm);
                    (seq[0], seq[swapIdx]) = (seq[swapIdx], seq[0]);
                }
            }

            return new string(seq);
        }

        #region Circular Peptide Support

        public static List<string> GetAllRotations(string peptide)
        {
            var rotations = new List<string>();
            for (int i = 0; i < peptide.Length; i++)
            {
                string rotation = peptide.Substring(i) + peptide.Substring(0, i);
                rotations.Add(rotation);
            }
            return rotations;
        }

        public static string GetRotation(string peptide, int rotationIndex)
        {
            int i = rotationIndex % peptide.Length;
            return peptide.Substring(i) + peptide.Substring(0, i);
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
            {
                var fitness = CalculateFitness(rotations[i], model.Predictions[i], generation);
                rotationFitnesses.Add(fitness);
            }

            double minFitness = rotationFitnesses.Min(f => f.OverallFitness);
            double avgFitness = rotationFitnesses.Average(f => f.OverallFitness);
            double maxFitness = rotationFitnesses.Max(f => f.OverallFitness);
            double fitnessStdDev = Math.Sqrt(rotationFitnesses
                .Select(f => Math.Pow(f.OverallFitness - avgFitness, 2)).Average());

            var worstRotation = rotationFitnesses.OrderBy(f => f.OverallFitness).First();
            var bestRotation = rotationFitnesses.OrderByDescending(f => f.OverallFitness).First();

            return new CircularFragmentationFitness
            {
                CanonicalSequence = peptide,
                Generation = generation,
                RotationFitnesses = rotationFitnesses,
                MinRotationFitness = minFitness,
                AvgRotationFitness = avgFitness,
                MaxRotationFitness = maxFitness,
                FitnessStdDev = fitnessStdDev,
                WorstRotation = worstRotation,
                BestRotation = bestRotation,
                OverallCircularFitness = 0.6 * minFitness + 0.3 * avgFitness + 0.1 * maxFitness
            };
        }

        public async Task<List<CircularFragmentationFitness>> EvaluateCircularPeptidesAsync(List<string> peptides, int generation)
        {
            var results = new List<CircularFragmentationFitness>();
            int batchSize = Math.Max(1, 100 / _peptideLength);

            for (int i = 0; i < peptides.Count; i += batchSize)
            {
                var batch = peptides.Skip(i).Take(batchSize).ToList();
                var batchTasks = batch.Select(p => CalculateCircularFitnessAsync(p, generation));
                var batchResults = await Task.WhenAll(batchTasks);
                results.AddRange(batchResults);
            }

            foreach (var result in results)
                AllCircularEvaluatedPeptides.Add(result);

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

                double diversity = (double)fitnessResults.Select(f => f.CanonicalSequence).Distinct().Count()
                                   / fitnessResults.Count;
                _lastGenerationDiversity = diversity;

                if (Math.Abs(bestThisGen.OverallCircularFitness - prevBestFitness) < 1e-5) stallCount++;
                else { stallCount = 0; prevBestFitness = bestThisGen.OverallCircularFitness; }

                double activeMutationRate = mutationRate;
                if (stallCount >= 3 || diversity < 0.5) activeMutationRate = Math.Min(mutationRate * 2.0, 0.6);
                if (stallCount >= 6 || diversity < 0.2) activeMutationRate = Math.Min(mutationRate * 4.0, 0.8);

                var linearResults = fitnessResults.Select(f => f.WorstRotation).ToList();
                LearnPatterns(linearResults);

                onGenerationComplete?.Invoke(gen, bestThisGen, diversity);

                var nextPopulation = new List<string>();

                int eliteCount = (int)(_populationSize * eliteRatio);
                foreach (var elite in fitnessResults.Take(eliteCount))
                    nextPopulation.Add(elite.CanonicalSequence);

                foreach (var champ in _championSequences)
                    if (!nextPopulation.Contains(champ))
                        nextPopulation.Add(champ);

                bool collapsed = fitnessResults.Take(10).Select(f => f.CanonicalSequence).Distinct().Count() == 1;
                int injectionCount = collapsed ? (int)(_populationSize * 0.20) : 0;
                for (int i = 0; i < injectionCount; i++)
                    nextPopulation.Add(GenerateRandomPeptide());

                while (nextPopulation.Count < _populationSize)
                {
                    var p1 = CircularTournamentSelect(fitnessResults, size: 5);
                    var p2 = CircularTournamentSelect(fitnessResults, size: 5);
                    var child = Mutate(Crossover(p1.CanonicalSequence, p2.CanonicalSequence), activeMutationRate);
                    if (_random.NextDouble() < 0.6)
                        child = ApplyLearnedPatterns(child);
                    nextPopulation.Add(child);
                }

                population = nextPopulation;
            }

            return bestEver!;
        }

        private CircularFragmentationFitness CircularTournamentSelect(List<CircularFragmentationFitness> pop, int size = 5)
        {
            CircularFragmentationFitness? best = null;
            for (int i = 0; i < size; i++)
            {
                var c = pop[_random.Next(pop.Count)];
                if (best == null || c.OverallCircularFitness > best.OverallCircularFitness) best = c;
            }
            return best!;
        }

        #endregion

        #region Reporting

        public string GetProgressReport()
        {
            var sb = new StringBuilder();
            sb.AppendLine($"=== OPTIMIZATION PROGRESS ===");
            sb.AppendLine($"Peptides evaluated: {AllEvaluatedPeptides.Count}");
            sb.AppendLine($"Generations: {BestPerGeneration.Count}");
            sb.AppendLine();
            sb.AppendLine("Best per generation:");
            for (int i = 0; i < BestPerGeneration.Count; i++)
                sb.AppendLine($"  Gen {i,3}: {BestPerGeneration[i].Sequence} - Fitness: {BestPerGeneration[i].OverallFitness:F4}");
            sb.AppendLine();
            sb.AppendLine("Top 5 overall:");
            foreach (var p in AllEvaluatedPeptides.OrderByDescending(f => f.OverallFitness).Take(5))
                sb.AppendLine($"  {p}");
            return sb.ToString();
        }

        public string GetLearnedPatternsSummary()
        {
            var sb = new StringBuilder();
            sb.AppendLine("=== LEARNED PATTERNS ===");
            sb.AppendLine();
            sb.AppendLine("Favorable dipeptides:");
            foreach (var p in LearnedPatterns.Where(kv => kv.Key.Length == 2).OrderByDescending(kv => kv.Value).Take(10))
                sb.AppendLine($"  {p.Key}: +{p.Value:F2}");
            sb.AppendLine();
            sb.AppendLine("Unfavorable dipeptides:");
            foreach (var p in LearnedPatterns.Where(kv => kv.Key.Length == 2).OrderBy(kv => kv.Value).Take(10))
                sb.AppendLine($"  {p.Key}: {p.Value:F2}");
            return sb.ToString();
        }

        public string GetCircularProgressReport()
        {
            var sb = new StringBuilder();
            sb.AppendLine("=== CIRCULAR PEPTIDE OPTIMIZATION PROGRESS ===");
            sb.AppendLine($"Peptides evaluated: {AllCircularEvaluatedPeptides.Count}");
            sb.AppendLine($"Generations: {BestCircularPerGeneration.Count}");
            sb.AppendLine($"Rotations per peptide: {_peptideLength}");
            sb.AppendLine();
            sb.AppendLine("Best circular fitness per generation:");
            for (int i = 0; i < BestCircularPerGeneration.Count; i++)
            {
                var best = BestCircularPerGeneration[i];
                sb.AppendLine($"  Gen {i,3}: {best.CanonicalSequence} - CircFit={best.OverallCircularFitness:F4} (Min={best.MinRotationFitness:F4}, Worst@{best.WorstRotation.Sequence})");
            }
            sb.AppendLine();
            sb.AppendLine("Top 5 circular peptides overall:");
            foreach (var p in AllCircularEvaluatedPeptides.OrderByDescending(f => f.OverallCircularFitness).Take(5))
            {
                sb.AppendLine($"  {p}");
                sb.AppendLine($"    Worst: {p.WorstRotation.Sequence} (Fit={p.WorstRotation.OverallFitness:F4})");
                sb.AppendLine($"    Best:  {p.BestRotation.Sequence} (Fit={p.BestRotation.OverallFitness:F4})");
            }
            return sb.ToString();
        }

        public string GenerateClaudePrompt()
        {
            var sb = new StringBuilder();
            sb.AppendLine("# PEPTIDE FRAGMENTATION OPTIMIZATION RESULTS");
            sb.AppendLine();
            sb.AppendLine("## GOAL");
            sb.AppendLine($"Design {_peptideLength}-AA peptides that fragment EVENLY across all backbone positions in HCD.");
            sb.AppendLine();
            sb.AppendLine("## CONSTRAINTS");
            sb.AppendLine($"- Exactly {_peptideLength} amino acids, each used only ONCE");
            sb.AppendLine($"- Available AAs: {string.Join(", ", AvailableAminoAcids)}");
            sb.AppendLine($"- Collision energy: {_collisionEnergy} eV, Charge: +{_precursorCharge}");
            sb.AppendLine();
            sb.AppendLine("## STATISTICS");
            sb.AppendLine($"- Peptides evaluated: {AllEvaluatedPeptides.Count:N0}");
            sb.AppendLine($"- Generations: {BestPerGeneration.Count}");
            sb.AppendLine($"- Fitness range: {AllEvaluatedPeptides.Min(p => p.OverallFitness):F4} to {AllEvaluatedPeptides.Max(p => p.OverallFitness):F4}");
            sb.AppendLine($"- Final diversity: {_lastGenerationDiversity:P1}");
            sb.AppendLine();
            sb.AppendLine("## TOP 10 PEPTIDES");
            sb.AppendLine("| Rank | Sequence | Fitness | B-CV | Y-CV | MinInt | Evenness |");
            sb.AppendLine("|------|----------|---------|------|------|--------|----------|");
            int rank = 1;
            foreach (var p in AllEvaluatedPeptides.OrderByDescending(f => f.OverallFitness).Take(10))
            {
                sb.AppendLine($"| {rank,4} | {p.Sequence} | {p.OverallFitness:F4} | {p.BIonCV:F4} | {p.YIonCV:F4} | {p.MinIntensity:F4} | {p.EvennessScore:F4} |");
                rank++;
            }
            sb.AppendLine();
            sb.AppendLine("## LEARNED PATTERNS");
            var dipeptides = LearnedPatterns.Where(kv => kv.Key.Length == 2).OrderByDescending(kv => kv.Value).ToList();
            sb.AppendLine("Favorable: " + string.Join(", ", dipeptides.Take(10).Select(p => $"{p.Key}({p.Value:+0.0;-0.0})")));
            sb.AppendLine("Unfavorable: " + string.Join(", ", dipeptides.TakeLast(10).Reverse().Select(p => $"{p.Key}({p.Value:+0.0;-0.0})")));
            sb.AppendLine();
            sb.AppendLine($"## SUGGEST 10 NEW {_peptideLength}-AA SEQUENCES");
            sb.AppendLine($"Use exactly {_peptideLength} UNIQUE amino acids from: {string.Join("", AvailableAminoAcids)}");
            return sb.ToString();
        }

        public string GenerateCircularClaudePrompt()
        {
            var sb = new StringBuilder();
            sb.AppendLine("# CIRCULAR PEPTIDE FRAGMENTATION OPTIMIZATION RESULTS");
            sb.AppendLine();
            sb.AppendLine("## GOAL");
            sb.AppendLine($"Design {_peptideLength}-AA CIRCULAR peptides that fragment evenly at ALL rotation points.");
            sb.AppendLine();
            sb.AppendLine("## CONSTRAINTS");
            sb.AppendLine($"- Exactly {_peptideLength} amino acids in a ring, each used only ONCE");
            sb.AppendLine($"- Available AAs: {string.Join(", ", AvailableAminoAcids)}");
            sb.AppendLine($"- Collision energy: {_collisionEnergy} eV, Charge: +{_precursorCharge}");
            sb.AppendLine();
            sb.AppendLine("## STATISTICS");
            sb.AppendLine($"- Circular peptides evaluated: {AllCircularEvaluatedPeptides.Count:N0}");
            sb.AppendLine($"- Total rotations tested: {AllCircularEvaluatedPeptides.Count * _peptideLength:N0}");
            sb.AppendLine($"- Generations: {BestCircularPerGeneration.Count}");
            sb.AppendLine($"- Final diversity: {_lastGenerationDiversity:P1}");
            sb.AppendLine();
            sb.AppendLine("## TOP 10 CIRCULAR PEPTIDES");
            sb.AppendLine("| Rank | Sequence | CircFit | MinRot | AvgRot | StdDev | Worst Rotation |");
            sb.AppendLine("|------|----------|---------|--------|--------|--------|----------------|");
            int rank = 1;
            foreach (var p in AllCircularEvaluatedPeptides.OrderByDescending(f => f.OverallCircularFitness).Take(10))
            {
                sb.AppendLine($"| {rank,4} | {p.CanonicalSequence} | {p.OverallCircularFitness:F4} | {p.MinRotationFitness:F4} | {p.AvgRotationFitness:F4} | {p.FitnessStdDev:F4} | {p.WorstRotation.Sequence} |");
                rank++;
            }
            sb.AppendLine();
            sb.AppendLine("## ROTATION BREAKDOWN (TOP 3)");
            foreach (var p in AllCircularEvaluatedPeptides.OrderByDescending(f => f.OverallCircularFitness).Take(3))
            {
                sb.AppendLine($"\n### {p.CanonicalSequence}");
                foreach (var rot in p.RotationFitnesses.OrderByDescending(r => r.OverallFitness))
                {
                    string marker = rot.Sequence == p.WorstRotation.Sequence ? " ← WORST" : "";
                    sb.AppendLine($"  {rot.Sequence} Fit={rot.OverallFitness:F4} B-CV={rot.BIonCV:F3} Y-CV={rot.YIonCV:F3}{marker}");
                }
            }
            sb.AppendLine();
            sb.AppendLine($"## SUGGEST 10 NEW {_peptideLength}-AA CIRCULAR SEQUENCES");
            sb.AppendLine($"Use exactly {_peptideLength} UNIQUE amino acids from: {string.Join("", AvailableAminoAcids)}");
            sb.AppendLine("Optimize for CONSISTENT fragmentation across ALL linearization points.");
            return sb.ToString();
        }

        #endregion
    }
}