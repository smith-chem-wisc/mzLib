using Chemistry;
using Microsoft.ML.OnnxRuntime;
using Microsoft.ML.OnnxRuntime.Tensors;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using PredictionClients.Koina.AbstractClasses;
using System.ComponentModel;
using System.Text.RegularExpressions;

namespace PredictionClients.LocalModels
{
    /// <summary>
    /// Local ONNX-based model for predicting internal fragment ion intensities in HCD spectra.
    /// Mirrors the Koina/Prosit API surface so it can be used as a drop-in replacement
    /// alongside remote models — same constructor pattern, same output types, same library
    /// generation pipeline.
    ///
    /// The model predicts TIC-normalized intensity (TicNI) for internal fragment ion candidates
    /// derived from double-backbone fragmentation events. It was trained on six proteases
    /// (LysC, Trypsin, ArgC, GluC, AspN, Chymotrypsin) and is protease-agnostic.
    ///
    /// USAGE
    /// -----
    /// var model = new InternalFragmentIntensityModel(
    ///     peptideSequences: sequences,      // bare or mzLib-modified sequences
    ///     precursorCharges: charges,
    ///     retentionTimes: rts,              // optional, used only for library metadata
    ///     out var warnings,
    ///     onnxModelPath: "path/to/model.onnx");
    ///
    /// await model.RunInferenceAsync();
    ///
    /// var spectra = model.PredictedSpectra; // List&lt;LibrarySpectrum&gt;, same type as Prosit
    ///
    /// DIFFERENCES FROM PROSIT
    /// -----------------------
    /// - No collision energy input (internal fragment model does not use NCE)
    /// - Inference is local ONNX Runtime, not HTTP to Koina
    /// - RunInferenceAsync() is CPU-bound but returns Task for API compatibility
    /// - Sequences may contain mzLib modification brackets; the bare AA sequence
    ///   is extracted automatically for feature computation
    /// - Internal fragment ions in PredictedSpectra use IsInternalFragment = true
    ///   and are annotated as bIb[start-end] (mzLib Product convention)
    ///
    /// MODEL DETAILS (v3, all-protease)
    /// ----------------------------------
    /// Input:  float32 [N, 18]  — 18 hand-crafted chemical features per candidate
    /// Output: float32 [N]      — predicted TicNI directly (no back-transform needed)
    /// Fragment lengths: 3–9 residues (beyond 9 model performance degrades)
    /// </summary>
    public class InternalFragmentIntensityModel : FragmentIntensityModel
    {
        #region Model Metadata

        public override string ModelName => "InternalFragmentScorer_v3_AllProteases";
        public override int MaxBatchSize => 50_000;
        public override int MaxPeptideLength => 100;
        public override int MinPeptideLength => 5; // need at least one internal fragment of length 3

        /// <summary>Default ONNX model filename.</summary>
        public const string DefaultModelFileName = "InternalFragmentScorer_v3_AllProteases.onnx";

        /// <summary>
        /// Resolves the default ONNX model path by searching common locations.
        /// </summary>
        public static string DefaultOnnxModelPath
        {
            get
            {
                // Check relative to executing assembly
                var assemblyDir = Path.GetDirectoryName(typeof(InternalFragmentIntensityModel).Assembly.Location) ?? "";
                
                var paths = new[]
                {
                    Path.Combine(assemblyDir, "LocalModels", DefaultModelFileName),
                    Path.Combine(assemblyDir, DefaultModelFileName),
                    Path.Combine("LocalModels", DefaultModelFileName),
                    DefaultModelFileName
                };

                foreach (var path in paths)
                {
                    if (File.Exists(path))
                        return Path.GetFullPath(path);
                }

                // Return default path even if not found (will throw at runtime)
                return Path.Combine(assemblyDir, "LocalModels", DefaultModelFileName);
            }
        }

        /// <summary>Charges 1–6 accepted; the charge feature is passed directly to the model.</summary>
        public override HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };

        /// <summary>
        /// All standard mzLib modifications are accepted — they are stripped from the
        /// sequence before feature computation and do not affect internal fragment mass predictions.
        /// Override ValidModificationUnimodMapping is unused here but present for base-class compat.
        /// </summary>
        public override Dictionary<string, string> ValidModificationUnimodMapping => new();

        #endregion

        #region Fragment Enumeration Parameters

        /// <summary>Minimum internal fragment length (residues). Below 3, specificity is low.</summary>
        public int MinFragmentLength { get; } = 3;

        /// <summary>Maximum internal fragment length (residues). Beyond 9, model performance degrades.</summary>
        public int MaxFragmentLength { get; } = 9;

        /// <summary>
        /// Minimum predicted TicNI for a fragment to be included in PredictedSpectra.
        /// Corresponds to approximately the 30th percentile of observed internal fragment intensities.
        /// </summary>
        public override double MinIntensityFilter { get; protected set; } = 0.002;

        /// <summary>
        /// Maximum internal fragment ions per peptide included in the spectral library.
        /// Applied after MinIntensityFilter, keeping the top N by predicted TicNI.
        /// </summary>
        public int MaxInternalIonsPerPeptide { get; private set; } = 20;

        #endregion

        #region Input/Output Data

        public override List<string> PeptideSequences { get; } = new();
        public override List<int> PrecursorCharges { get; } = new();
        public override List<double?> RetentionTimes { get; } = new();
        public override string? SpectralLibrarySavePath { get; }
        public override List<PeptideFragmentIntensityPrediction> Predictions { get; protected set; } = new();
        public override List<LibrarySpectrum> PredictedSpectra { get; protected set; } = new();

        #endregion

        #region Private State

        private readonly string _onnxModelPath;

        // Strips [modification brackets] to get the bare amino acid sequence.
        private static readonly Regex _modRegex = new(@"\[[^\]]+\]", RegexOptions.Compiled);

        #endregion

        #region Amino Acid Data

        // Monoisotopic residue masses (Da)
        private static readonly Dictionary<char, double> AaMass = new()
        {
            ['A'] = 71.03711,
            ['R'] = 156.10111,
            ['N'] = 114.04293,
            ['D'] = 115.02694,
            ['C'] = 103.00919,
            ['E'] = 129.04259,
            ['Q'] = 128.05858,
            ['G'] = 57.02146,
            ['H'] = 137.05891,
            ['I'] = 113.08406,
            ['L'] = 113.08406,
            ['K'] = 128.09496,
            ['M'] = 131.04049,
            ['F'] = 147.06841,
            ['P'] = 97.05276,
            ['S'] = 87.03203,
            ['T'] = 101.04768,
            ['W'] = 186.07931,
            ['Y'] = 163.06333,
            ['V'] = 99.06841,
        };

        private const double Proton = 1.007276;
        private const double Water = 18.010565;

        private static readonly HashSet<char> BasicResidues = new() { 'K', 'R', 'H' };
        private static readonly Dictionary<char, double> Hydrophobicity = new()
        {
            ['A'] = 1.8,
            ['V'] = 4.2,
            ['I'] = 4.5,
            ['L'] = 3.8,
            ['M'] = 1.9,
            ['F'] = 2.8,
            ['W'] = -0.9,
            ['P'] = -1.6,
            ['G'] = -0.4,
            ['S'] = -0.8,
            ['T'] = -0.7,
            ['C'] = 2.5,
            ['Y'] = -1.3,
            ['H'] = -3.2,
            ['K'] = -3.9,
            ['R'] = -4.5,
            ['E'] = -3.5,
            ['D'] = -3.5,
            ['N'] = -3.5,
            ['Q'] = -3.5,
        };

        #endregion

        #region Constructor

        /// <summary>
        /// Initializes the model using the default ONNX model path.
        /// </summary>
        public InternalFragmentIntensityModel(
            List<string> peptideSequences,
            List<int> precursorCharges,
            List<double?> retentionTimes,
            out WarningException? warnings,
            string? spectralLibrarySavePath = null,
            double minIntensityFilter = 0.002,
            int maxInternalIonsPerPeptide = 20)
            : this(peptideSequences, precursorCharges, retentionTimes, out warnings,
                   DefaultOnnxModelPath, spectralLibrarySavePath, minIntensityFilter, maxInternalIonsPerPeptide)
        {
        }

        /// <summary>
        /// Initializes the model and validates inputs. Invalid entries are filtered out
        /// with details recorded in the warnings output parameter.
        /// </summary>
        /// <param name="peptideSequences">
        /// Peptide sequences, optionally including mzLib modification brackets.
        /// All 20 canonical amino acids are accepted; modification brackets are stripped
        /// before feature computation. Example: "PEPTM[Common Variable:Oxidation on M]IDER"
        /// </param>
        /// <param name="precursorCharges">Precursor charge states (1–6). Must be same length as peptideSequences.</param>
        /// <param name="retentionTimes">
        /// Optional retention times for spectral library metadata. May contain nulls.
        /// Must be same length as peptideSequences.
        /// </param>
        /// <param name="warnings">Details of any filtered-out entries; null if all inputs are valid.</param>
        /// <param name="onnxModelPath">Path to InternalFragmentScorer_v3_AllProteases.onnx</param>
        /// <param name="spectralLibrarySavePath">If set, library is written to disk after inference.</param>
        /// <param name="minIntensityFilter">Minimum predicted TicNI to include a fragment. Default: 0.002.</param>
        /// <param name="maxInternalIonsPerPeptide">Top-N internal ions per peptide. Default: 20.</param>
        /// <exception cref="ArgumentException">Thrown when input list lengths differ.</exception>
        public InternalFragmentIntensityModel(
            List<string> peptideSequences,
            List<int> precursorCharges,
            List<double?> retentionTimes,
            out WarningException? warnings,
            string onnxModelPath,
            string? spectralLibrarySavePath = null,
            double minIntensityFilter = 0.002,
            int maxInternalIonsPerPeptide = 20)
        {
            if (peptideSequences.Count != precursorCharges.Count
                || precursorCharges.Count != retentionTimes.Count)
            {
                throw new ArgumentException("Input lists must have the same length.");
            }

            _onnxModelPath = onnxModelPath;
            SpectralLibrarySavePath = spectralLibrarySavePath;
            MinIntensityFilter = minIntensityFilter;
            MaxInternalIonsPerPeptide = maxInternalIonsPerPeptide;

            if (peptideSequences.Count == 0)
            {
                warnings = new WarningException("Inputs were empty. No predictions will be made.");
                return;
            }

            var invalidEntries = new List<string>();
            for (int i = 0; i < peptideSequences.Count; i++)
            {
                var seq = peptideSequences[i];
                var charge = precursorCharges[i];
                var bare = BareSequence(seq);

                if (!IsValidSequence(seq) || !AllowedPrecursorCharges.Contains(charge))
                {
                    invalidEntries.Add(
                        $"Index {i}: '{seq}' (bare length: {bare.Length}), charge: {charge}");
                }
                else
                {
                    PeptideSequences.Add(seq);
                    PrecursorCharges.Add(charge);
                    RetentionTimes.Add(retentionTimes[i]);
                }
            }

            warnings = null;
            if (invalidEntries.Count > 0)
            {
                warnings = new WarningException(
                    "The following entries are invalid and will be skipped:\n"
                    + string.Join("\n", invalidEntries)
                    + $"\nRequirements: bare sequence length {MinPeptideLength}–{MaxPeptideLength}, "
                    + $"charge in [{string.Join(", ", AllowedPrecursorCharges)}]");
            }
        }


        #endregion

        #region Inference

        /// <summary>
        /// Runs ONNX inference locally (CPU), then generates PredictedSpectra.
        /// Returns Task for API compatibility with the Koina async interface.
        /// </summary>
        public override async Task<WarningException?> RunInferenceAsync()
        {
            if (_disposed)
            {
                throw new ObjectDisposedException(
                    nameof(InternalFragmentIntensityModel),
                    "Cannot run inference on a disposed model instance. Results are still accessible.");
            }

            if (PeptideSequences.Count == 0)
            {
                Dispose();
                return null;
            }

            // Wrap CPU-bound work in Task.Run so callers can await without blocking their thread.
            await Task.Run(RunOnnxInference);

            WarningException? warning;
            if (SpectralLibrarySavePath is not null)
                SavePredictedSpectralLibrary(SpectralLibrarySavePath, out warning);
            else
                GenerateLibrarySpectraFromPredictions(out warning);

            Dispose();
            return warning;
        }

        /// <summary>
        /// Core ONNX inference loop. Enumerates all internal fragment candidates for every
        /// peptide, batches them through the ONNX session, and assembles Predictions.
        /// </summary>
        private void RunOnnxInference()
        {
            using var session = new InferenceSession(_onnxModelPath);

            for (int pepIdx = 0; pepIdx < PeptideSequences.Count; pepIdx++)
            {
                var fullSeq = PeptideSequences[pepIdx];
                var charge = PrecursorCharges[pepIdx];
                var bare = BareSequence(fullSeq);
                var n = bare.Length;

                // Enumerate all (i, j) internal fragment candidates
                var candidates = new List<(int i, int j, double neutralMass)>();
                for (int i = 0; i < n; i++)
                {
                    for (int j = i + MinFragmentLength - 1;
                         j < Math.Min(i + MaxFragmentLength, n); j++)
                    {
                        candidates.Add((i, j, InternalNeutralMass(bare, i, j)));
                    }
                }

                if (candidates.Count == 0)
                {
                    // Peptide too short to produce any internal fragments — emit empty prediction
                    Predictions.Add(new PeptideFragmentIntensityPrediction(
                        fullSeq, charge,
                        new List<string>(),
                        new List<double>(),
                        new List<double>()));
                    continue;
                }

                // Build feature matrix [nCandidates × 18]
                var featureData = new float[candidates.Count * 18];
                for (int c = 0; c < candidates.Count; c++)
                {
                    var (ci, cj, _) = candidates[c];
                    var features = ComputeFeatures(bare, ci, cj, charge);
                    for (int f = 0; f < 18; f++)
                        featureData[c * 18 + f] = features[f];
                }

                // ONNX inference
                var inputTensor = new DenseTensor<float>(featureData, new[] { candidates.Count, 18 });
                using var results = session.Run(new[]
                {
                    NamedOnnxValue.CreateFromTensor("X", inputTensor)
                });
                var predicted = results.First().AsEnumerable<float>().ToArray();

                // Convert to PeptideFragmentIntensityPrediction format
                // Annotation: "bIb[start-end]+1"  (mzLib internal ion convention, '+' for charge)
                var annotations = new List<string>(candidates.Count);
                var mzs = new List<double>(candidates.Count);
                var intensities = new List<double>(candidates.Count);

                for (int c = 0; c < candidates.Count; c++)
                {
                    var (ci, cj, neut) = candidates[c];
                    var predTicNI = (double)predicted[c];
                    annotations.Add($"bIb[{ci + 1}-{cj + 1}]+1");
                    mzs.Add((neut + Proton) / 1.0);   // charge 1 for internal ions
                    intensities.Add(predTicNI);
                }

                Predictions.Add(new PeptideFragmentIntensityPrediction(
                    fullSeq, charge, annotations, mzs, intensities));
            }
        }

        /// <summary>No batched HTTP requests — override is a no-op (ONNX batching handled internally).</summary>
        protected override List<Dictionary<string, object>> ToBatchedRequests()
            => new(); // unused; RunInferenceAsync calls RunOnnxInference directly

        /// <summary>No HTTP responses — override is a no-op.</summary>
        protected override void ResponseToPredictions(string[] response) { }

        #endregion

        #region Spectral Library Generation (override for internal ions)

        /// <summary>
        /// Converts Predictions into LibrarySpectrum objects.
        /// Overrides the base implementation to handle internal fragment ion annotations
        /// (bIb[start-end]+1) and to apply the top-N intensity filter.
        ///
        /// Compared to the base class:
        /// - Parses "bIb[3-6]+1" instead of "b5+1"
        /// - Constructs Product with SecondaryProductType set (IsInternalFragment = true)
        /// - Applies MaxInternalIonsPerPeptide limit per peptide
        /// - Normalizes intensities to max = 1.0 within each peptide
        /// </summary>
        protected override void GenerateLibrarySpectraFromPredictions(out WarningException? warning)
        {
            warning = null;
            if (Predictions.Count == 0) return;

            // Regex for "bIb[3-6]+1"
            var internalAnnotationRegex = new Regex(
                @"^([a-zA-Z]+)I([a-zA-Z]+)\[(\d+)-(\d+)\]\+(\d+)$", RegexOptions.Compiled);

            for (int predIdx = 0; predIdx < Predictions.Count; predIdx++)
            {
                var prediction = Predictions[predIdx];
                var bare = BareSequence(prediction.FullSequence);

                // Collect candidates that pass the minimum intensity threshold
                var passing = new List<(int annotIdx, double intensity)>();
                for (int f = 0; f < prediction.FragmentAnnotations.Count; f++)
                {
                    if (prediction.FragmentIntensities[f] >= MinIntensityFilter)
                        passing.Add((f, prediction.FragmentIntensities[f]));
                }

                // Apply top-N limit
                var topN = passing
                    .OrderByDescending(p => p.intensity)
                    .Take(MaxInternalIonsPerPeptide)
                    .ToList();

                if (topN.Count == 0)
                    continue;

                // Normalize intensities to max = 1.0 within this peptide
                double maxIntensity = topN.Max(p => p.intensity);

                var fragmentIons = new List<MatchedFragmentIon>(topN.Count);
                foreach (var (annotIdx, rawIntensity) in topN)
                {
                    var annotation = prediction.FragmentAnnotations[annotIdx];
                    var match = internalAnnotationRegex.Match(annotation);
                    if (!match.Success) continue;

                    var primaryType = Enum.Parse<ProductType>(match.Groups[1].Value, ignoreCase: true);
                    var secondaryType = Enum.Parse<ProductType>(match.Groups[2].Value, ignoreCase: true);
                    int startResidue = int.Parse(match.Groups[3].Value);
                    int endResidue = int.Parse(match.Groups[4].Value);
                    int fragCharge = int.Parse(match.Groups[5].Value);

                    var product = new Product(
                        productType: primaryType,
                        terminus: FragmentationTerminus.None,
                        neutralMass: prediction.FragmentMZs[annotIdx].ToMass(fragCharge),
                        fragmentNumber: startResidue,
                        residuePosition: startResidue,
                        neutralLoss: 0,
                        secondaryProductType: secondaryType,
                        secondaryFragmentNumber: endResidue);

                    fragmentIons.Add(new MatchedFragmentIon(
                        product,
                        experMz: prediction.FragmentMZs[annotIdx],
                        experIntensity: rawIntensity / maxIntensity,
                        charge: fragCharge));
                }

                if (fragmentIons.Count == 0) continue;

                // Compute precursor m/z from bare sequence
                double neutralPepMass = bare.Sum(aa => AaMass.GetValueOrDefault(aa, 111.0)) + Water;
                double precursorMz = (neutralPepMass + prediction.PrecursorCharge * Proton) / prediction.PrecursorCharge;

                PredictedSpectra.Add(new LibrarySpectrum(
                    sequence: prediction.FullSequence,
                    precursorMz: precursorMz,
                    chargeState: prediction.PrecursorCharge,
                    peaks: fragmentIons,
                    rt: RetentionTimes[predIdx]));
            }

            // Deduplicate by Name (Sequence/Charge)
            var unique = PredictedSpectra.DistinctBy(s => s.Name).ToList();
            if (unique.Count != PredictedSpectra.Count)
            {
                warning = new WarningException(
                    $"Duplicate spectra found. Reduced from {PredictedSpectra.Count} to {unique.Count} unique spectra.");
                PredictedSpectra = unique;
            }
        }

        #endregion

        #region Sequence and Feature Utilities

        /// <summary>Strips mzLib modification brackets to get the bare amino acid sequence.</summary>
        public static string BareSequence(string fullSeq)
            => _modRegex.Replace(fullSeq, string.Empty);

        /// <inheritdoc/>
        protected override bool IsValidSequence(string sequence)
        {
            var bare = BareSequence(sequence);
            return bare.Length >= MinPeptideLength
                && bare.Length <= MaxPeptideLength
                && bare.All(aa => AaMass.ContainsKey(aa));
        }

        /// <summary>
        /// Computes all 18 features for an internal fragment candidate.
        /// Mirrors InternalFragmentFeatureExtractor.cs exactly.
        ///
        /// In library generation mode (no observed spectrum), terminal ion intensity
        /// features are 0 and LocalIntensityRank = 5 (conservative mid-range).
        /// </summary>
        private static float[] ComputeFeatures(
            string seq, int i, int j, int charge,
            Dictionary<int, double>? bIntens = null,
            Dictionary<int, double>? yIntens = null,
            int rank = 5)
        {
            bIntens ??= new();
            yIntens ??= new();

            int n = seq.Length;
            int fragLen = j - i + 1;
            double relDist = (n - 1 - j) / (double)n;
            int basicY = seq[(j + 1)..].Count(aa => BasicResidues.Contains(aa));
            int basicB = seq[..i].Count(aa => BasicResidues.Contains(aa));
            int nBasic = seq.Count(aa => BasicResidues.Contains(aa));
            double bAtN = bIntens.GetValueOrDefault(i, 0.0);
            double yAtC = yIntens.GetValueOrDefault(n - 1 - j, 0.0);
            double maxT = bIntens.Values.Concat(yIntens.Values).DefaultIfEmpty(0.0).Max();
            double hydrophob = i > 0 ? Hydrophobicity.GetValueOrDefault(seq[i - 1], 0.0) : 0.0;
            char fn = seq[i];
            char ln = seq[j];

            return new float[]
            {
                (float)rank,                                    // [0]  LocalIntensityRank
                (float)charge,                                  // [1]  PrecursorCharge
                (float)n,                                       // [2]  PeptideLength
                (float)fragLen,                                 // [3]  FragmentLength
                (float)relDist,                                 // [4]  RelativeDistanceFromCTerm
                (float)basicY,                                  // [5]  BasicResiduesInYIonSpan
                (float)basicB,                                  // [6]  BasicResiduesInBIonSpan
                (float)nBasic,                                  // [7]  NumberOfBasicResidues
                (float)maxT,                                    // [8]  MaxTerminalIonIntensity
                (float)yAtC,                                    // [9]  YIonIntensityAtCTerm
                (float)bAtN,                                    // [10] BIonIntensityAtNTerm
                (float)hydrophob,                               // [11] NTerminalFlankingHydrophobicity
                (float)(bAtN * yAtC),                           // [12] BYProductScore
                fn == 'P' ? 1f : 0f,                           // [13] IsProlineAtInternalNTerminus
                (fn == 'P' || ln == 'P') ? 1f : 0f,           // [14] HasProlineAtEitherTerminus
                (fn == 'D' || ln == 'D') ? 1f : 0f,           // [15] HasAspartateAtEitherTerminus
                (float)(bAtN > 0 ? 1 : 0),                    // [16] HasMatchedBIonAtNTerm
                (float)(yAtC > 0 ? 1 : 0),                    // [17] HasMatchedYIonAtCTerm
            };
        }

        /// <summary>Neutral mass of internal fragment seq[i..j] (b-type, no water added).</summary>
        private static double InternalNeutralMass(string seq, int i, int j)
            => seq[i..(j + 1)].Sum(aa => AaMass.GetValueOrDefault(aa, 111.0));

        #endregion
    }
}
