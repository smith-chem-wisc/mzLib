using Chemistry;
using Microsoft.ML.OnnxRuntime;
using Microsoft.ML.OnnxRuntime.Tensors;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using PredictionClients.Koina.AbstractClasses;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

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
    ///     peptideSequences: sequences,
    ///     precursorCharges: charges,
    ///     retentionTimes: rts,
    ///     out var warnings,
    ///     onnxModelPath: "path/to/model.onnx");
    ///
    /// await model.RunInferenceAsync();
    /// var spectra = model.PredictedSpectra;
    ///
    /// MODEL DETAILS (v3, all-protease)
    /// ----------------------------------
    /// Input:  float32 [N, 18]
    /// Output: float32 [N]      — predicted TicNI
    /// Fragment lengths: 3–9 residues
    /// </summary>
    public class InternalFragmentIntensityModel : FragmentIntensityModel
    {
        #region Model Metadata

        public override string ModelName => "InternalFragmentScorer_v3_AllProteases";
        public override int MaxBatchSize => 50_000;
        public override int MaxPeptideLength => 100;
        public override int MinPeptideLength => 5;

        /// <summary>Default ONNX model filename.</summary>
        public const string DefaultModelFileName = "InternalFragmentScorer_v3_AllProteases.onnx";

        /// <summary>
        /// Resolves the default ONNX model path by searching common locations relative to the
        /// executing assembly.
        /// </summary>
        public static string DefaultOnnxModelPath
        {
            get
            {
                var assemblyDir = Path.GetDirectoryName(
                    typeof(InternalFragmentIntensityModel).Assembly.Location) ?? "";

                var candidates = new[]
                {
                    Path.Combine(assemblyDir, "LocalModels", DefaultModelFileName),
                    Path.Combine(assemblyDir, DefaultModelFileName),
                    Path.Combine("LocalModels", DefaultModelFileName),
                    DefaultModelFileName,
                };

                foreach (var p in candidates)
                    if (File.Exists(p)) return Path.GetFullPath(p);

                return Path.GetFullPath(candidates[0]);
            }
        }

        /// <summary>Charges 1–6 accepted.</summary>
        public override HashSet<int> AllowedPrecursorCharges => new() { 1, 2, 3, 4, 5, 6 };

        /// <summary>Modifications are stripped before feature computation.</summary>
        public override Dictionary<string, string> ValidModificationUnimodMapping => new();

        #endregion

        #region Fragment Enumeration Parameters

        public int MinFragmentLength { get; } = 3;
        public int MaxFragmentLength { get; } = 9;

        public override double MinIntensityFilter { get; protected set; } = 0.002;

        /// <summary>
        /// Maximum internal fragment ions per peptide in the spectral library.
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
        private static readonly Regex _modRegex = new(@"\[[^\]]+\]", RegexOptions.Compiled);

        #endregion

        #region Amino Acid Data

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

        #region Constructors

        /// <summary>
        /// Convenience constructor — resolves the ONNX model path automatically.
        /// Use when the model file is deployed alongside the assembly.
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
                   DefaultOnnxModelPath, spectralLibrarySavePath,
                   minIntensityFilter, maxInternalIonsPerPeptide)
        {
        }

        /// <summary>
        /// Full constructor with explicit ONNX model path.
        /// </summary>
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
                throw new ArgumentException("Input lists must have the same length.");

            _onnxModelPath = onnxModelPath;
            SpectralLibrarySavePath = spectralLibrarySavePath;
            MinIntensityFilter = minIntensityFilter;
            MaxInternalIonsPerPeptide = maxInternalIonsPerPeptide;

            if (peptideSequences.Count == 0)
            {
                warnings = new WarningException("Inputs were empty. No predictions will be made.");
                return;
            }

            var invalid = new List<string>();
            for (int i = 0; i < peptideSequences.Count; i++)
            {
                var seq = peptideSequences[i];
                var charge = precursorCharges[i];

                if (!IsValidSequence(seq) || !AllowedPrecursorCharges.Contains(charge))
                    invalid.Add($"Index {i}: '{seq}', charge: {charge}");
                else
                {
                    PeptideSequences.Add(seq);
                    PrecursorCharges.Add(charge);
                    RetentionTimes.Add(retentionTimes[i]);
                }
            }

            warnings = invalid.Count > 0
                ? new WarningException(
                    "The following entries are invalid and will be skipped:\n"
                    + string.Join("\n", invalid)
                    + $"\nRequirements: bare length {MinPeptideLength}–{MaxPeptideLength}, "
                    + $"charge in [{string.Join(", ", AllowedPrecursorCharges)}]")
                : null;
        }

        #endregion

        #region Inference

        /// <summary>
        /// Runs ONNX inference locally (CPU), then generates PredictedSpectra.
        /// </summary>
        public override async Task<WarningException?> RunInferenceAsync()
        {
            if (_disposed)
                throw new ObjectDisposedException(nameof(InternalFragmentIntensityModel),
                    "Cannot run inference on a disposed model instance.");

            if (PeptideSequences.Count == 0)
            {
                Dispose();
                return null;
            }

            await Task.Run(RunOnnxInference);

            WarningException? warning;
            if (SpectralLibrarySavePath is not null)
                SavePredictedSpectralLibrary(SpectralLibrarySavePath, out warning);
            else
                GenerateLibrarySpectraFromPredictions(out warning);

            Dispose();
            return warning;
        }

        private void RunOnnxInference()
        {
            using var session = new InferenceSession(_onnxModelPath);

            for (int pepIdx = 0; pepIdx < PeptideSequences.Count; pepIdx++)
            {
                var fullSeq = PeptideSequences[pepIdx];
                var charge = PrecursorCharges[pepIdx];
                var bare = BareSequence(fullSeq);
                int n = bare.Length;

                var candidates = new List<(int i, int j, double neutralMass)>();
                for (int i = 0; i < n; i++)
                    for (int j = i + MinFragmentLength - 1; j < Math.Min(i + MaxFragmentLength, n); j++)
                        candidates.Add((i, j, InternalNeutralMass(bare, i, j)));

                if (candidates.Count == 0)
                {
                    Predictions.Add(new PeptideFragmentIntensityPrediction(
                        fullSeq, charge, new(), new(), new()));
                    continue;
                }

                var featureData = new float[candidates.Count * 18];
                for (int c = 0; c < candidates.Count; c++)
                {
                    var features = ComputeFeatures(bare, candidates[c].i, candidates[c].j, charge);
                    for (int f = 0; f < 18; f++)
                        featureData[c * 18 + f] = features[f];
                }

                var inputTensor = new DenseTensor<float>(featureData, new[] { candidates.Count, 18 });
                using var results = session.Run(new[]
                {
                    NamedOnnxValue.CreateFromTensor("X", inputTensor)
                });
                var predicted = results.First().AsEnumerable<float>().ToArray();

                var annotations = new List<string>(candidates.Count);
                var mzs = new List<double>(candidates.Count);
                var intensities = new List<double>(candidates.Count);

                for (int c = 0; c < candidates.Count; c++)
                {
                    var (ci, cj, neut) = candidates[c];
                    annotations.Add($"bIb[{ci + 1}-{cj + 1}]+1");
                    mzs.Add(neut + Proton);
                    intensities.Add((double)predicted[c]);
                }

                Predictions.Add(new PeptideFragmentIntensityPrediction(
                    fullSeq, charge, annotations, mzs, intensities));
            }
        }

        protected override List<Dictionary<string, object>> ToBatchedRequests() => new();
        protected override void ResponseToPredictions(string[] response) { }

        #endregion

        #region Spectral Library Generation

        /// <summary>
        /// Overrides the base to parse internal fragment annotations (bIb[start-end]+1),
        /// construct Products with SecondaryProductType set (IsInternalFragment = true),
        /// apply the top-N filter, and normalise intensities to max = 1.0 per peptide.
        /// </summary>
        protected override void GenerateLibrarySpectraFromPredictions(out WarningException? warning)
        {
            warning = null;
            if (Predictions.Count == 0) return;

            var regex = new Regex(
                @"^([a-zA-Z]+)I([a-zA-Z]+)\[(\d+)-(\d+)\]\+(\d+)$",
                RegexOptions.Compiled);

            for (int predIdx = 0; predIdx < Predictions.Count; predIdx++)
            {
                var prediction = Predictions[predIdx];

                var passing = new List<(int idx, double intensity)>();
                for (int f = 0; f < prediction.FragmentAnnotations.Count; f++)
                    if (prediction.FragmentIntensities[f] >= MinIntensityFilter)
                        passing.Add((f, prediction.FragmentIntensities[f]));

                var topN = passing
                    .OrderByDescending(p => p.intensity)
                    .Take(MaxInternalIonsPerPeptide)
                    .ToList();

                if (topN.Count == 0) continue;

                double maxIntensity = topN.Max(p => p.intensity);
                var fragmentIons = new List<MatchedFragmentIon>(topN.Count);

                for (int t = 0; t < topN.Count; t++)
                {
                    int idx = topN[t].idx;
                    double rawIntensity = topN[t].intensity;

                    var annotation = prediction.FragmentAnnotations[idx];
                    var match = regex.Match(annotation);
                    if (!match.Success) continue;

                    var primaryType = Enum.Parse<ProductType>(match.Groups[1].Value, ignoreCase: true);
                    var secondaryType = Enum.Parse<ProductType>(match.Groups[2].Value, ignoreCase: true);
                    int startResidue = int.Parse(match.Groups[3].Value);
                    int endResidue = int.Parse(match.Groups[4].Value);
                    int fragCharge = int.Parse(match.Groups[5].Value);

                    var product = new Product(
                        productType: primaryType,
                        terminus: FragmentationTerminus.None,
                        neutralMass: prediction.FragmentMZs[idx].ToMass(fragCharge),
                        fragmentNumber: startResidue,
                        residuePosition: startResidue,
                        neutralLoss: 0,
                        secondaryProductType: secondaryType,
                        secondaryFragmentNumber: endResidue);

                    fragmentIons.Add(new MatchedFragmentIon(
                        product,
                        experMz: prediction.FragmentMZs[idx],
                        experIntensity: rawIntensity / maxIntensity,
                        charge: fragCharge));
                }

                if (fragmentIons.Count == 0) continue;

                var bare = BareSequence(prediction.FullSequence);
                double pepMass = bare.Sum(aa => AaMass.GetValueOrDefault(aa, 111.0)) + Water;
                double precursorMz = (pepMass + prediction.PrecursorCharge * Proton)
                                     / prediction.PrecursorCharge;

                PredictedSpectra.Add(new LibrarySpectrum(
                    sequence: prediction.FullSequence,
                    precursorMz: precursorMz,
                    chargeState: prediction.PrecursorCharge,
                    peaks: fragmentIons,
                    rt: RetentionTimes[predIdx]));
            }

            var unique = PredictedSpectra.DistinctBy(s => s.Name).ToList();
            if (unique.Count != PredictedSpectra.Count)
            {
                warning = new WarningException(
                    $"Duplicate spectra reduced from {PredictedSpectra.Count} to {unique.Count}.");
                PredictedSpectra = unique;
            }
        }

        #endregion

        #region Utilities

        public static string BareSequence(string fullSeq)
            => _modRegex.Replace(fullSeq, string.Empty);

        protected override bool IsValidSequence(string sequence)
        {
            var bare = BareSequence(sequence);
            return bare.Length >= MinPeptideLength
                && bare.Length <= MaxPeptideLength
                && bare.All(aa => AaMass.ContainsKey(aa));
        }

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
            double hydro = i > 0 ? Hydrophobicity.GetValueOrDefault(seq[i - 1], 0.0) : 0.0;
            char fn = seq[i];
            char ln = seq[j];

            return new float[]
            {
                (float)rank,                          // [0]  LocalIntensityRank
                (float)charge,                        // [1]  PrecursorCharge
                (float)n,                             // [2]  PeptideLength
                (float)fragLen,                       // [3]  FragmentLength
                (float)relDist,                       // [4]  RelativeDistanceFromCTerm
                (float)basicY,                        // [5]  BasicResiduesInYIonSpan
                (float)basicB,                        // [6]  BasicResiduesInBIonSpan
                (float)nBasic,                        // [7]  NumberOfBasicResidues
                (float)maxT,                          // [8]  MaxTerminalIonIntensity
                (float)yAtC,                          // [9]  YIonIntensityAtCTerm
                (float)bAtN,                          // [10] BIonIntensityAtNTerm
                (float)hydro,                         // [11] NTerminalFlankingHydrophobicity
                (float)(bAtN * yAtC),                 // [12] BYProductScore
                fn == 'P' ? 1f : 0f,                 // [13] IsProlineAtInternalNTerminus
                (fn == 'P' || ln == 'P') ? 1f : 0f, // [14] HasProlineAtEitherTerminus
                (fn == 'D' || ln == 'D') ? 1f : 0f, // [15] HasAspartateAtEitherTerminus
                bAtN > 0 ? 1f : 0f,                  // [16] HasMatchedBIonAtNTerm
                yAtC > 0 ? 1f : 0f,                  // [17] HasMatchedYIonAtCTerm
            };
        }

        private static double InternalNeutralMass(string seq, int i, int j)
            => seq[i..(j + 1)].Sum(aa => AaMass.GetValueOrDefault(aa, 111.0));

        #endregion
    }
}
