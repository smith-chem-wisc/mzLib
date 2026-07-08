#nullable enable
using MassSpectrometry;
using MzLibUtil;

namespace Readers
{
    /// <summary>
    /// DeconvolutionParameters that wrap a pre-computed list of MS1 features loaded
    /// from a feature file (FlashDeconv / TopFD <c>.ms1.feature</c>, Dinosaur
    /// <c>.feature.tsv</c>, ...). The format-specific reader is auto-detected from
    /// the file's extension via the Readers <see cref="FileReader"/> /
    /// <see cref="SupportedFileType"/> infrastructure.
    /// </summary>
    /// <remarks>
    /// Living in <c>Readers</c> (not <c>MassSpectrometry</c>) means this class can
    /// drive the format-detection machinery directly. The cost is that
    /// <c>MassSpectrometry</c>-side code can't reference this type by name —
    /// <see cref="Deconvoluter"/> discovers the algorithm polymorphically through
    /// <see cref="DeconvolutionParameters.CreateAlgorithm"/>.
    ///
    /// Features are stored sorted by m/z and indexed by a parallel <c>double[]</c>
    /// of m/z keys, allowing the algorithm to binary-search for the first feature
    /// at or above an m/z lower bound and then iterate forward until the upper
    /// bound — turning the per-MS2 query from O(N) to O(log N + k).
    ///
    /// Pairing the loaded features with MS2 scans is delegated to
    /// <see cref="FromFileDeconvolutionAlgorithm"/>.
    /// </remarks>
    public class FromFileDeconvolutionParameters : DeconvolutionParameters
    {
        public override DeconvolutionType DeconvolutionType { get; protected set; } = DeconvolutionType.FromFile;

        #region Lazy Loading Feature Files

        private sealed class FeatureCache(List<ISingleChargeMs1Feature> featuresMzAscending, double[] mzKeys)
        {
            public List<ISingleChargeMs1Feature> FeaturesMzAscending { get; } = featuresMzAscending;
            public double[] MzKeys { get; } = mzKeys;
        }

        private readonly object _featureCacheLock = new();

        private volatile FeatureCache? _featureCache;

        private FeatureCache GetOrLoadFeatureCache()
        {
            if (_featureCache is not null)
                return _featureCache;

            lock (_featureCacheLock)
            {
                if (_featureCache is not null)
                    return _featureCache;

                return _featureCache = LoadFeaturesFromFile();
            }
        }

        private FeatureCache LoadFeaturesFromFile()
        {
            if (string.IsNullOrWhiteSpace(FilePath))
                throw new InvalidOperationException(
                    "FromFileDeconvolutionParameters requires FilePath to be set before features can be loaded.");

            var resultFile = FileReader.ReadResultFile(FilePath);
            if (resultFile is not IMs1FeatureFile featureFile)
                throw new MzLibException(
                    $"File at '{FilePath}' is not a recognized MS1 feature file. " +
                    $"Detected type: {resultFile.GetType().Name}. " +
                    "Expected an Ms1Feature (.ms1.feature) or Dinosaur (.feature.tsv) file.");

            resultFile.LoadResults();
            return BuildFeatureCache(NormaliseRetentionTimes(featureFile.GetMs1Features()));
        }

        private static FeatureCache BuildFeatureCache(IEnumerable<ISingleChargeMs1Feature> features)
        {
            var sorted = features.OrderBy(f => f.Mz).ToList();
            var keys = new double[sorted.Count];
            for (int i = 0; i < sorted.Count; i++) keys[i] = sorted[i].Mz;
            return new FeatureCache(sorted, keys);
        }

        #endregion

        /// <summary>
        /// The pre-loaded per-charge features the algorithm filters against, sorted
        /// by ascending <see cref="ISingleChargeMs1Feature.Mz"/>.
        /// </summary>
        public IReadOnlyList<ISingleChargeMs1Feature> Features => GetOrLoadFeatureCache().FeaturesMzAscending;

        private string? _filePath;
        public string? FilePath
        {
            get => _filePath;
            set
            {
                if (string.Equals(_filePath, value, StringComparison.Ordinal))
                    return;

                _filePath = value;
                _featureCache = null;
            }
        }

        /// <summary>
        /// True if the feature file's retention times were detected as seconds at load
        /// (max RetentionTimeEnd &gt; 500) and divided by 60 to minutes. Lets callers
        /// surface that a unit conversion happened instead of it being silent.
        /// </summary>
        public bool RetentionTimeNormalizedFromSeconds { get; private set; }

        /// <summary>
        /// Constructs from a feature-file path. Reader is auto-detected from the
        /// extension. The file's per-charge expansion is materialized lazily on
        /// first use; this object is intended to be reused across many MS2 scans
        /// rather than reconstructed per call.
        /// </summary>
        public FromFileDeconvolutionParameters(string filePath, int minCharge, int maxCharge,
            Polarity polarity = Polarity.Positive)
            : base(minCharge, maxCharge, polarity)
        {
            if (filePath is null) throw new ArgumentNullException(nameof(filePath));
            FilePath = filePath;
        }

        /// <summary>
        /// Detects RT-in-seconds inputs (FlashDeconv / TopFD canonical convention) and
        /// converts them to minutes so the algorithm's RT comparisons against
        /// <see cref="MsDataScan.RetentionTime"/> (always minutes in mzML / Thermo's API)
        /// line up. Heuristic: if the largest <c>RetentionTimeEnd</c> in the file
        /// exceeds 500, it can't be minutes -- no realistic LC run is over 8 hours
        /// long -- so the file is in seconds and every RT value is divided by 60.
        /// Files already in minutes pass through unchanged.
        ///
        /// Without this normalisation, FlashDeconv / TopFD files that emit RT in
        /// seconds (e.g. <c>Time_begin = 2787</c> for a 46-minute LC run) silently
        /// produce zero matching features at search time, because the algorithm
        /// compares 2787 (seconds-as-loaded) against ~46 (the scan's RT in minutes)
        /// and never sees overlap on any scan.
        /// </summary>
        private IEnumerable<ISingleChargeMs1Feature> NormaliseRetentionTimes(
            IEnumerable<ISingleChargeMs1Feature> features)
        {
            var materialised = features.ToList();
            if (materialised.Count == 0) return materialised;

            // Ignore non-finite RetentionTimeEnd values when sniffing units: a single
            // NaN/Infinity would otherwise make Max(...) non-finite and (NaN <= 500) false,
            // spuriously flipping a minutes file into a /60 conversion. If no finite value
            // exists at all, we can't tell the units apart, so leave them unchanged.
            double maxRt = materialised
                .Select(f => f.RetentionTimeEnd)
                .Where(double.IsFinite)
                .DefaultIfEmpty(double.NaN)
                .Max();
            if (double.IsNaN(maxRt) || maxRt <= 500.0) return materialised;

            // Looks like seconds (no realistic LC run exceeds 8 h = 500 min). Convert to
            // minutes, and record + announce it -- a unit guess on scientific data should be
            // visible to the caller (via the property below) and the console, not silent.
            RetentionTimeNormalizedFromSeconds = true;
            Console.Error.WriteLine(
                "[FromFileDeconvolution] Feature-file retention times look like seconds " +
                $"(max RetentionTimeEnd {maxRt:F0} > 500 min); normalising to minutes (/60).");

            return materialised.Select(f => (ISingleChargeMs1Feature)new SingleChargeMs1Feature(
                f.Mz,
                f.Charge,
                f.RetentionTimeStart / 60.0,
                f.RetentionTimeEnd / 60.0,
                f.Intensity,
                f.NumberOfIsotopes));
        }

        /// <summary>
        /// Constructs from a pre-loaded sequence of features. Internal because public
        /// callers should go through the file-path constructor; this overload exists
        /// for unit tests that want to seed the parameters synthetically without
        /// writing a temporary feature file to disk.
        /// </summary>
        internal FromFileDeconvolutionParameters(IEnumerable<ISingleChargeMs1Feature> features,
            int minCharge, int maxCharge, Polarity polarity = Polarity.Positive)
            : base(minCharge, maxCharge, polarity)
        {
            if (features is null) throw new ArgumentNullException(nameof(features));
            _featureCache = BuildFeatureCache(features);
        }

        /// <summary>
        /// Returns the index of the first feature whose <c>Mz</c> is at or above
        /// <paramref name="minMz"/>. If no such feature exists, returns
        /// <c>Features.Count</c>. Callers should iterate forward from this index
        /// and break once <c>Features[i].Mz</c> exceeds their upper bound.
        /// </summary>
        internal int FindFirstIndexAtOrAbove(double minMz)
        {
            int idx = Array.BinarySearch(GetOrLoadFeatureCache().MzKeys, minMz);
            return idx < 0 ? ~idx : idx;
        }

        #region IEquatable<FromFileDeconvolutionParameters>

        protected override bool EqualProperties(DeconvolutionParameters other)
        {
            var o = (FromFileDeconvolutionParameters)other;
            if (Features.Count != o.Features.Count) return false;
            if (FilePath != o.FilePath) return false;
            return true;
        }

        protected override void AddHashCodes(HashCode hash)
        {
            hash.Add(Features.Count);
            hash.Add(FilePath);
        }

        #endregion

        public override FromFileDeconvolutionParameters Clone()
        {
            var clone = !string.IsNullOrWhiteSpace(FilePath) 
                ? new FromFileDeconvolutionParameters(FilePath!, MinAssumedChargeState, MaxAssumedChargeState, Polarity) 
                : new FromFileDeconvolutionParameters(Features.ToList(), MinAssumedChargeState, MaxAssumedChargeState, Polarity);

            clone.UseGenericScore = UseGenericScore;
            clone.ExpectedIsotopeSpacing = ExpectedIsotopeSpacing;
            clone.AverageResidueModel = AverageResidueModel;
            return clone;
        }

        /// <summary>
        /// Polymorphic factory hook: returns a <see cref="FromFileDeconvolutionAlgorithm"/>
        /// bound to these parameters. <see cref="Deconvoluter"/> calls this first,
        /// before falling back to its enum-based switch.
        /// </summary>
        public override DeconvolutionAlgorithm CreateAlgorithm() => new FromFileDeconvolutionAlgorithm(this);

        /// <summary>No decoy support — masses are taken as authoritative.</summary>
        public override DeconvolutionParameters? ToDecoyParameters() => null;
    }
}
