#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
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

        // Sorted by Mz; FeaturesMzAscending and _mzKeys are kept in lockstep.
        private readonly List<ISingleChargeMs1Feature> _featuresMzAscending;
        private readonly double[] _mzKeys;

        /// <summary>
        /// The pre-loaded per-charge features the algorithm filters against, sorted
        /// by ascending <see cref="ISingleChargeMs1Feature.Mz"/>.
        /// </summary>
        public IReadOnlyList<ISingleChargeMs1Feature> Features => _featuresMzAscending;

        /// <summary>
        /// Constructs from a feature-file path. Reader is auto-detected from the
        /// extension. The file's per-charge expansion is materialized in-memory
        /// once; this object is intended to be reused across many MS2 scans rather
        /// than reconstructed per call.
        /// </summary>
        public FromFileDeconvolutionParameters(string filePath, int minCharge, int maxCharge,
            Polarity polarity = Polarity.Positive)
            : base(minCharge, maxCharge, polarity)
        {
            if (filePath is null) throw new ArgumentNullException(nameof(filePath));

            var resultFile = FileReader.ReadResultFile(filePath);
            if (resultFile is not IMs1FeatureFile featureFile)
                throw new MzLibException(
                    $"File at '{filePath}' is not a recognized MS1 feature file. " +
                    $"Detected type: {resultFile.GetType().Name}. " +
                    "Expected an Ms1Feature (.ms1.feature) or Dinosaur (.feature.tsv) file.");

            resultFile.LoadResults();
            (_featuresMzAscending, _mzKeys) = SortAndIndex(featureFile.GetMs1Features());
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
            (_featuresMzAscending, _mzKeys) = SortAndIndex(features);
        }

        private static (List<ISingleChargeMs1Feature> sorted, double[] keys)
            SortAndIndex(IEnumerable<ISingleChargeMs1Feature> features)
        {
            var sorted = features.OrderBy(f => f.Mz).ToList();
            var keys = new double[sorted.Count];
            for (int i = 0; i < sorted.Count; i++) keys[i] = sorted[i].Mz;
            return (sorted, keys);
        }

        /// <summary>
        /// Returns the index of the first feature whose <c>Mz</c> is at or above
        /// <paramref name="minMz"/>. If no such feature exists, returns
        /// <c>Features.Count</c>. Callers should iterate forward from this index
        /// and break once <c>Features[i].Mz</c> exceeds their upper bound.
        /// </summary>
        internal int FindFirstIndexAtOrAbove(double minMz)
        {
            int idx = Array.BinarySearch(_mzKeys, minMz);
            return idx < 0 ? ~idx : idx;
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
