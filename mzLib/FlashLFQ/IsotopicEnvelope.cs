using System;
using Easy.Common.Extensions;
using MassSpectrometry.MzSpectra;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    /// <summary>
    /// Contains the summed intensities of all isotope peaks detected in a single MS1 scan for a given species.
    /// </summary>
    public class IsotopicEnvelope
    {
        /// <summary>
        /// The most abundant isotopic peak used for peak finding.
        /// </summary>
        public readonly IIndexedMzPeak IndexedPeak;
        public readonly int ChargeState;
        /// <summary>
        /// The summed intensity of all isotope peaks detected in one MS1 scan. This sum may contain 
        /// imputed intensity values for expected isotopes that weren't observed, but only if the observed 
        /// isotopic distribution was otherwise similar to the expected isotopic distribution.
        /// </summary>
        public double Intensity { get; private set; }

        /// <summary>
        /// The Pearson correlation coefficient between the observed isotopic distribution and the theoretical isotopic distribution.
        /// </summary>
        public double PearsonCorrelation { get; init; }

        /// <summary>
        /// Contains the summed intensities of all isotope peaks detected in a single MS1 scan for a given species.
        /// </summary>
        /// <param name="mostAbundantIsotopologuePeak">The IndexedMzPeak for the theoretical most abundant isotopologue of a species used during peakfinding</param>
        public IsotopicEnvelope(IIndexedMzPeak mostAbundantIsotopologuePeak, int chargeState, double intensity, double pearsonCorrelation)
        {
            IndexedPeak = mostAbundantIsotopologuePeak;
            ChargeState = chargeState;
            Intensity = intensity / chargeState;
            PearsonCorrelation = pearsonCorrelation;
        }

        public IsotopicEnvelope(IndexedMassSpectralPeak monoisotopicPeak, int chargeState, double intensity, double pearsonCorrelation) :
            this((IIndexedMzPeak)monoisotopicPeak, chargeState, intensity, pearsonCorrelation) { }


        public List<IIndexedMzPeak> IsotopologuePeaks { get; private set; }
        public IsotopicEnvelope(IIndexedMzPeak mostAbundantIsotopologuePeak, int chargeState, double intensity, double pearsonCorrelation, List<IIndexedMzPeak> isotopePeaks) :
            this(mostAbundantIsotopologuePeak, chargeState, intensity, pearsonCorrelation)
        {
            IsotopologuePeaks = isotopePeaks.OrderBy(p => p.Mz).ToList();
        }

        public HashSet<IIndexedMzPeak> PeakSet
        {
            get
            {
                _peakSet ??= IsotopologuePeaks.ToHashSet();
                return _peakSet;
            }
        }

        private HashSet<IIndexedMzPeak> _peakSet;
        public double[] MzArray => IsotopologuePeaks.IsNotNullOrEmpty() ? IsotopologuePeaks.Select(p => p.Mz).ToArray() : [];
        public double[] IntensityArray => IsotopologuePeaks.IsNotNullOrEmpty() ? IsotopologuePeaks.Select(p => p.Intensity).ToArray() : [];

        public double CheckSimilarity(IsotopicEnvelope other)
        {
            double? cosineSimilarity = new SpectralSimilarity(MzArray, IntensityArray, other.MzArray, other.IntensityArray, SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak,
                toleranceInPpm: 5, allPeaks: true, filterOutBelowThisMz: 100).CosineSimilarity();
            return cosineSimilarity ?? -1;
        }

        public void Normalize(double normalizationFactor)
        {
            Intensity *= normalizationFactor;
        }

        public override string ToString()
        {
            return "+" + ChargeState + "|" + Intensity.ToString("F0") + "|" + IndexedPeak.RetentionTime.ToString("F3") + "|" + IndexedPeak.ZeroBasedScanIndex;
        }

        public override bool Equals(object obj)
        {
            var otherEnv = (IsotopicEnvelope)obj;

            return otherEnv != null
                && otherEnv.ChargeState == this.ChargeState
                && otherEnv.IndexedPeak.Equals(this.IndexedPeak);
        }

        public override int GetHashCode()
        {
            return ChargeState.GetHashCode() + IndexedPeak.GetHashCode();
        }
    }
}