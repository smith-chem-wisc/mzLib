using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;

namespace MassSpectrometry
{
    public class IsotopicEnvelope : IHasMass, IEquatable<IsotopicEnvelope>
    {
        public readonly List<(double mz, double intensity)> Peaks;
        public double MonoisotopicMass { get; private set; }

        /// <summary>
        /// Mass of most abundant observed isotopic peak, not accounting for addition or subtraction or protons due to ESI charge state induction
        /// </summary>
        internal double MostAbundantObservedIsotopicMass { get; private set; }
        public readonly int Charge;
        public readonly double TotalIntensity;
        public readonly int PrecursorId;

        /// <summary>
        /// Algorithm-specific score. Set at construction time and is the value used by each
        /// deconvolution algorithm internally (e.g. Classic uses an empirical heuristic, IsoDec
        /// passes its DLL-computed cosine score). Different algorithms place values on
        /// different scales — do not compare directly across algorithms.
        /// </summary>
        public double Score { get; private set; }

        /// <summary>
        /// Optional, post-construction generic deconvolution score in [0, 1] computed by
        /// <see cref="DeconvolutionScorer.ScoreEnvelope"/>. Null until <see cref="SetGenericScore"/>
        /// is called. Designed to be comparable across deconvolution algorithms because it is
        /// recomputed from the envelope's peak list and an Averagine model — independent of
        /// which algorithm produced the envelope.
        /// </summary>
        public double? GenericScore { get; private set; }

        /// <summary>
        /// Used for an isotopic envelope that mzLib deconvoluted (e.g., from a mass spectrum)
        /// </summary>
        public IsotopicEnvelope(List<(double mz, double intensity)> bestListOfPeaks, double bestMonoisotopicMass, int bestChargeState, double bestTotalIntensity, double bestStDev)
        {
            Peaks = bestListOfPeaks;
            MonoisotopicMass = bestMonoisotopicMass;
            MostAbundantObservedIsotopicMass = bestListOfPeaks.MaxBy(p => p.intensity).mz * Math.Abs(bestChargeState);
            Charge = bestChargeState;
            TotalIntensity = bestTotalIntensity;
            Score = ScoreIsotopeEnvelope(bestStDev);
        }

        /// <summary>
        /// Used for a neutral mass read in from a deconvoluted file
        /// Assumes the mass is correct: score is max value
        /// </summary>
        public IsotopicEnvelope(double monoisotopicMass, double intensity, int charge)
        {
            MonoisotopicMass = monoisotopicMass;
            Charge = charge;
            TotalIntensity = intensity;
            Score = double.MaxValue;
            Peaks = [(monoisotopicMass.ToMz(charge), intensity)];
        }

        /// <summary>
        /// Used for A deconvolution method that calculates its own score. 
        /// </summary>
        /// <param name="id">All missed mono products of the same peak will share an ID if enabled in IsoDec</param>
        /// <param name="peaks"></param>
        /// <param name="monoisotopicmass"></param>
        /// <param name="chargestate"></param>
        /// <param name="intensity"></param>
        /// <param name="score"></param>
        public IsotopicEnvelope(int id, List<(double mz, double intensity)> peaks, double monoisotopicmass, int chargestate, double intensity, double score)
        {
            PrecursorId = id;
            Peaks = peaks;
            MonoisotopicMass = monoisotopicmass;
            Charge = chargestate;
            TotalIntensity = intensity;
            Score = score;
            MostAbundantObservedIsotopicMass = peaks.MaxBy(p => p.intensity).mz * Math.Abs(chargestate);
        }

        /// <summary>
        /// Sets the optional generic deconvolution score for this envelope. The generic score is
        /// produced by <see cref="DeconvolutionScorer"/> and is comparable across deconvolution
        /// algorithms. Calling this method does not change <see cref="Score"/>.
        /// </summary>
        /// <param name="genericScore">Generic deconvolution score, expected in [0, 1].</param>
        public void SetGenericScore(double genericScore)
        {
            GenericScore = genericScore;
        }

        /// <summary>
        /// True if a generic deconvolution score has been computed and stored on this envelope.
        /// Equivalent to <c>GenericScore.HasValue</c>; provided for readability at call sites.
        /// </summary>
        public bool HasGenericScore => GenericScore.HasValue;

        /// <summary>
        /// Returns <see cref="GenericScore"/> if set, otherwise <see cref="Score"/>. Callers that
        /// want a single number per envelope and don't care which scoring system produced it should
        /// prefer <see cref="GenericScore"/> directly, since it is the only score that is comparable
        /// across deconvolution algorithms.
        /// </summary>
        /// <remarks>
        /// This fallback exists for callers that have not yet wired up generic scoring — it lets
        /// them write threshold or ranking code that works whether or not the score has been computed.
        /// New code should not rely on this fallback; it should compute the generic score explicitly
        /// (e.g. via <see cref="IsotopicEnvelopeExtensions.GetOrComputeGenericScore(IsotopicEnvelope, AverageResidue)"/>).
        /// </remarks>
        public double GenericOrFallbackScore => GenericScore ?? Score;

        public override string ToString()
        {
            return Charge + "\t" + Peaks[0].mz.ToString("G8") + "\t" + Peaks.Count + "\t" + TotalIntensity;
        }

        private double ScoreIsotopeEnvelope(double stDev) //likely created by Stefan Solntsev using peptide data
        {
            return Peaks.Count >= 2 ?
                TotalIntensity / Math.Pow(stDev, 0.13) * Math.Pow(Peaks.Count, 0.4) / Math.Pow(Math.Abs(Charge), 0.06) :
                0;
        }

        public void AggregateChargeStateScore(IsotopicEnvelope chargeStateEnvelope)
        {
            Score += chargeStateEnvelope.Score;
        }

        public void SetMedianMonoisotopicMass(List<double> monoisotopicMassPredictions)
        {
            MonoisotopicMass = monoisotopicMassPredictions.Median();
        }

        public bool Equals(IsotopicEnvelope other)
        {
            if (other is null) return false;
            if (ReferenceEquals(this, other)) return true;
            if (Charge != other.Charge || Peaks.Count != other.Peaks.Count) return false;
            if (Math.Abs(TotalIntensity - other.TotalIntensity) >= 0.001) return false;
            if (Math.Abs(MonoisotopicMass - other.MonoisotopicMass) >= 0.001) return false;
            if (Math.Abs(MostAbundantObservedIsotopicMass - other.MostAbundantObservedIsotopicMass) >= 0.001) return false;

            for (int i = 0; i < Peaks.Count; i++)
            {
                var p1 = Peaks[i];
                var p2 = other.Peaks[i];
                if (Math.Abs(p1.mz - p2.mz) >= 0.001 || Math.Abs(p1.intensity - p2.intensity) >= 0.001)
                    return false;
            }
            return true;
        }

        public override bool Equals(object obj)
        {
            if (obj is null) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != GetType()) return false;
            return Equals((IsotopicEnvelope)obj);
        }

        public override int GetHashCode()
        {
            unchecked
            {
                int hash = 17;
                foreach (var peak in Peaks)
                {
                    hash = hash * 23 + peak.mz.GetHashCode();
                    hash = hash * 23 + peak.intensity.GetHashCode();
                }
                hash = hash * 23 + Charge.GetHashCode();
                hash = hash * 23 + TotalIntensity.GetHashCode();
                hash = hash * 23 + MonoisotopicMass.GetHashCode();
                hash = hash * 23 + MostAbundantObservedIsotopicMass.GetHashCode();
                return hash;
            }
        }
    }
}