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
        /// Most abundant observed isotopic peak as m/z × |charge|, i.e. <b>not</b> proton-corrected to a
        /// neutral mass (it does not subtract the charge-carrier protons added during ESI). This is a
        /// charge-scaled m/z, not a neutral mass; for the proton-corrected neutral mass use
        /// <see cref="MostAbundantObservedNeutralMass"/>. Retained for existing deconvolution-quality tests.
        ///
        /// <b>Sentinel:</b> this value is only meaningful for an envelope that carries observed isotopic
        /// peaks. Constructors that have no observed envelope (e.g. a neutral mass read from a
        /// pre-deconvoluted file) set it to <c>-1</c> to mark "no most-abundant peak available", rather
        /// than a misleading zero or a synthetic value.
        /// </summary>
        public double MostAbundantObservedIsotopicMass { get; private set; }

        /// <summary>
        /// The <b>most-abundant observed neutral mass</b>: the neutral mass of the single most intense
        /// (tallest) observed isotopic peak of this charge-state envelope — the highest-signal, directly
        /// measured point of the envelope. The most intense peak's m/z is converted to a neutral mass via
        /// the envelope charge (i.e. proton-corrected, <c>mz.ToMass(Charge)</c>), so it is directly
        /// comparable to a theoretical proteoform neutral mass. It is the precursor mass used for candidate
        /// selection in most-abundant mode. Computed from <see cref="Peaks"/> and <see cref="Charge"/>
        /// (both readonly).
        ///
        /// Terminology (the three masses are distinct):
        ///  • <b>most-abundant (observed neutral) mass</b> — this property; the tallest single isotopic
        ///    peak, proton-corrected to a neutral mass.
        ///  • <b>monoisotopic mass</b> (<see cref="MonoisotopicMass"/>) — the all-light-isotope mass; for
        ///    large proteoforms this isotopologue is rare and often undetectable.
        ///  • <b>average / centroid mass</b> — the intensity-weighted mean over the whole envelope, used
        ///    for isotopically unresolved (high-mass) species. That property is added separately, with the
        ///    unresolved-envelope work, and is not part of this resolved most-abundant feature.
        /// This is the neutral-mass form of <see cref="MostAbundantObservedIsotopicMass"/> (which is the
        /// SAME tallest peak expressed as m/z × |charge|, NOT proton-corrected): the two differ by exactly
        /// <c>|charge| × ProtonMass</c>. Deriving it from that single source means the two values can never
        /// disagree, and it carries the same <c>-1</c> "no observed peak" sentinel.
        /// </summary>
        public double MostAbundantObservedNeutralMass =>
            MostAbundantObservedIsotopicMass < 0 ? -1 : MostAbundantObservedIsotopicMass - Charge * Constants.ProtonMass;

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
            // A neutral mass read from a pre-deconvoluted file has no observed isotopic envelope, so
            // there is no most-abundant observed peak to report — use the -1 sentinel.
            MostAbundantObservedIsotopicMass = -1;
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