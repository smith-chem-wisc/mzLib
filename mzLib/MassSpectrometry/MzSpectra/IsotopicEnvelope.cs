using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;

namespace MassSpectrometry
{
    public class IsotopicEnvelope
    {
        public readonly List<(double mz, double intensity)> Peaks;
        public double MonoisotopicMass { get; private set; }
        public readonly int Charge;

        // Legacy fields used in the ClassicDeconvolutionAlgorithm
        public readonly double TotalIntensity;
        public readonly double StDev;
        public readonly int MassIndex;

        public double[] MzArray => Peaks.OrderBy(p => p.mz).Select(p => p.mz).ToArray();
        public double[] IntensityArray => Peaks.OrderBy(p => p.mz).Select(p => p.intensity).ToArray();

        public double MostAbundantObservedIsotopicMass => _mostAbundantObservedIsotopicMass ?? 0;
        public double SecondMostAbundantObservedIsotopicMass => _secondMostAbundantObservedIsotopicMass ?? 0;
        private double? _mostAbundantObservedIsotopicMass;
        private double? _secondMostAbundantObservedIsotopicMass;
        public double AmbiguityRatioMinimum;
        public double Score { get; private set; }

        public IsotopicEnvelope(List<(double mz, double intensity)> bestListOfPeaks, double bestMonoisotopicMass,
            int bestChargeState, double bestTotalIntensity, double bestStDev, int bestMassIndex)
        {
            Peaks = bestListOfPeaks;
            MonoisotopicMass = bestMonoisotopicMass;
            Charge = bestChargeState;
            FindMostAbundantObservedIsotopicMass();

            TotalIntensity = bestTotalIntensity;
            StDev = bestStDev;
            MassIndex = bestMassIndex;
            Score = ScoreIsotopeEnvelope();
        }

        /// <summary>
        /// Takes in an Isotopic Distribution and a given charge state and converts it to an IsotopicEnvelope object
        /// TODO: Test this function specifically
        /// </summary>
        /// <param name="theoreticalDistribution"> An IsotopicDistribution generated from a ChemicalFormula</param>
        /// <param name="charge"> The charge state (corresponding to the z value of m/z) </param>
        public IsotopicEnvelope(IsotopicDistribution theoreticalDistribution, int charge, double ambiguityRatioMinimum = 0.9)
        {
            Peaks = theoreticalDistribution.Masses.Zip(theoreticalDistribution.Intensities,
                (first, second) => (first.ToMz(charge), (double)second)).ToList();
            MonoisotopicMass = theoreticalDistribution.MonoIsotopicMass; // I think this is right, need to test it tho
            Charge = charge;
            AmbiguityRatioMinimum = ambiguityRatioMinimum;
            FindMostAbundantObservedIsotopicMass();

        }

        /// <summary>
        /// Finds the m/z value of the greatest intensity peak. If the second most intense peak
        /// is within 90% of the most intense peak, the m/z value of that peak is stored 
        /// in the _secondMostAbundantObservedIsotopicMass field
        /// </summary>
        /// <returns></returns>
        public void FindMostAbundantObservedIsotopicMass()
        {
            if (!_mostAbundantObservedIsotopicMass.HasValue)
            {
                List<(double mz, double intensity)> intensityOrderedPeaks = Peaks.OrderByDescending(p => p.intensity).ToList();
                _mostAbundantObservedIsotopicMass = intensityOrderedPeaks.Select(p => p.mz).First().ToMz(Charge);
                if (intensityOrderedPeaks[1].intensity / intensityOrderedPeaks[0].intensity >= AmbiguityRatioMinimum &&
                    AmbiguityRatioMinimum > 0)
                {
                    _secondMostAbundantObservedIsotopicMass = intensityOrderedPeaks[1].mz;
                }
            }
        }

        public override string ToString()
        {
            return Charge + "\t" + Peaks[0].mz.ToString("G8") + "\t" + Peaks.Count + "\t" + TotalIntensity;
        }

        // This should be done using a Strategy pattern
        private double ScoreIsotopeEnvelope() //likely created by Stefan Solntsev using peptide data
        {
            return Peaks.Count >= 2 ?
                TotalIntensity / Math.Pow(StDev, 0.13) * Math.Pow(Peaks.Count, 0.4) / Math.Pow(Charge, 0.06) :
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

    }
}