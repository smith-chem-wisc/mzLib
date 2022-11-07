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
        public double MostAbundantObservedIsotopicMass { get; }
        private double? _mostAbundantObservedIsotopicMass;
        public readonly int Charge;
        public readonly double TotalIntensity;
        public readonly double StDev;
        public readonly int MassIndex;

        public double Score { get; private set; }

        public IsotopicEnvelope(List<(double mz, double intensity)> bestListOfPeaks, double bestMonoisotopicMass,
            int bestChargeState, double bestTotalIntensity, double bestStDev, int bestMassIndex)
        {
            Peaks = bestListOfPeaks;
            MonoisotopicMass = bestMonoisotopicMass;
            MostAbundantObservedIsotopicMass = GetMostAbundantObservedIsotopicMass(bestListOfPeaks, bestChargeState);
            Charge = bestChargeState;
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
        public IsotopicEnvelope(IsotopicDistribution theoreticalDistribution, int charge)
        {
            Peaks = theoreticalDistribution.Masses.Zip(theoreticalDistribution.Intensities,
                (first, second) => (first.ToMz(charge), (double)second)).ToList();
            MonoisotopicMass = theoreticalDistribution.Masses.Min(); // I think this is right, need to test it tho
            MostAbundantObservedIsotopicMass = GetMostAbundantObservedIsotopicMass(Peaks, charge);
            Charge = charge;

        }

        // This is terrifying. It sure looks like the most abundant observed isotopic mass was calculated by multiplying by charge, without
        // any correction for the protons present
        public double GetMostAbundantObservedIsotopicMass(List<(double mz, double intensity)> peaks, int charge)
        {
            if (!_mostAbundantObservedIsotopicMass.HasValue)
            {
                _mostAbundantObservedIsotopicMass = (peaks.OrderByDescending(p => p.intensity).ToList()[0].Item1) * charge;
            }

            return (double)_mostAbundantObservedIsotopicMass;
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