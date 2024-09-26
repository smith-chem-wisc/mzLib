using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;

namespace MassSpectrometry
{
    public class IsotopicEnvelope : IHasMass
    {
        public readonly List<(double mz, double intensity)> Peaks;
        public double MonoisotopicMass { get; private set; }

        /// <summary>
        /// Mass of most abundant observed isotopic peak, not accounting for addition or subtraction or protons due to ESI charge state induction
        /// </summary>
        public double MostAbundantObservedIsotopicMass { get; private set; }
        public readonly int Charge;
        public readonly double TotalIntensity;
        public readonly double StDev;
        public readonly int MassIndex;

        public double Score { get; private set; }

        public IsotopicEnvelope(List<(double mz, double intensity)> bestListOfPeaks, double bestMonoisotopicMass, int bestChargeState, double bestTotalIntensity, double bestStDev, int bestMassIndex)
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

        public double GetMostAbundantObservedIsotopicMass(List<(double mz, double intensity)> peaks, int charge)
        {
            return peaks.MaxBy(p => p.intensity).mz * Math.Abs(charge);
        }

        public override string ToString()
        {
            return Charge + "\t" + Peaks[0].mz.ToString("G8") + "\t" + Peaks.Count + "\t" + TotalIntensity;
        }

        private double ScoreIsotopeEnvelope() //likely created by Stefan Solntsev using peptide data
        {
            return Peaks.Count >= 2 ?
                TotalIntensity / Math.Pow(StDev, 0.13) * Math.Pow(Peaks.Count, 0.4) / Math.Pow(Math.Abs(Charge), 0.06) :
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