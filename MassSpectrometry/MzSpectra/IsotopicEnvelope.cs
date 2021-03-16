using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;

namespace MassSpectrometry
{
    public class IsotopicEnvelope
    {
        public readonly List<(double mz, double intensity)> Peaks;
        public double MonoisotopicMass { get; private set; }
        public readonly int Charge;
        public readonly double TotalIntensity;
        public readonly double StDev;
        public readonly int MassIndex;
        public double Score { get; set; }
        public double SN { get; set; }
        public double PearsonCorrelation { get; set; }
        public double FracIntensityObserved { get; set; }
        public double Noise { get; set; }

        public IsotopicEnvelope(List<(double mz, double intensity)> bestListOfPeaks, double bestMonoisotopicMass, int bestChargeState, double bestTotalIntensity, double bestStDev, int bestMassIndex)
        {
            Peaks = bestListOfPeaks;
            MonoisotopicMass = bestMonoisotopicMass;
            Charge = bestChargeState;
            TotalIntensity = bestTotalIntensity;
            StDev = bestStDev;
            MassIndex = bestMassIndex;
            Score = ScoreIsotopeEnvelope();
        }

        public IsotopicEnvelope(List<(double mz, double intensity)> peaks, double monoMass, int charge, double intensity, int massIndex, double score)
        {
            Peaks = peaks;
            MonoisotopicMass = monoMass;
            Charge = charge;
            TotalIntensity = intensity;
            StDev = double.NaN;
            MassIndex = massIndex;
            Score = score;
        }

        public override string ToString()
        {
            //return Charge + "\t" + Peaks[0].mz.ToString("G8") + "\t" + Peaks.Count + "\t" + TotalIntensity;
            return MonoisotopicMass.ToString("F1") + "; Peaks: " + Peaks.Count + "; z=" + Charge;
        }

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