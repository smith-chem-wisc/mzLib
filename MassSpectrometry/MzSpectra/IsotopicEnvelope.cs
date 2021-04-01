using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    public class IsotopicEnvelope
    {
        public List<(double mz, double intensity)> Peaks { get; private set; }
        public double MonoisotopicMass { get; private set; }
        public readonly int Charge;
        public double TotalIntensity { get; private set; }
        public readonly double StDev;
        public readonly int MassIndex;
        public double Score { get; set; }
        public double SN { get; set; }
        public double PearsonCorrelation { get; private set; }
        public double FracIntensityObserved { get; private set; }
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

        public IsotopicEnvelope(List<(double mz, double intensity)> peaks, double monoMass, int charge, double intensity, int massIndex, double pearsonCorr, double fractionIntensityObserved)
        {
            Peaks = peaks;
            MonoisotopicMass = monoMass;
            Charge = charge;
            TotalIntensity = intensity;
            StDev = double.NaN;
            MassIndex = massIndex;
            PearsonCorrelation = pearsonCorr;
            FracIntensityObserved = fractionIntensityObserved;
            Score = fractionIntensityObserved * pearsonCorr;//pearsonCorr * peaks.Count;
        }

        public override string ToString()
        {
            //return Charge + "\t" + Peaks[0].mz.ToString("G8") + "\t" + Peaks.Count + "\t" + TotalIntensity;
            return MonoisotopicMass.ToString("F1") + "; Peaks: " + Peaks.Count + "; z=" + Charge;

            //env.MonoisotopicMass + "\t" +
            //                   string.Join(";", env.Peaks.Select(p => p.mz.ToString("F3"))) + "\t" +
            //                   string.Join(";", env.Peaks.Select(p => p.intensity.ToString("F1"))) + "\t" +
            //                   env.Charge + "\t" +
            //                   scan.MsnOrder + "\t" +
            //                   env.PearsonCorrelation + "\t" +
            //                   env.FracIntensityObserved + "\t" +
            //                   env.SN + "\t" +
            //                   Math.Log(env.TotalIntensity, 2) + "\t" +
            //                   env.Noise + "\t" +
            //                   envNum;

            //return MonoisotopicMass + "\t" +
            //    string.Join(";", Peaks.Select(p => p.mz.ToString("F3"))) + "\t" +
            //                   string.Join(";", Peaks.Select(p => p.intensity.ToString("F1"))) + "\t" +
            //                   Charge + "\t" +
            //                   PearsonCorrelation + "\t" +
            //                   FracIntensityObserved + "\t" +
            //                   Math.Log(TotalIntensity, 2) + "\t";

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