using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    public class IsotopicEnvelope
    {
        public int EnvelopeIdentifier { get; set; }
        public MsDataScan Scan { get; set; }

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
            Score = Peaks.Count * pearsonCorr;
        }

        public override string ToString()
        {
            return MonoisotopicMass.ToString("F1") + "; Peaks: " + Peaks.Count + "; z=" + Charge;
        }

        public static string OutputHeader()
        {
            return "Scan Number" + "\t" +
                "Monoisotopic Mass" + "\t" +
                "Peaks m/z list" + "\t" +
                "Peaks intensity list" + "\t" +
                "Charge" + "\t" +
                "MS Order" + "\t" +
                "Pearson Correlation to Averagine" + "\t" +
                "Fraction of Intensity Observed" + "\t" +
                "S/N" + "\t" +
                "Total Intensity" + "\t" +
                "Noise" + "\t" +
                "ID" + "\t";
        }

        public string ToOutputString()
        {
            return Scan.OneBasedScanNumber + "\t" +
                MonoisotopicMass + "\t" +
                string.Join(";", Peaks.Select(p => p.mz.ToString("F3"))) + "\t" +
                string.Join(";", Peaks.Select(p => p.intensity.ToString("F1"))) + "\t" +
                Charge + "\t" +
                Scan.MsnOrder + "\t" +
                PearsonCorrelation + "\t" +
                FracIntensityObserved + "\t" +
                SN + "\t" +
                Math.Log(TotalIntensity, 2) + "\t" +
                Noise + "\t" +
                EnvelopeIdentifier;
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