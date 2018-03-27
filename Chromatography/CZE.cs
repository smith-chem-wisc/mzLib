using System;
using System.Collections.Generic;
using System.Text;

namespace Chromatography
{
    public class CZE
    {
        public double ExperimentalElectrophoreticMobility(double length, double timeMin, double voltsPerMeter)
        {
            if (length >= 0 && timeMin >= 0)
                return (length / (60 * timeMin * voltsPerMeter) * 1e9);
            else
                return (-1);
        }

        public double PredictedElectrophoreticMobility(string peptideSequence, double observedMass)
        {
            double predictedMu = 0;

            predictedMu = 3.069 + 386 * Math.Log(1d + 0.35 * PredictedChargeCorrected(peptideSequence)) / (Math.Pow(observedMass, 0.411) + Offset(PredictedChargeCorrected(peptideSequence), peptideSequence.Length));
            return predictedMu;
        }

        private double PredictedCharge(string peptideSequence)
        {
            string substitutedString = peptideSequence.Replace("R","").Replace("K","").Replace("H","").ToString();
            return (1d + (double)(peptideSequence.Length - substitutedString.Length));
        }

        private double PredictedChargeCorrected(string peptideSequence)
        {
            double runningSum = 0;
            string internalString = peptideSequence.Substring(3, peptideSequence.Length - 5);


            switch (peptideSequence.Substring(0,1))//first AA
            {
                case "D":
                    runningSum -=0.26741;
                    break;
                case "E":
                    runningSum -= 0.06852;
                    break;
                case "N":
                    runningSum += 0.011699;
                    break;
                case "Q":
                    runningSum += 0;
                    break;
                default:
                    runningSum += 0;
                    break;
            }
            switch (peptideSequence.Substring(1, 1))//second AA
            {
                case "D":
                    runningSum -= 0.10947;
                    break;
                case "E":
                    runningSum -= 0.04011;
                    break;
                case "N":
                    runningSum += 0.012535;
                    break;
                case "Q":
                    runningSum += 0.011699;
                    break;
                default:
                    runningSum += 0;
                    break;
            }
            switch (peptideSequence.Substring(2, 1))//third AA
            {
                case "D":
                    runningSum -= 0.08022;
                    break;
                case "E":
                    runningSum -= 0.03426;
                    break;
                case "N":
                    runningSum += 0.016713;
                    break;
                case "Q":
                    runningSum += 0.00585;
                    break;
                default:
                    runningSum += 0;
                    break;
            }
            switch (peptideSequence.Substring(peptideSequence.Length-2, 1))//second to last AA
            {
                case "D":
                    runningSum -= 0.03844;
                    break;
                case "E":
                    runningSum -= 0.01337;
                    break;
                case "N":
                    runningSum += 0.026741;
                    break;
                case "Q":
                    runningSum -= 0.00084;
                    break;
                default:
                    runningSum += 0;
                    break;
            }
            switch (peptideSequence.Substring(peptideSequence.Length - 1, 1))//last AA
            {
                case "D":
                    runningSum -= 0.02256;
                    break;
                case "E":
                    runningSum -= 0.00418;
                    break;
                case "N":
                    runningSum += 0.010864;
                    break;
                case "Q":
                    runningSum -= 0.0117;
                    break;
                default:
                    runningSum += 0;
                    break;
            }

            if (internalString.Contains("D"))
                runningSum -= 0.05014;
            if (internalString.Contains("E"))
                runningSum -= 0.01922;
            if (internalString.Contains("N"))
                runningSum += 0.012535;
            if (internalString.Contains("Q"))
                runningSum -= 0.000251;

            runningSum += PredictedCharge(peptideSequence);

            return runningSum;
        }

        private double Offset(double correctedCharge, int length)
        {
            return 0;
            //should fit 5th order polynomical to plot of (ExperimentalElectrophoreticMobility - PredictedElectrophoreticMobility) vs. (Zc/N) where N is peptidelength.
        }
    }
}
