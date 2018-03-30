using System;
using System.Collections.Generic;
using System.Text;

namespace Chromatography
{
    /// <summary>
    /// This class will return theoretical retention times, hydrobphobicites, electrophoretic mobilities and etc. for peptides.
    /// These values would be useful for comparision with experimentally observed retention times. This information might be
    /// informative for evaluation of false positives and also for discerning the prescence of certain PTMs that would
    /// alter the experimental chromatographic behavior.
    /// </summary>




    /// <summary>
    /// The public methods of this class are limited to returning hydrophobicity calculated using an algorithm such as SSRCALC.
    /// </summary>
    public class HPLC
    {
        
        //public void PredictedHydrophobicity()
        //{

        //}  
    }


    /// <summary>
    /// The public methods of this class are limited to electrophoretic mobilities of peptides detected in a CZE-MS/MS experiment.
    /// </summary>
    public class CZE
    {
        /// <summary>
        /// This class returns calculated electrophoretic mobility for an observed peptide. The calculation requires use of an
        /// observed retention time(min), the total capillary length(m) and the applied voltage (V/m)
        /// </summary>
        /// <param name="length"></param>
        /// <param name="timeMin"></param>
        /// <param name="voltsPerMeter"></param>
        /// <returns></returns>
        public double ExperimentalElectrophoreticMobility(double length, double timeMin, double voltsPerMeter)
        {
            if (length >= 0 && timeMin >= 0)
                return (length / (60 * timeMin * voltsPerMeter) * 1e9);
            else
                return (-1);
        }

        /// <summary>
        /// This calculated the predicted electrophoretic mobility of a peptide.
        ///
        /// See for reference 
        /// Anal Chem. 2017 Feb 7;89(3):2000-2008. doi: 10.1021/acs.analchem.6b04544. Epub 2017 Jan 19.
        /// Predicting Electrophoretic Mobility of Tryptic Peptides for High-Throughput CZE-MS Analysis.
        /// Krokhin OV, Anderson G, Spicer V, Sun L1, Dovichi NJ2.
        /// https://www.ncbi.nlm.nih.gov/pubmed/28208305
        /// 
        /// </summary>
        /// <param name="peptideSequence"></param>
        /// <param name="observedMass"></param>
        /// <returns></returns>
        public double PredictedElectrophoreticMobility(string peptideSequence, double observedMass)
        {
            double predictedMu = 0;

            //calculation described in Anal Chem. 2017 Feb 7;89(3):2000-2008
            //3.069 and 386 are coefficients applied to align output with experimentally measured values(slope 1 and intercept 0 in). I think we may need to reset these.
            //other values from best fit model of Cifuentes and Poppe (J. Chromatogr. A 1994, 680, 321−340) used as described in the AC paper.
            predictedMu = 3.069 + 386 * Math.Log(1d + 0.35 * PredictedChargeCorrected(peptideSequence)) / 
                (Math.Pow(observedMass, 0.411) + Offset(PredictedChargeCorrected(peptideSequence), peptideSequence.Length));

            return predictedMu;
        }

        /// <summary>
        /// The predicted charge is plus 1 for the N-terminal and plus for the count of lysine(K), arginine(R) and histidine(H).
        /// </summary>
        /// <param name="peptideSequence"></param>
        /// <returns></returns>
        private double PredictedCharge(string peptideSequence)
        {
            string substitutedString = peptideSequence.Replace("R","").Replace("K","").Replace("H","").ToString();
            return (1d + (double)(peptideSequence.Length - substitutedString.Length));
        }

        /// <summary>
        /// minimal charge correction is position dependenat and predominantly at the peptide termini. Adjustments are made for presence of D, E, N and Q
        /// at the ends and in the middle. 
        /// 
        /// In the future, I would like to use linear algebra to estimate these more accurately for each dataset separately. Currently
        /// these numbers are from a table in Anal Chem. 2017 Feb 7;89(3):2000-2008. doi: 10.1021/acs.analchem.6b04544. Epub 2017 Jan 19.
        /// 
        /// </summary>
        /// <param name="peptideSequence"></param>
        /// <returns></returns>
        private double PredictedChargeCorrected(string peptideSequence)
        {
            double runningSum = 0;
            string internalString = peptideSequence.Substring(3, peptideSequence.Length - 5);


            switch (peptideSequence[0])//first AA
            {
                case 'D':
                    runningSum -=0.26741;
                    break;
                case 'E':
                    runningSum -= 0.06852;
                    break;
                case 'N':
                    runningSum += 0.011699;
                    break;
                case 'Q':
                    runningSum += 0;
                    break;
                default:
                    runningSum += 0;
                    break;
            }
            switch (peptideSequence[1])//second AA
            {
                case 'D':
                    runningSum -= 0.10947;
                    break;
                case 'E':
                    runningSum -= 0.04011;
                    break;
                case 'N':
                    runningSum += 0.012535;
                    break;
                case 'Q':
                    runningSum += 0.011699;
                    break;
                default:
                    runningSum += 0;
                    break;
            }
            switch (peptideSequence[2])//third AA
            {
                case 'D':
                    runningSum -= 0.08022;
                    break;
                case 'E':
                    runningSum -= 0.03426;
                    break;
                case 'N':
                    runningSum += 0.016713;
                    break;
                case 'Q':
                    runningSum += 0.00585;
                    break;
                default:
                    runningSum += 0;
                    break;
            }
            switch (peptideSequence[peptideSequence.Length-2])//second to last AA
            {
                case 'D':
                    runningSum -= 0.03844;
                    break;
                case 'E':
                    runningSum -= 0.01337;
                    break;
                case 'N':
                    runningSum += 0.026741;
                    break;
                case 'Q':
                    runningSum -= 0.00084;
                    break;
                default:
                    runningSum += 0;
                    break;
            }
            switch (peptideSequence[peptideSequence.Length - 1])//last AA
            {
                case 'D':
                    runningSum -= 0.02256;
                    break;
                case 'E':
                    runningSum -= 0.00418;
                    break;
                case 'N':
                    runningSum += 0.010864;
                    break;
                case 'Q':
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

        /// <summary>
        /// 
        /// The offset in the AC paper is a 5th order polynomial best fit to a plot of Zc/N versus the difference between experimental and predicted electrophoretic mobility. 
        /// This smells of dead fish. I'm leaving it out for not but it might need to be used as some point.
        /// 
        /// </summary>
        /// <param name="correctedCharge"></param>
        /// <param name="length"></param>
        /// <returns></returns>
        private double Offset(double correctedCharge, int length)
        {
            return 0;
            //should fit 5th order polynomical to plot of (ExperimentalElectrophoreticMobility - PredictedElectrophoreticMobility) vs. (Zc/N) where N is peptidelength.
        }
    }
}
