using System;

namespace Proteomics.RetentionTimePrediction
{
    /// <summary>
    /// This class will return theoretical retention times, hydrobphobicites, electrophoretic mobilities and etc. for peptides.
    /// These values would be useful for comparision with experimentally observed retention times. This information might be
    /// informative for evaluation of false positives and also for discerning the prescence of certain PTMs that would
    /// alter the experimental chromatographic behavior.
    ///
    /// This class returns calculated electrophoretic mobility for an observed peptide. The calculation requires use of an
    /// observed retention time(min), the total capillary length(m) and the applied voltage (V/m)
    ///
    /// The public methods of this class are limited to electrophoretic mobilities of peptides detected in a CZE-MS/MS experiment.
    /// </summary>
    public class CZE
    {
        private readonly double ColumnLength; //in meters
        private readonly double VoltsPerMeter; //in volts/meter
        public CZE(double columnLength, double voltsPerMeter)
        {
            ColumnLength = columnLength;
            VoltsPerMeter = voltsPerMeter;
        }

        /// <summary>
        /// This method returns calculated electrophoretic mobility for an observed peptide. The calculation requires use of an
        /// observed retention time(min), the total capillary length(m) and the applied voltage (V/m)
        /// </summary>
        /// <param name="timeMin"></param>
        /// <returns></returns>
        public double ExperimentalElectrophoreticMobility(double timeMin)
        {
            if (ColumnLength >= 0 && timeMin >= 0)
            {
                return ColumnLength / (60 * timeMin * VoltsPerMeter) * 1e9;
            }
            else
            {
                return -1;
            }
        }

        /// <summary>
        /// This method returns an expected retention time for a given electrophoretic mobility and experiment. The calculation requires use of an
        /// electrophoretic mobility, the total capillary length(m) and the applied voltage (V/m)
        /// </summary>
        /// <param name="electrophoreticMobility"></param>
        /// <returns></returns>
        public double TheoreticalElutionTime(double electrophoreticMobility)
        {
            if (ColumnLength >= 0)
            {
                return (ColumnLength * 1e9) / (60 * VoltsPerMeter *electrophoreticMobility);
            }
            else
            {
                return -1;
            }
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
        public static double PredictedElectrophoreticMobility(string peptideSequence, double observedMass)
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
        private static double PredictedCharge(string peptideSequence)
        {
            string substitutedString = peptideSequence.Replace("R", "").Replace("K", "").Replace("H", "").ToString();
            return (1d + (peptideSequence.Length - substitutedString.Length));
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
        private static double PredictedChargeCorrected(string peptideSequence)
        {
            double runningSum = 0;
            string internalString = peptideSequence.Substring(3, peptideSequence.Length - 5);

            char firstAA = peptideSequence[0];
            if (firstAA == 'D')
            {
                runningSum -= 0.26741;
            }
            else if (firstAA == 'E')
            {
                runningSum -= 0.06852;
            }
            else if (firstAA == 'N')
            {
                runningSum += 0.011699;
            }
            else
            {
                //change nothing
            }

            char secondAA = peptideSequence[1];
            if (secondAA == 'D')
            {
                runningSum -= 0.10947;
            }
            else if (secondAA == 'E')
            {
                runningSum -= 0.04011;
            }
            else if (secondAA == 'N')
            {
                runningSum += 0.012535;
            }
            else if (secondAA == 'Q')
            {
                runningSum += 0.011699;
            }
            else
            {
                //change nothing
            }

            char thirdAA = peptideSequence[2];
            if (thirdAA == 'D')
            {
                runningSum -= 0.08022;
            }
            else if (thirdAA == 'E')
            {
                runningSum -= 0.03426;
            }
            else if (thirdAA == 'N')
            {
                runningSum += 0.016713;
            }
            else if (thirdAA == 'Q')
            {
                runningSum += 0.00585;
            }
            else
            {
                //change nothing
            }

            char secondToLastAA = peptideSequence[peptideSequence.Length - 2];
            if (secondToLastAA == 'D')
            {
                runningSum -= 0.03844;
            }
            else if (secondToLastAA == 'E')
            {
                runningSum -= 0.01337;
            }
            else if (secondToLastAA == 'N')
            {
                runningSum += 0.026741;
            }
            else if (secondToLastAA == 'Q')
            {
                runningSum -= 0.00084;
            }
            else
            {
                //change nothing
            }

            char lastAA = peptideSequence[peptideSequence.Length - 1];
            if (lastAA == 'D')
            {
                runningSum -= 0.02256;
            }
            else if (lastAA == 'E')
            {
                runningSum -= 0.00418;
            }
            else if (lastAA == 'N')
            {
                runningSum += 0.010864;
            }
            else if (lastAA == 'Q')
            {
                runningSum -= 0.0117;
            }
            else
            {
                //change nothing
            }

            //consider internal residues
            if (internalString.Contains("D"))
            {
                runningSum -= 0.05014;
            }
            if (internalString.Contains("E"))
            {
                runningSum -= 0.01922;
            }
            if (internalString.Contains("N"))
            {
                runningSum += 0.012535;
            }
            if (internalString.Contains("Q"))
            {
                runningSum -= 0.000251;
            }

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
        private static double Offset(double correctedCharge, int length)
        {
            return 0;
            //should fit 5th order polynomical to plot of (ExperimentalElectrophoreticMobility - PredictedElectrophoreticMobility) vs. (Zc/N) where N is peptidelength.
        }
    }
}