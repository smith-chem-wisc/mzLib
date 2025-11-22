namespace Chromatography.RetentionTimePrediction.CZE;

/// <summary>
/// Core calculations for CZE electrophoretic mobility predictions.
/// Based on: Krokhin et al., Anal Chem. 2017 Feb 7;89(3):2000-2008.
/// "Predicting Electrophoretic Mobility of Tryptic Peptides for High-Throughput CZE-MS Analysis"
/// https://www.ncbi.nlm.nih.gov/pubmed/28208305
/// </summary>
public static class CZECalculations
{
    /// <summary>
    /// Predicts electrophoretic mobility of a peptide based on sequence and mass.
    /// 
    /// Calculation described in Anal Chem. 2017 Feb 7;89(3):2000-2008.
    /// Uses Cifuentes and Poppe model (J. Chromatogr. A 1994, 680, 321−340)
    /// with coefficients aligned to experimental values (slope 1, intercept 0).
    /// </summary>
    /// <param name="peptideSequence">Base amino acid sequence</param>
    /// <param name="observedMass">Observed monoisotopic mass (Da)</param>
    /// <returns>Predicted electrophoretic mobility (10^-9 m^2/(V·s))</returns>
    public static double PredictedElectrophoreticMobility(string peptideSequence, double observedMass)
    {
        if (string.IsNullOrEmpty(peptideSequence) || observedMass <= 0)
            return -1;

        // Coefficients 3.069 and 386 align output with experimentally measured values
        // Other values from best fit model of Cifuentes and Poppe
        double correctedCharge = PredictedChargeCorrected(peptideSequence);
        double predictedMu = 3.069 + 386 * Math.Log(1.0 + 0.35 * correctedCharge) /
            (Math.Pow(observedMass, 0.411) + Offset(correctedCharge, peptideSequence.Length));

        return predictedMu;
    }

    /// <summary>
    /// Calculates base charge: +1 for N-terminal and +1 for each K, R, H
    /// </summary>
    private static double PredictedCharge(string peptideSequence)
    {
        int basicResidues = 0;
        foreach (char aa in peptideSequence)
        {
            if (aa == 'R' || aa == 'K' || aa == 'H')
                basicResidues++;
        }
        return 1.0 + basicResidues;
    }

    /// <summary>
    /// Applies position-dependent charge corrections for D, E, N, Q residues.
    /// Corrections are more significant at peptide termini than in the middle.
    /// 
    /// Correction values from table in Anal Chem. 2017 Feb 7;89(3):2000-2008.
    /// Future: could use linear algebra to estimate these per-dataset.
    /// </summary>
    private static double PredictedChargeCorrected(string peptideSequence)
    {
        if (peptideSequence.Length < 6)
            return PredictedCharge(peptideSequence); // Not enough residues for correction

        double correction = 0;

        // N-terminal corrections (positions 0, 1, 2)
        correction += GetNTerminalCorrection(peptideSequence[0], 0);
        correction += GetNTerminalCorrection(peptideSequence[1], 1);
        correction += GetNTerminalCorrection(peptideSequence[2], 2);

        // C-terminal corrections (last two positions)
        int lastIdx = peptideSequence.Length - 1;
        correction += GetCTerminalCorrection(peptideSequence[lastIdx - 1], 1);
        correction += GetCTerminalCorrection(peptideSequence[lastIdx], 0);

        // Internal residue corrections (if any exist between position 3 and length-3)
        if (peptideSequence.Length > 5)
        {
            string internalString = peptideSequence.Substring(3, peptideSequence.Length - 5);
            correction += GetInternalCorrection(internalString);
        }

        return PredictedCharge(peptideSequence) + correction;
    }

    private static double GetNTerminalCorrection(char aa, int positionFromNTerm)
    {
        // Position 0 (first residue)
        if (positionFromNTerm == 0)
        {
            return aa switch
            {
                'D' => -0.26741,
                'E' => -0.06852,
                'N' => +0.011699,
                _ => 0
            };
        }
        // Position 1 (second residue)
        else if (positionFromNTerm == 1)
        {
            return aa switch
            {
                'D' => -0.10947,
                'E' => -0.04011,
                'N' => +0.012535,
                'Q' => +0.011699,
                _ => 0
            };
        }
        // Position 2 (third residue)
        else if (positionFromNTerm == 2)
        {
            return aa switch
            {
                'D' => -0.08022,
                'E' => -0.03426,
                'N' => +0.016713,
                'Q' => +0.00585,
                _ => 0
            };
        }

        return 0;
    }

    private static double GetCTerminalCorrection(char aa, int positionFromCTerm)
    {
        // Position 0 (last residue)
        if (positionFromCTerm == 0)
        {
            return aa switch
            {
                'D' => -0.02256,
                'E' => -0.00418,
                'N' => +0.010864,
                'Q' => -0.0117,
                _ => 0
            };
        }
        // Position 1 (second to last)
        else if (positionFromCTerm == 1)
        {
            return aa switch
            {
                'D' => -0.03844,
                'E' => -0.01337,
                'N' => +0.026741,
                'Q' => -0.00084,
                _ => 0
            };
        }

        return 0;
    }

    private static double GetInternalCorrection(string internalSequence)
    {
        double correction = 0;

        // Apply correction if ANY of these residues are present (not per occurrence)
        if (internalSequence.Contains('D')) correction -= 0.05014;
        if (internalSequence.Contains('E')) correction -= 0.01922;
        if (internalSequence.Contains('N')) correction += 0.012535;
        if (internalSequence.Contains('Q')) correction -= 0.000251;

        return correction;
    }

    /// <summary>
    /// Offset adjustment (currently disabled).
    /// 
    /// The offset in the AC paper is a 5th order polynomial best fit to a plot of Zc/N 
    /// versus the difference between experimental and predicted electrophoretic mobility.
    /// This could be implemented as a per-dataset calibration in the future.
    /// </summary>
    private static double Offset(double correctedCharge, int length)
    {
        return 0;
        // Future: fit 5th order polynomial to plot of 
        // (ExperimentalMobility - PredictedMobility) vs. (Zc/N) where N is peptide length
    }
}