namespace Chromatography.RetentionTimePrediction.CZE;

/// <summary>
/// CZE-based retention time predictor using electrophoretic mobility calculations.
/// Predicts migration times for peptides in capillary zone electrophoresis.
/// 
/// Based on: Krokhin et al., Anal Chem. 2017 Feb 7;89(3):2000-2008.
/// "Predicting Electrophoretic Mobility of Tryptic Peptides for High-Throughput CZE-MS Analysis"
/// </summary>
public class CZERetentionTimePredictor : RetentionTimePredictor
{
    private readonly double _columnLengthMeters;
    private readonly double _voltsPerMeter;

    public override string PredictorName => "CZE";
    public override SeparationType SeparationType => SeparationType.CZE;

    /// <summary>
    /// Initializes a new CZE predictor with custom instrument parameters
    /// </summary>
    /// <param name="modHandlingMode">How to handle modifications (CZE uses mass, so this mainly affects which sequence to use)</param>
    /// <param name="columnLengthMeters">Total capillary length in meters (default: 1.0m)</param>
    /// <param name="voltsPerMeter">Applied voltage gradient in V/m (default: 300,000 V/m)</param>
    public CZERetentionTimePredictor(
        IncompatibleModHandlingMode modHandlingMode = IncompatibleModHandlingMode.UsePrimarySequence,
        double columnLengthMeters = 1.0,
        double voltsPerMeter = 300000)
        : base(modHandlingMode)
    {
        if (columnLengthMeters <= 0)
            throw new ArgumentException("Column length must be positive", nameof(columnLengthMeters));
        
        if (voltsPerMeter <= 0)
            throw new ArgumentException("Voltage must be positive", nameof(voltsPerMeter));

        _columnLengthMeters = columnLengthMeters;
        _voltsPerMeter = voltsPerMeter;
    }

    public override string? GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
    {
        failureReason = null;
        return peptide.BaseSequence;
    }

    protected override bool ValidateBasicConstraints(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
    {
        if (peptide.MonoisotopicMass <= 0)
        {
            failureReason = RetentionTimeFailureReason.InvalidMass;
            return false;
        }

        return base.ValidateBasicConstraints(peptide, out failureReason);
    }

    /// <summary>
    /// Predicts migration time (minutes) for a peptide based on its sequence and mass
    /// </summary>
    protected override double? PredictCore(IRetentionPredictable peptide, string? formattedSequence = null)
    {            
        // Calculate electrophoretic mobility
        double mobility = CZECalculations.PredictedElectrophoreticMobility(peptide.BaseSequence, peptide.MonoisotopicMass);

        if (mobility <= 0)
            return null;

        // Convert mobility to migration time
        return TheoreticalElutionTime(mobility);
    }

    /// <summary>
    /// Calculates theoretical elution time from electrophoretic mobility
    /// </summary>
    public double TheoreticalElutionTime(double electrophoreticMobility)
    {
        if (_columnLengthMeters <= 0 || electrophoreticMobility <= 0)
            return -1;

        return _columnLengthMeters * 1e9 / (60 * _voltsPerMeter * electrophoreticMobility);
    }

    /// <summary>
    /// Calculates experimental electrophoretic mobility from an observed migration time.
    /// This is useful for calibration or comparison with predictions.
    /// </summary>
    public double ExperimentalElectrophoreticMobility(double migrationTimeMinutes)
    {
        if (_columnLengthMeters <= 0 || migrationTimeMinutes <= 0)
            return -1;

        return _columnLengthMeters / (60 * migrationTimeMinutes * _voltsPerMeter) * 1e9;
    }
}