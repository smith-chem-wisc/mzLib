using Chromatography.RetentionTimePrediction.Util;

namespace Chromatography.RetentionTimePrediction.CZE
{
    /// <summary>
    /// CZE-based retention time predictor using electrophoretic mobility calculations.
    /// Predicts migration times for peptides in capillary zone electrophoresis.
    /// 
    /// Based on: Krokhin et al., Anal Chem. 2017 Feb 7;89(3):2000-2008.
    /// "Predicting Electrophoretic Mobility of Tryptic Peptides for High-Throughput CZE-MS Analysis"
    /// </summary>
    public class CZERetentionTimePredictor : RetentionTimePredictorBase
    {
        private readonly double _columnLengthMeters;
        private readonly double _voltsPerMeter;

        public override string PredictorName => "CZE";
        public override SeparationType SeparationType => SeparationType.CZE;
        
        // CZE predictor works with base sequence and mass - no complex modification handling needed
        public override bool RequiresModificationChecking => false;

        /// <summary>
        /// Initializes a new CZE predictor with default parameters (100cm column, 300kV/m)
        /// </summary>
        public CZERetentionTimePredictor(
            IncompatibleModHandlingMode modHandlingMode = IncompatibleModHandlingMode.UsePrimarySequence)
            : this(modHandlingMode, columnLengthMeters: 1.0, voltsPerMeter: 300000)
        {
        }

        /// <summary>
        /// Initializes a new CZE predictor with custom instrument parameters
        /// </summary>
        /// <param name="modHandlingMode">How to handle modifications (CZE uses mass, so this mainly affects which sequence to use)</param>
        /// <param name="columnLengthMeters">Total capillary length in meters (default: 1.0m)</param>
        /// <param name="voltsPerMeter">Applied voltage gradient in V/m (default: 300,000 V/m)</param>
        public CZERetentionTimePredictor(
            IncompatibleModHandlingMode modHandlingMode,
            double columnLengthMeters,
            double voltsPerMeter)
            : base(modHandlingMode)
        {
            if (columnLengthMeters <= 0)
                throw new ArgumentException("Column length must be positive", nameof(columnLengthMeters));
            
            if (voltsPerMeter <= 0)
                throw new ArgumentException("Voltage must be positive", nameof(voltsPerMeter));

            _columnLengthMeters = columnLengthMeters;
            _voltsPerMeter = voltsPerMeter;
        }

        protected override bool ValidateBasicConstraints(IRetentionPredictable peptide, out string? failureReason)
        {
            if (string.IsNullOrEmpty(peptide.BaseSequence))
            {
                failureReason = "Empty sequence";
                return false;
            }

            // CZE predictor requires minimum length for charge correction calculations
            if (peptide.BaseSequence.Length < 6)
            {
                failureReason = "CZE predictor requires peptides of at least 6 amino acids for accurate charge correction";
                return false;
            }

            if (peptide.MonoisotopicMass <= 0)
            {
                failureReason = "Invalid monoisotopic mass";
                return false;
            }

            failureReason = null;
            return true;
        }

        protected override double? PredictCore(IRetentionPredictable peptide)
        {
            // CZE prediction uses base sequence and total mass (includes modifications)
            return PredictMigrationTime(peptide.BaseSequence, peptide.MonoisotopicMass);
        }

        protected override double? PredictWithBaseSequence(string baseSequence)
        {
            // Without mass information, we can't predict - this shouldn't be called for CZE
            // but we implement it for interface compliance
            throw new InvalidOperationException(
                "CZE prediction requires peptide mass. Use PredictCore with full IRetentionPredictable.");
        }

        /// <summary>
        /// Predicts migration time (minutes) for a peptide based on its sequence and mass
        /// </summary>
        private double? PredictMigrationTime(string baseSequence, double monoisotopicMass)
        {
            // Calculate electrophoretic mobility
            double mobility = CZECalculations.PredictedElectrophoreticMobility(baseSequence, monoisotopicMass);

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
}