namespace Quantification.Strategies
{
    /// <summary>
    /// Rolls up by taking the largest observed value in each column. Zero means "not observed" and is
    /// excluded, so a column with no observations rolls up to zero.
    ///
    /// This is how FlashLFQ sets a peptide's intensity in a file: the most intense of its chromatographic
    /// peaks, not their sum. Rolling peaks (as spectral matches) up to peptides with it reproduces
    /// FlashLFQ's peptide intensities.
    /// </summary>
    public class MaxRollUp : AggregatingRollUp
    {
        public MaxRollUp() : base(new MaxAggregation())
        {
        }
    }
}
