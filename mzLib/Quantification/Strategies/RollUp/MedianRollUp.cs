namespace Quantification.Strategies
{
    /// <summary>
    /// Rolls up by taking the median of the observed values in each column. Zero means "not observed"
    /// and is excluded, so a column with no observations rolls up to zero.
    ///
    /// Compared to <see cref="SumRollUp"/>, median roll-up is more robust to outliers and is not
    /// biased toward entities with more identified lower-level features.
    /// </summary>
    public class MedianRollUp : AggregatingRollUp
    {
        public MedianRollUp() : base(new MedianAggregation())
        {
        }
    }
}
