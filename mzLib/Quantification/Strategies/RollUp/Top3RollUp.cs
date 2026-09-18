namespace Quantification.Strategies
{
    /// <summary>
    /// Rolls peptides up to proteins by averaging the three most intense observed peptides in each
    /// sample -- the "Top3" estimator.
    ///
    /// Unlike <see cref="SumRollUp"/>, the result does not scale with how many peptides the search
    /// happened to identify, which is what makes it usable for comparing abundance between different
    /// proteins rather than one protein across samples.
    /// </summary>
    public class Top3RollUp : AggregatingRollUp
    {
        public Top3RollUp() : base(new Top3Aggregation())
        {
        }
    }
}
