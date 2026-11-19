namespace Quantification.Strategies
{
    /// <summary>
    /// Rolls up by summing. The natural choice for reporter-ion intensities, where a peptide's value
    /// is the total signal its spectral matches carried.
    ///
    /// Biased toward entities with more identified lower-level features -- a protein with twice the
    /// peptides reports roughly twice the intensity -- which is intended for isobaric data and is the
    /// reason <see cref="MedianRollUp"/> exists as an alternative.
    /// </summary>
    public class SumRollUp : AggregatingRollUp
    {
        public SumRollUp() : base(new SumAggregation())
        {
        }
    }
}
