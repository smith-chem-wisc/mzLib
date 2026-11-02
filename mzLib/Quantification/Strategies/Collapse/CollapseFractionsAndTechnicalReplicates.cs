namespace Quantification.Strategies
{
    /// <summary>
    /// Combines fractions and technical replicates in one step, leaving one column per
    /// condition and biological replicate.
    ///
    /// This is the grouping the former SumCollapse and MeanCollapse both performed; pair it with
    /// <see cref="SumAggregation"/> or <see cref="MeanAggregation"/> to get their behaviour.
    /// </summary>
    public class CollapseFractionsAndTechnicalReplicates : CollapseStrategyBase
    {
        public override string Name => "Collapse Fractions and Technical Replicates";

        protected override bool CollapsesBiologicalReplicate => false;
        protected override bool CollapsesFraction => true;
        protected override bool CollapsesTechnicalReplicate => true;
    }
}
