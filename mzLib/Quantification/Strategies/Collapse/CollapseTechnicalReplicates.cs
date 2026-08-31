namespace Quantification.Strategies
{
    /// <summary>
    /// Combines the technical replicates of a sample, keeping condition, biological replicate and
    /// fraction distinct.
    /// </summary>
    public class CollapseTechnicalReplicates : CollapseStrategyBase
    {
        public override string Name => "Collapse Technical Replicates";

        protected override bool CollapsesBiologicalReplicate => false;
        protected override bool CollapsesFraction => false;
        protected override bool CollapsesTechnicalReplicate => true;
    }
}
