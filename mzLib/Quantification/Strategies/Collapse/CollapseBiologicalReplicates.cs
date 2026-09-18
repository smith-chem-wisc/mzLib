namespace Quantification.Strategies
{
    /// <summary>
    /// Combines the biological replicates of a condition, keeping fraction and technical replicate
    /// distinct. Less common, but useful for condition-level summaries.
    /// </summary>
    public class CollapseBiologicalReplicates : CollapseStrategyBase
    {
        public override string Name => "Collapse Biological Replicates";

        protected override bool CollapsesBiologicalReplicate => true;
        protected override bool CollapsesFraction => false;
        protected override bool CollapsesTechnicalReplicate => false;
    }
}
