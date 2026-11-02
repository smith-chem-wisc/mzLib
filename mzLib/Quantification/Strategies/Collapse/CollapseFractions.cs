namespace Quantification.Strategies
{
    /// <summary>
    /// Combines the fractions of a sample, keeping condition, biological replicate and technical
    /// replicate distinct. The usual first step for a fractionated experiment.
    /// </summary>
    public class CollapseFractions : CollapseStrategyBase
    {
        public override string Name => "Collapse Fractions";

        protected override bool CollapsesBiologicalReplicate => false;
        protected override bool CollapsesFraction => true;
        protected override bool CollapsesTechnicalReplicate => false;
    }
}
