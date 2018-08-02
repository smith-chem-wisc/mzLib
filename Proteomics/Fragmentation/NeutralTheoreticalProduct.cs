using Chemistry;

namespace Proteomics.Fragmentation
{
    public class NeutralTheoreticalProduct
    {
        protected static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        
        public readonly double NeutralMass;
        public readonly ProductType ProductType;
        public readonly NeutralTerminusFragment TerminusFragment;
        public readonly string NeutralLossLabel;

        /// <summary>
        ///
        /// </summary>
        public NeutralTheoreticalProduct(ProductType productType, NeutralTerminusFragment terminusFragment, string neutralLossLabel)
        {
            TerminusFragment = terminusFragment;
            ProductType = productType;
            this.NeutralLossLabel = neutralLossLabel;
            NeutralMass = DissociationTypeCollection.ProductTypeSpecificFragmentNeutralMass(terminusFragment.NeutralMass, productType);
        }

        /// <summary>
        /// Summarizes a TheoreticalFragmentIon into a string for debug purposes
        /// </summary>
        public override string ToString()
        {
            return ProductType + "" + TerminusFragment.AminoAcidPositionNumber + ";" + NeutralMass + "-" + NeutralLossLabel;
        }
    }
}