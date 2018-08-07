namespace Proteomics.Fragmentation
{
    public class Product
    {
        public readonly double NeutralMass;
        public readonly ProductType ProductType;
        public readonly NeutralTerminusFragment TerminusFragment;
        public readonly double NeutralLoss;

        /// <summary>
        ///
        /// </summary>
        public Product(ProductType productType, NeutralTerminusFragment terminusFragment, double neutralLoss)
        {
            TerminusFragment = terminusFragment;
            ProductType = productType;
            this.NeutralLoss = neutralLoss;
            NeutralMass = DissociationTypeCollection.ProductTypeSpecificFragmentNeutralMass(terminusFragment.NeutralMass, productType) - neutralLoss;
        }

        /// <summary>
        /// Summarizes a TheoreticalFragmentIon into a string for debug purposes
        /// </summary>
        public override string ToString()
        {
            return ProductType + "" + TerminusFragment.FragmentNumber + ";" + NeutralMass + "-" + NeutralLoss;
        }
    }
}