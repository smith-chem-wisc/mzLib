using System.Text;

namespace Proteomics.Fragmentation
{
    public class Product
    {
        public readonly double NeutralMass;
        public readonly ProductType ProductType;
        public readonly NeutralTerminusFragment TerminusFragment;
        public readonly double NeutralLoss;

        /// <summary>
        /// A product is the individual neutral fragment from an MS dissociation. A fragmentation product here contains one of the two termini (N- or C-). 
        /// The ProductType describes where along the backbone the fragmentaiton occurred (e.g. b-, y-, c-, zdot-). The neutral loss mass (if any) that 
        /// occurred from a mod on the fragment is listed as a mass. Finally the neutral mass of the whole fragment is provided.
        /// </summary>
        public Product(ProductType productType, NeutralTerminusFragment terminusFragment, double neutralLoss)
        {
            TerminusFragment = terminusFragment;
            ProductType = productType;
            NeutralLoss = neutralLoss;
            NeutralMass = DissociationTypeCollection.ProductTypeSpecificFragmentNeutralMass(terminusFragment.NeutralMass, productType) - neutralLoss;
        }

        public string Annotation
        {
            get
            {
                StringBuilder sb = new StringBuilder();
                bool containsNeutralLoss = false;

                if (NeutralLoss != 0)
                {
                    containsNeutralLoss = true;
                }

                if (containsNeutralLoss)
                {
                    sb.Append("(");
                }

                sb.Append(ProductType);
                sb.Append(TerminusFragment.FragmentNumber);

                if (containsNeutralLoss)
                {
                    sb.Append("-");
                    sb.Append(NeutralLoss.ToString("F2"));
                    sb.Append(")");
                }

                return sb.ToString();
            }
        }

        /// <summary>
        /// Summarizes a Product into a string for debug purposes
        /// </summary>
        public override string ToString()
        {
            return ProductType + "" + TerminusFragment.FragmentNumber + ";" + NeutralMass + "-" + NeutralLoss;
        }

        public override bool Equals(object obj)
        {
            Product other = (Product)obj;

            return this.ProductType == other.ProductType 
                && this.TerminusFragment == other.TerminusFragment 
                && this.NeutralLoss == other.NeutralLoss;
        }

        public override int GetHashCode()
        {
            return NeutralMass.GetHashCode();
        }
    }
}