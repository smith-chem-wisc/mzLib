using System.Text;
using Chemistry;

namespace Omics.Fragmentation
{
    public class Product : IHasMass, IEquatable<Product>
    {
        public double NeutralMass { get; }
        public ProductType ProductType { get; }
        public double NeutralLoss { get; }
        public FragmentationTerminus Terminus { get; }
        public int FragmentNumber { get; }
        public int ResiduePosition { get; }
        public int AminoAcidPosition => ResiduePosition;
        public ProductType? SecondaryProductType { get; } //used for internal fragment ions
        public int SecondaryFragmentNumber { get; } //used for internal fragment ions
        public double MonoisotopicMass => NeutralMass;

        /// <summary>
        /// A product is the individual neutral fragment from an MS dissociation. A fragmentation product here contains one of the two termini (N- or C-). 
        /// The ProductType describes where along the backbone the fragmentaiton occurred (e.g. b-, y-, c-, zdot-). The neutral loss mass (if any) that 
        /// occurred from a mod on the fragment is listed as a mass. Finally the neutral mass of the whole fragment is provided.
        /// </summary>
        public Product(ProductType productType, FragmentationTerminus terminus, double neutralMass,
            int fragmentNumber, int residuePosition, double neutralLoss, ProductType? secondaryProductType = null,
            int secondaryFragmentNumber = 0)
        {
            NeutralMass = neutralMass;
            ProductType = productType;
            NeutralLoss = neutralLoss;
            Terminus = terminus;
            FragmentNumber = fragmentNumber;
            ResiduePosition = residuePosition;
            SecondaryProductType = secondaryProductType;
            SecondaryFragmentNumber = secondaryFragmentNumber;
        }

        public string Annotation
        {
            get
            {
                StringBuilder sb = new StringBuilder();

                if (SecondaryProductType == null)
                {
                    sb.Append(ProductType);

                    // for "normal" fragments this is just the fragment number (e.g., the 3 in the b3 ion)
                    // for diagnostic ions, it's the m/z assuming z=1
                    // (e.g., a diagnostic ion with neutral mass 100 Da will be reported as the D101 fragment)
                    sb.Append(FragmentNumber);
                }
                else
                {
                    //internal fragment ion, annotation used here: 10.1007/s13361-015-1078-1
                    //example: yIb[18-36]
                    sb.Append(ProductType + "I" + SecondaryProductType.Value + "[" + FragmentNumber + "-" + SecondaryFragmentNumber + "]");
                }
                if (NeutralLoss != 0)
                {
                    sb.Append("-");
                    sb.Append(NeutralLoss.ToString("F2"));
                }

                return sb.ToString();
            }
        }

        /// <summary>
        /// Summarizes a Product into a string for debug purposes
        /// </summary>
        public override string ToString()
        {
            if (SecondaryProductType == null)
            {
                return ProductType + "" + FragmentNumber + ";" + NeutralMass.ToString("F5") + "-" +
                       string.Format("{0:0.##}", NeutralLoss);
            }
            else
            {
                return ProductType + "I" + SecondaryProductType.Value + "[" + FragmentNumber + "-" +
                       SecondaryFragmentNumber + "]" + ";" + NeutralMass.ToString("F5") + "-" +
                       string.Format("{0:0.##}", NeutralLoss);
            }
        }

        public override bool Equals(object obj)
        {
            return obj is Product other && Equals(other);
        }

        public bool Equals(Product? product)
        {
            return product != null
                   && NeutralMass.Equals(product.NeutralMass)
                   && ProductType == product.ProductType
                   && NeutralLoss.Equals(product.NeutralLoss)
                   && Terminus == product.Terminus
                   && FragmentNumber == product.FragmentNumber
                   && ResiduePosition == product.ResiduePosition
                   && SecondaryProductType == product.SecondaryProductType
                   && SecondaryFragmentNumber == product.SecondaryFragmentNumber
                   && MonoisotopicMass.Equals(product.MonoisotopicMass);
        }

        public override int GetHashCode()
        {
            return NeutralMass.GetHashCode();
        }
    }
}
