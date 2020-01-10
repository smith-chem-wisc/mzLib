using System;
using System.Text;

namespace Proteomics.Fragmentation
{
    public struct Product
    {
        public readonly double NeutralMass;
        public readonly ProductType ProductType;
        public readonly double NeutralLoss;
        public readonly FragmentationTerminus Terminus;
        public readonly int FragmentNumber;
        public readonly int AminoAcidPosition;

        /// <summary>
        /// A product is the individual neutral fragment from an MS dissociation. A fragmentation product here contains one of the two termini (N- or C-). 
        /// The ProductType describes where along the backbone the fragmentaiton occurred (e.g. b-, y-, c-, zdot-). The neutral loss mass (if any) that 
        /// occurred from a mod on the fragment is listed as a mass. Finally the neutral mass of the whole fragment is provided.
        /// </summary>
        public Product(ProductType productType, FragmentationTerminus terminus, double neutralMass,
            int fragmentNumber, int aminoAcidPosition, double neutralLoss)
        {
            NeutralMass = neutralMass;
            ProductType = productType;
            NeutralLoss = neutralLoss;
            Terminus = terminus;
            FragmentNumber = fragmentNumber;
            AminoAcidPosition = aminoAcidPosition;
        }

        public string Annotation
        {
            get
            {
                StringBuilder sb = new StringBuilder();

                sb.Append(ProductType);

                // for "normal" fragments this is just the fragment number (e.g., the 3 in the b3 ion)
                // for diagnostic ions, it's the m/z assuming z=1
                // (e.g., a diagnostic ion with neutral mass 100 Da will be reported as the D101 fragment)
                sb.Append(FragmentNumber);

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
            return ProductType + "" + FragmentNumber + ";" + NeutralMass.ToString("F5") + "-" + string.Format("{0:0.##}", NeutralLoss);
        }

        public override bool Equals(object obj)
        {
            Product other = (Product)obj;

            return this.ProductType == other.ProductType
                && this.NeutralMass == other.NeutralMass
                && this.FragmentNumber == other.FragmentNumber
                && this.NeutralLoss == other.NeutralLoss;
        }

        public override int GetHashCode()
        {
            return NeutralMass.GetHashCode();
        }
    }
}