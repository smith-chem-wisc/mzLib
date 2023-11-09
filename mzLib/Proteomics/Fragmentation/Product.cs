using System;
using System.Text;
using Omics.Fragmentation;

namespace Proteomics.Fragmentation
{
    public class Product : IProduct
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
            int fragmentNumber, int aminoAcidPosition, double neutralLoss, ProductType? secondaryProductType = null, int secondaryFragmentNumber = 0)
        {
            NeutralMass = neutralMass;
            ProductType = productType;
            NeutralLoss = neutralLoss;
            Terminus = terminus;
            FragmentNumber = fragmentNumber;
            ResiduePosition = aminoAcidPosition;
            SecondaryProductType = secondaryProductType;
            SecondaryFragmentNumber = secondaryFragmentNumber;
        }
        public string Annotation => (this as IProduct).GetAnnotation();
        

        /// <summary>
        /// Summarizes a Product into a string for debug purposes
        /// </summary>
        public override string ToString()
        {
            if (SecondaryProductType == null)
            {
                return ProductType + "" + FragmentNumber + ";" + NeutralMass.ToString("F5") + "-" + string.Format("{0:0.##}", NeutralLoss);
            }
            else
            {
                return ProductType + "I" + SecondaryProductType.Value + "[" + FragmentNumber + "-" + SecondaryFragmentNumber + "]" + ";" + NeutralMass.ToString("F5") + "-" + string.Format("{0:0.##}", NeutralLoss);
            }
        }

        public override bool Equals(object obj)
        {
            return obj is Product other && Equals(other);
        }

        public bool Equals(IProduct product)
        {
            return this.ProductType.Equals(product.ProductType)
                   && this.NeutralMass.Equals(product.NeutralMass)
                   && this.FragmentNumber == product.FragmentNumber
                   && this.NeutralLoss.Equals(product.NeutralLoss)
                   && this.SecondaryFragmentNumber == product.SecondaryFragmentNumber
                   && this.SecondaryProductType == product.SecondaryProductType;
        }

        public override int GetHashCode()
        {
            return NeutralMass.GetHashCode();
        }
    }
}