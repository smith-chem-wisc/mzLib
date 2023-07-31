using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Xml;
using Chemistry;
using MassSpectrometry;

namespace Proteomics.Fragmentation
{
    // 230731 - Changed this from a struct to a class, may need to revise in the future if large performance losses occur
    public class Product : IProduct
    {
        public double NeutralMass { get; }
        public ProductType ProductType { get; }
        public double NeutralLoss { get; }
        public FragmentationTerminus Terminus { get; }
        public int FragmentNumber { get; }
        public int AminoAcidPosition { get; }
        public int ResiduePosition => AminoAcidPosition; // added for interface compatibility 
        public ProductType? SecondaryProductType { get; } //used for internal fragment ions
        public int SecondaryFragmentNumber { get; } //used for internal fragment ions
        public double MonoisotopicMass => NeutralMass;
        public string Annotation => (this as IProduct).GetAnnotation();

        /// <summary>
        /// A product is the individual neutral fragment from an MS dissociation. A fragmentation product here contains one of the two termini (N- or C-). 
        /// The ProductType describes where along the backbone the fragmentaiton occurred (e.g. b-, y-, c-, zdot-). The neutral loss mass (if any) that 
        /// occurred from a mod on the fragment is listed as a mass. Finally the neutral mass of the whole fragment is provided.
        /// </summary>
        public Product(ProductType productType, FragmentationTerminus terminus, double neutralMass,
            int fragmentNumber, int aminoAcidPosition, double neutralLoss, ProductType? secondaryProductType = null, 
            int secondaryFragmentNumber = 0)
        {
            NeutralMass = neutralMass;
            ProductType = productType;
            NeutralLoss = neutralLoss;
            Terminus = terminus;
            FragmentNumber = fragmentNumber;
            AminoAcidPosition = aminoAcidPosition;
            SecondaryProductType = secondaryProductType;
            SecondaryFragmentNumber = secondaryFragmentNumber;
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

        public bool Equals(IProduct other)
        {
            return NeutralMass.Equals(other.NeutralMass) 
                   && ProductType == other.ProductType
                   && NeutralLoss.Equals(other.NeutralLoss) 
                   && Terminus == other.Terminus 
                   && FragmentNumber == other.FragmentNumber 
                   && ResiduePosition == other.ResiduePosition
                   && SecondaryProductType == other.SecondaryProductType
                   && SecondaryFragmentNumber == other.SecondaryFragmentNumber 
                   && MonoisotopicMass.Equals(other.MonoisotopicMass);
        }

        public override bool Equals(object obj)
        {
            return obj is Product other && Equals(other);
        }

        public override int GetHashCode()
        {
            var hashCode = new HashCode();
            hashCode.Add(NeutralMass);
            hashCode.Add((int)ProductType);
            hashCode.Add(NeutralLoss);
            hashCode.Add((int)Terminus);
            hashCode.Add(FragmentNumber);
            hashCode.Add(AminoAcidPosition);
            hashCode.Add(SecondaryProductType);
            hashCode.Add(SecondaryFragmentNumber);
            hashCode.Add(MonoisotopicMass);
            return hashCode.ToHashCode();
        }
    }
}