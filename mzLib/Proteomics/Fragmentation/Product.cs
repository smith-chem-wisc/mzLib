using System;
using System.Text;
using Chemistry;
using MassSpectrometry;

namespace Proteomics.Fragmentation
{
    public readonly struct Product : IFragment<ProductType, FragmentationTerminus>
    {
        public double NeutralMass { get; }
        public  ProductType ProductType { get; }
        public double NeutralLoss { get; }
        public  FragmentationTerminus Terminus { get; }
        public int FragmentNumber { get; }
        public int AminoAcidPosition { get; } 
        public int ResiduePosition => AminoAcidPosition; // added for interface compatibility 
        public ProductType? SecondaryProductType { get; } //used for internal fragment ions
        public int? SecondaryFragmentNumber { get; } //used for internal fragment ions
        public double MonoisotopicMass => NeutralMass;

        public string Annotation => (this as IFragment<ProductType, FragmentationTerminus>).Annotation;

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
            AminoAcidPosition = aminoAcidPosition;
            SecondaryProductType = secondaryProductType;
            SecondaryFragmentNumber = secondaryFragmentNumber;
        }

        public override string ToString() => (this as IFragment<ProductType, FragmentationTerminus>).GetStringRepresentation();

        public bool Equals(IFragment<ProductType, FragmentationTerminus> other)
        {
            return NeutralMass.Equals(other.NeutralMass) && ProductType == other.ProductType &&
                   NeutralLoss.Equals(other.NeutralLoss) && Terminus == other.Terminus &&
                   FragmentNumber == other.FragmentNumber && Math.Abs(ResiduePosition - other.ResiduePosition) < 0.0001 &&
                   SecondaryProductType == other.SecondaryProductType &&
                   SecondaryFragmentNumber == other.SecondaryFragmentNumber &&
                   MonoisotopicMass.Equals(other.MonoisotopicMass);
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