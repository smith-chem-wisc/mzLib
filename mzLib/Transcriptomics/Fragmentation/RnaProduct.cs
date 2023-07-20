using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;

namespace Transcriptomics
{
    public class RnaProduct : IProduct
    {
        public double NeutralMass { get; }
        public ProductType ProductType { get; }
        public double NeutralLoss { get; }
        public FragmentationTerminus Terminus { get; }
        public int FragmentNumber { get; }
        public int ResiduePosition { get; }
        public ProductType? SecondaryProductType { get; }
        public int SecondaryFragmentNumber { get; }
        public double MonoisotopicMass => NeutralMass;
        public string Annotation => (this as IProduct).GetAnnotation();
        public NucleicAcid? Parent { get; }

        public RnaProduct(ProductType productType, FragmentationTerminus terminus, double neutralMass,
            int fragmentNumber, int residuePosition, double neutralLoss, ProductType? secondaryProductType = null, 
            int secondaryFragmentNumber = 0, NucleicAcid? parent = null)
        {
            NeutralMass = neutralMass;
            ProductType = productType;
            NeutralLoss = neutralLoss;
            Terminus = terminus;
            FragmentNumber = fragmentNumber;
            ResiduePosition = residuePosition;
            SecondaryProductType = secondaryProductType;
            SecondaryFragmentNumber = secondaryFragmentNumber;
            Parent = parent;
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
            if (other is null) return false;
            return NeutralMass.Equals(other.NeutralMass) && ProductType == other.ProductType &&
                   NeutralLoss.Equals(other.NeutralLoss) && Terminus == other.Terminus &&
                   FragmentNumber == other.FragmentNumber && ResiduePosition == other.ResiduePosition &&
                   SecondaryProductType == other.SecondaryProductType &&
                   SecondaryFragmentNumber == other.SecondaryFragmentNumber;
        }

        public override bool Equals(object? obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;
            return Equals((RnaProduct)obj);
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(NeutralMass, (int)ProductType, NeutralLoss, (int)Terminus, FragmentNumber, ResiduePosition, SecondaryProductType, SecondaryFragmentNumber);
        }
    }
}
