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
            return NeutralMass.Equals(other.NeutralMass) && ProductType == other.ProductType &&
                   NeutralLoss.Equals(other.NeutralLoss) && Terminus == other.Terminus &&
                   FragmentNumber == other.FragmentNumber && Math.Abs(ResiduePosition - other.ResiduePosition) < 0.0001 &&
                   SecondaryProductType == other.SecondaryProductType &&
                   SecondaryFragmentNumber == other.SecondaryFragmentNumber &&
                   MonoisotopicMass.Equals(other.MonoisotopicMass);
        }

        public override bool Equals(object obj)
        {
            return obj is RnaProduct other && Equals(other);
        }

        public override int GetHashCode()
        {
            var hashCode = new HashCode();
            hashCode.Add(NeutralMass);
            hashCode.Add((int)ProductType);
            hashCode.Add(NeutralLoss);
            hashCode.Add((int)Terminus);
            hashCode.Add(FragmentNumber);
            hashCode.Add(ResiduePosition);
            hashCode.Add(SecondaryProductType);
            hashCode.Add(SecondaryFragmentNumber);
            hashCode.Add(MonoisotopicMass);
            return hashCode.ToHashCode();
        }
    }
}
