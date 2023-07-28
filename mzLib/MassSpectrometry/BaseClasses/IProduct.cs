using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;

namespace MassSpectrometry
{
    public interface IProduct : IHasMass, IEquatable<IProduct>
    {
        public double NeutralMass { get; }
        public ProductType ProductType { get; }
        public double NeutralLoss { get; }
        public FragmentationTerminus Terminus { get; }
        public int FragmentNumber { get; }
        public int ResiduePosition { get; } // added for interface compatibility 
        public int AminoAcidPosition => ResiduePosition;
        public ProductType? SecondaryProductType { get; } //used for internal fragment ions
        public int SecondaryFragmentNumber { get; } //used for internal fragment ions
        public double MonoisotopicMass => NeutralMass;
        public string Annotation => GetAnnotation();

        public string GetAnnotation()
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


}
