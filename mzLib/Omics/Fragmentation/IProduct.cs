using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Chemistry;

namespace Omics.Fragmentation
{
    public interface IProduct : IHasMass, IEquatable<IProduct>
    {
        double NeutralMass { get; }
        ProductType ProductType { get; }
        double NeutralLoss { get; }
        FragmentationTerminus Terminus { get; }
        int FragmentNumber { get; }
        int ResiduePosition { get; }
        int AminoAcidPosition => ResiduePosition;
        ProductType? SecondaryProductType { get; } //used for internal fragments
        int SecondaryFragmentNumber { get; } //used for internal fragment ions
        string Annotation => GetAnnotation();
        string GetAnnotation()
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

        string ToString();
        int GetHashCode();
    }
}
