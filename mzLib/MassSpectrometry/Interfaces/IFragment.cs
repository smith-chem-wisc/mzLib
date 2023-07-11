using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;

namespace MassSpectrometry
{
    public interface IFragment<TProductType, TTerminusType> : IHasMass, IEquatable<IFragment<TProductType, TTerminusType>>
        where TProductType : struct, System.Enum
        where TTerminusType : struct, System.Enum
    {
        public double NeutralMass { get; }
        public int FragmentNumber { get; }
        public int ResiduePosition { get; }
        public TProductType ProductType { get; }
        public TTerminusType Terminus { get; }
        public double NeutralLoss { get; }
        public TProductType? SecondaryProductType { get; }
        public int? SecondaryFragmentNumber { get; }

        public virtual string Annotation
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
        public string ToString()
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

    }
}
