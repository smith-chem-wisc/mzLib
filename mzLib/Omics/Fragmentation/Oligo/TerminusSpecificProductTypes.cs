using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Omics.Fragmentation.Oligo
{
    public static class TerminusSpecificProductTypes
    {
        public static List<ProductType> GetRnaTerminusSpecificProductTypes(
            this FragmentationTerminus fragmentationTerminus)
        {
            return ProductIonTypesFromSpecifiedTerminus[fragmentationTerminus];
        }

        /// <summary>
        /// The types of ions that can be generated from an oligo fragment, based on the terminus of the fragment
        /// </summary>
        public static Dictionary<FragmentationTerminus, List<ProductType>> ProductIonTypesFromSpecifiedTerminus = new Dictionary<FragmentationTerminus, List<ProductType>>
        {
            { 
                FragmentationTerminus.FivePrime, new List<ProductType>
                {
                    ProductType.a, ProductType.aWaterLoss, ProductType.aBaseLoss,
                    ProductType.b, ProductType.bWaterLoss, ProductType.bBaseLoss,
                    ProductType.c, ProductType.cWaterLoss, ProductType.cBaseLoss,
                    ProductType.d, ProductType.dWaterLoss, ProductType.dBaseLoss,
                }
            },
            { 
                FragmentationTerminus.ThreePrime, new List<ProductType>
                {
                    ProductType.w, ProductType.wWaterLoss, ProductType.wBaseLoss,
                    ProductType.x, ProductType.xWaterLoss, ProductType.xBaseLoss,
                    ProductType.y, ProductType.yWaterLoss, ProductType.yBaseLoss,
                    ProductType.z, ProductType.zWaterLoss, ProductType.zBaseLoss,
                }
            },
            { 
                FragmentationTerminus.Both, new List<ProductType>
                {

                    ProductType.a, ProductType.aWaterLoss, ProductType.aBaseLoss,
                    ProductType.b, ProductType.bWaterLoss, ProductType.bBaseLoss,
                    ProductType.c, ProductType.cWaterLoss, ProductType.cBaseLoss,
                    ProductType.d, ProductType.dWaterLoss, ProductType.dBaseLoss,
                    ProductType.w, ProductType.wWaterLoss, ProductType.wBaseLoss,
                    ProductType.x, ProductType.xWaterLoss, ProductType.xBaseLoss,
                    ProductType.y, ProductType.yWaterLoss, ProductType.yBaseLoss,
                    ProductType.z, ProductType.zWaterLoss, ProductType.zBaseLoss,
                    ProductType.M
                }

            },
            { 
                FragmentationTerminus.None, new List<ProductType>()
            }
        };


        public static FragmentationTerminus GetRnaTerminusType(this ProductType fragmentType)
        {
            switch (fragmentType)
            {
                case ProductType.a:
                case ProductType.aWaterLoss:
                case ProductType.aBaseLoss:
                case ProductType.b:
                case ProductType.bWaterLoss:
                case ProductType.bBaseLoss:
                case ProductType.c:
                case ProductType.cWaterLoss:
                case ProductType.cBaseLoss:
                case ProductType.d:
                case ProductType.dWaterLoss:
                case ProductType.dBaseLoss:
                case ProductType.w:
                case ProductType.wWaterLoss:
                case ProductType.wBaseLoss:
                case ProductType.x:
                case ProductType.xWaterLoss:
                case ProductType.xBaseLoss:
                case ProductType.y:
                case ProductType.yWaterLoss:
                case ProductType.yBaseLoss:
                case ProductType.z:
                case ProductType.zWaterLoss:
                case ProductType.zBaseLoss:
                case ProductType.M:
                    return ProductTypeToFragmentationTerminus[fragmentType];

                case ProductType.aStar:
                case ProductType.aDegree:
                case ProductType.bAmmoniaLoss:
                case ProductType.yAmmoniaLoss:
                case ProductType.zPlusOne:
                case ProductType.D:
                case ProductType.Ycore:
                case ProductType.Y:
                default:
                    throw new ArgumentOutOfRangeException(nameof(fragmentType), fragmentType, null);
            }
        }


        /// <summary>
        /// The terminus of the oligo fragment that the product ion is generated from
        /// </summary>
        public static Dictionary<ProductType, FragmentationTerminus> ProductTypeToFragmentationTerminus = new Dictionary<ProductType, FragmentationTerminus>
        {
            { ProductType.a, FragmentationTerminus.FivePrime },
            { ProductType.aWaterLoss, FragmentationTerminus.FivePrime },
            { ProductType.aBaseLoss, FragmentationTerminus.FivePrime },
            { ProductType.b, FragmentationTerminus.FivePrime },
            { ProductType.bWaterLoss, FragmentationTerminus.FivePrime },
            { ProductType.bBaseLoss, FragmentationTerminus.FivePrime },
            { ProductType.c, FragmentationTerminus.FivePrime },
            { ProductType.cWaterLoss, FragmentationTerminus.FivePrime },
            { ProductType.cBaseLoss, FragmentationTerminus.FivePrime },
            { ProductType.d, FragmentationTerminus.FivePrime },
            { ProductType.dWaterLoss, FragmentationTerminus.FivePrime },
            { ProductType.dBaseLoss, FragmentationTerminus.FivePrime },

            { ProductType.w, FragmentationTerminus.ThreePrime },
            { ProductType.wWaterLoss, FragmentationTerminus.ThreePrime },
            { ProductType.wBaseLoss, FragmentationTerminus.ThreePrime },
            { ProductType.x, FragmentationTerminus.ThreePrime },
            { ProductType.xWaterLoss, FragmentationTerminus.ThreePrime },
            { ProductType.xBaseLoss, FragmentationTerminus.ThreePrime },
            { ProductType.y, FragmentationTerminus.ThreePrime },
            { ProductType.yWaterLoss, FragmentationTerminus.ThreePrime },
            { ProductType.yBaseLoss, FragmentationTerminus.ThreePrime },
            { ProductType.z, FragmentationTerminus.ThreePrime },
            { ProductType.zWaterLoss, FragmentationTerminus.ThreePrime },
            { ProductType.zBaseLoss, FragmentationTerminus.ThreePrime },

            { ProductType.M, FragmentationTerminus.Both }
        };
    }
}
