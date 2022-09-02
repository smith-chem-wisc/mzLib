using System;
using System.Collections.Generic;
using System.Text;

namespace Proteomics.Fragmentation
{
    public class TerminusSpecificProductTypes
    {
        public static Dictionary<FragmentationTerminus, List<ProductType>> ProductIonTypesFromSpecifiedTerminus = new Dictionary<FragmentationTerminus, List<ProductType>>
        {
            {FragmentationTerminus.N, new List<ProductType>{ ProductType.a, ProductType.a_H2O, ProductType.a_NH3, ProductType.b, ProductType.b_H2O, ProductType.b_NH3, ProductType.c } }, //all ion types that include the N-terminus
            {FragmentationTerminus.C, new List<ProductType>{ ProductType.x, ProductType.y, ProductType.YH2O, ProductType.y_NH3, ProductType.zDot, ProductType.zPlusOne } }, //all ion types that include the C-terminus
            {FragmentationTerminus.Both, new List<ProductType>{ ProductType.a, ProductType.a_H2O, ProductType.a_NH3, ProductType.b, ProductType.b_H2O, ProductType.b_NH3, ProductType.c, ProductType.x, ProductType.y, ProductType.YH2O, ProductType.y_NH3, ProductType.zDot, ProductType.zPlusOne} },
            {FragmentationTerminus.None, new List<ProductType>() }
        };

        public static Dictionary<ProductType, FragmentationTerminus> ProductTypeToFragmentationTerminus = new Dictionary<ProductType, FragmentationTerminus>
        {
            { ProductType.a, FragmentationTerminus.N },
            { ProductType.a_H2O, FragmentationTerminus.N },
            { ProductType.a_NH3, FragmentationTerminus.N },
            { ProductType.b, FragmentationTerminus.N },
            { ProductType.b_H2O, FragmentationTerminus.N },
            { ProductType.b_NH3, FragmentationTerminus.N },
            { ProductType.c, FragmentationTerminus.N },
            { ProductType.x, FragmentationTerminus.C },
            { ProductType.y, FragmentationTerminus.C },
            { ProductType.YH2O, FragmentationTerminus.C },
            { ProductType.y_NH3, FragmentationTerminus.C },
            { ProductType.zDot, FragmentationTerminus.C },
            { ProductType.zPlusOne, FragmentationTerminus.C },
        };

    }
}
