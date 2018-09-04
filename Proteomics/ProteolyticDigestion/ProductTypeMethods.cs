using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;

namespace Proteomics.ProteolyticDigestion
{
    public static class ProductTypeMethods
    {
        public static FragmentationTerminus IdentifyTerminusType(List<ProductType> productTypes)
        {
            if ((productTypes.Contains(ProductType.b) || productTypes.Contains(ProductType.c) || productTypes.Contains(ProductType.aDegree))
                && (productTypes.Contains(ProductType.y) || productTypes.Contains(ProductType.zPlusOne) || productTypes.Contains(ProductType.x)))
            {
                return FragmentationTerminus.Both;
            }
            else if (productTypes.Contains(ProductType.y) || productTypes.Contains(ProductType.zPlusOne) || productTypes.Contains(ProductType.x))
            {
                return FragmentationTerminus.C;
            }
            else //"lp.Contains(ProductType.B) || lp.Contains(ProductType.BnoB1ions) || lp.Contains(ProductType.C) || lp.Contains(ProductType.Adot))"
            {
                return FragmentationTerminus.N;
            }
        }

        public static List<List<ProductType>> SeparateIonsByTerminus(List<ProductType> ionTypes)
        {
            List<ProductType> nIons = new List<ProductType>();
            List<ProductType> cIons = new List<ProductType>();
            foreach (ProductType productType in ionTypes)
            {
                if (productType == ProductType.b || productType == ProductType.c)
                {
                    nIons.Add(productType);
                }
                else // Y and Z
                {
                    cIons.Add(productType);
                }
            }
            if (nIons.Count != 0 && cIons.Count != 0)
            {
                return new List<List<ProductType>> { nIons, cIons };
            }
            else if (nIons.Count != 0)
            {
                return new List<List<ProductType>> { nIons };
            }
            else if (cIons.Count != 0)
            {
                return new List<List<ProductType>> { cIons };
            }
            else
            {
                throw new ArgumentException("No ions types were selected.");
            }
        }
    }
}