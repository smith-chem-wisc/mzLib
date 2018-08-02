using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;

namespace Proteomics.ProteolyticDigestion
{
    public static class ProductTypeMethods
    {
        public static FragmentationTerminus IdentifyTerminusType(List<ProductType> productTypes)
        {
            if ((productTypes.Contains(ProductType.B) || productTypes.Contains(ProductType.C) || productTypes.Contains(ProductType.Adot))
                && (productTypes.Contains(ProductType.Y) || productTypes.Contains(ProductType.Z) || productTypes.Contains(ProductType.X)))
            {
                return FragmentationTerminus.Both;
            }
            else if (productTypes.Contains(ProductType.Y) || productTypes.Contains(ProductType.Z) || productTypes.Contains(ProductType.X))
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
                if (productType == ProductType.B || productType == ProductType.C)
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