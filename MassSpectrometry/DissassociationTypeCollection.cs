using System.Collections.Generic;

namespace MassSpectrometry
{
    public enum DissassociationType
    {
        HCD,
        ECD,
        ETD,
        CID,
        EThCD,
        Any,
        Custom
    }

    public class DissassociationTypeCollection
    {
        public List<ProductType> ProductTypes { get; private set; }

        public DissassociationTypeCollection(DissassociationType dissassociationType, List<ProductType> customPtypes = null)
        {
            switch (dissassociationType)
            {
                case DissassociationType.HCD:
                    this.ProductTypes = HCD_ProductTypes();
                    break;

                case DissassociationType.ECD:
                    this.ProductTypes = ECD_ProductTypes();
                    break;

                case DissassociationType.ETD:
                    this.ProductTypes = ETD_ProductTypes();
                    break;

                case DissassociationType.CID:
                    this.ProductTypes = CID_ProductTypes();
                    break;

                case DissassociationType.EThCD:
                    this.ProductTypes = EThCD_ProductTypes();
                    break;

                case DissassociationType.Any:
                    this.ProductTypes = Any_ProductTypes();
                    break;

                case DissassociationType.Custom:
                    this.ProductTypes = Custom_ProductTypes(customPtypes);
                    break;

                default:
                    break;
            }
        }

        private List<ProductType> HCD_ProductTypes()
        {
            List<ProductType> pTypes = new List<ProductType>
            {
                ProductType.BnoB1ions,
                ProductType.Y
            };

            return pTypes;
        }

        private List<ProductType> ECD_ProductTypes()
        {
            List<ProductType> pTypes = new List<ProductType>
            {
                ProductType.C,
                ProductType.Y,
                ProductType.Z
            };
            return pTypes;
        }

        private List<ProductType> ETD_ProductTypes()
        {
            List<ProductType> pTypes = new List<ProductType>
            {
                ProductType.C,
                ProductType.Y,
                ProductType.Z
            };
            return pTypes;
        }

        private List<ProductType> CID_ProductTypes()
        {
            List<ProductType> pTypes = new List<ProductType>
            {
                ProductType.B,
                ProductType.Y
            };
            return pTypes;
        }

        private List<ProductType> EThCD_ProductTypes()
        {
            List<ProductType> pTypes = new List<ProductType>
            {
                ProductType.BnoB1ions,
                ProductType.Y
            };
            return pTypes;
        }

        private List<ProductType> Any_ProductTypes()
        {
            List<ProductType> pTypes = new List<ProductType>
            {
                ProductType.B,
                ProductType.Y
            };
            return pTypes;
        }

        private List<ProductType> Custom_ProductTypes(List<ProductType> customPtypes)
        {
            List<ProductType> pTypes = new List<ProductType>();
            pTypes.AddRange(customPtypes);
            return pTypes;
        }
    }
}