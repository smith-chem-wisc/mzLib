using System.Collections.Generic;

namespace MassSpectrometry
{

    public class DissociationTypeCollection
    {
        public List<ProductType> ProductTypes { get; private set; }

        public static Dictionary<DissociationType, List<ProductType>> ProductsFromDissociationType = new Dictionary<DissociationType, List<ProductType>>
        {
            {DissociationType.HCD, new List<ProductType>{ ProductType.BnoB1ions, ProductType.Y } },
            {DissociationType.ECD, new List<ProductType>{ ProductType.C, ProductType.Y, ProductType.Z } },
            {DissociationType.ETD, new List<ProductType>{ ProductType.C, ProductType.Y, ProductType.Z } },
            {DissociationType.CID, new List<ProductType>{ ProductType.B, ProductType.Y } },
            {DissociationType.EThCD, new List<ProductType>{ ProductType.BnoB1ions, ProductType.Y } },
            {DissociationType.AnyActivationType, new List<ProductType>{ ProductType.B, ProductType.Y } },
            {DissociationType.Custom, new List<ProductType>() }
        };

    }
}