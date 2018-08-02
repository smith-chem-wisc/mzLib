using System;
using System.Collections.Generic;
using System.Text;

namespace Proteomics.Fragmentation
{
    public class TerminusSpecificProductTypes
    {

        public static Dictionary<FragmentationTerminus, List<ProductType>> ProductIonTypesFromSpecifiedTerminus = new Dictionary<FragmentationTerminus, List<ProductType>>
        {
            {FragmentationTerminus.N, new List<ProductType>{ ProductType.A, ProductType.Adot, ProductType.Astar, ProductType.B, ProductType.Bdot, ProductType.Bstar, ProductType.C } }, //all ion types that include the N-terminus
            {FragmentationTerminus.C, new List<ProductType>{ ProductType.X, ProductType.Y, ProductType.Ydot, ProductType.Ystar, ProductType.Z } }, //all ion types that include the C-terminus
            {FragmentationTerminus.Both, new List<ProductType>{ ProductType.A, ProductType.Adot, ProductType.Astar, ProductType.B, ProductType.Bdot, ProductType.Bstar, ProductType.C, ProductType.X, ProductType.Y, ProductType.Ydot, ProductType.Ystar, ProductType.Z} },
            {FragmentationTerminus.None, new List<ProductType>() }
        };

    }
}
