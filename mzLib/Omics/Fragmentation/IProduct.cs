using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Omics.Fragmentation
{
    public interface IProduct
    {
        double NeutralMass { get; }
        ProductType ProductType { get; }
        double NeutralLoss { get; }
        FragmentationTerminus Terminus { get; }
        int FragmentNumber { get; }
        int Position { get; }
        ProductType? SecondaryProductType { get; } //used for internal fragments
        int SecondaryFragmentNumber { get; } //used for internal fragment ions
        string Annotation { get; }

    }
}
