using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public enum FragmentationTerminus
    {
        Both, //N- and C-terminus
        N, //N-terminus only
        C, //C-terminus only
        None, //used for internal fragments, could be used for top down intact mass?
        FivePrime, // 5' for NucleicAcids
        ThreePrime, // 3' for NucleicAcids
    }
}
