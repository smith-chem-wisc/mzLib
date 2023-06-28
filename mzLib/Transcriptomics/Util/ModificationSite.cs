using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Transcriptomics
{
    [Flags]
    public enum ModificationSite
    {
        // bit shift operations to quickly flag relevant modifications
        None = 0,
        A = 1 << 0,
        C = 1 << 1,
        G = 1 << 2,
        U = 1 << 3,
        FivePrimeTerminus = 1 << 4,
        ThreePrimeTerminus = 1 << 5,
        All = (1 << 6) - 1, // Handy way of setting all below the 5th bit
        Any = 1 << 7 // Acts like none, but is equal to all
    }

    public static class ModificationSiteExtensions
    {

    }
}
