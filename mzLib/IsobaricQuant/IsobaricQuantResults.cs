using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IsobaricQuant
{
    internal class IsobaricQuantResults
    {
        public readonly List<TmtFileInfo> TmtSpectraFiles;


        public IsobaricQuantResults(List<TmtFileInfo> tmtSpectraFiles)
        {
            TmtSpectraFiles = tmtSpectraFiles;
        }
    }
}
