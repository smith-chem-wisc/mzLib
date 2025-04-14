using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    internal interface ISeparable
    {
        // Separation domain can be time (RT) or ion mobility (IMS)
        double SeparationDomainValue { get; }
        double Intensity { get; }
        int ZeroBasedScanIndex { get; }
    }
}
