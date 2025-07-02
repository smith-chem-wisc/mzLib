using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class PeakFindingParameters
    {
        public int BinSize { get; set; }
        public double MinMass { get; set; } // for neutral mass indexing only
        public int MinCharge { get; set; } // for neutral mass indexing only
        public Tolerance PeakFindingTolerance { get; set; }
        public double MaxPeakHalfWidth { get; set; }
        public int MinPeakCount { get; set; }
        public int MaxMissedScansAllowed { get; set; }

        public PeakFindingParameters(Tolerance peakFindingTolerance, int binSize = 100, double minMass = 0, int minCharge = 1, double maxPeakHalfWidth = 1, int minPeakCount = 3, int maxMissedScansAllowed = 2)
        {
            BinSize = binSize;
            MinMass = minMass;
            MinCharge = minCharge;
            PeakFindingTolerance = peakFindingTolerance;
            MaxPeakHalfWidth = maxPeakHalfWidth;
            MinPeakCount = minPeakCount;
            MaxMissedScansAllowed = maxMissedScansAllowed;
        }
    }
}
