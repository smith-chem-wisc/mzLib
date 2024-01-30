using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    internal class RtInfo
    {
        public double PredictedRt { get; }
        public double Width { get; }
        public double? RtSd { get; }
        public double? RtInterquartileRange { get; }

        internal RtInfo(double predictedRt, double width, double? rtSd, double? rtInterquartileRange)
        {
            PredictedRt = predictedRt;
            Width = width;
            RtSd = rtSd;
            RtInterquartileRange = rtInterquartileRange;
        }
    }
}
