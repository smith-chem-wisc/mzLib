using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    public class RtInfo
    {
        public double PredictedRt { get; }
        public double Width { get; set; }
        public double RtStartHypothesis => PredictedRt - (Width / 2.0);
        public double RtEndHypothesis => PredictedRt + (Width / 2.0);

        public RtInfo(double predictedRt, double width)
        {
            PredictedRt = predictedRt;
            Width = width;
        }
    }
}
