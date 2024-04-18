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
        public double Width { get; }
        // the Math.Max components ensure that the width of an RT Window is at least _minimumWindowWidth wide
        private double _minimumWindowWidth = 0.5;
        public double RtStartHypothesis => PredictedRt - Math.Max((Width / 2.0), _minimumWindowWidth/2);
        public double RtEndHypothesis => PredictedRt + Math.Max((Width / 2.0), _minimumWindowWidth/2);

        public RtInfo(double predictedRt, double width)
        {
            PredictedRt = predictedRt;
            Width = width;
        }
    }
}
