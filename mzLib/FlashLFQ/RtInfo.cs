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
        public double RtStartHypothesis => PredictedRt - (Width / 2.0);
        public double RtEndHypothesis => PredictedRt + (Width / 2.0);

        // These will be introduced in a later PR. For now, we're sticking with the classic version
        //private double _minimumWindowWidth = 0.5;
        //public double RtStartHypothesis => PredictedRt - Math.Max((Width / 2.0), _minimumWindowWidth/2); // the Math.Max components ensure that the width of an RT Window is at least _minimumWindowWidth wide
        //public double RtEndHypothesis => PredictedRt + Math.Max((Width / 2.0), _minimumWindowWidth/2);

        public RtInfo(double predictedRt, double width)
        {
            PredictedRt = predictedRt;
            Width = width;
        }
    }
}
