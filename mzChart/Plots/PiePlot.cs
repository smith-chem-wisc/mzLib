using BayesianEstimation;
using OxyPlot.Wpf;
using System;
using System.Collections.Generic;
using System.Text;
using MzLibUtil;

namespace mzPlot
{
    public class PiePlot : Plot
    {
        public PiePlot(PlotView oxyPlotView, List<Datum> data) : base(oxyPlotView)
        {
            AddPiePlot(data);
        }
    }
}
