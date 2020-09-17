using BayesianEstimation;
using OxyPlot.Wpf;
using System;
using System.Collections.Generic;
using System.Text;
using MzLibUtil;

namespace mzPlot
{
    public class LinePlot : Plot
    {
        public LinePlot(PlotView oxyPlotView, List<Datum> data) : base(oxyPlotView)
        {
            AddLinePlot(data);
        }
    }
}
