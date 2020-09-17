using BayesianEstimation;
using OxyPlot.Wpf;
using System;
using System.Collections.Generic;
using System.Text;
using MzLibUtil;

namespace mzPlot
{
    public class ScatterPlot : Plot
    {
        public ScatterPlot(PlotView oxyPlotView, List<Datum> data) : base(oxyPlotView)
        {
            AddScatterPlot(data);
        }
    }
}
