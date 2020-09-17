using BayesianEstimation;
using OxyPlot.Wpf;
using System;
using System.Collections.Generic;
using System.Text;
using MzLibUtil;

namespace mzPlot
{
    public class BarPlot : Plot
    {
        public BarPlot(PlotView oxyPlotView, List<Datum> data) : base(oxyPlotView)
        {
            AddBarPlot(data);
        }
    }
}
