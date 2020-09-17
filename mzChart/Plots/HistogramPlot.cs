using BayesianEstimation;
using OxyPlot.Wpf;
using System;
using System.Collections.Generic;
using System.Text;
using MzLibUtil;

namespace mzPlot
{
    public class HistogramPlot : Plot
    {
        public HistogramPlot(PlotView oxyPlotView, List<Datum> data) : base(oxyPlotView)
        {
            AddHistogram(data);
        }
    }
}
