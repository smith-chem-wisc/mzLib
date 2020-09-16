using BayesianEstimation;
using OxyPlot.Wpf;
using System;
using System.Collections.Generic;
using System.Text;

namespace mzPlot
{
    public class ScatterPlot : Plot
    {
        public ScatterPlot(PlotView oxyPlotView, List<Datum> data) : base(oxyPlotView)
        {
            AddData(PlotType.Scatter, data);
        }
    }
}
