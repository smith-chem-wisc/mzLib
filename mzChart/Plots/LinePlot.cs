using BayesianEstimation;
using OxyPlot.Wpf;
using System;
using System.Collections.Generic;
using System.Text;

namespace mzPlot
{
    public class LinePlot : Plot
    {
        public LinePlot(PlotView oxyPlotView, List<Datum> data) : base(oxyPlotView)
        {
            AddData(PlotType.Line, data);
        }
    }
}
