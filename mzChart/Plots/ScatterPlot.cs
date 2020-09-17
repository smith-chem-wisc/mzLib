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
        /// <summary>
        /// Creates a scatter plot. The data X value is the X-coordinate, the data Y value is the Y-coordinate.
        /// </summary>
        public ScatterPlot(PlotView oxyPlotView, IEnumerable<Datum> data) : base(oxyPlotView)
        {
            AddScatterPlot(data);
        }
    }
}
