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
        /// <summary>
        /// Creates a line plot. The data X value is the X-coordinate, the data Y value is the Y-coordinate.
        /// </summary>
        public LinePlot(PlotView oxyPlotView, IEnumerable<Datum> data) : base(oxyPlotView)
        {
            AddLinePlot(data);
        }
    }
}
