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
        /// <summary>
        /// Creates a bar plot. The data X value is the height of the bar, and the data label is the X-axis label under the bar.
        /// </summary>
        public BarPlot(PlotView oxyPlotView, IEnumerable<Datum> data) : base(oxyPlotView)
        {
            AddBarPlot(data);
        }
    }
}
