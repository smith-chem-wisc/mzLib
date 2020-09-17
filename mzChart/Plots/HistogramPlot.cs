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
        /// <summary>
        /// Creates a histogram plot. The data X value is used to bin the data.
        /// </summary>
        public HistogramPlot(PlotView oxyPlotView, IEnumerable<Datum> data) : base(oxyPlotView)
        {
            AddHistogram(data);
        }
    }
}
