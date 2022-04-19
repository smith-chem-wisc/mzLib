using BayesianEstimation;
using System;
using System.Collections.Generic;
using System.Text;
using MzLibUtil;
using OxyPlot;
using System.Linq;
using OxyPlot.Series;
using OxyPlot.Axes;

namespace mzPlot
{
    public class HistogramPlot : Plot
    {
        /// <summary>
        /// Creates a histogram plot. The data X value is used to bin the data.
        /// </summary>
        public HistogramPlot(OxyPlot.Wpf.PlotView oxyPlotView, IEnumerable<Datum> data, int numBins, OxyColor? borderColor = null, OxyColor? fillColor = null,
            double borderThickness = 1, string xAxisLabel = null, string chartTitle = null, string chartSubtitle = null,
            bool refreshAfterAddingData = true) : base(oxyPlotView)
        {
            AddHistogram(data, numBins, borderColor, fillColor, borderThickness, xAxisLabel, chartTitle, chartSubtitle, refreshAfterAddingData);
        }
    }
}
