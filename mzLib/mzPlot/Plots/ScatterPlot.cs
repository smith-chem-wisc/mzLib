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
    public class ScatterPlot : Plot
    {
        /// <summary>
        /// Creates a scatter plot. The data X value is the X-coordinate, the data Y value is the Y-coordinate.
        /// </summary>
        public ScatterPlot(OxyPlot.Wpf.PlotView oxyPlotView, IEnumerable<Datum> data, OxyColor? markerColor = null, double markerSize = 3,
            string xAxisLabel = null, string yAxisLabel = null, string chartTitle = null, string chartSubtitle = null,
            bool addToLegend = true, string seriesTitle = null, bool refreshAfterAddingData = true) : base(oxyPlotView)
        {
            AddScatterPlot(data, markerColor, markerSize, xAxisLabel, yAxisLabel, chartTitle, chartSubtitle, addToLegend,
                seriesTitle, refreshAfterAddingData);
        }
    }
}
