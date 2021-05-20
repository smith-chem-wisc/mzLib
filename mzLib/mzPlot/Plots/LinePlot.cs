using BayesianEstimation;
using System;
using System.Collections.Generic;
using System.Text;
using MzLibUtil;
using OxyPlot;
using OxyPlot.Series;
using System.Linq;
using OxyPlot.Axes;

namespace mzPlot
{
    public class LinePlot : Plot
    {
        /// <summary>
        /// Creates a line plot. The data X value is the X-coordinate, the data Y value is the Y-coordinate.
        /// </summary>
        public LinePlot(OxyPlot.Wpf.PlotView oxyPlotView, IEnumerable<Datum> data, OxyColor? lineColor = null, double lineThickness = 2,
            string xAxisLabel = null, string yAxisLabel = null, string chartTitle = null, string chartSubtitle = null,
            bool addToLegend = true, string seriesTitle = null, bool refreshAfterAddingData = true) : base(oxyPlotView)
        {
            AddLinePlot(data, lineColor, lineThickness, xAxisLabel, yAxisLabel, chartTitle, chartSubtitle, 
                addToLegend, seriesTitle, refreshAfterAddingData);
        }
    }
}
