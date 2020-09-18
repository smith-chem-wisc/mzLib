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
    public class BarPlot : Plot
    {
        /// <summary>
        /// Creates a bar plot. The data X value is the height of the bar, and the data label is the X-axis label under the bar.
        /// </summary>
        public BarPlot(OxyPlot.Wpf.PlotView oxyPlotView, IEnumerable<Datum> data, OxyColor? borderColor = null, OxyColor? fillColor = null,
            double borderThickness = 1, string xAxisLabel = null, string yAxisLabel = null, string chartTitle = null,
            string chartSubtitle = null, bool addToLegend = true, string seriesTitle = null, bool refreshAfterAddingData = true) 
            : base(oxyPlotView)
        {
            AddBarPlot(data, borderColor, fillColor, borderThickness, xAxisLabel, yAxisLabel, chartTitle, chartSubtitle,
                addToLegend, seriesTitle, refreshAfterAddingData);
        }
    }
}
