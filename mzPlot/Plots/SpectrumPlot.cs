using BayesianEstimation;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Text;
using MzLibUtil;
using System.Linq;
using OxyPlot;
using OxyPlot.Axes;

namespace mzPlot
{
    public class SpectrumPlot : Plot
    {
        /// <summary>
        /// Adds a spectrum plot. The data X value is the X-coordinate of the spectral line, the data Y value is the height of 
        /// the spectral line at X.
        /// </summary>
        public SpectrumPlot(OxyPlot.Wpf.PlotView oxyPlotView, IEnumerable<Datum> data, OxyColor? lineColor = null, double lineThickness = 0.5,
            string xAxisLabel = null, string yAxisLabel = null, string chartTitle = null, string chartSubtitle = null,
            bool addToLegend = true, string seriesTitle = null, bool refreshAfterAddingData = true) : base(oxyPlotView)
        {
            AddSpectrumPlot(data, lineColor, lineThickness, xAxisLabel, yAxisLabel, chartTitle, chartSubtitle,
                addToLegend, seriesTitle, refreshAfterAddingData);
        }
    }
}
