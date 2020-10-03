using BayesianEstimation;
using System;
using System.Collections.Generic;
using System.Text;
using MzLibUtil;
using OxyPlot;
using System.Linq;
using OxyPlot.Series;
using OxyPlot.Axes;
using System.Drawing;

namespace mzPlot
{
    public class ScatterPlot : Plot
    {
        public string Name;
        public MarkerType MarkerType;
        public Color? MarkerColor { get; private set; }
        public double MarkerSize { get; private set; }

        /// <summary>
        /// Creates a scatter plot. The data X value is the X-coordinate, the data Y value is the Y-coordinate.
        /// </summary>
        public ScatterPlot(OxyPlot.Wpf.PlotView oxyPlotView, IEnumerable<Datum> data, Color? markerColor = null, double markerSize = 2, MarkerType markerType = MarkerType.Circle,
            string xAxisLabel = null, string yAxisLabel = null, string chartTitle = null, string chartSubtitle = null,
            bool addToLegend = true, string seriesTitle = null, bool refreshAfterAddingData = true, string name = null) : base(oxyPlotView)
        {
            Name = name;
            MarkerColor = markerColor;
            MarkerSize = markerSize;

            AddScatterSeries(
                data: data,
                markerColor: MarkerColor,
                markerSize: MarkerSize,
                markerType: MarkerType,
                xAxisLabel: xAxisLabel,
                yAxisLabel: yAxisLabel,
                chartTitle: chartTitle,
                chartSubtitle: chartSubtitle,
                addToLegend: addToLegend,
                seriesTitle: seriesTitle,
                refreshAfterAddingData: refreshAfterAddingData);
        }

        /// <summary>
        /// Copies the settings (color, axes, etc.) from another plot, with new data.
        /// </summary>
        public ScatterPlot(OxyPlot.Wpf.PlotView oxyPlotView, IEnumerable<Datum> data, ScatterPlot existingScatterPlot) : base(oxyPlotView)
        {
            this.Name = existingScatterPlot.Name;
            this.MarkerColor = existingScatterPlot.MarkerColor;
            this.MarkerSize = existingScatterPlot.MarkerSize;
        }

        public override string ToString()
        {
            return Name;
        }
    }
}
