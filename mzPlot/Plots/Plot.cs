using BayesianEstimation;
using OxyPlot;
using OxyPlot.Axes;
using OxyPlot.Series;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Text;
using MzLibUtil;
using MassSpectrometry;
using System.Drawing;
using Nett;

namespace mzPlot
{
    public abstract class Plot : INotifyPropertyChanged
    {
        private PlotModel privateModel;

        /// <summary>
        /// The OxyPlot model for the chart.
        /// </summary>
        [TomlIgnore]
        public PlotModel Model
        {
            get
            {
                return this.privateModel;
            }
            set
            {
                this.privateModel = value;
                RefreshChart();
            }
        }

        /// <summary>
        /// Adds a scatter plot. The data X value is the X-coordinate, the data Y value is the Y-coordinate.
        /// </summary>
        public void AddScatterSeries(IEnumerable<Datum> data, Color? markerColor = null, double markerSize = 2,
            MarkerType markerType = MarkerType.Circle, string xAxisLabel = null, string yAxisLabel = null,
            bool addToLegend = true, string seriesTitle = null, bool refreshAfterAddingData = true)
        {
            if (!data.Any())
            {
                return;
            }

            // add data series
            var scatter = new ScatterSeries
            {
                Title = seriesTitle,
                RenderInLegend = addToLegend,
                MarkerSize = markerSize,
                MarkerType = markerType
            };

            if (markerColor != null)
            {
                scatter.MarkerFill = OxyColor.FromArgb(markerColor.Value.A, markerColor.Value.R, markerColor.Value.G, markerColor.Value.B);
                scatter.MarkerStroke = scatter.MarkerFill;
            }

            foreach (Datum datum in data)
            {
                scatter.Points.Add(new ScatterPoint(datum.X, datum.Y.Value));
            }

            Model.Series.Add(scatter);

            // add axes
            if (!Model.Axes.Any())
            {
                var xAxis = GetLinearAxis(xAxisLabel, AxisPosition.Bottom);
                var yAxis = GetLinearAxis(yAxisLabel, AxisPosition.Left);
                AddAxes(xAxis, yAxis);
            }

            // refresh chart
            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

        /// <summary>
        /// Adds a line plot. The data X value is the X-coordinate, the data Y value is the Y-coordinate.
        /// </summary>
        public void AddLineSeries(IEnumerable<Datum> data, OxyColor? lineColor = null, double lineThickness = 2,
            string xAxisLabel = null, string yAxisLabel = null, bool addToLegend = true, string seriesTitle = null, 
            bool refreshAfterAddingData = true)
        {
            if (!data.Any())
            {
                return;
            }

            // add data series
            var line = new LineSeries
            {
                Title = seriesTitle,
                RenderInLegend = addToLegend,
                StrokeThickness = lineThickness
            };

            if (lineColor != null)
            {
                line.Color = lineColor.Value;
            }

            foreach (Datum datum in data)
            {
                line.Points.Add(new DataPoint(datum.X, datum.Y.Value));
            }

            Model.Series.Add(line);

            // add axes
            if (!Model.Axes.Any())
            {
                var xAxis = GetLinearAxis(xAxisLabel, AxisPosition.Bottom);
                var yAxis = GetLinearAxis(yAxisLabel, AxisPosition.Left);
                AddAxes(xAxis, yAxis);
            }

            // refresh chart
            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

        /// <summary>
        /// Adds a histogram. The data X value is used to bin the data.
        /// </summary>
        public void AddHistogramSeries(IEnumerable<Datum> data, int numBins, OxyColor? borderColor = null, OxyColor? fillColor = null,
            double borderThickness = 0, string xAxisLabel = null, bool refreshAfterAddingData = true)
        {
            if (!data.Any())
            {
                return;
            }

            // add data series
            double min = data.Min(p => p.X);
            double max = data.Max(p => p.X);

            double binWidth = (max - min) / numBins;
            var histogramSeries = new RectangleBarSeries();

            double binFloor = min;
            for (int i = 0; i < numBins; i++)
            {
                var binData = data.Where(p => p.X >= binFloor && p.X < binFloor + binWidth);

                var bin = new RectangleBarItem(binFloor, 0, binFloor + binWidth, binData.Count());
                histogramSeries.Items.Add(bin);

                binFloor += binWidth;
            }

            histogramSeries.StrokeThickness = borderThickness;

            if (fillColor.HasValue)
            {
                histogramSeries.FillColor = fillColor.Value;
            }
            if (borderColor.HasValue)
            {
                histogramSeries.StrokeColor = borderColor.Value;
            }

            Model.Series.Add(histogramSeries);

            // add axes
            if (!Model.Axes.Any())
            {
                var xAxis = new LinearAxis()
                {
                    MajorStep = (max - min) / 5,
                    Position = AxisPosition.Bottom,
                    StringFormat = "F2",
                    MinorTickSize = 0,
                    Title = xAxisLabel
                };

                var yAxis = GetLinearAxis("Count", AxisPosition.Left);
                yAxis.AbsoluteMinimum = 0;

                AddAxes(xAxis, yAxis);
            }

            // refresh chart
            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

        /// <summary>
        /// Adds a spectrum plot. The data X value is the X-coordinate of the spectral line, the data Y value is the height of 
        /// the spectral line at X.
        /// </summary>
        public void AddSpectrumSeries(IEnumerable<Datum> data, OxyColor? lineColor = null, double lineThickness = 0.5,
            string xAxisLabel = null, string yAxisLabel = null, bool addToLegend = true, string seriesTitle = null, 
            bool refreshAfterAddingData = true)
        {
            if (!data.Any())
            {
                return;
            }

            // add data series
            if (!lineColor.HasValue)
            {
                lineColor = OxyColors.DimGray;
            }

            foreach (Datum datum in data)
            {
                AddLineSeries(
                    data: new List<Datum> { new Datum(datum.X, 0), datum },
                    lineColor: lineColor,
                    lineThickness: lineThickness,
                    addToLegend: addToLegend,
                    seriesTitle: seriesTitle,
                    refreshAfterAddingData: false);
            }

            // add axes
            if (!Model.Axes.Any())
            {
                var xAxis = GetLinearAxis(xAxisLabel, AxisPosition.Bottom);
                var yAxis = GetLinearAxis(yAxisLabel, AxisPosition.Left);
                yAxis.AbsoluteMinimum = 0;
                AddAxes(xAxis, yAxis);
            }

            // refresh chart
            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

        /// <summary>
        /// Adds a bar plot. The data X value is the height of the bar, and the data label is the X-axis label under the bar.
        /// </summary>
        public void AddBarSeries(IEnumerable<Datum> data, OxyColor? borderColor = null, OxyColor? fillColor = null,
            double borderThickness = 1, string xAxisLabel = null, string yAxisLabel = null, string chartTitle = null,
            string chartSubtitle = null, bool addToLegend = true, string seriesTitle = null, bool refreshAfterAddingData = true)
        {
            // set chart title and subtitle
            SetCommonChartProperties(chartTitle, chartSubtitle);

            if (!data.Any())
            {
                return;
            }

            // add data series
            var barSeries = new ColumnSeries
            {
                Title = seriesTitle,
                RenderInLegend = addToLegend,
            };

            if (borderColor != null)
            {
                barSeries.StrokeColor = borderColor.Value;
            }

            barSeries.StrokeThickness = borderThickness;

            if (fillColor != null)
            {
                barSeries.FillColor = fillColor.Value;
            }

            foreach (Datum datum in data)
            {
                barSeries.Items.Add(new ColumnItem(datum.X));
            }

            Model.Series.Add(barSeries);

            // add axes
            if (!Model.Axes.Any())
            {
                var xAxis = GetCategoryAxis(xAxisLabel, AxisPosition.Bottom);
                xAxis.Labels.AddRange(data.Select(p => p.Label));

                var yAxis = GetLinearAxis(yAxisLabel, AxisPosition.Left);
                AddAxes(xAxis, yAxis);
            }

            // refresh chart
            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

        /// <summary>
        /// Adds text to the plot. The x and y coordinates refer to the location on the chart itself, not on the data x and y axes.
        /// </summary>
        public void AddAnnotation(mzPlotAnnotation annotation)
        {
            Model.Annotations.Add(annotation);
            RefreshChart();
        }

        /// <summary>
        /// Clears the current plot.
        /// </summary>
        public void ClearChart()
        {
            Model = new PlotModel() { Title = "Data Annotation", Subtitle = "using OxyPlot" };
            RefreshChart();
        }

        /// <summary>
        /// Refreshes the chart.
        /// </summary>
        public void RefreshChart()
        {
            Model.InvalidatePlot(true);
            NotifyPropertyChanged("Model");
        }

        /// <summary>
        /// Exports the plot to a .pdf file.
        /// </summary>
        public void ExportToPdf(string path, double width = 800, double height = 600)
        {
            using (var s = File.Create(path))
            {
                PdfExporter.Export(Model, s, width, height);
            }
        }

        /// <summary>
        /// Exports the plot to a .png file.
        /// </summary>
        public void ExportToPng(string path, int width = 800, int height = 600)
        {
            using (var s = File.Create(path))
            {
                var pngExporter = new OxyPlot.Wpf.PngExporter { Width = width, Height = height, Background = OxyColors.White };
                pngExporter.Export(Model, s);
            }
        }

        /// <summary>
        /// Exports the plot to an .svg file. The resulting .svg files seem to not render properly in Google Chrome, 
        /// but work in Mozilla Firefox, and Microsoft Internet Explorer/Edge.
        /// </summary>
        public void ExportToSvg(string path, int width = 800, int height = 600)
        {
            using (var s = File.Create(path))
            {
                var svgExporter = new SvgExporter { Width = width, Height = height, UseVerticalTextAlignmentWorkaround = true };
                svgExporter.Export(Model, s);
            }
        }

        public void SetDefaultColors(List<Color> colors)
        {
            Model.DefaultColors.Clear();

            foreach (var color in colors)
            {
                Model.DefaultColors.Add(OxyColor.FromArgb(color.A, color.R, color.G, color.B));
            }
        }

        public event PropertyChangedEventHandler PropertyChanged;

        protected Plot(OxyPlot.Wpf.PlotView plotView, string title, string subtitle)
        {
            ClearChart();

            if (plotView != null)
            {
                plotView.DataContext = this;
                plotView.Model = Model;
            }

            SetCommonChartProperties(title, subtitle);
        }

        protected void AddAxes(Axis xAxis, Axis yAxis)
        {
            Model.Axes.Add(xAxis);
            Model.Axes.Add(yAxis);
        }

        protected LinearAxis GetLinearAxis(string axisLabel, AxisPosition position)
        {
            var axis = new LinearAxis() { Title = axisLabel, Position = position, MinorTickSize = 0 };
            return axis;
        }

        protected LogarithmicAxis GetLogarithmicAxis(string axisLabel, AxisPosition position)
        {
            var axis = new LogarithmicAxis() { Title = axisLabel, Position = position, MinorTickSize = 0 };
            return axis;
        }

        protected CategoryAxis GetCategoryAxis(string axisLabel, AxisPosition position)
        {
            var axis = new CategoryAxis() { Title = axisLabel, Position = position, MinorTickSize = 0 };
            return axis;
        }

        protected void NotifyPropertyChanged(string propertyName)
        {
            PropertyChangedEventHandler handler = PropertyChanged;
            if (handler != null)
            {
                handler(this, new PropertyChangedEventArgs(propertyName));
            }
        }

        private void SetCommonChartProperties(string chartTitle, string chartSubtitle)
        {
            if (chartTitle != null)
            {
                Model.Title = chartTitle;
            }

            if (chartSubtitle != null)
            {
                Model.Subtitle = chartSubtitle;
            }

            Model.LegendBorderThickness = 1;
            Model.LegendBorder = OxyColors.Black;
            Model.LegendBackground = OxyColors.White;
        }
    }
}
