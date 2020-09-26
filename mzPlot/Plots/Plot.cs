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

namespace mzPlot
{
    public abstract class Plot : INotifyPropertyChanged
    {
        private PlotModel privateModel;

        /// <summary>
        /// The OxyPlot model for the chart.
        /// </summary>
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
        public void AddScatterPlot(IEnumerable<Datum> data, OxyColor? markerColor = null, double markerSize = 3,
            string xAxisLabel = null, string yAxisLabel = null, string chartTitle = null, string chartSubtitle = null,
            bool addToLegend = true, string seriesTitle = null, bool refreshAfterAddingData = true)
        {
            SetCommonChartProperties(chartTitle, chartSubtitle);

            if (!data.Any())
            {
                return;
            }

            var scatter = new ScatterSeries
            {
                Title = seriesTitle,
                RenderInLegend = addToLegend,
                MarkerSize = markerSize,
            };

            if (markerColor != null)
            {
                scatter.MarkerFill = markerColor.Value;
            }

            foreach (Datum datum in data)
            {
                scatter.Points.Add(new ScatterPoint(datum.X, datum.Y.Value));
            }

            if (!Model.Axes.Any())
            {
                AddLinearAxis(xAxisLabel, AxisPosition.Bottom);
                AddLinearAxis(yAxisLabel, AxisPosition.Left);
            }

            Model.Series.Add(scatter);

            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

        /// <summary>
        /// Adds a line plot. The data X value is the X-coordinate, the data Y value is the Y-coordinate.
        /// </summary>
        public void AddLinePlot(IEnumerable<Datum> data, OxyColor? lineColor = null, double lineThickness = 2,
            string xAxisLabel = null, string yAxisLabel = null, string chartTitle = null, string chartSubtitle = null,
            bool addToLegend = true, string seriesTitle = null, bool refreshAfterAddingData = true)
        {
            SetCommonChartProperties(chartTitle, chartSubtitle);

            if (!data.Any())
            {
                return;
            }

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

            if (!Model.Axes.Any())
            {
                AddLinearAxis(xAxisLabel, AxisPosition.Bottom);
                AddLinearAxis(yAxisLabel, AxisPosition.Left);
            }

            Model.Series.Add(line);

            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

        /// <summary>
        /// Adds a histogram. The data X value is used to bin the data.
        /// </summary>
        public void AddHistogram(IEnumerable<Datum> data, int numBins, OxyColor? borderColor = null, OxyColor? fillColor = null,
            double borderThickness = 1, string xAxisLabel = null, string chartTitle = null, string chartSubtitle = null,
            bool refreshAfterAddingData = true)
        {
            SetCommonChartProperties(chartTitle, chartSubtitle);

            if (!data.Any())
            {
                return;
            }

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

            var xAxis = new LinearAxis()
            {
                MajorStep = binWidth,
                Position = AxisPosition.Bottom,
                StringFormat = "F2",
                MinorTickSize = 0,
                Title = xAxisLabel
            };

            var yAxis = new LinearAxis() { Position = AxisPosition.Left, MinorTickSize = 0, Title = "Count" };

            if (!Model.Axes.Any())
            {
                Model.Axes.Add(xAxis);
                Model.Axes.Add(yAxis);
            }

            Model.Series.Add(histogramSeries);

            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

        /// <summary>
        /// Adds a spectrum plot. The data X value is the X-coordinate of the spectral line, the data Y value is the height of 
        /// the spectral line at X.
        /// </summary>
        public void AddSpectrumPlot(IEnumerable<Datum> data, OxyColor? lineColor = null, double lineThickness = 0.5,
            string xAxisLabel = null, string yAxisLabel = null, string chartTitle = null, string chartSubtitle = null,
            bool addToLegend = true, string seriesTitle = null, bool refreshAfterAddingData = true)
        {
            SetCommonChartProperties(chartTitle, chartSubtitle);

            if (!data.Any())
            {
                return;
            }

            if (!lineColor.HasValue)
            {
                lineColor = OxyColors.DimGray;
            }

            foreach (Datum datum in data)
            {
                AddLinePlot(
                    data: new List<Datum> { new Datum(datum.X, 0), datum },
                    lineColor: lineColor,
                    lineThickness: lineThickness,
                    addToLegend: addToLegend,
                    seriesTitle: seriesTitle,
                    refreshAfterAddingData: false);
            }

            if (!Model.Axes.Any())
            {
                AddLinearAxis(xAxisLabel, AxisPosition.Bottom);
                AddLinearAxis(yAxisLabel, AxisPosition.Left);
            }

            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

        /// <summary>
        /// Adds a bar plot. The data X value is the height of the bar, and the data label is the X-axis label under the bar.
        /// </summary>
        public void AddBarPlot(IEnumerable<Datum> data, OxyColor? borderColor = null, OxyColor? fillColor = null,
            double borderThickness = 1, string xAxisLabel = null, string yAxisLabel = null, string chartTitle = null,
            string chartSubtitle = null, bool addToLegend = true, string seriesTitle = null, bool refreshAfterAddingData = true)
        {
            SetCommonChartProperties(chartTitle, chartSubtitle);

            if (!data.Any())
            {
                return;
            }

            var barSeries = new ColumnSeries
            {
                Title = seriesTitle,
                RenderInLegend = addToLegend
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

            var xAxis = new CategoryAxis { Title = xAxisLabel, Position = AxisPosition.Bottom };
            xAxis.Labels.AddRange(data.Select(p => p.Label));

            foreach (Datum datum in data)
            {
                barSeries.Items.Add(new ColumnItem(datum.X));
            }

            if (!Model.Axes.Any())
            {
                Model.Axes.Add(xAxis);
                AddLinearAxis(yAxisLabel, AxisPosition.Left);
            }

            Model.Series.Add(barSeries);

            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

        /// <summary>
        /// Adds text to the plot. The x and y coordinates refer to the location on the chart itself, not on the data x and y axes.
        /// </summary>
        public PlotTextAnnotation AddTextAnnotationToPlotArea(string text, double x, double y, OxyColor? textColor = null,
            string font = "Times New Roman", double fontSize = 20)
        {
            var annotation = new PlotTextAnnotation()
            {
                Text = text,
                FontSize = fontSize,
                Font = font,
                TextColor = OxyColors.Black,
                X = x,
                Y = y
            };

            if (textColor.HasValue)
            {
                annotation.TextColor = textColor.Value;
            }

            Model.Annotations.Add(annotation);

            RefreshChart();

            return annotation;
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

        public event PropertyChangedEventHandler PropertyChanged;

        protected Plot(OxyPlot.Wpf.PlotView plotView)
        {
            ClearChart();
            plotView.DataContext = this;
            plotView.Model = Model;
        }

        protected void AddLinearAxis(string axisLabel, AxisPosition position)
        {
            var axis = new LinearAxis() { Title = axisLabel, Position = position, MinorTickSize = 0 };

            Model.Axes.Add(axis);
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
        }
    }
}
