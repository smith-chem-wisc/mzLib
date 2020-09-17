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
        public void AddScatterPlot(IEnumerable<Datum> data, OxyColor? borderColor = null, OxyColor? fillColor = null, double strokeThickness = 2,
            bool addToLegend = true, string seriesTitle = "", bool refreshAfterAddingData = true)
        {
            if (!data.Any())
            {
                return;
            }

            var scatter = new ScatterSeries
            {
                Title = seriesTitle,
                RenderInLegend = addToLegend,
                MarkerStrokeThickness = strokeThickness,
            };

            if (borderColor != null)
            {
                scatter.MarkerFill = borderColor.Value;
            }

            foreach (Datum datum in data)
            {
                scatter.Points.Add(new ScatterPoint(datum.X, datum.Y.Value));
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
        public void AddLinePlot(IEnumerable<Datum> data, OxyColor? borderColor = null, OxyColor? fillColor = null, double strokeThickness = 2,
            bool addToLegend = true, string seriesTitle = "", bool refreshAfterAddingData = true)
        {
            if (!data.Any())
            {
                return;
            }

            var line = new LineSeries
            {
                Title = seriesTitle,
                RenderInLegend = addToLegend,
                StrokeThickness = strokeThickness
            };

            if (borderColor != null)
            {
                line.Color = borderColor.Value;
            }

            foreach (Datum datum in data)
            {
                line.Points.Add(new DataPoint(datum.X, datum.Y.Value));
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
        public void AddHistogram(IEnumerable<Datum> data, int numBins, OxyColor? borderColor = null, OxyColor? fillColor = null, double strokeThickness = 2,
            bool addToLegend = true, string seriesTitle = "", bool refreshAfterAddingData = true)
        {
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

            Model.Series.Add(histogramSeries);

            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

        /// <summary>
        /// Adds a spectrum plot. The data X value is the X-coordinate of the spectral line, the data Y value is the height of the spectral line at X.
        /// </summary>
        public void AddSpectrumPlot(IEnumerable<Datum> data, OxyColor? borderColor = null, OxyColor? fillColor = null, double strokeThickness = 2,
            bool addToLegend = true, string seriesTitle = "", bool refreshAfterAddingData = true)
        {
            if (!data.Any())
            {
                return;
            }

            foreach (Datum datum in data)
            {
                AddLinePlot(new List<Datum> { new Datum(datum.X, 0), datum }, refreshAfterAddingData: false);
            }

            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

        /// <summary>
        /// Adds a pie plot. Not implemented yet.
        /// </summary>
        public void AddPiePlot(IEnumerable<Datum> data, OxyColor? borderColor = null, OxyColor? fillColor = null, double strokeThickness = 2,
            bool addToLegend = true, string seriesTitle = "", bool refreshAfterAddingData = true)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Adds a bar plot. The data X value is the height of the bar, and the data label is the X-axis label under the bar.
        /// </summary>
        public void AddBarPlot(IEnumerable<Datum> data, OxyColor? borderColor = null, OxyColor? fillColor = null, double strokeThickness = 2,
            bool addToLegend = true, string seriesTitle = "", bool refreshAfterAddingData = true)
        {
            if (!data.Any())
            {
                return;
            }

            var barSeries = new ColumnSeries
            {
                Title = seriesTitle,
            };

            if (borderColor != null)
            {
                barSeries.StrokeColor = borderColor.Value;
            }

            if (fillColor != null)
            {
                barSeries.FillColor = fillColor.Value;
            }

            var xAxis = new CategoryAxis { Position = AxisPosition.Bottom };
            xAxis.Labels.AddRange(data.Select(p => p.Label));
            Model.Axes.Add(xAxis);

            foreach (Datum datum in data)
            {
                barSeries.Items.Add(new ColumnItem(datum.X));
            }

            Model.Series.Add(barSeries);

            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

        /// <summary>
        /// Adds a heatmap. Not implemented yet.
        /// </summary>
        public void AddHeatMap(IEnumerable<Datum> data, OxyColor? borderColor = null, OxyColor? fillColor = null, double strokeThickness = 2,
            bool addToLegend = true, string seriesTitle = "", bool refreshAfterAddingData = true)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Adds a boxplot. Not implemented yet.
        /// </summary>
        public void AddBoxPlot(IEnumerable<Datum> data, OxyColor? borderColor = null, OxyColor? fillColor = null, double strokeThickness = 2,
            bool addToLegend = true, string seriesTitle = "", bool refreshAfterAddingData = true)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Adds an extracted ion chromatogram.
        /// </summary>
        public void AddXicPlot(ExtractedIonChromatogram xic, OxyColor? borderColor = null, OxyColor? fillColor = null, double strokeThickness = 2,
            bool addToLegend = true, string seriesTitle = "", bool refreshAfterAddingData = true)
        {
            if (!xic.Data.Any())
            {
                return;
            }

            AddLinePlot(xic.Data, refreshAfterAddingData: refreshAfterAddingData);
        }

        /// <summary>
        /// Adds text to the plot. The x and y coordinates refer to the location on the chart itself, not on the data x and y axes.
        /// </summary>
        public void AddTextAnnotationToPlot(string text, double x, double y, OxyColor? textColor = null, double fontSize = 20)
        {
            var annotation = new PlotTextAnnotation()
            {
                Text = text,
                FontSize = fontSize,
                Font = "Times New Roman",
                TextColor = OxyColors.Black,
                X = x,
                Y = y
            };

            this.Model.Annotations.Add(annotation);

            RefreshChart();
        }

        /// <summary>
        /// Adds text to the plot. The x and y coordinates refer to the location on the data x and y axes, not the chart area on the screen.
        /// Not implemented yet.
        /// </summary>
        public void AddTextAnnotationToData(string text, double x, double y, OxyColor? textColor = null, double fontSize = 20)
        {
            throw new NotImplementedException();
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
            using (var s = File.OpenWrite(path))
            {
                PdfExporter.Export(Model, s, width, height);
            }
        }

        /// <summary>
        /// Exports the plot to a .png file.
        /// </summary>
        public void ExportToPng(string path, int width = 800, int height = 600)
        {
            using (var s = File.OpenWrite(path))
            {
                var pngExporter = new OxyPlot.Wpf.PngExporter { Width = width, Height = height, Background = OxyColors.White };
                pngExporter.Export(Model, s);
            }
        }

        /// <summary>
        /// Exports the plot to an .svg file. Not implemented yet.
        /// </summary>
        public void ExportToSvg(string path)
        {
            throw new NotImplementedException();
        }

        public event PropertyChangedEventHandler PropertyChanged;

        protected Plot(OxyPlot.Wpf.PlotView plotView)
        {
            ClearChart();
            plotView.DataContext = this;
        }

        protected void NotifyPropertyChanged(string propertyName)
        {
            PropertyChangedEventHandler handler = PropertyChanged;
            if (handler != null)
            {
                handler(this, new PropertyChangedEventArgs(propertyName));
            }
        }
    }
}
