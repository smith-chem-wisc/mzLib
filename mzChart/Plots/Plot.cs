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

        public void AddHistogram(IEnumerable<Datum> data, OxyColor? borderColor = null, OxyColor? fillColor = null, double strokeThickness = 2,
            bool addToLegend = true, string seriesTitle = "", bool refreshAfterAddingData = true)
        {
            if (!data.Any())
            {
                return;
            }

            int numBins = data.Count() / 10;
            numBins = Math.Max(numBins, 1);
            numBins = Math.Min(numBins, 100);

            double min = data.Min(p => p.X);
            double max = data.Max(p => p.X);

            double binWidth = (max - min) / numBins;

            double binFloor = min;

            var histXAxis = new CategoryAxis { Position = AxisPosition.Bottom };

            var histogramSeries = new ColumnSeries
            {
                Title = seriesTitle,
            };

            for (int i = 0; i < numBins; i++)
            {
                var binData = data.Where(p => p.X >= binFloor && p.X < binFloor + binWidth);

                histXAxis.Labels.Add(binFloor.ToString("F2"));

                histogramSeries.Items.Add(new ColumnItem(binData.Count()));

                binFloor += binWidth;
            }

            Model.Axes.Add(histXAxis);
            Model.Series.Add(histogramSeries);

            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

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

        public void AddPiePlot(IEnumerable<Datum> data, OxyColor? borderColor = null, OxyColor? fillColor = null, double strokeThickness = 2,
            bool addToLegend = true, string seriesTitle = "", bool refreshAfterAddingData = true)
        {
            throw new NotImplementedException();

            if (!data.Any())
            {
                return;
            }

            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

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

        public void AddHeatMap(IEnumerable<Datum> data, OxyColor? borderColor = null, OxyColor? fillColor = null, double strokeThickness = 2,
            bool addToLegend = true, string seriesTitle = "", bool refreshAfterAddingData = true)
        {
            throw new NotImplementedException();

            if (!data.Any())
            {
                return;
            }

            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

        public void AddBoxPlot(IEnumerable<Datum> data, OxyColor? borderColor = null, OxyColor? fillColor = null, double strokeThickness = 2,
            bool addToLegend = true, string seriesTitle = "", bool refreshAfterAddingData = true)
        {
            throw new NotImplementedException();

            if (!data.Any())
            {
                return;
            }

            if (refreshAfterAddingData)
            {
                RefreshChart();
            }
        }

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
        /// Exports the plot to a .png file. Not implemented yet.
        /// </summary>
        public void ExportToPng(string path)
        {
            throw new NotImplementedException();
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
