using BayesianEstimation;
using OxyPlot;
using OxyPlot.Series;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Text;

namespace mzPlot
{
    public enum PlotType { Scatter, Line, Bar, Histogram, Pie };

    public abstract class Plot : INotifyPropertyChanged
    {
        private PlotModel privateModel;

        /// <summary>
        /// Clears the current plot.
        /// </summary>
        public void ClearChart()
        {
            Model = new PlotModel() { Title = "Data Annotation", Subtitle = "using OxyPlot" };
            RefreshChart();
        }

        /// <summary>
        /// Refreshes the chart. Call this method after adding a plot to the chart.
        /// </summary>
        public void RefreshChart()
        {
            Model.InvalidatePlot(true);
            NotifyPropertyChanged("Model");
        }

        public void AddData(PlotType plotType, List<Datum> data, OxyColor? borderColor = null, OxyColor? fillColor = null, double strokeThickness = 2,
            bool addToLegend = true, string seriesTitle = "", bool refreshAfterAddingData = true)
        {
            switch (plotType)
            {
                // scatter plot
                case PlotType.Scatter:
                    var scatter = new ScatterSeries
                    {
                        Title = seriesTitle,
                        RenderInLegend = addToLegend,
                        //MarkerStrokeThickness = strokeThickness,
                    };

                    if (borderColor != null)
                    {
                        scatter.MarkerFill = borderColor.Value;
                    }

                    foreach (Datum datum in data)
                    {
                        scatter.Points.Add(new ScatterPoint(datum.Dimension1, datum.Dimension2.Value));
                    }

                    Model.Series.Add(scatter);

                    break;

                // line chart
                case PlotType.Line:
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
                        line.Points.Add(new DataPoint(datum.Dimension1, datum.Dimension2.Value));
                    }

                    Model.Series.Add(line);

                    break;

                // bar chart
                case PlotType.Bar:
                    break;
            }

            RefreshChart();
        }

        public event PropertyChangedEventHandler PropertyChanged;

        protected Plot(OxyPlot.Wpf.PlotView plotView)
        {
            ClearChart();
            plotView.DataContext = this;
        }

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
