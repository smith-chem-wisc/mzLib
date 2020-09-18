using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using mzPlot;
using BayesianEstimation;
using System.Threading;
using OxyPlot.Series;
using MzLibUtil;

namespace Test
{
    [TestFixture, Apartment(ApartmentState.STA)]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestMzPlot
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public static void TestScatterPlot()
        {
            // the PlotView is a WPF control that's created in the .xaml code
            OxyPlot.Wpf.PlotView examplePlotView = new OxyPlot.Wpf.PlotView();

            // just some example data to plot
            var datapoint1 = new Datum(0, 1);
            var datapoint2 = new Datum(2, 3);

            // create the plot
            Plot plot = new ScatterPlot(examplePlotView, new List<Datum> { datapoint1, datapoint2 });

            // check to make sure the data was plotted
            Assert.That(plot.Model.Series.Count == 1);
            var series = plot.Model.Series[0];

            var points = ((ScatterSeries)series).ActualPoints;
            Assert.That(points.Count == 2);
            Assert.That(points[0].X == 0);
            Assert.That(points[0].Y == 1);
            Assert.That(points[1].X == 2);
            Assert.That(points[1].Y == 3);
        }

        [Test]
        public static void TestLinePlot()
        {
            // the PlotView is a WPF control that's created in the .xaml code
            OxyPlot.Wpf.PlotView examplePlotView = new OxyPlot.Wpf.PlotView();

            // just some example data to plot
            var datapoint1 = new Datum(0, 1);
            var datapoint2 = new Datum(2, 3);

            // create the plot
            Plot plot = new LinePlot(examplePlotView, new List<Datum> { datapoint1, datapoint2 });

            // check to make sure the data was plotted
            Assert.That(plot.Model.Series.Count == 1);
            var series = plot.Model.Series[0];

            var points = ((LineSeries)series).Points;
            Assert.That(points.Count == 2);
            Assert.That(points[0].X == 0);
            Assert.That(points[0].Y == 1);
            Assert.That(points[1].X == 2);
            Assert.That(points[1].Y == 3);

            plot.ExportToPdf(@"C:\Data\LVS_TD_Yeast\MSConvertMzml\TestPdfExport.pdf");
        }

        [Test]
        public static void TestBarPlot()
        {
            // the PlotView is a WPF control that's created in the .xaml code
            OxyPlot.Wpf.PlotView examplePlotView = new OxyPlot.Wpf.PlotView();

            // just some example data to plot
            var datapoint1 = new Datum(0, 1);
            var datapoint2 = new Datum(2, 3);

            // create the plot
            Plot plot = new BarPlot(examplePlotView, new List<Datum> { datapoint1, datapoint2 });

            // check to make sure the data was plotted
            Assert.That(plot.Model.Series.Count == 1);
            var series = plot.Model.Series[0];

            var points = ((LineSeries)series).Points;
            Assert.That(points.Count == 2);
            Assert.That(points[0].X == 0);
            Assert.That(points[0].Y == 1);
            Assert.That(points[1].X == 2);
            Assert.That(points[1].Y == 3);
        }

        [Test]
        public static void TestHistogramPlot()
        {
            // the PlotView is a WPF control that's created in the .xaml code
            OxyPlot.Wpf.PlotView examplePlotView = new OxyPlot.Wpf.PlotView();

            // just some example data to plot
            var datapoint1 = new Datum(0, 1);
            var datapoint2 = new Datum(2, 3);

            // create the plot
            Plot plot = new HistogramPlot(examplePlotView, new List<Datum> { datapoint1, datapoint2 }, 5);

            // check to make sure the data was plotted
            Assert.That(plot.Model.Series.Count == 1);
            var series = plot.Model.Series[0];

            var points = ((LineSeries)series).Points;
            Assert.That(points.Count == 2);
            Assert.That(points[0].X == 0);
            Assert.That(points[0].Y == 1);
            Assert.That(points[1].X == 2);
            Assert.That(points[1].Y == 3);
        }

        [Test]
        public static void TestSpectrumPlot()
        {
            // the PlotView is a WPF control that's created in the .xaml code
            OxyPlot.Wpf.PlotView examplePlotView = new OxyPlot.Wpf.PlotView();

            // just some example data to plot
            var datapoint1 = new Datum(0, 1);
            var datapoint2 = new Datum(2, 3);

            // create the plot
            Plot plot = new SpectrumPlot(examplePlotView, new List<Datum> { datapoint1, datapoint2 });

            // check to make sure the data was plotted
            Assert.That(plot.Model.Series.Count == 1);
            var series = plot.Model.Series[0];

            var points = ((LineSeries)series).Points;
            Assert.That(points.Count == 2);
            Assert.That(points[0].X == 0);
            Assert.That(points[0].Y == 1);
            Assert.That(points[1].X == 2);
            Assert.That(points[1].Y == 3);
        }

        [Test]
        public static void TestDoublePlot()
        {
            // the PlotView is a WPF control that's created in the .xaml code
            OxyPlot.Wpf.PlotView examplePlotView = new OxyPlot.Wpf.PlotView();

            // just some example data to plot
            var datapoint1 = new Datum(0, 1);
            var datapoint2 = new Datum(2, 3);

            // create the plot
            LinePlot plot = new LinePlot(examplePlotView, new List<Datum> { datapoint1, datapoint2 });

            var datapoint3 = new Datum(4, 5);
            var datapoint4 = new Datum(4, 6);
            plot.AddScatterPlot(new List<Datum> { datapoint3, datapoint4 });

            // check to make sure the data was plotted
            // the chart should have a line plot and a scatter plot on the same chart
            Assert.That(plot.Model.Series.Count == 2);

            var lineSeries = plot.Model.Series[0];
            var lineSeriesPoints = ((LineSeries)lineSeries).Points;
            Assert.That(lineSeriesPoints.Count == 2);
            Assert.That(lineSeriesPoints[0].X == 0);
            Assert.That(lineSeriesPoints[0].Y == 1);
            Assert.That(lineSeriesPoints[1].X == 2);
            Assert.That(lineSeriesPoints[1].Y == 3);

            var scatterSeries = plot.Model.Series[1];
            var scatterPoints = ((ScatterSeries)scatterSeries).Points;
            Assert.That(scatterPoints.Count == 2);
            Assert.That(scatterPoints[0].X == 4);
            Assert.That(scatterPoints[0].Y == 5);
            Assert.That(scatterPoints[1].X == 4);
            Assert.That(scatterPoints[1].Y == 6);
        }
    }
}
