﻿using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using mzPlot;
using BayesianEstimation;
using System.Threading;
using OxyPlot.Series;
using MzLibUtil;
using System.IO;
using OxyPlot;

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
            Plot plot = new ScatterPlot(examplePlotView, new List<Datum> { datapoint1, datapoint2 }, markerColor: OxyColors.Blue,
                xAxisLabel: "xAxis", yAxisLabel: "yAxis", chartTitle: "title", chartSubtitle: "subtitle");

            // check to make sure the data was plotted
            Assert.That(plot.Model.Series.Count == 1);
            var series = plot.Model.Series[0];

            var points = ((ScatterSeries)series).ActualPoints;
            Assert.That(points.Count == 2);
            Assert.That(points[0].X == 0);
            Assert.That(points[0].Y == 1);
            Assert.That(points[1].X == 2);
            Assert.That(points[1].Y == 3);

            Assert.That(((ScatterSeries)series).ActualMarkerFillColor == OxyColors.Blue);

            Assert.That(plot.Model.Title == "title");
            Assert.That(plot.Model.Subtitle == "subtitle");
            Assert.That(plot.Model.Axes[0].Title == "xAxis");
            Assert.That(plot.Model.Axes[1].Title == "yAxis");
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
            Plot plot = new LinePlot(examplePlotView, new List<Datum> { datapoint1, datapoint2 }, lineColor: OxyColors.Blue);

            // check to make sure the data was plotted
            Assert.That(plot.Model.Series.Count == 1);
            var series = plot.Model.Series[0];

            var points = ((LineSeries)series).Points;
            Assert.That(points.Count == 2);
            Assert.That(points[0].X == 0);
            Assert.That(points[0].Y == 1);
            Assert.That(points[1].X == 2);
            Assert.That(points[1].Y == 3);

            Assert.That(((LineSeries)series).ActualColor == OxyColors.Blue);
        }

        [Test]
        public static void TestBarPlot()
        {
            // the PlotView is a WPF control that's created in the .xaml code
            OxyPlot.Wpf.PlotView examplePlotView = new OxyPlot.Wpf.PlotView();

            // just some example data to plot
            var datapoint1 = new Datum(0, label: "category1");
            var datapoint2 = new Datum(2, label: "category2");

            // create the plot
            Plot plot = new BarPlot(examplePlotView, new List<Datum> { datapoint1, datapoint2 }, fillColor: OxyColors.Blue);

            // check to make sure the data was plotted
            Assert.That(plot.Model.Series.Count == 1);
            var series = plot.Model.Series[0];

            var points = ((ColumnSeries)series).Items;
            Assert.That(points.Count == 2);
            Assert.That(points[0].Value == 0);
            Assert.That(points[1].Value == 2);

            Assert.That(((ColumnSeries)series).ActualFillColor == OxyColors.Blue);
        }

        [Test]
        public static void TestHistogramPlot()
        {
            // the PlotView is a WPF control that's created in the .xaml code
            OxyPlot.Wpf.PlotView examplePlotView = new OxyPlot.Wpf.PlotView();

            int numBins = 10;

            // just some example data to plot
            var data = new List<Datum>()
            {
                new Datum(0),
                new Datum(0.75),
                new Datum(0),
                new Datum(0.5),
                new Datum(3),
                new Datum(3.5),
                new Datum(3),
                new Datum(3),
                new Datum(10),
                new Datum(10),
            };

            // create the plot
            Plot plot = new HistogramPlot(examplePlotView, data, numBins, fillColor: OxyColors.Blue);

            // check to make sure the data was plotted
            Assert.That(plot.Model.Series.Count == 1);
            var series = plot.Model.Series[0];

            var points = ((RectangleBarSeries)series).Items;
            Assert.That(points.Count == numBins);
            Assert.That(points[0].X0 == 0); // bin starts at 0
            Assert.That(points[0].X1 == 1); // bin ends at 1
            Assert.That(points[0].Y1 == 4); // 4 data points between 0 and 1

            Assert.That(((RectangleBarSeries)series).ActualFillColor == OxyColors.Blue);
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
            Plot plot = new SpectrumPlot(examplePlotView, new List<Datum> { datapoint1, datapoint2 }, lineColor: OxyColors.Blue);

            // check to make sure the data was plotted
            Assert.That(plot.Model.Series.Count == 2);
            var series = plot.Model.Series[0];

            var points = ((LineSeries)series).Points;
            Assert.That(points.Count == 2);
            Assert.That(points[0].X == 0);
            Assert.That(points[0].Y == 0);
            Assert.That(points[1].X == 0);
            Assert.That(points[1].Y == 1);

            series = plot.Model.Series[1];
            points = ((LineSeries)series).Points;
            Assert.That(points.Count == 2);
            Assert.That(points[0].X == 2);
            Assert.That(points[0].Y == 0);
            Assert.That(points[1].X == 2);
            Assert.That(points[1].Y == 3);

            Assert.That(((LineSeries)series).ActualColor == OxyColors.Blue);
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

        [Test]
        public static void TestPdfExport()
        {
            // the PlotView is a WPF control that's created in the .xaml code
            OxyPlot.Wpf.PlotView examplePlotView = new OxyPlot.Wpf.PlotView();

            // create the plot
            Plot plot = new ScatterPlot(examplePlotView, new List<Datum> { new Datum(0, 1), new Datum(2, 3) });

            string exportPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "testPdfExport.pdf");
            plot.ExportToPdf(exportPath);
            Assert.That(File.Exists(exportPath));

            File.Delete(exportPath);
        }

        [Test]
        public static void TestPngExport()
        {
            // the PlotView is a WPF control that's created in the .xaml code
            OxyPlot.Wpf.PlotView examplePlotView = new OxyPlot.Wpf.PlotView();

            // create the plot
            Plot plot = new ScatterPlot(examplePlotView, new List<Datum> { new Datum(0, 1), new Datum(2, 3) });

            string exportPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "testPngExport.png");
            plot.ExportToPng(exportPath);
            Assert.That(File.Exists(exportPath));

            File.Delete(exportPath);
        }

        [Test]
        public static void TestDataToString()
        {
            Datum datum = new Datum(0, 1, 2, "test", 0.5);
            Assert.That(datum.ToString() == "0, 1, 2; test; 0.5");
        }

        [Test]
        public static void TestTextAnnotation()
        {
            // the PlotView is a WPF control that's created in the .xaml code
            OxyPlot.Wpf.PlotView examplePlotView = new OxyPlot.Wpf.PlotView();

            // create the plot
            Plot plot = new ScatterPlot(examplePlotView, new List<Datum> { new Datum(0, 1), new Datum(2, 3) });

            plot.AddTextAnnotationToPlotArea("PEPTIDESEQUENCE", 100, -10, OxyColors.Blue);

            Assert.That(plot.Model.Annotations.Count == 1);
            Assert.That(((PlotTextAnnotation)plot.Model.Annotations[0]).Text == "PEPTIDESEQUENCE");
            Assert.That(((PlotTextAnnotation)plot.Model.Annotations[0]).TextColor == OxyColors.Blue);

            string exportPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "testTextAnnotatedPdfExport.pdf");
            plot.ExportToPdf(exportPath);
            Assert.That(File.Exists(exportPath));

            File.Delete(exportPath);
        }

        [Test]
        public static void TestEmptyDataCharts()
        {
            // the PlotView is a WPF control that's created in the .xaml code
            OxyPlot.Wpf.PlotView examplePlotView = new OxyPlot.Wpf.PlotView();

            List<Datum> data = new List<Datum>();

            Plot plot = new ScatterPlot(examplePlotView, data);
            Assert.That(plot.Model.Series.Count == 0);

            plot = new LinePlot(examplePlotView, data);
            Assert.That(plot.Model.Series.Count == 0);

            plot = new BarPlot(examplePlotView, data);
            Assert.That(plot.Model.Series.Count == 0);

            plot = new HistogramPlot(examplePlotView, data, 10);
            Assert.That(plot.Model.Series.Count == 0);

            plot = new SpectrumPlot(examplePlotView, data);
            Assert.That(plot.Model.Series.Count == 0);
        }
    }
}
