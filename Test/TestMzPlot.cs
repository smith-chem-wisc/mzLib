using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using mzPlot;
using BayesianEstimation;
using System.Threading;
using OxyPlot.Series;

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
            ScatterPlot plot = new ScatterPlot(examplePlotView, new List<Datum> { datapoint1, datapoint2 });

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
            LinePlot plot = new LinePlot(examplePlotView, new List<Datum> { datapoint1, datapoint2 });

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
    }
}
