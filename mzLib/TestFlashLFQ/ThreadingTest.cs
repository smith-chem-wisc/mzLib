using FlashLFQ;
using FlashLFQ.Alex_project;
using MathNet.Numerics.Distributions;
using NUnit.Framework;
using Plotly.NET;
using Plotly.NET.CSharp;
using Chart = Plotly.NET.CSharp.Chart;
using GenericChartExtensions = Plotly.NET.CSharp.GenericChartExtensions;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading;
using System.Threading.Tasks;

namespace TestFlashLFQ
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal class ThreadingTest
    {
        [Test]
        public void TestThreading()
        {
            Random random = new Random();
            Dictionary<int, XIC> xicDictionary = new Dictionary<int, XIC>();
            XIC xic = null;
            for (int i = 0; i < 230; i++)
            {
                
                int randomIensitry = random.Next(1,3);
                int randomRetention = random.Next(-3, 3);

                List<double> timesPoints = Enumerable.Range(0,200).Select(t => 10 + (double)t/10.0).ToList();
                List<double> intensities = new List<double>();
                var peak_1 = new Normal(14.5, 0.6);    // first peak
                var peak_2 = new Normal(16.5, 0.6);    // second peak

                for (int k = 0; k < 200; k++) //creat the intensity point
                {
                    intensities.Add( (peak_1.Density(timesPoints[k]) + peak_2.Density(timesPoints[k]) + new Normal(0, 0.02).Sample()) * randomIensitry);
                }

                List<IndexedMassSpectralPeak> peakList = new List<IndexedMassSpectralPeak>(); //creat the timepoint
                for (int j = 0; j < 200; j++)
                {
                   peakList.Add(new IndexedMassSpectralPeak(mz: 500, intensity: intensities[j], j, timesPoints[j] + randomRetention ));                                   
                }

                if (i % 10 == 0)
                {
                    xic = new XIC(peakList, 500, new SpectraFileInfo("", "", 1, 1, 1),true);
                }

                else
                {
                    xic = new XIC(peakList, 500, new SpectraFileInfo("", "", 1, 1, 1), false);
                }                
                xicDictionary.Add(i, xic);
            }

            ParallelSearch ps = new ParallelSearch(xicDictionary);
            ps.run();


        }



        [Test]
        public void TestSingleThread()
        {
            Random random = new Random();
            Dictionary<int, XIC> xicDictionary = new Dictionary<int, XIC>();
            XIC xic = null;
            for (int i = 0; i < 230; i++)
            {

                int randomIensitry = random.Next(1, 3);
                int randomRetention = random.Next(-3, 3);

                List<double> timesPoints = Enumerable.Range(0, 200).Select(t => 10 + (double)t / 10.0).ToList();
                List<double> intensities = new List<double>();
                var peak_1 = new Normal(14.5, 0.6);    // first peak
                var peak_2 = new Normal(16.5, 0.6);    // second peak

                for (int k = 0; k < 200; k++) //creat the intensity point
                {
                    intensities.Add((peak_1.Density(timesPoints[k]) + peak_2.Density(timesPoints[k]) + new Normal(0, 0.02).Sample()) * randomIensitry);
                }

                List<IndexedMassSpectralPeak> peakList = new List<IndexedMassSpectralPeak>(); //creat the timepoint
                for (int j = 0; j < 200; j++)
                {
                    peakList.Add(new IndexedMassSpectralPeak(mz: 500, intensity: intensities[j], j, timesPoints[j] + randomRetention));
                }

                if (i % 10 == 0)
                {
                    xic = new XIC(peakList, 500, new SpectraFileInfo("", "", 1, 1, 1), true);
                }
                else
                {
                    xic = new XIC(peakList, 500, new SpectraFileInfo("", "", 1, 1, 1), false);
                }
                xicDictionary.Add(i, xic);
            }



        }
    }
}
