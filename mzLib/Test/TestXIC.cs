using FlashLFQ.Alex_project;
using FlashLFQ;
using MathNet.Numerics.Distributions;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test
{
    internal class TestXIC
    {
        [Test]
        public void ExtremunTesting()
        {
            Random random = new Random();

            List<double> timesPoints = Enumerable.Range(0, 200).Select(t => 10 + (double)t / 10.0).ToList(); //The timeLine from 10 to 30 mins


            List<double> intensities_P1 = new List<double>();
            List<double> intensities_P2 = new List<double>();
            List<double> intensities_P3 = new List<double>();

            var peak_1 = new Normal(20, 0.6);    // first peak

            for (int k = 0; k < 200; k++) //creat the intensity point
            {
                intensities_P1.Add((peak_1.Density(timesPoints[k])));
                intensities_P2.Add((peak_1.Density(timesPoints[k] - 3)));
                intensities_P3.Add((peak_1.Density(timesPoints[k] + 3)));
            }

            List<IndexedMassSpectralPeak> p1List = new List<IndexedMassSpectralPeak>(); //creat the timepoint
            List<IndexedMassSpectralPeak> p2List = new List<IndexedMassSpectralPeak>(); 
            List<IndexedMassSpectralPeak> p3List = new List<IndexedMassSpectralPeak>(); 

            for (int j = 0; j < 200; j++)
            {
                p1List.Add(new IndexedMassSpectralPeak(mz: 500, intensity: intensities_P1[j], j, timesPoints[j] ));
                p2List.Add(new IndexedMassSpectralPeak(mz: 500, intensity: intensities_P2[j], j, timesPoints[j] ));
                p3List.Add(new IndexedMassSpectralPeak(mz: 500, intensity: intensities_P3[j], j, timesPoints[j] ));
            }

            XIC xic = new XIC(p1List, 500, new SpectraFileInfo("", "", 1, 1, 1), true);
            XIC xic_2 = new XIC(p2List, 500, new SpectraFileInfo("", "", 1, 1, 1), false);
            XIC xic_3 = new XIC(p3List, 500, new SpectraFileInfo("", "", 1, 1, 1), false);

            XICGroups xicGroups = new XICGroups(new List<XIC> { xic, xic_2, xic_3 },0.7);

        }
    }
}
