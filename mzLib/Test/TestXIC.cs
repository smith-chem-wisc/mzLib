using FlashLFQ.Alex_project;
using FlashLFQ;
using MathNet.Numerics.Distributions;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Interpolation;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;

namespace Test
{
    internal class TestXIC
    {
        [Test]
        public static void TestXICConstructor()
        {
            // Arrange
            var peaks = new List<IndexedMassSpectralPeak>
            {
                new IndexedMassSpectralPeak(100, 200, 0, 1.0),
                new IndexedMassSpectralPeak(100, 210, 1, 2.0),
                new IndexedMassSpectralPeak(100, 220, 2, 3.0),
                new IndexedMassSpectralPeak(100, 230, 3, 4.0),
                new IndexedMassSpectralPeak(100, 240, 4, 5.0)
            };
            double peakFindingMass = 100.0;
            var spectraFile = new SpectraFileInfo("path/to/file.raw", "Condition", 1, 1, 1);
            bool isReference = true;
            var ids = new List<Identification>();

            // Act
            var xic = new XIC(peaks, peakFindingMass, spectraFile, isReference, ids);

            // Assert
            Assert.AreEqual(peaks, xic.Ms1Peaks);
            Assert.AreEqual(peakFindingMass, xic.PeakFindingMz);
            Assert.AreEqual(spectraFile, xic.SpectraFile);
            Assert.AreEqual(isReference, xic.Reference);
            Assert.AreEqual(ids, xic.Ids);
            Assert.IsNotNull(xic.LinearSpline);
            Assert.IsNotNull(xic.SmoothedCubicSpline);
        }

        [Test]
        public static void TestLinearSpline()
        {
            // Arrange
            var peaks = new List<IndexedMassSpectralPeak>
            {
                new IndexedMassSpectralPeak(100, 100, 0, 1.0),
                new IndexedMassSpectralPeak(100, 200, 1, 2.0),
                new IndexedMassSpectralPeak(100, 300, 2, 3.0),
                new IndexedMassSpectralPeak(100, 400, 3, 4.0),
                new IndexedMassSpectralPeak(100, 500, 4, 5.0)
            };
            var xic = new XIC(peaks, 100, new SpectraFileInfo("file.raw", "Condition", 1, 1, 1));


            // Assert
            Assert.IsNotNull(xic.LinearSpline);
            Assert.AreEqual(xic.LinearSpline.Interpolate(1.5),150);
            Assert.AreEqual(xic.LinearSpline.Interpolate(2.5),250);
            Assert.AreEqual(xic.LinearSpline.Interpolate(3.5), 350);
            Assert.AreEqual(xic.LinearSpline.Interpolate(4.5), 450);

            for (double timePoint = 1.0; timePoint < 5.0; timePoint = timePoint + 1.0)
            {
                Assert.AreEqual(xic.LinearSpline.Differentiate(timePoint), 100); // The slope of the linear spline should be 100
                Assert.AreEqual(xic.LinearSpline.Differentiate2(timePoint), 0); // The second derivative of the linear spline should be 0
            }
        }

        [Test]
        public static void TestPeakAlignment()
        {
            // Arrange

            //The Apex of the peak1 is at 3.0
            var peaks1 = new List<IndexedMassSpectralPeak>
            {
                new IndexedMassSpectralPeak(100, 10, 0, 1),
                new IndexedMassSpectralPeak(100, 20, 1, 2),
                new IndexedMassSpectralPeak(100, 30, 2, 3),
                new IndexedMassSpectralPeak(100, 20, 3, 4),
                new IndexedMassSpectralPeak(100, 10, 4, 5)
            };

            //The Apex of the peak2 is at 3.1
            var peaks2 = new List<IndexedMassSpectralPeak>
            {
                new IndexedMassSpectralPeak(100, 10, 0, 1.1),
                new IndexedMassSpectralPeak(100, 20, 1, 2.1),
                new IndexedMassSpectralPeak(100, 30, 2, 3.1),
                new IndexedMassSpectralPeak(100, 20, 3, 4.1),
                new IndexedMassSpectralPeak(100, 10, 4, 5.1)
            };

            //The Apex of the peak is at 2.9
            var peaks3 = new List<IndexedMassSpectralPeak>
            {
                new IndexedMassSpectralPeak(100, 10, 0, 0.9),
                new IndexedMassSpectralPeak(100, 20, 1, 1.9),
                new IndexedMassSpectralPeak(100, 30, 2, 2.9),
                new IndexedMassSpectralPeak(100, 20, 3, 3.9),
                new IndexedMassSpectralPeak(100, 10, 4, 4.9)
            };

            var spectraFile = new SpectraFileInfo("test.raw", "Condition", 1, 1, 1);
            var xic1 = new XIC(peaks1, 100, spectraFile, true);
            var xic2 = new XIC(peaks2, 100, spectraFile, false);
            var xic3 = new XIC(peaks3, 100, spectraFile, false);

            double rtTolerance = 0.0011;
            // Assert
            Assert.AreEqual(-0.1, xic2.AlignXICs(xic1), rtTolerance); // Peak 2 should be shifted to the left by 0.1
            Assert.AreEqual(0.1, xic3.AlignXICs(xic1), rtTolerance); // Peak 3 should be shifted to the right by 0.1
        }

        [Test]
        public static void TestBuildSmoothedCubicSpline_LessPoint()
        {
            // Arrange
            var peaks = new List<IndexedMassSpectralPeak>
            {
                new IndexedMassSpectralPeak(100, 100, 0, 1.0),
                new IndexedMassSpectralPeak(100, 200, 1, 2.0),
                new IndexedMassSpectralPeak(100, 300, 2, 3.0),
                new IndexedMassSpectralPeak(100, 500, 4, 5.0)
            };
            var xic = new XIC(peaks, 100, new SpectraFileInfo("file.raw", "Condition", 1, 1, 1));

            // If the number of points is less than 5, the cubic spline should not be built
            Assert.IsNull(xic.SmoothedCubicSpline);
        }

        [Test]
        public static void TestBuildSmoothedCubicSpline()
        {
            //

            //Create a normal distribution peak with the retention time of 20 min

            //Peak setting
            var peak = new Normal(20, 0.6);
            double peakFindingMass = 100.0;
            List<double> timesPoints = Enumerable.Range(0, 200).Select(t => 10 + (double)t / 10.0).ToList(); //The timeLine from 10 to 30 mins

            List<IndexedMassSpectralPeak> peaks = new List<IndexedMassSpectralPeak>();
            for (int k = 0; k < 200; k++)
            {
                peaks.Add(new IndexedMassSpectralPeak(mz: peakFindingMass, intensity: peak.Density(timesPoints[k]), 0, retentionTime: timesPoints[k]));
            }

            var spectraFile = new SpectraFileInfo("test.raw", "TestCondition", 1, 1, 1);
            var xic = new XIC(peaks, peakFindingMass, spectraFile);

            // Assert
            Assert.IsNotNull(xic.SmoothedCubicSpline);
            Assert.AreEqual(xic.SmoothedCubicSpline.GetType(), typeof(CubicSpline));
        }

        [Test]
        public static void TestExtrema()
        {
            var peak = new Normal(20, 0.6);
            double peakFindingMass = 100.0;
            List<double> timesPoints = Enumerable.Range(0, 200).Select(t => 10 + (double)t / 10.0).ToList(); //The timeLine from 10 to 30 mins

            List<IndexedMassSpectralPeak> peaks = new List<IndexedMassSpectralPeak>();
            for (int k = 0; k < 200; k++)
            {
                peaks.Add(new IndexedMassSpectralPeak(mz: peakFindingMass, intensity: peak.Density(timesPoints[k]), 0, retentionTime: timesPoints[k]));
            }

            var spectraFile = new SpectraFileInfo("test.raw", "TestCondition", 1, 1, 1);
            var xic = new XIC(peaks, peakFindingMass, spectraFile);
            xic.FindExtrema();
            Assert.IsNotNull(xic.Extrema);
        }

    }


    internal class TestXICGroup
    {
        [Test]
        public static void TestXICGroupConstructor()
        {
            // Arrange
            var peaks = new List<IndexedMassSpectralPeak>
            {
                new IndexedMassSpectralPeak(100, 200, 0, 1.0),
                new IndexedMassSpectralPeak(100, 210, 1, 2.0),
                new IndexedMassSpectralPeak(100, 220, 2, 3.0),
                new IndexedMassSpectralPeak(100, 230, 3, 4.0),
                new IndexedMassSpectralPeak(100, 240, 4, 5.0)
            };

            var peaks2 = new List<IndexedMassSpectralPeak>
            {
                new IndexedMassSpectralPeak(100, 200, 0, 1.0),
                new IndexedMassSpectralPeak(100, 210, 1, 2.0),
                new IndexedMassSpectralPeak(100, 220, 2, 3.0),
                new IndexedMassSpectralPeak(100, 230, 3, 4.0),
                new IndexedMassSpectralPeak(100, 240, 4, 5.0)
            };

            var spectraFile = new SpectraFileInfo("path/to/file.raw", "Condition", 1, 1, 1);

            var id1 = new Identification(spectraFile, "BaseSequence1", "ModifiedSequence1", 100.0, 1.0, 1, new List<ProteinGroup>(), null, true, 0, 0, false);
            var id2 = new Identification(spectraFile, "BaseSequence2", "ModifiedSequence2", 200.0, 2.0, 2, new List<ProteinGroup>(), null, true, 0, 0, false);
            
            var xic = new XIC(peaks, 100.0, spectraFile, true, new List<Identification>(){ id1 });
            var xic2 = new XIC(peaks2, 100.0, spectraFile, false, new List<Identification>(){ id2 });
            var xicGroup = new XICGroups(new List<XIC>(){ xic, xic2 });
            // Assert
            Assert.IsNotNull(xicGroup);
            Assert.AreEqual(xicGroup.ReferenceXIC, xic);
            Assert.IsNotNull(xicGroup.XICs);
            Assert.AreEqual(xicGroup.XICs.Count, 2);
            Assert.AreEqual(xic, xicGroup.XICs[0]);
            Assert.AreEqual(xic2, xicGroup.XICs[1]);


            Assert.IsNotNull(xicGroup.RTDict);
            Assert.AreEqual(xicGroup.RTDict.Count, 2);
            Assert.AreEqual(xicGroup.RTDict[0], 0.0);  //reference XIC Rt is 0
            Assert.AreEqual(xicGroup.RTDict[1], 0.0);  //non-reference XIC Rt is 0

            Assert.IsNotNull(xicGroup.IdList);
            Assert.AreEqual(xicGroup.IdList.Count, 2);
            var idList = new List<Identification>() { id1, id2};
            var movedIdList = new List<Identification>() { id1, id2 };
            CollectionAssert.AreEqual(idList.Select(p=>p.ModifiedSequence), xicGroup.IdList.Select(p => p.ModifiedSequence));
            // Because the XICGroup will create a new MovedId base on the original Id, so the IdList should be the same as the original IdList. But cannot use the Assert to test.

        }

        [Test]
        public static void TestXICGroup_RtDict()
        {
            // Arrange

            //The Apex of the peak1 is at 3.0
            var peaks1 = new List<IndexedMassSpectralPeak>
            {
                new IndexedMassSpectralPeak(100, 10, 0, 1),
                new IndexedMassSpectralPeak(100, 20, 1, 2),
                new IndexedMassSpectralPeak(100, 30, 2, 3),
                new IndexedMassSpectralPeak(100, 20, 3, 4),
                new IndexedMassSpectralPeak(100, 10, 4, 5)
            };

            //The Apex of the peak2 is at 3.1
            var peaks2 = new List<IndexedMassSpectralPeak>
            {
                new IndexedMassSpectralPeak(100, 10, 0, 1.1),
                new IndexedMassSpectralPeak(100, 20, 1, 2.1),
                new IndexedMassSpectralPeak(100, 30, 2, 3.1),
                new IndexedMassSpectralPeak(100, 20, 3, 4.1),
                new IndexedMassSpectralPeak(100, 10, 4, 5.1)
            };

            //The Apex of the peak is at 2.9
            var peaks3 = new List<IndexedMassSpectralPeak>
            {
                new IndexedMassSpectralPeak(100, 10, 0, 0.9),
                new IndexedMassSpectralPeak(100, 20, 1, 1.9),
                new IndexedMassSpectralPeak(100, 30, 2, 2.9),
                new IndexedMassSpectralPeak(100, 20, 3, 3.9),
                new IndexedMassSpectralPeak(100, 10, 4, 4.9)
            };

            var spectraFile = new SpectraFileInfo("test.raw", "Condition", 1, 1, 1);
            var xic1 = new XIC(peaks1, 100, spectraFile, true);
            var xic2 = new XIC(peaks2, 100, spectraFile, false);
            var xic3 = new XIC(peaks3, 100, spectraFile, false);

            var xicGroup = new XICGroups(new List<XIC> { xic1, xic2, xic3 });

            // Assert
            Assert.AreEqual(xicGroup.RTDict[0], xic1.AlignXICs(xic1)); 
            Assert.AreEqual(xicGroup.RTDict[1], xic2.AlignXICs(xic1)); 
            Assert.AreEqual(xicGroup.RTDict[2], xic3.AlignXICs(xic1)); 
        }

        [Test]
        public static void TestXICGroup_IdList()
        {
            // Arrange
            var peaks = new List<IndexedMassSpectralPeak>
            {
                new IndexedMassSpectralPeak(100, 10, 0, 1),
                new IndexedMassSpectralPeak(100, 20, 1, 2),
                new IndexedMassSpectralPeak(100, 30, 2, 3),
                new IndexedMassSpectralPeak(100, 20, 3, 4),
                new IndexedMassSpectralPeak(100, 10, 4, 5)
            };

            var file1 = new SpectraFileInfo("file1", "Condition", 1, 1, 1);
            var file2 = new SpectraFileInfo("file2", "Condition", 1, 1, 1);

            var id1 = new Identification(file1, "BaseSequence1", "ModifiedSequence1", 100.0, 1.0, 1, new List<ProteinGroup>(), null, true, 0, 0, false);
            var id2 = new Identification(file1, "BaseSequence2", "ModifiedSequence2", 200.0, 2.0, 2, new List<ProteinGroup>(), null, true, 0, 0, false);
            var id3 = new Identification(file2, "BaseSequence3", "ModifiedSequence3", 200.0, 2.0, 2, new List<ProteinGroup>(), null, true, 0, 0, false);

            var xic = new XIC(peaks, 100.0, file1, true, new List<Identification>() { id1 });
            var xic2 = new XIC(peaks, 100.0, file2, false, new List<Identification>() { id2 }); // The XIC borrowed the id from the first XIC, so the id will not be add into the IdList
            var xic3 = new XIC(peaks, 100.0, file2, false, new List<Identification>() { id3 });
            var xicGroup = new XICGroups(new List<XIC>() { xic, xic2, xic3 });

            // Assert
            Assert.IsNotNull(xicGroup.IdList);
            Assert.AreEqual(xicGroup.IdList.Count, 2);
            var idList = new List<Identification>() { id1, id3 };
        }

        [Test]
        public void TestXICGroup_Tracking()
        {
            List<double> timesPoints = Enumerable.Range(0, 200).Select(t => 10 + (double)t / 10.0).ToList(); //The timeLine from 10 to 30 mins

            List<double> intensities_P1 = new List<double>();
            List<double> intensities_P2 = new List<double>();
            List<double> intensities_P3 = new List<double>();

            var peak_1 = new Normal(20, 0.6);    // first peak
            var intesity = 1E6; // In order to make the peak intensity high enough to be detected
            // Create three peaks with different retention times
            for (int k = 0; k < 200; k++)
            {
                intensities_P1.Add(( peak_1.Density(timesPoints[k])));
                intensities_P2.Add((peak_1.Density(timesPoints[k] - 3)));
                intensities_P3.Add((peak_1.Density(timesPoints[k] + 3)));
            }

            List<IndexedMassSpectralPeak> p1List = new List<IndexedMassSpectralPeak>(); // create the timepoints
            List<IndexedMassSpectralPeak> p2List = new List<IndexedMassSpectralPeak>();
            List<IndexedMassSpectralPeak> p3List = new List<IndexedMassSpectralPeak>();

            for (int j = 0; j < 200; j++)
            {
                p1List.Add(new IndexedMassSpectralPeak(mz: 500, intensity: intensities_P1[j], j, timesPoints[j]));
                p2List.Add(new IndexedMassSpectralPeak(mz: 500, intensity: intensities_P2[j], j, timesPoints[j]));
                p3List.Add(new IndexedMassSpectralPeak(mz: 500, intensity: intensities_P3[j], j, timesPoints[j]));
            }

            XIC xic = new XIC(p1List, 500, new SpectraFileInfo("", "", 1, 1, 1), true);
            XIC xic_2 = new XIC(p2List, 500, new SpectraFileInfo("", "", 1, 1, 1), false);
            XIC xic_3 = new XIC(p3List, 500, new SpectraFileInfo("", "", 1, 1, 1), false);

            var xicGroup = new XICGroups(new List<XIC> { xic, xic_2, xic_3 }, cutoff: 0);

            // Assert RT shift
            Assert.AreEqual(xicGroup.RTDict[0], 0, 0.001);
            Assert.AreEqual(xicGroup.RTDict[1], -3, 0.001);
            Assert.AreEqual(xicGroup.RTDict[2], 3, 0.001);

            // Assert the sharedPeak projection, the time point should be the same as the reference XIC
            Assert.AreEqual(xicGroup.SharedExtrema.Count, xicGroup.ExtremaInRef.Count);
            Assert.AreEqual(xicGroup.SharedExtrema.Select(p=>p.RetentionTime), xicGroup.ExtremaInRef.Select(p=>p.Key));

            // Assert the Apex in this case. The apex should be the maximum peak in the shared peaks
            var ApexExetremum = xicGroup.SharedExtrema.Where(p => p.RetentionTime < 20.1 && p.RetentionTime > 19.9); // The apex should be around 20
            Assert.IsNotNull(ApexExetremum);
            Assert.AreEqual(ApexExetremum.Count(), 1);
            Extremum Apex = ApexExetremum.First();
            Assert.AreEqual(Apex.Type, ExtremumType.Maximum);

            // Assert shared peaks
            var sharedPeaks = xicGroup.SharedPeaks;
            Assert.IsNotNull(sharedPeaks);
            Assert.IsNotNull(sharedPeaks[Apex.RetentionTime]); // The shared peak should contain the apex peak




        }

    }

}
