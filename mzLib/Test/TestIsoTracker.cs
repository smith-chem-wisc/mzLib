using FlashLFQ;
using FlashLFQ.IsoTracker;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Interpolation;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Chemistry;
using Easy.Common.Extensions;
using MzLibUtil;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestIsoTracker
    {
        // Test the IsobaricPeptideGroup class
        [Test]
        public static void TestIsobaricPeptideGroup()
        {
            // Description: Test the IsobaricPeptideGroup class
            // In this testing, we will create a new IsobaricPeptideGroup and check the properties
            List<Identification> ids = new List<Identification>
            {
                new Identification(null, "PEPTIDE", "PEPTIDE", 500.0001, 5.0, 2, new List<ProteinGroup>(), null, true, 0.9, 0, false),
                new Identification(null, "PEPTIDE", "PEPTIDE", 500.0007, 5.0, 2, new List<ProteinGroup>(), null, true, 0.9, 0, false)
            };
            IsobaricPeptideGroup group = new IsobaricPeptideGroup("PEPTIDE", 500.0, ids);
            Assert.IsNotNull(group);
            Assert.AreEqual(group.BaseSequence, "PEPTIDE");
            Assert.AreEqual(group.MonoisotopicMass, 500.0);
            Assert.IsNotNull(group.Identifications);
            Assert.AreEqual(group.Identifications.Count, 2);

            IsobaricPeptideGroup group2 = new IsobaricPeptideGroup("PEPTIDE", 500.0, null);
            Assert.IsNotNull(group2); 
            Assert.IsNotNull(group2.Identifications); // The identifications should be initialized to an empty list
            Assert.AreEqual(group2, group); // Because the two groups have the same BaseSequence and MonoisotopicMass, they should be equal
            Assert.AreEqual(group.GetHashCode(), group2.GetHashCode());

            IsobaricPeptideGroup group3 = new IsobaricPeptideGroup("PEPTIDE", 501.0, ids);
            Assert.IsNotNull(group3);
            Assert.IsNotNull(group3.Identifications);
            Assert.AreNotEqual(group3, group); // Because the two groups have difference MonoisotopicMass, they should be not equal
            Assert.AreNotEqual(group3.GetHashCode(), group.GetHashCode());

        }

        //Test the basic function in IsoTracker
        [Test]
        public static void TestIsoTrackerIdFilter_constructor()
        {
            // Test the constructor of IsoTrackerIdFilter
            IsoTrackerIdFilter target = new IsoTrackerIdFilter(new List<char>() { 'A', 'B', 'C' });
            Assert.IsNotNull(target.TargetMotifs);
            Assert.AreEqual(target.TargetMotifs.Count, 3);
            Assert.IsNotNull(target.TargetMotifPattern);
            Assert.AreEqual(target.TargetMotifPattern.Count, 3);
            // TargetMotif is null
            IsoTrackerIdFilter target_null = new IsoTrackerIdFilter();
            Assert.IsNotNull(target_null.TargetMotifs);
            Assert.AreEqual(target_null.TargetMotifs.Count, 0);
            Assert.IsNull(target_null.TargetMotifPattern);
            // TargetMotif is invalid char
            IsoTrackerIdFilter target_invalid = new IsoTrackerIdFilter(new List<char>() { '8', '#' });
            Assert.IsNotNull(target_invalid.TargetMotifs);
            Assert.AreEqual(target_invalid.TargetMotifs.Count, 0);
            Assert.IsNull(target_invalid.TargetMotifPattern);
        }

        [Test]
        public static void TestIsoTrackerIdFilter_Pattern()
        {
            List<char> targetMotifs = new List<char>() { 'S', 'T', 'Y', 'N' }; //Motif: S, T, Y, N
            IsoTrackerIdFilter target = new IsoTrackerIdFilter(targetMotifs);
            Assert.IsNotNull(target.TargetMotifs);
            Assert.AreEqual(target.TargetMotifs.Count, 4);
            Assert.IsNotNull(target.TargetMotifPattern);
            Assert.AreEqual(target.TargetMotifPattern.Count, 4);
            Assert.AreEqual(target.ContainsAcceptableModifiedResidue(null), false);
            Assert.AreEqual(target.ContainsAcceptableModifiedResidue("PEPTIEKNY"), false);
            Assert.AreEqual(target.ContainsAcceptableModifiedResidue("PEPTIEKN[mod]Y"), true);
            Assert.AreEqual(target.ContainsAcceptableModifiedResidue("PEP[mod]TIEKNY"), false);
            Assert.AreEqual(target.ContainsAcceptableModifiedResidue("PEP[mod]TIEKN[mod]Y"), true);
            Assert.AreEqual(target.ContainsAcceptableModifiedResidue("P[mod]E[mod]P[mod]TI[mod]E[mod]K[mod]NY"), false); 
        }

        [Test]
        public static void TestIsoTrackerIdFilter_FilterPeptide()
        {
            List<char> modificationOption = new List<char>() { 'S', 'T', 'Y', 'N' }; //Motif: S, T, Y, N
            IsoTrackerIdFilter searchingTarget = new IsoTrackerIdFilter(modificationOption);
            Assert.IsNotNull(searchingTarget.TargetMotifs);

            HashSet<ProteinGroup> proteinGroups = new HashSet<ProteinGroup>() {new ProteinGroup("ProteinA", "gene", "ass")};
            List<Peptide> peptides = new List<Peptide>()
            {
                new Peptide(null,null,false,proteinGroups), // We can handle the null peptide
                new Peptide("PEPTIEKNY", "PEPTIEKNY", false,  proteinGroups),
              new Peptide("PEPT[Modification]IEKNY", "PEPTIEKNY", false,  proteinGroups), // motif is T  
              new Peptide("PEPTIEKN[modification]Y", "PEPTIEKNY", false, proteinGroups), // motif is N  
              new Peptide("PEPTIEKNY[modificaiton]", "PEPTIEKNY", false, proteinGroups), // motif is Y  
              new Peptide("PES[modification]TIEKNY", "PESTIEKNY", false, proteinGroups)  // motif is S
            };
            List<Peptide> filteredPeptides = peptides.Where(p=> searchingTarget.ContainsAcceptableModifiedResidue(p.Sequence)).ToList();
            Assert.AreEqual(filteredPeptides.Count, 4);
            Assert.IsTrue(filteredPeptides[0].Sequence != "PEPTIEKNY");

            List<char> modificationOption_2 = new List<char>() { 'S', 'T'};
            IsoTrackerIdFilter searchingTarget_2 = new IsoTrackerIdFilter(modificationOption_2);
            Assert.IsNotNull(searchingTarget_2.TargetMotifs);
            List<Peptide> filteredPeptides_2 = peptides.Where(p => searchingTarget_2.ContainsAcceptableModifiedResidue(p.Sequence)).ToList();
            Assert.AreEqual(filteredPeptides_2.Count, 2);
            Assert.IsTrue(filteredPeptides_2[0].Sequence == "PEPT[Modification]IEKNY");
            Assert.IsTrue(filteredPeptides_2[1].Sequence == "PES[modification]TIEKNY");

        }

        [Test]
        public static void TestGetTargeMz_case1()
        {
            // Description: Test the GetTargetMz function in FlashLfqEngine
            // In this testing, we will check the isobaricPeptideGroup and targetMzs output
            // All three ids are isobaric peptides with the same monoisotopic mass, so they should be grouped together and generate only 5 target m/z values

            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "XICData");
            string file1 = "20100604_Velos1_TaGe_SA_A549_3_first_noRt";
            string file2 = "20100604_Velos1_TaGe_SA_A549_3_second_noRt";
            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(testDataDirectory, file1 + ".mzML"), "one", 1, 1, 1);
            SpectraFileInfo f1r2 = new SpectraFileInfo(Path.Combine(testDataDirectory, file2 + ".mzML"), "two", 1, 1, 1);

            List<Identification> ids = new List<Identification>();
            ids.Add(new Identification(f1r1, "PEPTIDE", "PE[A]PTIDE", 500.000, 5.0, 2, new List<ProteinGroup>(), null, true, 0.9, 0, false));
            ids.Add(new Identification(f1r1, "PEPTIDE", "PEPT[A]IDE", 500.000, 5.0, 2, new List<ProteinGroup>(), null, true, 0.9, 0, false));
            ids.Add(new Identification(f1r1, "PEPTIDE", "PEPTID[A]E", 500.000, 5.0, 2, new List<ProteinGroup>(), null, true, 0.9, 0, false));

            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: false,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: false,
                isoTracker: true,
                requireMultipleIdsInOneFiles: false,
                maxThreads: 1);
            var results = engine.Run();

            // Check the grouping function by manually create a new isobaric peptide group,
            // then with the engine's peptide groups
            var isobaricGroup = ids.Where(p => p.BaseSequence != p.ModifiedSequence && !p.IsDecoy)
                    .GroupBy(p => (p.BaseSequence, MonoisotopicMassGroup: Math.Round(p.MonoisotopicMass / 0.001)))
                    .Select(g => new IsobaricPeptideGroup(g.Key.BaseSequence, g.Key.MonoisotopicMassGroup, g.ToList()))
                    .ToList();
            Assert.IsTrue(isobaricGroup.Count == engine.PeptideGroupsForIsoTracker.Count); // The group count should be the same
            Assert.AreEqual(isobaricGroup, engine.PeptideGroupsForIsoTracker); // The monoMass and the baseSeq should be the same

            // Manually calculate the targetMzs
            var peptideIsotopicDistribution = engine.ModifiedSequenceToIsotopicDistribution;
            var interestedMass = new HashSet<float>();

            // Randomly select one sequence from the isobaric group.
            // The isobaric group should have the same monoisotopic mass, so we can use any sequence to calculate the target m/z values
            var randomOneSeq = isobaricGroup.SelectMany(p => p.Identifications)
                .First().ModifiedSequence;
            var monoMass = ids.FirstOrDefault(id => id.ModifiedSequence == randomOneSeq)?.MonoisotopicMass;
            var distribution = peptideIsotopicDistribution[randomOneSeq];
            List<float> massesOfInterest = distribution.Select(p => (float)(p.massShift + monoMass)).ToList();
            var minMass = massesOfInterest.Min();
            var maxMass = massesOfInterest.Max();
            massesOfInterest.Add(minMass - (float)Constants.C13MinusC12); // Except the isotopic distribution, we also add three mass (one before the lowest, two over the biggest)
            massesOfInterest.Add(maxMass + (float)Constants.C13MinusC12);
            massesOfInterest.Add(maxMass + 2f * (float)Constants.C13MinusC12);
            interestedMass.AddRange(massesOfInterest);
            var manuallyValues = interestedMass.Select(p => p.ToMz(2)).ToList();
            manuallyValues.Sort();

            // Check the target m/z values in experimental values and manually calculated values
            var experimentalValues = engine.GetTargetMz();
            Assert.AreEqual(experimentalValues.Count, 5); // Each peptide has 5 isotopic peaks, but they all have the same monoMass, so the total should be 5
            Assert.IsNotNull(experimentalValues);
            CollectionAssert.AreEqual(experimentalValues, manuallyValues);
        }

        [Test]
        public static void TestGetTargeMz_case2()
        {
            // Description: Test the GetTargetMz function in FlashLfqEngine
            // In this testing, we will check the isobaricPeptideGroup and targetMzs output
            // All three ids are isobaric peptides with the different monoisotopic mass, so they should be grouped together and generate 15(3*5) target m/z values

            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "XICData");
            string file1 = "20100604_Velos1_TaGe_SA_A549_3_first_noRt";
            string file2 = "20100604_Velos1_TaGe_SA_A549_3_second_noRt";
            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(testDataDirectory, file1 + ".mzML"), "one", 1, 1, 1);

            List<Identification> ids = new List<Identification>();
            ids.Add(new Identification(f1r1, "PEPTIDE", "PE[A]PTIDE", 500.0096, 5.0, 2, new List<ProteinGroup>(), null, true, 0.9, 0, false));
            ids.Add(new Identification(f1r1, "PEPTIDE", "PEPT[A]IDE", 500.0104, 5.0, 2, new List<ProteinGroup>(), null, true, 0.9, 0, false));
            ids.Add(new Identification(f1r1, "PEPTIDE", "PEPTID[A]E", 500.01, 5.0, 2, new List<ProteinGroup>(), null, true, 0.9, 0, false));
            
            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: false,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: false,
                isoTracker: true,
                requireMultipleIdsInOneFiles: false,
                maxThreads: 1);
            var results = engine.Run();

            // Check the grouping function by manually create a new isobaric peptide group,
            // then with the engine's peptide groups
            var isobaricGroup = ids.Where(p => p.BaseSequence != p.ModifiedSequence && !p.IsDecoy)
                    .GroupBy(p => (p.BaseSequence, MonoisotopicMassGroup: Math.Round(p.MonoisotopicMass / 0.001)))
                    .Select(g => new IsobaricPeptideGroup(g.Key.BaseSequence, g.Key.MonoisotopicMassGroup, g.ToList()))
                    .ToList();
            Assert.IsTrue(isobaricGroup.Count == engine.PeptideGroupsForIsoTracker.Count); // The group count should be the same
            Assert.AreEqual(isobaricGroup, engine.PeptideGroupsForIsoTracker); // The monoMass and the baseSeq should be the same

            // Manually calculate the targetMzs
            var peptideIsotopicDistribution = engine.ModifiedSequenceToIsotopicDistribution;
            var interestedMass = new HashSet<float>();

            IEnumerable<string> uniqueFullSeq = isobaricGroup.SelectMany(p => p.Identifications)
                .Select(p => p.ModifiedSequence)
                .Distinct();
            foreach (var seq in uniqueFullSeq)
            {
                var monoMass = ids.FirstOrDefault(id => id.ModifiedSequence == seq)?.MonoisotopicMass;
                var distribution = peptideIsotopicDistribution[seq];
                List<float> massesOfInterest = distribution.Select(p => (float)(p.massShift + monoMass)).ToList();
                var minMass = massesOfInterest.Min();
                var maxMass = massesOfInterest.Max();
                massesOfInterest.Add(minMass - (float)Constants.C13MinusC12); // Except the isotopic distribution, we also add three mass (one before the lowest, two over the biggest)
                massesOfInterest.Add(maxMass + (float)Constants.C13MinusC12);
                massesOfInterest.Add(maxMass + 2f * (float)Constants.C13MinusC12);
                interestedMass.AddRange(massesOfInterest);
            }
            var manuallyValues = interestedMass.Select(p => p.ToMz(2)).ToList();
            manuallyValues.Sort();

            // Check the target m/z values in experimental values and manually calculated values
            var experimentalValues = engine.GetTargetMz();
            Assert.AreEqual(experimentalValues.Count, 15); // Each peptide has 5 isotopic peaks, and we have 3 peptides with difference monoMass, so the total should be 15
            Assert.IsNotNull(experimentalValues);
            CollectionAssert.AreEqual(experimentalValues, manuallyValues);
        }

        [Test]
        public static void TestIndexPeakPrune()
        {
            // Description: Test the peak indexing engine pruning function
            // In this test, we will create the targetMzs from the ids to prune the indexPeaks.
            // After pruning, the index engine should only keep the peaks with the target m/z values.
            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "XICData");
            string file1 = "20100604_Velos1_TaGe_SA_A549_3_first_noRt";
            string file2 = "20100604_Velos1_TaGe_SA_A549_3_second_noRt";
            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(testDataDirectory, file1 + ".mzML"), "one", 1, 1, 1);
            SpectraFileInfo f1r2 = new SpectraFileInfo(Path.Combine(testDataDirectory, file2 + ".mzML"), "two", 1, 1, 1);
            List<SpectraFileInfo> fileInfos = new List<SpectraFileInfo>{ f1r1, f1r2 };
            List<Identification> ids = new List<Identification>();
            ids.Add(new Identification(f1r1, "PEPTIDE", "PE[A]PTIDE", 1201.53952, 5.0, 2, new List<ProteinGroup>(), null, true, 0.9, 0, false));
            ids.Add(new Identification(f1r2, "PEPTIDE", "PEPT[A]IDE", 1201.53952, 5.0, 2, new List<ProteinGroup>(), null, true, 0.9, 0, false));

            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: false,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: false,
                isoTracker: true,
                requireMultipleIdsInOneFiles: false,
                maxThreads: 1);
            var results = engine.Run();

            Tolerance ppTolerance = new PpmTolerance(10);

            // Before pruning, the index engine should contain all the peaks
            var indexEngine = PeakIndexingEngine.InitializeIndexingEngine(f1r1);
            var xic = indexEngine.GetXic(866.9613, 108.02801, ppTolerance, 2);
            var xic_2 = indexEngine.GetXic(601.77704, 108.72372, ppTolerance, 2);
            Assert.IsTrue(xic.Any()); // Those two peaks are found
            Assert.IsTrue(xic_2.Any());

            // Prune the index engine that only keep the mz with 601.77704
            var targetMzs = engine.GetTargetMz(); // Use the Ids to generate the target m/z values
            indexEngine.PruneIndex(targetMzs);
            xic = indexEngine.GetXic(866.9613, 108.02801, ppTolerance, 2);
            xic_2 = indexEngine.GetXic(601.77704, 108.72372, ppTolerance, 2);
            Assert.IsTrue(!xic.Any()); // Then there is not peak found with mz 866.9613
            Assert.IsTrue(xic_2.Any()); // But the peak with mz 601.77704 should still be found
        }


        // Test the XIC class

        [Test]
        public static void TestXICConstructor()
        {
            // Arrange
            var peaks = new List<IIndexedPeak>
            {
                new IndexedMassSpectralPeak(100, 200, 0, 1.0f),
                new IndexedMassSpectralPeak(100, 210, 1, 2.0f),
                new IndexedMassSpectralPeak(100, 220, 2, 3.0f),
                new IndexedMassSpectralPeak(100, 230, 3, 4.0f),
                new IndexedMassSpectralPeak(100, 240, 4, 5.0f)
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
            //Description: Test the linear spline interpolation and differentiation
            //The testing model is a linear function y = 100x, where x is the time point and y is the intensity
            //The slope will be 100 and the second derivative will be 0

            var peaks = new List<IIndexedPeak>
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
            //Description: Test the peak alignment function
            //The testing model is a triangle peak with the Apex.
            //The Apex of three peaks are 3, 3.1, 2.9 min
            //The time shift should be 0.1 min for the peak2 and -0.1 min for the peak3

            //The Apex of the peak1 is at 3.0
            var peaks1 = new List<IIndexedPeak>
            {
                new IndexedMassSpectralPeak(100, 10, 0, 1),
                new IndexedMassSpectralPeak(100, 20, 1, 2),
                new IndexedMassSpectralPeak(100, 30, 2, 3),
                new IndexedMassSpectralPeak(100, 20, 3, 4),
                new IndexedMassSpectralPeak(100, 10, 4, 5)
            };

            //The Apex of the peak2 is at 3.1
            var peaks2 = new List<IIndexedPeak>
            {
                new IndexedMassSpectralPeak(100, 10, 0, 1.1),
                new IndexedMassSpectralPeak(100, 20, 1, 2.1),
                new IndexedMassSpectralPeak(100, 30, 2, 3.1),
                new IndexedMassSpectralPeak(100, 20, 3, 4.1),
                new IndexedMassSpectralPeak(100, 10, 4, 5.1)
            };

            //The Apex of the peak is at 2.9
            var peaks3 = new List<IIndexedPeak>
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
            //Description: Test the cubic spline interpolation
            //The testing model has less than 5 points that cannot build the cubic spline
            //The cubic spline should be null

            // Arrange
            var peaks = new List<IIndexedPeak>
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
            //Peak setting: normal disstriubtion with mean 20 and standard deviation 0.6
            var peak = new Normal(20, 0.6);
            double peakFindingMass = 100.0;
            List<double> timesPoints = Enumerable.Range(0, 200).Select(t => 10 + (double)t / 10.0).ToList(); //The timeLine from 10 to 30 mins

            List<IIndexedPeak> peaks = new List<IIndexedPeak>();
            for (int k = 0; k < 200; k++)
            {
                peaks.Add(new IndexedMassSpectralPeak(mz:peakFindingMass, intensity: peak.Density(timesPoints[k]), 0, retentionTime: timesPoints[k]));
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
            //Peak setting: normal disstriubtion with mean 20 and standard deviation 0.6
            var peak = new Normal(20, 0.6);
            double peakFindingMass = 100.0;
            List<double> timesPoints = Enumerable.Range(0, 200).Select(t => 10 + (double)t / 10.0).ToList(); //The timeLine from 10 to 30 mins

            List<IIndexedPeak> peaks = new List<IIndexedPeak>();
            for (int k = 0; k < 200; k++)
            {
                peaks.Add(new IndexedMassSpectralPeak(mz: peakFindingMass, intensity: peak.Density(timesPoints[k]), 0, retentionTime: timesPoints[k]));
            }

            var spectraFile = new SpectraFileInfo("test.raw", "TestCondition", 1, 1, 1);
            var xic = new XIC(peaks, peakFindingMass, spectraFile);
            xic.FindExtrema();
            Assert.IsNotNull(xic.Extrema);

            var Apex = xic.Extrema.Where(p => p.Type == ExtremumType.Maximum).OrderBy(p => p.Intensity).Last(); // The apex should be located at 20 min
            Assert.AreEqual(Apex.RetentionTime, 20, 0.1);
            
        }


        // Test the XICGroup class
        [Test]
        public static void TestXICGroupConstructor()
        {
            // XIC building
            var peaks = new List<IIndexedPeak>
            {
                new IndexedMassSpectralPeak(100, 200, 0, 1.0),
                new IndexedMassSpectralPeak(100, 210, 1, 2.0),
                new IndexedMassSpectralPeak(100, 220, 2, 3.0),
                new IndexedMassSpectralPeak(100, 230, 3, 4.0),
                new IndexedMassSpectralPeak(100, 240, 4, 5.0)
            };

            var peaks2 = new List<IIndexedPeak>
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

            var xic = new XIC(peaks, 100.0, spectraFile, true, new List<Identification>() { id1 });
            var xic2 = new XIC(peaks2, 100.0, spectraFile, false, new List<Identification>() { id2 });
            var xicGroup = new XICGroups(new List<XIC>() { xic, xic2 });
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
            var idList = new List<Identification>() { id1, id2 };
            CollectionAssert.AreEqual(idList.Select(p => p.ModifiedSequence), xicGroup.IdList.Select(p => p.ModifiedSequence));
            // Because the XICGroup will create a new MovedId base on the original Id, so the IdList should be the same as the original IdList. But cannot use the Assert to test.

        }

        [Test]
        public static void TestXICGroup_RtDict()
        {
            //Description: Test the peakAlignment function in the XICGroup
            //The testing has three normal distribution XIC peaks.
            //The Apex of three peaks are 3, 3.1, 2.9 min
            //The time shift should be 0.1 min for the peak2 and -0.1 min for the peak3

            //The Apex of the peak1 is at 3.0
            var peaks1 = new List<IIndexedPeak>
            {
                new IndexedMassSpectralPeak(100, 10, 0, 1),
                new IndexedMassSpectralPeak(100, 20, 1, 2),
                new IndexedMassSpectralPeak(100, 30, 2, 3),
                new IndexedMassSpectralPeak(100, 20, 3, 4),
                new IndexedMassSpectralPeak(100, 10, 4, 5)
            };

            //The Apex of the peak2 is at 3.1
            var peaks2 = new List<IIndexedPeak>
            {
                new IndexedMassSpectralPeak(100, 10, 0, 1.1),
                new IndexedMassSpectralPeak(100, 20, 1, 2.1),
                new IndexedMassSpectralPeak(100, 30, 2, 3.1),
                new IndexedMassSpectralPeak(100, 20, 3, 4.1),
                new IndexedMassSpectralPeak(100, 10, 4, 5.1)
            };

            //The Apex of the peak is at 2.9
            var peaks3 = new List<IIndexedPeak>
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
            //Description: Test the IdList in the XICGroup
            //The testing model has three XICs, one of this XIC has no Id, then it borrows one Id from the first XIC
            //If the Id is borrowed, the Id will not be added into the IdList
            //The IdList should contain the Ids from the first and the third XIC

            var peaks = new List<IIndexedPeak>
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
        public static void TestXICGroup_Tracking()
        {
            //Description: Test the peak tracking function in the XICGroup
            //The testing has three normal distribution XIC peaks.
            //The Apex of three peaks are 20, 23, 17 min
            //The time shift should be +3 min for the peak2 and -3 min for the peak3
            //One shared peak should contain the Apex at 20 min.

            //The timeLine from 10 to 30 mins
            List<double> timesPoints = Enumerable.Range(0, 200).Select(t => 10 + (double)t / 10.0)
                        .ToList();

            List<double> intensities_P1 = new List<double>();
            List<double> intensities_P2 = new List<double>();
            List<double> intensities_P3 = new List<double>();

            var peak = new Normal(20, 0.6); // first peak

            // Create three peaks with different retention times
            for (int k = 0; k < 200; k++)
            {
                intensities_P1.Add((peak.Density(timesPoints[k])));
                intensities_P2.Add((peak.Density(timesPoints[k] - 3)));
                intensities_P3.Add((peak.Density(timesPoints[k] + 3)));
            }

            List<IIndexedPeak> p1List = new List<IIndexedPeak>(); // create the timepoints
            List<IIndexedPeak> p2List = new List<IIndexedPeak>();
            List<IIndexedPeak> p3List = new List<IIndexedPeak>();

            for (int j = 0; j < 200; j++)
            {
                p1List.Add(new IndexedMassSpectralPeak(mz: 500, intensity: intensities_P1[j], j, timesPoints[j]));
                p2List.Add(new IndexedMassSpectralPeak(mz: 500, intensity: intensities_P2[j], j, timesPoints[j]));
                p3List.Add(new IndexedMassSpectralPeak(mz: 500, intensity: intensities_P3[j], j, timesPoints[j]));
            }

            XIC xic = new XIC(p1List, 500, new SpectraFileInfo("", "", 1, 1, 1), true);
            XIC xic_2 = new XIC(p2List, 500, new SpectraFileInfo("", "", 1, 1, 1), false);
            XIC xic_3 = new XIC(p3List, 500, new SpectraFileInfo("", "", 1, 1, 1), false);

            var xicGroup = new XICGroups(new List<XIC> { xic, xic_2, xic_3 }, cutOff: 0);


            // Assert RT shift
            Assert.AreEqual(xicGroup.RTDict[0], 0, 0.001);
            Assert.AreEqual(xicGroup.RTDict[1], -3, 0.001);
            Assert.AreEqual(xicGroup.RTDict[2], 3, 0.001);

            // Assert the sharedPeak projection, the time point should be the same as the reference XIC
            Assert.AreEqual(xicGroup.SharedExtrema.Count, xicGroup.ExtremaInRef.Count);
            Assert.AreEqual(xicGroup.SharedExtrema.Select(p => p.RetentionTime), xicGroup.ExtremaInRef.Select(p => p.Key));

            // Assert the Apex in this case. The apex should be the at 20 min and extremum type is maximum
            var Apex = xicGroup.SharedExtrema.Where(p => p.Type == ExtremumType.Maximum).OrderBy(p => p.Intensity).Last();
            Assert.IsNotNull(Apex);
            Assert.AreEqual(Apex.RetentionTime, 20, 0.1);

            // Assert shared peaks
            var sharedPeaks = xicGroup.SharedPeaks;
            Assert.IsNotNull(sharedPeaks);

            //There should be only one peak that contains the apex
            var ApexPeaks = sharedPeaks
                .Where(p => p.StartRT <= Apex.RetentionTime && p.EndRT >= Apex.RetentionTime).ToList();
            Assert.IsNotNull(ApexPeaks);
            Assert.AreEqual(ApexPeaks.Count, 1);
        }


        [Test]
        public static void TestCombinedSearching()
        {
            //Description: Test the IsoTracker in the FlashLFQ, checking items include the peak tracking and the peak output
            //There are three XIC included isobaric peaks that with 3 min gap.

            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "XICData");
            string outputDirectory = Path.Combine(testDataDirectory, "testFlash");
            Directory.CreateDirectory(outputDirectory);

            string psmFile = Path.Combine(testDataDirectory, "AllPSMs.psmtsv");
            string file1 = "20100604_Velos1_TaGe_SA_A549_3_first_noRt";
            string file2 = "20100604_Velos1_TaGe_SA_A549_3_second_noRt";
            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(testDataDirectory, file1 + ".mzML"), "one", 1, 1, 1);
            SpectraFileInfo f1r2 = new SpectraFileInfo(Path.Combine(testDataDirectory, file2 + ".mzML"), "two", 1, 1, 1);

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });
                //skip the header
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains(file1))
                {
                    file = f1r1;
                }
                else if (split[0].Contains(file2))
                {
                    file = f1r2;
                }


                var decoy = split[33];
                var contaminant = split[32];
                var qvalue = double.Parse(split[51]);
                var qvalueNotch = double.Parse(split[54]);
                string baseSequence = split[13];
                string fullSequence = split[14];

                if (baseSequence.Contains("|") || fullSequence.Contains("|"))
                {
                    continue;
                }
                if (decoy.Contains("Y") || contaminant.Contains("Y") || qvalue > 0.01 || qvalueNotch > 0.01)
                {
                    continue;
                }

                double monoMass = double.Parse(split[23].Split(new char[] { '|' }).First());
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[26].Split(new char[] { '|' });
                List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups);
                ids.Add(id);

            }

            // Test there is no crush while turning on the whole searching options
            Assert.DoesNotThrow(() =>
            {
                var engine = new FlashLfqEngine(ids,
                    matchBetweenRuns: true,
                    requireMsmsIdInCondition: true,
                    useSharedPeptidesForProteinQuant: true,
                    isoTracker: true,
                    requireMultipleIdsInOneFiles: false,
                    maxThreads: 1);
                var results = engine.Run();
            });

            // Test there is no crush while turning on the IsoTracker and MBR
            Assert.DoesNotThrow(() =>
            {
                var engine = new FlashLfqEngine(ids,
                    matchBetweenRuns: true,
                    requireMsmsIdInCondition: false,
                    useSharedPeptidesForProteinQuant: false,
                    isoTracker: true,
                    requireMultipleIdsInOneFiles: false,
                    maxThreads: 1);
                var results = engine.Run();
            });


        }

        // Test the IsoTracker search function
        [Test]
        public static void TestPeakOutput()
        {
            //Description: Test the IsoTracker in the FlashLFQ, checking items include the peak tracking and the peak output
            //There are three XIC included isobaric peaks that with 3 min gap.

            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "XICData");
            string outputDirectory = Path.Combine(testDataDirectory, "testFlash");
            Directory.CreateDirectory(outputDirectory);

            string psmFile = Path.Combine(testDataDirectory, "AllPSMs.psmtsv");
            string file1 = "20100604_Velos1_TaGe_SA_A549_3_first_noRt";
            string file2 = "20100604_Velos1_TaGe_SA_A549_3_second_noRt";
            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(testDataDirectory, file1 + ".mzML"), "one", 1, 1, 1);
            SpectraFileInfo f1r2 = new SpectraFileInfo(Path.Combine(testDataDirectory, file2 + ".mzML"), "two", 1, 1, 1);

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });
                //skip the header
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains(file1))
                {
                    file = f1r1;
                }
                else if (split[0].Contains(file2))
                {
                    file = f1r2;
                }


                var decoy = split[33];
                var contaminant = split[32];
                var qvalue = double.Parse(split[51]);
                var qvalueNotch = double.Parse(split[54]);
                string baseSequence = split[13];
                string fullSequence = split[14];

                if (baseSequence.Contains("|") || fullSequence.Contains("|"))
                {
                    continue;
                }
                if (decoy.Contains("Y") || contaminant.Contains("Y") || qvalue > 0.01 || qvalueNotch > 0.01)
                {
                    continue;
                }

                double monoMass = double.Parse(split[23].Split(new char[] { '|' }).First());
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[26].Split(new char[] { '|' });
                List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups);
                ids.Add(id);

            }

            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: false,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: false,
                isoTracker: true,
                requireMultipleIdsInOneFiles: false,
                maxThreads: 1);
            var results = engine.Run();

            results.WriteResults(Path.Combine(outputDirectory, "peaks.tsv"), Path.Combine(outputDirectory, "peptides.tsv"), Path.Combine(outputDirectory, "proteins.tsv"), Path.Combine(outputDirectory, "bayesian.tsv"), true);


            List<string> peaksList = File.ReadAllLines(Path.Combine(outputDirectory, "peaks.tsv")).ToList();
            List<string> peptidesList = File.ReadAllLines(Path.Combine(outputDirectory, "peptides.tsv")).ToList();
            List<string> proteinsList = File.ReadAllLines(Path.Combine(outputDirectory, "proteins.tsv")).ToList();

            // Assert new peptide headers for the Isotracker
            var peptideHeader = File.ReadLines(Path.Combine(outputDirectory, "peptides.tsv")).ToList();
            List<string> peptideHeaderList = peptideHeader[0].Split('\t').ToList();
            List<string> expectedPeptideHeaders = new List<string> { "Sequence", "Base Sequence", "Peak Order","Protein Groups", "Gene Names", "Organism", "Intensity_" + file1, "Intensity_" + file2,
                "RetentionTime (min)_"+file1, "RetentionTime (min)_" + file2, "Detection Type_"+file1, "Detection Type_" + file2 };
            Assert.AreEqual(12, peptideHeaderList.Count);
            CollectionAssert.AreEqual(expectedPeptideHeaders, peptideHeaderList);

            //Assert the isopeaks data
            //Only modified sequences can have isobaric peaks
            //In this case the intensity and retention time of the isobaric peaks from two runs should be the same
            //The detection type should be different, one from MSMSScan and the other from MBR
            foreach (var line in File.ReadLines(Path.Combine(outputDirectory, "peptides.tsv")).Skip(1))
            {
                var fullSeq = line.Split('\t')[0];
                var baseSeq = line.Split('\t')[1];
                if (fullSeq.Contains("Isopeptide", StringComparison.OrdinalIgnoreCase)) // The full sequence for Isobaric case should contain the word "peak"
                {
                    CollectionAssert.AreNotEqual(baseSeq, fullSeq);
                    var intensityRun1 = line.Split('\t')[6];
                    var intensityRun2 = line.Split('\t')[7];
                    var retentionTimeRun1 = line.Split('\t')[8];
                    var retentionTimeRun2 = line.Split('\t')[9];
                    var detectionTypeRun1 = line.Split('\t')[10];
                    var detectionTypeRun2 = line.Split('\t')[11];
                    CollectionAssert.AreEqual(intensityRun1, intensityRun2);
                    CollectionAssert.AreEqual(retentionTimeRun1, retentionTimeRun2);
                    CollectionAssert.AreNotEqual(detectionTypeRun1, detectionTypeRun2);
                }

            }

            //check that all rows including header have the same number of elements
            Assert.AreEqual(1, peaksList.Select(l => l.Split('\t').Length).Distinct().ToList().Count);
            Assert.AreEqual(1, peptidesList.Select(l => l.Split('\t').Length).Distinct().ToList().Count);
            Assert.AreEqual(1, proteinsList.Select(l => l.Split('\t').Length).Distinct().ToList().Count);

            Directory.Delete(outputDirectory, true);
        }


        [Test]
        public static void TestIsoSequence_Ambiguous()
        {
            //Description: Test the IsoTracker in the FlashLFQ, checking the algorithm can correctly recognize the IsoID
            //IsoID: DIVENY[Common Variable:Oxidation on M]FMR   should be the same as DIVENYFM[Common Variable:Oxidation on M]R

            //Try to turn on the MBR and Isotracker at the same time
            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "XICData");
            string outputDirectory = Path.Combine(testDataDirectory, "testFlash");
            Directory.CreateDirectory(outputDirectory);

            string psmFile = Path.Combine(testDataDirectory, "AllPSMs_IsoID.psmtsv");
            string file1 = "20100604_Velos1_TaGe_SA_A549_3_first_noRt";
            string file2 = "20100604_Velos1_TaGe_SA_A549_3_second_noRt";
            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(testDataDirectory, file1 + ".mzML"), "one", 1, 1, 1);
            SpectraFileInfo f1r2 = new SpectraFileInfo(Path.Combine(testDataDirectory, file2 + ".mzML"), "two", 1, 1, 1);

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });
                //skip the header
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains(file1))
                {
                    file = f1r1;
                }
                else if (split[0].Contains(file2))
                {
                    file = f1r2;
                }


                var decoy = split[33];
                var contaminant = split[32];
                var qvalue = double.Parse(split[51]);
                var qvalueNotch = double.Parse(split[54]);
                string baseSequence = split[13];
                string fullSequence = split[14];

                if (baseSequence.Contains("|") || fullSequence.Contains("|"))
                {
                    continue;
                }
                if (decoy.Contains("Y") || contaminant.Contains("Y") || qvalue > 0.01 || qvalueNotch > 0.01)
                {
                    continue;
                }

                double monoMass = double.Parse(split[23].Split(new char[] { '|' }).First());
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[26].Split(new char[] { '|' });
                List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups);
                ids.Add(id);

            }

            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: false,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: false,
                isoTracker: true,
                requireMultipleIdsInOneFiles: false,
                maxThreads: 1);
            var results = engine.Run();

            results.WriteResults(Path.Combine(outputDirectory, "peaks.tsv"), Path.Combine(outputDirectory, "peptides.tsv"), Path.Combine(outputDirectory, "proteins.tsv"), null, true);


            List<string> peaksList = File.ReadAllLines(Path.Combine(outputDirectory, "peaks.tsv")).Skip(1).ToList();
            List<string> peptidesList = File.ReadAllLines(Path.Combine(outputDirectory, "peptides.tsv")).Skip(1).ToList();
            List<string> proteinsList = File.ReadAllLines(Path.Combine(outputDirectory, "proteins.tsv")).Skip(1).ToList();


            // Check the output: there are one kind of IsoID(modifiedPeptide), then two isobaricPeptide, and four isobaricPeaks
            int modifiedPeptideNum = ids.GroupBy(p => new { p.BaseSequence, p.MonoisotopicMass }).Count();
            int isobaricPeptideNum = modifiedPeptideNum * 2;
            int isobaricPeakNum = isobaricPeptideNum * 2;

            Assert.AreEqual(proteinsList.Count, modifiedPeptideNum);
            Assert.AreEqual(peptidesList.Count, isobaricPeptideNum);
            Assert.AreEqual(peaksList.Count, isobaricPeakNum);
            List<string> expectedSequence = new List<string> { "DIVENYFM[Common Variable:Oxidation on M]R", "DIVENYF[Common Variable:Oxidation on M]MR", "DIVENY[Common Variable:Oxidation on M]FMR" };

            //Check the detectionType of each peak, in this case all  peaks is IsoTrack_Ambiguous
            //The output sequence should be the same as the expected sequence
            foreach (var peak in peaksList)
            {
                var peakSeq = peak.Split('\t')[2].Split('|').ToList();
                var detectionType = peak.Split('\t')[16];
                Assert.AreEqual(peakSeq, expectedSequence);
                CollectionAssert.AreEqual(detectionType, "IsoTrack_Ambiguous");
            }

            //Check the detectionType of each peptide, in this case all peptides are IsoTrack_Ambiguous
            //The output sequence should be the same as the expected sequence
            foreach (var peptide in peptidesList)
            {
                string sequences = peptide.Split('\t')[0];
                var detectionType_run1 = peptide.Split('\t')[10];
                var detectionType_run2 = peptide.Split('\t')[11];
                Assert.AreEqual(detectionType_run1, "IsoTrack_Ambiguous");
                Assert.AreEqual(detectionType_run2, "IsoTrack_Ambiguous");

                foreach (var seq in expectedSequence)
                {
                    Assert.IsTrue(sequences.Contains(seq));
                }
            }

            Directory.Delete(outputDirectory, true);
        }

        [Test]
        public static void TestIsoSequence_MonoIsotopicMassTolerance()
        {
            //Description: Test the IsoTracker in the FlashLFQ, checking the algorithm can correctly recognize the IsoID
            //IsoID: DIVENY[Common Variable:Oxidation on M]FMR   should be the same as DIVENYFM[Common Variable:Oxidation on M]R
            //The Monoisotopic mass are 1201.5436, 1201.5437, 1201.5438, they should be recognized as the same IsoID

            //Try to turn on the MBR and Isotracker at the same time
            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "XICData");
            string outputDirectory = Path.Combine(testDataDirectory, "testFlash");
            Directory.CreateDirectory(outputDirectory);

            string psmFile = Path.Combine(testDataDirectory, "AllPSMs_IsoID_MonoIsotopicmassTolerance.psmtsv");
            string file1 = "20100604_Velos1_TaGe_SA_A549_3_first_noRt";
            string file2 = "20100604_Velos1_TaGe_SA_A549_3_second_noRt";
            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(testDataDirectory, file1 + ".mzML"), "one", 1, 1, 1);
            SpectraFileInfo f1r2 = new SpectraFileInfo(Path.Combine(testDataDirectory, file2 + ".mzML"), "two", 1, 1, 1);

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });
                //skip the header
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains(file1))
                {
                    file = f1r1;
                }
                else if (split[0].Contains(file2))
                {
                    file = f1r2;
                }


                var decoy = split[33];
                var contaminant = split[32];
                var qvalue = double.Parse(split[51]);
                var qvalueNotch = double.Parse(split[54]);
                string baseSequence = split[13];
                string fullSequence = split[14];

                if (baseSequence.Contains("|") || fullSequence.Contains("|"))
                {
                    continue;
                }
                if (decoy.Contains("Y") || contaminant.Contains("Y") || qvalue > 0.01 || qvalueNotch > 0.01)
                {
                    continue;
                }

                double monoMass = double.Parse(split[23].Split(new char[] { '|' }).First());
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[26].Split(new char[] { '|' });
                List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups);
                ids.Add(id);

            }

            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: false,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: false,
                isoTracker: true,
                requireMultipleIdsInOneFiles: false,
                maxThreads: 1);
            var results = engine.Run();

            results.WriteResults(Path.Combine(outputDirectory, "peaks.tsv"), Path.Combine(outputDirectory, "peptides.tsv"), Path.Combine(outputDirectory, "proteins.tsv"), null, true);


            List<string> peaksList = File.ReadAllLines(Path.Combine(outputDirectory, "peaks.tsv")).Skip(1).ToList();
            List<string> peptidesList = File.ReadAllLines(Path.Combine(outputDirectory, "peptides.tsv")).Skip(1).ToList();
            List<string> proteinsList = File.ReadAllLines(Path.Combine(outputDirectory, "proteins.tsv")).Skip(1).ToList();


            // Check the output: there are one kind of IsoID(modifiedPeptide), then two isobaricPeptide, and four isobaricPeaks
            int modifiedPeptideNum = 1;
            int isobaricPeptideNum = modifiedPeptideNum * 2;
            int isobaricPeakNum = isobaricPeptideNum * 2;

            Assert.AreEqual(proteinsList.Count, modifiedPeptideNum);
            Assert.AreEqual(peptidesList.Count, isobaricPeptideNum);
            Assert.AreEqual(peaksList.Count, isobaricPeakNum);


            Directory.Delete(outputDirectory, true);
        }

        [Test]
        public static void TestIsoSequence_Mixture()
        {
            //Test the IsoTracker in the FlashLFQ, checking the algorithm can correctly recognize the IsoID
            //There is one ID for the peak1, and two IsoID for the peak2
            //The quantified peak1 will be MSMS and Isotrack_MBR, and the quantified peak2 will be IsoTrack_Ambiguous

            //Try to turn on the MBR and Isotracker at the same time
            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "XICData");
            string outputDirectory = Path.Combine(testDataDirectory, "testFlash");
            Directory.CreateDirectory(outputDirectory);

            string psmFile = Path.Combine(testDataDirectory, "AllPSMs_IsoID_Mixture.psmtsv");
            string file1 = "20100604_Velos1_TaGe_SA_A549_3_first_noRt";
            string file2 = "20100604_Velos1_TaGe_SA_A549_3_second_noRt";
            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(testDataDirectory, file1 + ".mzML"), "one", 1, 1, 1);
            SpectraFileInfo f1r2 = new SpectraFileInfo(Path.Combine(testDataDirectory, file2 + ".mzML"), "two", 1, 1, 1);

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });
                //skip the header
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains(file1))
                {
                    file = f1r1;
                }
                else if (split[0].Contains(file2))
                {
                    file = f1r2;
                }


                var decoy = split[33];
                var contaminant = split[32];
                var qvalue = double.Parse(split[51]);
                var qvalueNotch = double.Parse(split[54]);
                string baseSequence = split[13];
                string fullSequence = split[14];

                if (baseSequence.Contains("|") || fullSequence.Contains("|"))
                {
                    continue;
                }
                if (decoy.Contains("Y") || contaminant.Contains("Y") || qvalue > 0.01 || qvalueNotch > 0.01)
                {
                    continue;
                }

                double monoMass = double.Parse(split[23].Split(new char[] { '|' }).First());
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[26].Split(new char[] { '|' });
                List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups);
                ids.Add(id);

            }

            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: false,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: false,
                isoTracker: true,
                requireMultipleIdsInOneFiles: false,
                maxThreads: 1);
            var results = engine.Run();

            results.WriteResults(Path.Combine(outputDirectory, "peaks.tsv"), Path.Combine(outputDirectory, "peptides.tsv"), Path.Combine(outputDirectory, "proteins.tsv"), null, true);


            List<string> peaksList = File.ReadAllLines(Path.Combine(outputDirectory, "peaks.tsv")).Skip(1).ToList();
            List<string> peptidesList = File.ReadAllLines(Path.Combine(outputDirectory, "peptides.tsv")).Skip(1).ToList();
            List<string> proteinsList = File.ReadAllLines(Path.Combine(outputDirectory, "proteins.tsv")).Skip(1).ToList();

            // Check the output: there are one kind of IsoID(modifiedPeptide), then two isobaricPeptide, and four isobaricPeaks
            int modifiedPeptideNum = ids.GroupBy(p => new { p.BaseSequence, p.MonoisotopicMass }).Count();
            int isobaricPeptideNum = modifiedPeptideNum * 2;
            int isobaricPeakNum = isobaricPeptideNum * 2;

            Assert.AreEqual(proteinsList.Count, modifiedPeptideNum);
            Assert.AreEqual(peptidesList.Count, isobaricPeptideNum);
            Assert.AreEqual(peaksList.Count, isobaricPeakNum);
            List<string> expectedSequence_Peak1 = new List<string> { "DIVENYFM[Common Variable:Oxidation on M]R" };
            List<string> expectedSequence_Peak2 = new List<string> { "DIVENYF[Common Variable:Oxidation on M]MR", "DIVENY[Common Variable:Oxidation on M]FMR" };


            //Check the detectionType of each peptide, in this case all peptides are IsoTrack_Ambiguous
            //The output sequence should be the same as the expected sequence
            foreach (var peptide in peptidesList)
            {
                string outputSequence = peptide.Split('\t')[0];
                var peakOrder = int.Parse(peptide.Split('\t')[2]);

                if (peakOrder == 1)
                {
                    // Check the detectionType of each peptide, in this case all peptides are IsoTrack_Ambiguous
                    var detectionType_run1 = peptide.Split('\t')[10];
                    var detectionType_run2 = peptide.Split('\t')[11];
                    Assert.AreEqual(detectionType_run1, "IsoTrack_MBR");
                    Assert.AreEqual(detectionType_run2, "IsoTrack_MSMS");

                    //The output sequence should be the same as the expected sequence
                    foreach (var seq in expectedSequence_Peak1)
                    {
                        Assert.IsTrue(outputSequence.Contains(seq));
                    }
                }
                else if (peakOrder == 2)
                {
                    // Check the detectionType of each peptide, in this case all peptides are IsoTrack_Ambiguous
                    var detectionType_run1 = peptide.Split('\t')[10];
                    var detectionType_run2 = peptide.Split('\t')[11];
                    Assert.AreEqual(detectionType_run1, "IsoTrack_Ambiguous");
                    Assert.AreEqual(detectionType_run2, "IsoTrack_Ambiguous");

                    //The output sequence should be the same as the expected sequence
                    foreach (var seq in expectedSequence_Peak2)
                    {
                        Assert.IsTrue(outputSequence.Contains(seq));
                    }

                }

            }

            Directory.Delete(outputDirectory, true);
        }

        [Test]
        public static void TestIsoSequence_CombinedTesting()
        {
            //In this test, there are two unmodifiedPeptide, one ambiguous isoPeptide set, and two normal isobaric peptides.


            //Try to turn on the MBR and Isotracker at the same time
            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "XICData");
            string outputDirectory = Path.Combine(testDataDirectory, "testFlash");
            Directory.CreateDirectory(outputDirectory);

            string psmFile = Path.Combine(testDataDirectory, "AllPSMs_IsoID_Combined.psmtsv");
            string file1 = "20100604_Velos1_TaGe_SA_A549_3_first_noRt";
            string file2 = "20100604_Velos1_TaGe_SA_A549_3_second_noRt";
            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(testDataDirectory, file1 + ".mzML"), "one", 1, 1, 1);
            SpectraFileInfo f1r2 = new SpectraFileInfo(Path.Combine(testDataDirectory, file2 + ".mzML"), "two", 1, 1, 1);

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });
                //skip the header
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains(file1))
                {
                    file = f1r1;
                }
                else if (split[0].Contains(file2))
                {
                    file = f1r2;
                }


                var decoy = split[33];
                var contaminant = split[32];
                var qvalue = double.Parse(split[51]);
                var qvalueNotch = double.Parse(split[54]);
                string baseSequence = split[13];
                string fullSequence = split[14];

                if (baseSequence.Contains("|") || fullSequence.Contains("|"))
                {
                    continue;
                }
                if (decoy.Contains("Y") || contaminant.Contains("Y") || qvalue > 0.01 || qvalueNotch > 0.01)
                {
                    continue;
                }

                double monoMass = double.Parse(split[23].Split(new char[] { '|' }).First());
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[26].Split(new char[] { '|' });
                List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups);
                ids.Add(id);

            }

            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: false,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: false,
                isoTracker: true,
                requireMultipleIdsInOneFiles: false,
                maxThreads: 1);
            var results = engine.Run();

            results.WriteResults(Path.Combine(outputDirectory, "peaks.tsv"), Path.Combine(outputDirectory, "peptides.tsv"), Path.Combine(outputDirectory, "proteins.tsv"), null, true);


            List<string> peaksList = File.ReadAllLines(Path.Combine(outputDirectory, "peaks.tsv")).Skip(1).ToList();
            List<string> peptidesList = File.ReadAllLines(Path.Combine(outputDirectory, "peptides.tsv")).Skip(1).ToList();
            List<string> proteinsList = File.ReadAllLines(Path.Combine(outputDirectory, "proteins.tsv")).Skip(1).ToList();

            // Check the output: there are one kind of IsoID(modifiedPeptide), then two isobaricPeptide, and four isobaricPeaks
            int ambiguityPeptideNum = ids.GroupBy(p => new { p.BaseSequence, p.MonoisotopicMass })
                .Where(p => p.Count() > 2).Count();
            int normalIsoPeptideNum = ids.GroupBy(p => new { p.BaseSequence, p.MonoisotopicMass })
                .Where(p => p.Count() == 2).Count();
            int unModifiedPeptideNum = ids.Where(p => p.BaseSequence == p.ModifiedSequence).Count();

            Assert.AreEqual(ambiguityPeptideNum, 1);
            Assert.AreEqual(normalIsoPeptideNum, 2);
            Assert.AreEqual(unModifiedPeptideNum, 2);

            int QPeptideNum = ambiguityPeptideNum * 2 + normalIsoPeptideNum * 2 + unModifiedPeptideNum;
            Assert.IsTrue(QPeptideNum == 8 && QPeptideNum == peptidesList.Count);


            //Check the detectionType of each Isobaric peptide, in this case
            //For isobaric peptide: one is MSMS and the other is IsoTrack_MBR
            //For the ambiguous peptide: both are IsoTrack_Ambiguous
            foreach (var peptide in peptidesList.Where(p => p.Split('\t')[0].Contains("Isopeptide")))
            {
                string baseSequence = peptide.Split('\t')[1];
                var peakOrder = int.Parse(peptide.Split('\t')[2]);
                var detectionType_run1 = peptide.Split('\t')[10];
                var detectionType_run2 = peptide.Split('\t')[11];

                if (baseSequence == "GICEMELLDGK")
                {
                    if (peakOrder == 1)
                    {
                        Assert.AreEqual(detectionType_run1, "IsoTrack_MSMS");
                        Assert.AreEqual(detectionType_run2, "IsoTrack_MBR");
                    }
                    if (peakOrder == 2)
                    {
                        Assert.AreEqual(detectionType_run1, "IsoTrack_MBR");
                        Assert.AreEqual(detectionType_run2, "IsoTrack_MSMS");
                    }
                }
                if (baseSequence == "LLQFYAETCPAPER")
                {
                    if (peakOrder == 1)
                    {
                        Assert.AreEqual(detectionType_run1, "IsoTrack_MSMS");
                        Assert.AreEqual(detectionType_run2, "IsoTrack_MBR");
                    }
                    if (peakOrder == 2)
                    {
                        Assert.AreEqual(detectionType_run1, "IsoTrack_MBR");
                        Assert.AreEqual(detectionType_run2, "IsoTrack_MSMS");
                    }
                }
                if (baseSequence == "DIVENYFMR")
                {
                    if (peakOrder == 1)
                    {
                        Assert.AreEqual(detectionType_run1, "IsoTrack_MBR");
                        Assert.AreEqual(detectionType_run2, "IsoTrack_MSMS");
                    }
                    if (peakOrder == 2)
                    {
                        Assert.AreEqual(detectionType_run1, "IsoTrack_Ambiguous");
                        Assert.AreEqual(detectionType_run2, "IsoTrack_Ambiguous");
                    }
                }
            }

            Directory.Delete(outputDirectory, true);
        }

        // Test the IsoTracker searchTarget
        [Test]
        public static void TestRun_SearchingTarget()
        {
            //Description: we will upload a motifList for IsoTracker
            //Only peptide with motif on N can be searched
            //In this case, only one kind of peptide can be searched: baseSequence PEPNINEN -> PEPN[Mod]INEN, PEPNIN[Mod]EN, PEPNINEN[Mod]
            // Run 1 with PEPNIN[Mod]EN, PEPNINEN[Mod]
            // Run 2 with PEPN[Mod]INEN
            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "XICData");
            string outputDirectory = Path.Combine(testDataDirectory, "testFlash");
            Directory.CreateDirectory(outputDirectory);

            string psmFile = Path.Combine(testDataDirectory, "AllPSMs_SearchTarget.psmtsv");
            string file1 = "20100604_Velos1_TaGe_SA_A549_3_first_noRt";
            string file2 = "20100604_Velos1_TaGe_SA_A549_3_second_noRt";
            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(testDataDirectory, file1 + ".mzML"), "one", 1, 1, 1);
            SpectraFileInfo f1r2 = new SpectraFileInfo(Path.Combine(testDataDirectory, file2 + ".mzML"), "two", 1, 1, 1);

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });
                //skip the header
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains(file1))
                {
                    file = f1r1;
                }
                else if (split[0].Contains(file2))
                {
                    file = f1r2;
                }


                var decoy = split[33];
                var contaminant = split[32];
                var qvalue = double.Parse(split[51]);
                var qvalueNotch = double.Parse(split[54]);
                string baseSequence = split[13];
                string fullSequence = split[14];

                if (baseSequence.Contains("|") || fullSequence.Contains("|"))
                {
                    continue;
                }
                if (decoy.Contains("Y") || contaminant.Contains("Y") || qvalue > 0.01 || qvalueNotch > 0.01)
                {
                    continue;
                }

                double monoMass = double.Parse(split[23].Split(new char[] { '|' }).First());
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[26].Split(new char[] { '|' });
                List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups);
                ids.Add(id);

            }

            var motifList = new List<char> { 'N' };
            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: false,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: false,
                isoTracker: true,
                requireMultipleIdsInOneFiles: false,
                maxThreads: 1,
                motifsList: motifList);
            var results = engine.Run();

            results.WriteResults(Path.Combine(outputDirectory, "peaks.tsv"), Path.Combine(outputDirectory, "peptides.tsv"), Path.Combine(outputDirectory, "proteins.tsv"), null, true);
            List<string> peptidesList = File.ReadAllLines(Path.Combine(outputDirectory, "peptides.tsv")).Skip(1).ToList();

            // Every isobaric peptide should contain the motif "N" (asparagine)
            foreach (var peptide in peptidesList)
            {
                var IsIsobaric = peptide.Split('\t')[0].Contains("Isopeptide");
                if (IsIsobaric)
                {
                    var fullSequence = peptide.Split('\t')[0];
                    var baseSequence = peptide.Split('\t')[1];
                    Assert.IsTrue(fullSequence.Contains("N["));
                    Assert.IsTrue(fullSequence.Contains("]"));
                    Assert.IsTrue(baseSequence.Contains("N"));
                }
            }
            Directory.Delete(outputDirectory, true);

            // In second test, we don't include the motif list. Every isobaric peptide should be searched
            var engine_noMotif = new FlashLfqEngine(ids,
                matchBetweenRuns: false,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: false,
                isoTracker: true,
                requireMultipleIdsInOneFiles: false,
                maxThreads: 1,
                motifsList: null);
            var results_noMotif = engine_noMotif.Run();
            var isobaricPeptide = results_noMotif.IsobaricPeptideDict;
            // Total 3 isobaric peptides
            Assert.IsTrue(isobaricPeptide.Count == 3);
        }

        [Test]
        public static void TestRun_IDChecking()
        {
            //Description: we will turn on the IDchecking for IsoTracker
            //Only when one XIC with more than one id, we do the searching
            //In this case, run 1 has 4 ids (pepA_1, pepA_2, pepB_1, pepC_1)
            //run 2 has 3 ids (pepA_1, pepB_1, pepC_1)
            //Only pepA can be searched.
            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "XICData");
            string outputDirectory = Path.Combine(testDataDirectory, "testFlash");
            Directory.CreateDirectory(outputDirectory);

            string psmFile = Path.Combine(testDataDirectory, "AllPSMs_IdChecking.psmtsv");
            string file1 = "20100604_Velos1_TaGe_SA_A549_3_first_noRt";
            string file2 = "20100604_Velos1_TaGe_SA_A549_3_second_noRt";
            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(testDataDirectory, file1 + ".mzML"), "one", 1, 1, 1);
            SpectraFileInfo f1r2 = new SpectraFileInfo(Path.Combine(testDataDirectory, file2 + ".mzML"), "two", 1, 1, 1);

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });
                //skip the header
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;
                if (split[0].Contains(file1))
                {
                    file = f1r1;
                }
                else if (split[0].Contains(file2))
                {
                    file = f1r2;
                }

                var decoy = split[33];
                var contaminant = split[32];
                var qvalue = double.Parse(split[51]);
                var qvalueNotch = double.Parse(split[54]);
                string baseSequence = split[13];
                string fullSequence = split[14];

                if (baseSequence.Contains("|") || fullSequence.Contains("|"))
                {
                    continue;
                }
                if (decoy.Contains("Y") || contaminant.Contains("Y") || qvalue > 0.01 || qvalueNotch > 0.01)
                {
                    continue;
                }

                double monoMass = double.Parse(split[23].Split(new char[] { '|' }).First());
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[26].Split(new char[] { '|' });
                List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups);
                ids.Add(id);

            }

            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: false,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: false,
                isoTracker: true,
                maxThreads: 1);
            var results = engine.Run();

            results.WriteResults(Path.Combine(outputDirectory, "peaks.tsv"), Path.Combine(outputDirectory, "peptides.tsv"), Path.Combine(outputDirectory, "proteins.tsv"), null, true);
            // only one isobaric peptide should be searched
            Assert.IsTrue(results.IsobaricPeptideDict.Count == 1);

            // through the IDchecking, we can only get the one isobaric peptide (pepA)
            List<string> peptidesList = File.ReadAllLines(Path.Combine(outputDirectory, "peptides.tsv")).Skip(1).ToList();
            foreach (var peptide in peptidesList)
            {
                var IsIsobaric = peptide.Split('\t')[0].Contains("Isopeptide");
                if (IsIsobaric)
                {
                    var baseSequence = peptide.Split('\t')[1];
                    // only peptdieA can be isopeptide
                    Assert.AreEqual(baseSequence, "PEPTIDEA");
                    Assert.AreNotEqual(baseSequence, "PEPTIDEB");
                    Assert.AreNotEqual(baseSequence, "PEPTIDEC");
                }
            }

            // I remove the first id in the list, therefore both run only have one Id for each peptide.
            // There is no any isobaric peptide
            ids.RemoveAt(1);
            var new_engine = new FlashLfqEngine(ids,
                matchBetweenRuns: false,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: false,
                isoTracker: true,
                maxThreads: 1);
            var results_2 = new_engine.Run();
            Assert.IsTrue(results_2.IsobaricPeptideDict.Count == 0);
            Directory.Delete(outputDirectory, true);
        }

        [Test]
        public static void Test_mergeConflict()
        {
            //In this test, we want to ensure there is no conflict when peaks were merged in RunErrorCheck.

            //There are two isobaric peptide, Iso_A and Iso_B (with similar mass)
            //RunErrorCheck will merge the two peaks, and we hope there is no crush and the detectionType will be set as MSMSAmbiguousPeakfinding
            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "XICData");
            string outputDirectory = Path.Combine(testDataDirectory, "testFlash");
            Directory.CreateDirectory(outputDirectory);

            string psmFile = Path.Combine(testDataDirectory, "AllPSMs_IDMergeConflict.psmtsv");
            string file1 = "20100604_Velos1_TaGe_SA_A549_3_first_noRt";
            string file2 = "20100604_Velos1_TaGe_SA_A549_3_second_noRt";
            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(testDataDirectory, file1 + ".mzML"), "one", 0, 0, 1);
            SpectraFileInfo f1r2 = new SpectraFileInfo(Path.Combine(testDataDirectory, file2 + ".mzML"), "one", 1, 0, 1);

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });
                //skip the header
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains(file1))
                {
                    file = f1r1;
                }
                else if (split[0].Contains(file2))
                {
                    file = f1r2;
                }


                var decoy = split[33];
                var contaminant = split[32];
                var qvalue = double.Parse(split[51]);
                var qvalueNotch = double.Parse(split[54]);
                string baseSequence = split[13];
                string fullSequence = split[14];

                if (baseSequence.Contains("|") || fullSequence.Contains("|"))
                {
                    continue;
                }
                if (decoy.Contains("Y") || contaminant.Contains("Y") || qvalue > 0.01 || qvalueNotch > 0.01)
                {
                    continue;
                }

                double monoMass = double.Parse(split[23].Split(new char[] { '|' }).First());
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[26].Split(new char[] { '|' });
                List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups);
                ids.Add(id);

            }

            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: true,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: false,
                isoTracker: true,
                requireMultipleIdsInOneFiles: false,
                normalize: false,
                maxThreads: 1);
            var results = engine.Run();

            results.WriteResults(Path.Combine(outputDirectory, "peaks.tsv"), Path.Combine(outputDirectory, "peptides.tsv"), Path.Combine(outputDirectory, "proteins.tsv"), null, true);

            List<string> peaksList = File.ReadAllLines(Path.Combine(outputDirectory, "peaks.tsv")).Skip(1).ToList();
            List<string> peptidesList = File.ReadAllLines(Path.Combine(outputDirectory, "peptides.tsv")).Skip(1).ToList();
            foreach (var peak in peaksList)
            {
                var fullSeq = peak.Split('\t')[2];
                if (fullSeq == "PEPTIDEA|PEPTIDEB")
                {
                    Assert.AreEqual(peak.Split('\t')[16], "MSMSAmbiguousPeakfinding");
                }
            }

            foreach (var pep in peptidesList)
            {
                var baseSeq = pep.Split('\t')[1];
                var peakOrder = pep.Split('\t')[2];
                var detectionType_File1 = pep.Split('\t')[10];
                var detectionType_File2 = pep.Split('\t')[11];
                if (baseSeq == "PEPTIDEA")
                {
                    if (peakOrder == "1")
                    {
                        Assert.AreEqual(detectionType_File1, "IsoTrack_MSMS");
                        Assert.AreEqual(detectionType_File2, "IsoTrack_MBR");
                    }

                    if (peakOrder == "2")
                    {
                        Assert.AreEqual(detectionType_File1, "MSMSAmbiguousPeakfinding");
                        Assert.AreEqual(detectionType_File2, "MSMSAmbiguousPeakfinding");
                    }
                }

                if (baseSeq == "PEPTIDEB")
                {
                    Assert.AreEqual(peakOrder, "");
                    Assert.AreEqual(detectionType_File1, "MSMSAmbiguousPeakfinding");
                    Assert.AreEqual(detectionType_File2, "MSMSAmbiguousPeakfinding");
                }
            }

        }

        [Test]
        public static void Test_mergeConflict_test2()
        {
            //In this test, we want to ensure there is no conflict when peaks were merged in RunErrorCheck.

            //There are two isobaric peptide, Iso_A and Iso_B (with similar mass)
            //RunErrorCheck will merge the two peaks, total four merge isoPeaks
            //According to the definition of IsoTracker, the detectionType should be set as MSMSAmbiguousPeakfinding, and their intensity should be set as 0
            string testDataDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "XICData");
            string outputDirectory = Path.Combine(testDataDirectory, "testFlash");
            Directory.CreateDirectory(outputDirectory);

            string psmFile = Path.Combine(testDataDirectory, "AllPSMs_IDMergeConflict_2.psmtsv");
            string file1 = "20100604_Velos1_TaGe_SA_A549_3_first_noRt";
            string file2 = "20100604_Velos1_TaGe_SA_A549_3_second_noRt";
            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(testDataDirectory, file1 + ".mzML"), "one", 0, 0, 1);
            SpectraFileInfo f1r2 = new SpectraFileInfo(Path.Combine(testDataDirectory, file2 + ".mzML"), "one", 1, 0, 1);

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });
                //skip the header
                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains(file1))
                {
                    file = f1r1;
                }
                else if (split[0].Contains(file2))
                {
                    file = f1r2;
                }


                var decoy = split[33];
                var contaminant = split[32];
                var qvalue = double.Parse(split[51]);
                var qvalueNotch = double.Parse(split[54]);
                string baseSequence = split[13];
                string fullSequence = split[14];

                if (baseSequence.Contains("|") || fullSequence.Contains("|"))
                {
                    continue;
                }
                if (decoy.Contains("Y") || contaminant.Contains("Y") || qvalue > 0.01 || qvalueNotch > 0.01)
                {
                    continue;
                }

                double monoMass = double.Parse(split[23].Split(new char[] { '|' }).First());
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[26].Split(new char[] { '|' });
                List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups);
                ids.Add(id);

            }

            var engine = new FlashLfqEngine(ids,
                matchBetweenRuns: true,
                requireMsmsIdInCondition: false,
                useSharedPeptidesForProteinQuant: false,
                isoTracker: true,
                requireMultipleIdsInOneFiles: false,
                normalize: false,
                maxThreads: 1);
            var results = engine.Run();

            results.WriteResults(Path.Combine(outputDirectory, "peaks.tsv"), Path.Combine(outputDirectory, "peptides.tsv"), Path.Combine(outputDirectory, "proteins.tsv"), null, true);

            List<string> peaksList = File.ReadAllLines(Path.Combine(outputDirectory, "peaks.tsv")).Skip(1).ToList();
            List<string> peptidesList = File.ReadAllLines(Path.Combine(outputDirectory, "peptides.tsv")).Skip(1).ToList();
            foreach (var peak in peaksList)
            {
                var fullSeq = peak.Split('\t')[2];
                var retentionTime = peak.Split('\t')[6];
                double.TryParse(peak.Split('\t')[9], out double intensity);
                Assert.AreNotEqual(intensity, 0);
                if (fullSeq == "PEPTIDEA_1|PEPTIDEB_1")
                {
                    if (retentionTime == "")
                    {
                        Assert.AreEqual(peak.Split('\t')[16], "MSMSAmbiguousPeakfinding");
                    }
                }
            }
            // In IsoTracker, each peak are label as MSMSAmbiguousPeakfinding by RunErrorCheck
            foreach (var pep in peptidesList)
            {
                double.TryParse(pep.Split('\t')[6], out double intensity_File1);
                double.TryParse(pep.Split('\t')[7], out double intensity_File2);
                var detectionType_File1 = pep.Split('\t')[10];
                var detectionType_File2 = pep.Split('\t')[11];
                Assert.AreEqual(intensity_File1, 0);
                Assert.AreEqual(intensity_File2, 0);
                Assert.AreEqual(detectionType_File1, "MSMSAmbiguousPeakfinding");
                Assert.AreEqual(detectionType_File2, "MSMSAmbiguousPeakfinding");
            }

        }
    }

}
