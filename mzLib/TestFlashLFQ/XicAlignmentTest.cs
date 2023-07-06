using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Easy.Common.Extensions;
using FlashLFQ.PeakPicking;
using MassSpectrometry.MzSpectra;
using Test.FileReadingTests;
using UsefulProteomicsDatabases;
using ChromatographicPeak = FlashLFQ.ChromatographicPeak;
using Stopwatch = System.Diagnostics.Stopwatch;
using MathNet.Numerics.Interpolation;
using SharpLearning.Containers.Extensions;

namespace TestFlashLFQ
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal class XicAlignmentTest
    {
        [Test]
        public void TestSignalAlignment()
        {
            // get the raw file paths
            SpectraFileInfo inflix = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_Inflix_Tryp_60s_3-calib.mzML"),
                "inflix", 0, 0, 0);
            SpectraFileInfo nist = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_NIST_Tryp_60s_3-calib.mzML"),
                "nist", 1, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            double monoisotopicMass = 2005.980136305;

            // In the actual data, these are identical but have mods on different positions.
            double peakFindingMass = 670.00;

            string sequence1 = "PWEPIYYLGGVFQLEK"; // CF3 on W
            double ms2RetentionTime1Inflix = 32.464;
            double ms2RetentionTime1Nist = 32.488;

            string sequence2 = "PWYYEPILGGVFQLEK"; // CF3 on Y
            double ms2RetentionTime2Inflix = 33.393; // No Ms2 was actual present here, this is the peak apex RT from MBR
            double ms2RetentionTime2Nist = 33.224;

            Identification pep1Inflix = new Identification(inflix, sequence1, sequence1, monoisotopicMass,
                ms2RetentionTime1Inflix, 3, new List<ProteinGroup> { pg });
            Identification pep1Nist = new Identification(nist, sequence1, sequence1, monoisotopicMass,
                ms2RetentionTime1Nist, 3, new List<ProteinGroup> { pg });
            
            Identification pep2Inflix = new Identification(inflix, sequence2, sequence2, monoisotopicMass,
                ms2RetentionTime2Inflix, 3, new List<ProteinGroup> { pg });
            Identification pep2Nist = new Identification(nist, sequence2, sequence2, monoisotopicMass,
                ms2RetentionTime2Nist, 3, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { pep1Inflix, pep1Nist, pep2Inflix, pep2Nist },
                normalize: true, maxThreads: 1, matchBetweenRuns: true); // peaks are only serialized if match between runs = true

            // run the engine and grab XICs
            engine.Run();
            var indexingEngine = engine.GetIndexingEngine();
            var inflixXic = indexingEngine.GetXIC(peakFindingMass, inflix);
            var nistXic = indexingEngine.GetXIC(peakFindingMass, nist);

            // generate a fake XIC with a known shift and test that AlignPeaks returns the correct shift
            double rtShift = 0.15;
            var shiftedInflixXic = inflixXic.Select(p =>
                    new IndexedMassSpectralPeak(p.Mz, p.Intensity, p.ZeroBasedMs1ScanIndex, p.RetentionTime + rtShift))
                .ToList();
            double shiftedOffset = XicProcessing.AlignPeaks(inflixXic, shiftedInflixXic, 100);
            Assert.That(shiftedOffset, Is.EqualTo(-1*rtShift).Within(0.01));

            // Calculate similarity of the two real XICs with and w/o alignment.
            SpectralSimilarity xicSimilarity = new SpectralSimilarity(
                P_XArray: inflixXic.Select(p => p.RetentionTime).ToArray(),
                P_YArray: inflixXic.Select(p => p.Intensity).ToArray(),
                Q_XArray: nistXic.Select(p => p.RetentionTime).ToArray(),
                Q_YArray: nistXic.Select(p => p.Intensity).ToArray(),
                SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak,
                toleranceInPpm: 250, // Not sure what the ideal ppmTolerance is here. 100 ppm resulted in 34 not matched points, 250 ppm matches all points
                allPeaks: true,
                filterOutBelowThisMz: 0);
            double? rawXicAngle = xicSimilarity.SpectralContrastAngle();

            double nistOffset = XicProcessing.AlignPeaks(inflixXic, nistXic);
            
            SpectralSimilarity alignedXicSimilarity = new SpectralSimilarity(
                P_XArray: inflixXic.Select(p => p.RetentionTime).ToArray(),
                P_YArray: inflixXic.Select(p => p.Intensity).ToArray(),
                Q_XArray: nistXic.Select(p => p.RetentionTime + nistOffset).ToArray(),
                Q_YArray: nistXic.Select(p => p.Intensity).ToArray(),
                SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak,
                toleranceInPpm: 250, // Not sure what the ideal ppmTolerance is here. 100 ppm resulted in 34 not matched points, 250 ppm matches all points
                allPeaks: true,
                filterOutBelowThisMz: 0);
            double? alignedXicAngle = alignedXicSimilarity.SpectralContrastAngle();

            Assert.That(alignedXicAngle, Is.GreaterThan(rawXicAngle));
        }

        [Test]
        public static void TestReconcileExtrema()
        {
            SpectraFileInfo nist = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_NIST_Tryp_60s_3-calib.mzML"),
                "nist", 1, 0, 0);
            SpectraFileInfo inflix = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_Inflix_Tryp_60s_3-calib.mzML"),
                "inflix", 0, 0, 0);

            // create IDs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            double monoisotopicMass = 2005.980136305;
            double peakFindingMass = 670.00;

            string baseSequence = "PWYEPIYLGGVFQLEK";
            string modSeq = "PWYEPIY[CF3:CF3 on Y]LGGVFQLEK";
            double rt = 33.22371;

            Identification y7ModId = new Identification(nist, baseSequence, modSeq, monoisotopicMass,
                rt, 3, new List<ProteinGroup> { pg });
            Identification y7ModIdInflix = new Identification(inflix, baseSequence, modSeq, monoisotopicMass,
                rt, 3, new List<ProteinGroup> { pg });
            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { y7ModId, y7ModIdInflix },
                normalize: false, maxThreads: 1, matchBetweenRuns: true, quantifyAmbiguousPeptides: true); // peaks are only serialized if match between runs = true

            // run the engine and grab XICs
            var results = engine.Run();

            var indexingEngine = engine.GetIndexingEngine();
            var nistPeaks = indexingEngine.GetXIC(peakFindingMass, nist);
            var inflixPeaks = indexingEngine.GetXIC(peakFindingMass, inflix);
            double inflixAdjustment = XicProcessing.AlignPeaks(nistPeaks, inflixPeaks, 100);

            Xic nistXic = new Xic(nistPeaks, peakFindingMass, nist, 0, referenceXic: true);
            Xic inflixXic = new Xic(inflixPeaks, peakFindingMass, inflix, inflixAdjustment, referenceXic: false);

            var matchedExtrema = XicProcessing.ReconcileExtrema(nistXic.Extrema,
                new List<List<Extremum>> { inflixXic.Extrema });

            var refType = matchedExtrema[0].Select(e => e.Type).ToArray();

            Assert.AreEqual(matchedExtrema[0].Length, 16);
            Assert.AreEqual(matchedExtrema[1].Length, 16);
        }

        [Test]
        public static void GetRtPairsTest()
        {
            var ex1 = new Extremum(1, 0, ExtremumType.Maximum);
            var ex2 = new Extremum(2, 0, ExtremumType.Maximum);
            var ex3 = new Extremum(3, 0, ExtremumType.Maximum);
            var ex4 = new Extremum(4, 0, ExtremumType.Maximum);
            var ex5 = new Extremum(5, 0, ExtremumType.Maximum);
            Extremum[] refArray = { ex1, ex2, ex3, ex4, ex5 };
            Extremum[] expArray = { ex1, ex3, ex5 };

            var pairs = XicProcessing.GetRetentionTimePairs(refArray, expArray);
            Assert.AreEqual(pairs, new List<(Extremum, Extremum)>
            {
                (ex1, ex1),
                (ex3, ex3),
                (ex5, ex5)
            });

            var ex32 = new Extremum(3.2, 0, ExtremumType.Maximum);
            var ex47 = new Extremum(4.7, 0, ExtremumType.Maximum);
            var ex51 = new Extremum(5.1, 0, ExtremumType.Maximum);
            var ex53 = new Extremum(5.3, 0, ExtremumType.Maximum);
            var ex55 = new Extremum(5.5, 0, ExtremumType.Maximum);
            
            Extremum[] refArray2 = { ex1, ex3, ex51, ex53, ex55 };
            Extremum[] expArray2 = { ex1, ex1, ex3, ex32, ex4, ex47, ex5, ex5 };
            pairs = XicProcessing.GetRetentionTimePairs(refArray2, expArray2);
            Assert.AreEqual(pairs, new List<(Extremum, Extremum)>
            {
                (ex1, ex1),
                (ex3, ex3),
                (ex51, ex5)
            });
        }

        [Test]
        public static void SplineScratchPad()
        {
            // get the raw file paths
            SpectraFileInfo inflix = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_Inflix_Tryp_60s_3-calib.mzML"),
                "inflix", 0, 0, 0);
            SpectraFileInfo nist = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_NIST_Tryp_60s_3-calib.mzML"),
                "nist", 1, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            double monoisotopicMass = 2005.980136305;

            // In the actual data, these are identical but have mods on different positions.
            double peakFindingMass = 670.00;

            string sequence1 = "PWEPIYYLGGVFQLEK"; // CF3 on W
            double ms2RetentionTime1Inflix = 32.464;
            double ms2RetentionTime1Nist = 32.488;

            string sequence2 = "PWYYEPILGGVFQLEK"; // CF3 on Y
            double ms2RetentionTime2Inflix = 33.393; // No Ms2 was actual present here, this is the peak apex RT from MBR
            double ms2RetentionTime2Nist = 33.224;

            Identification pep1Inflix = new Identification(inflix, sequence1, sequence1, monoisotopicMass,
                ms2RetentionTime1Inflix, 3, new List<ProteinGroup> { pg });
            Identification pep1Nist = new Identification(nist, sequence1, sequence1, monoisotopicMass,
                ms2RetentionTime1Nist, 3, new List<ProteinGroup> { pg });

            Identification pep2Inflix = new Identification(inflix, sequence2, sequence2, monoisotopicMass,
                ms2RetentionTime2Inflix, 3, new List<ProteinGroup> { pg });
            Identification pep2Nist = new Identification(nist, sequence2, sequence2, monoisotopicMass,
                ms2RetentionTime2Nist, 3, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { pep1Inflix, pep1Nist, pep2Inflix, pep2Nist },
                normalize: true, maxThreads: 1, matchBetweenRuns: true); // peaks are only serialized if match between runs = true

            // run the engine and grab XICs
            engine.Run();
            var indexingEngine = engine.GetIndexingEngine();
            var inflixXic = indexingEngine.GetXIC(peakFindingMass, inflix);
            var nistXic = indexingEngine.GetXIC(peakFindingMass, nist);

            double nistOffset = XicProcessing.AlignPeaks(inflixXic, nistXic);

            // Creates a spline, or a smoothed representation of the XIC which allows for up-sampling of the data
            // Unsure about whether to use cubic or linear. Can maybe be a delegate argument?
            var inflixSpline = CubicSpline.InterpolateAkimaSorted(
                inflixXic.Select(p => p.RetentionTime).ToArray(), inflixXic.Select(p => p.Intensity).ToArray());
            var nistSpline = CubicSpline.InterpolateAkimaSorted(
                nistXic.Select(p => p.RetentionTime + nistOffset).ToArray(), nistXic.Select(p => p.Intensity).ToArray());

            var inflixStationaryPoint = inflixSpline.StationaryPoints();
            Array.Sort(inflixStationaryPoint);
            // This is probably unneccesary. Just check the first 2nd derivate and make sure it's positive (local minima)
            // Then, the sorted array is going to alternate between positive and negative (local min, local max)
            var inflixSecondDerivatives = inflixStationaryPoint.Select(zero => inflixSpline.Differentiate2(zero)).ToArray();

            var nistStationaryPoint = nistSpline.StationaryPoints();
            Array.Sort(nistStationaryPoint);
            var nistSecondDerivates = nistStationaryPoint.Select(zero => nistSpline.Differentiate2(zero)).ToArray();

            //var test = XicProcessing.GetRetentionTimePairs(inflixStationaryPoint, nistStationaryPoint);

            var placeholder = 0;

            


        }

        [Test]
        public static void TestModPeptideWithMultiplePeaks()
        {
            // In FPOP and PLIMB data, we observe cases where a modification on an aromatic residue can result in 
            // multiple chromatographic peaks belonging to the same ID ( think F modified at ortho, meta, and para positions)
            // get the raw file paths

            // This test is designed to simulate a situtation where the same peptides elutes twice as two well behaved peaks.
            // The first peak is MS/MS id'd in one file (inflix), and the second peak is MS/MS id'd in the second (nist)
            SpectraFileInfo inflix = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_Inflix_Tryp_60s_3-calib.mzML"),
                "inflix", 0, 0, 0);
            SpectraFileInfo nist = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_NIST_Tryp_60s_3-calib.mzML"),
                "nist", 1, 0, 0);

            // create IDs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            double monoisotopicMass = 2005.980136305;
            double peakFindingMass = 670.00;

            string baseSequence = "PWYEPIYLGGVFQLEK"; // 
            string fullSequence = "PW[CF3:CF3 on W]YEPIYLGGVFQLEK";
            double ms2RetentionTimeInflix = 32.464; // MS2 from first, largest peak
            double ms2RetentionTimeNist = 33.064; // MS2 from second, smaller peak

            // Intensities were found by examining the XICs at the RT flashLFQ defines as the peak apex
            double firstPeakInflixIntensity = 22710240; // MSMS Peak
            double secondPeakInflixIntensity = 771447; // No MSMS

            double firstPeakNistIntensity = 26988416; // No MSMS
            double secondPeakNistIntensity = 955824.06; // MSMS Peak

            Identification pepInflix = new Identification(inflix, baseSequence, fullSequence, monoisotopicMass,
                ms2RetentionTimeInflix, 3, new List<ProteinGroup> { pg });
            Identification pepNist = new Identification(nist, baseSequence, fullSequence, monoisotopicMass,
                ms2RetentionTimeNist, 3, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { pepInflix, pepNist },
                normalize: false, maxThreads: 1, matchBetweenRuns: false); // peaks are only serialized if match between runs = true

            // run the engine and grab XICs
            var results = engine.Run();

            Assert.That(results.PeptideModifiedSequences[fullSequence].GetIntensity(inflix), 
                Is.EqualTo(firstPeakInflixIntensity + secondPeakInflixIntensity).Within(1));

            Assert.That(results.PeptideModifiedSequences[fullSequence].GetIntensity(nist),
                Is.EqualTo(firstPeakNistIntensity + secondPeakNistIntensity).Within(1));
        }

        [Test]
        public static void TestPeakRecognitionForCoelutingIsobars()
        {
            // In FPOP and PLIMB data, we observe cases where a modification on an aromatic residue can result in 
            // multiple chromatographic peaks belonging to the same ID ( think F modified at ortho, meta, and para positions)
            // get the raw file paths

            // This test simulates a situation where three peaks elute very close together - close enough that FlashLFQ considers it 
            // one peak.
            // Peak 1 - One broad peak that is most likely caused by two co-eluting (only small rt shift) peptides,
            //      "PWYEPIY[CF3:CF3 on Y]LGGVFQLEK", score 16 and "PWY[CF3:CF3 on Y]EPIYLGGVFQLEK", score 6
            // Peak 2 - One broad peak w/o a good MS2 scan. Based on other runs, this peak most likely belongs to
            //      "PW[CF3:CF3 on W]YEPIYLGGVFQLEK". However, in this run, only one ambiguous MS2 scan was collected at the start
            //      of the peak elution. Ideally, this peak would be identified by MBR. Another test will be written for that case
            // Peak 3 - Undetermined, most likely two co-eluting species. Unimportant for the purposes of this test

            SpectraFileInfo nist = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_NIST_Tryp_60s_3-calib.mzML"),
                "nist", 1, 0, 0);

            // create IDs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            double monoisotopicMass = 2005.980136305;
            double peakFindingMass = 670.00;

            string baseSequence = "PWYEPIYLGGVFQLEK";
            string y7ModSeq = "PWYEPIY[CF3:CF3 on Y]LGGVFQLEK";
            string y3ModSeq = "PWY[CF3:CF3 on Y]EPIYLGGVFQLE";
            string ambiguousModSeq = "PW[CF3:CF3 on W]YEPIYLGGVFQLEK|PWY[CF3:CF3 on Y]EPIYLGGVFQLEK|PWYE[CF3:CF3 on E]PIYLGGVFQLEK"; // score of 5

            double y7ModRt = 33.22371; 
            double y3ModRt = 33.27658;
            double ambiguousModRt = 33.2862;

            // Intensities were found by examining the XICs and finding a local maxima
            double firstPeakIntensity = 533157;
            double secondPeakIntensity = 1778368;

            Identification y7ModId = new Identification(nist, baseSequence, y7ModSeq, monoisotopicMass,
                y7ModRt, 3, new List<ProteinGroup> { pg });
            Identification y3ModId = new Identification(nist, baseSequence, y3ModSeq, monoisotopicMass,
                y3ModRt, 3, new List<ProteinGroup> { pg });
            Identification ambiguousModId = new Identification(nist, baseSequence, ambiguousModSeq, monoisotopicMass,
                ambiguousModRt, 3, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { y7ModId, y3ModId, ambiguousModId },
                normalize: false, maxThreads: 1, matchBetweenRuns: false, quantifyAmbiguousPeptides: true); // peaks are only serialized if match between runs = true

            // run the engine and grab XICs
            var results = engine.Run();

            // Remeber that intensity is the sum of the isotopic envelope, so you have to do something else here.
            Assert.That(results.PeptideModifiedSequences.ContainsKey(y7ModSeq + "|" + y3ModSeq));
            Assert.That(results.PeptideModifiedSequences[y7ModSeq + "|" + y3ModSeq].GetIntensity(nist),
                Is.EqualTo(firstPeakIntensity).Within(1));

            Assert.That(results.PeptideModifiedSequences.ContainsKey("PW[CF3:CF3 on W]YEPIYLGGVFQLEK|PWYE[CF3:CF3 on E]PIYLGGVFQLEK"));
            // Questionable about whether or not we want to assign any intensity to this one. Probably not
            Assert.That(results.PeptideModifiedSequences["PW[CF3:CF3 on W]YEPIYLGGVFQLEK|PWYE[CF3:CF3 on E]PIYLGGVFQLEK"].GetIntensity(nist),
                Is.EqualTo(null));

            // Currently, y3 and the ambiguous are getting assigned to the same peak (the second one)
            // We need a better peak finding algorithm. I'm thinking splines

        }

        [Test]
        public static void TestXicSizeEqualization()
        {
            double[] smallArray = { 0, 1, 2, 3, 4, 5, 6 };
            double[] bigArray = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
            double[] expectedArrayPostTrim =  { 3, 4, 5, 6, 7, 8, 9, 10 };
            List<IndexedMassSpectralPeak> smallPeakList =
                smallArray.Select(d => new IndexedMassSpectralPeak(mz: d, 0, 0, 0)).ToList();
            List<IndexedMassSpectralPeak> bigPeakList =
                bigArray.Select(d => new IndexedMassSpectralPeak(mz: d, 0, 0, 0)).ToList();

            XicProcessing.EqualizeListLength(ref smallPeakList, ref bigPeakList);
            // Removing an even number of peaks
            Assert.AreEqual(bigPeakList.Select(p => p.Mz).ToArray(),
                expectedArrayPostTrim);

            bigArray = new double[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13 };
            bigPeakList =
                bigArray.Select(d => new IndexedMassSpectralPeak(mz: d, 0, 0, 0)).ToList();
            XicProcessing.EqualizeListLength(ref bigPeakList, ref smallPeakList);
            // Removing an odd number of peaks. EqualizeListSize should remove more from the end
            Assert.AreEqual(bigPeakList.Select(p => p.Mz).ToArray(),
                expectedArrayPostTrim);
        }
    }
}
