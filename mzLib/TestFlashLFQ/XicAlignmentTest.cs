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
using static Nett.TomlObjectFactory;
using System.Windows.Documents;
using System.Reflection;
using NUnit.Framework.Internal;

namespace TestFlashLFQ
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal static class XicAlignmentTest
    {
        [Test]
        public static void TestSignalAlignment()
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
            var inflixXic = indexingEngine.ExtractPeaks(peakFindingMass, inflix);
            var nistXic = indexingEngine.ExtractPeaks(peakFindingMass, nist);

            CollectionAssert.AreEqual(inflixXic, inflixXic.OrderBy(p => p.RetentionTime).ToList());

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
        public static void MatchExtremaTest()
        {
            var ex1 = new Extremum(1, 0, ExtremumType.Maximum);
            var ex2 = new Extremum(2, 0, ExtremumType.Maximum);
            var ex3 = new Extremum(3, 0, ExtremumType.Maximum);
            var ex4 = new Extremum(4, 0, ExtremumType.Maximum);
            var ex5 = new Extremum(5, 0, ExtremumType.Maximum);

            Extremum[] refArray = { ex1, ex2, ex3, ex4, ex5 };
            Extremum[] expArray = { ex1, ex3, ex5 };
            var pairs = XicProcessing.MatchExtrema(refArray, expArray);
            Assert.AreEqual(pairs, new Extremum[] { ex1, null, ex3, null, ex5 });

            var ex32 = new Extremum(3.2, 0, ExtremumType.Maximum);
            var ex47 = new Extremum(4.7, 0, ExtremumType.Maximum);
            var ex51 = new Extremum(5.1, 0, ExtremumType.Maximum);
            var ex53 = new Extremum(5.3, 0, ExtremumType.Maximum);
            var ex55 = new Extremum(5.5, 0, ExtremumType.Maximum);
            
            Extremum[] refArray2 = { ex1, ex3, ex51, ex53, ex55 };
            Extremum[] expArray2 = { ex1, ex1, ex3, ex32, ex4, ex47, ex5, ex5 };
            pairs = XicProcessing.MatchExtrema(refArray2, expArray2);
            Assert.AreEqual(pairs, new Extremum[]
                { ex1, ex3, ex5, null, null });

            var ex21 = new Extremum(2.1, 0, ExtremumType.Maximum);
            var ex22 = new Extremum(2.2, 0, ExtremumType.Maximum);
            var ex23 = new Extremum(2.3, 0, ExtremumType.Maximum);

            Extremum[] refArray3 = { ex1, ex2, ex3, ex4};
            Extremum[] expArray3 = { ex2, ex21, ex22, ex23 };
            pairs = XicProcessing.MatchExtrema(refArray3, expArray3);
            Assert.AreEqual(pairs, new Extremum[]
                { null, ex2, ex23, null });
        }

        [Test]
        public static void ReconcileExtremaTest()
        {
            // Test simple case, non-reference array all missing the same entries
            var ex1 = new Extremum(1, 0, ExtremumType.Maximum);
            var ex2 = new Extremum(2, 0, ExtremumType.Maximum);
            var ex3 = new Extremum(3, 0, ExtremumType.Maximum);
            var ex4 = new Extremum(4, 0, ExtremumType.Maximum);
            var ex5 = new Extremum(5, 0, ExtremumType.Maximum);

            Extremum[] refArray = { ex1, ex2, ex3, ex4, ex5 };
            Extremum[] expArray = { ex1, null, ex3, null, ex5 };
            //Extremum[] expArray2 = { ex1, null, ex3, null, ex5 };
            //Extremum[] expArray3 = { ex1, null, ex3, null, ex5 };
            var matrix = XicProcessing.ReconcileExtrema(
                refArray.ToList(),
                new List<List<Extremum>>
                {
                    expArray.ToList(),
                    expArray.ToList(),
                    expArray.ToList()
                });

            var firstRow = Enumerable.Range(0, 5)
                .Select(col => matrix[0, col])
                .ToArray();
            CollectionAssert.AreEqual(firstRow, expArray);

            // Only 2/3 non-reference arrays have missing entries. 2/4 (ref+non-ref) are defined at each index.
            // every row should have every entry defined in the matrix
            matrix = XicProcessing.ReconcileExtrema(
                refArray.ToList(),
                new List<List<Extremum>>
                {
                    refArray.ToList(),
                    expArray.ToList(),
                    expArray.ToList()
                });

            firstRow = Enumerable.Range(0, 5)
                .Select(col => matrix[0, col])
                .ToArray();
            CollectionAssert.AreEqual(firstRow, refArray);

            var fourthRow = Enumerable.Range(0, 5)
                .Select(col => matrix[3, col])
                .ToArray();
            CollectionAssert.AreNotEqual(fourthRow, expArray);

            var imputed2 = new Extremum(2, -1, ExtremumType.Maximum);
            var imputed4 = new Extremum(4, -1, ExtremumType.Maximum);
            CollectionAssert.AreEqual(fourthRow, new Extremum[] { ex1, imputed2, ex3, imputed4, ex5 });


            // 3/4 non-reference arrays have missing entries. 2/5 (ref+non-ref) are defined at each index.
            // All rows should be null for the second and fourth columns
            matrix = XicProcessing.ReconcileExtrema(
                refArray.ToList(),
                new List<List<Extremum>>
                {
                    refArray.ToList(),
                    expArray.ToList(),
                    expArray.ToList(),
                    expArray.ToList()
                });

            firstRow = Enumerable.Range(0, 5)
                .Select(col => matrix[0, col])
                .ToArray();
            CollectionAssert.AreEqual(firstRow, expArray);

            var secondRow = Enumerable.Range(0, 5)
                .Select(col => matrix[1, col])
                .ToArray();
            CollectionAssert.AreEqual(secondRow, expArray);
        }

        [Test]
        public static void TestVVRealData()
        {
            // Laptop

            //SpectraFileInfo nist2 = new SpectraFileInfo(
            //    @"C:\Users\asolivai\Documents\Immuto\JD020823_TNFa_Tryp_60s_AllSamples\Task1-CalibrateTask\JD020823_TNFa_NIST_Tryp_60s_2-calib.mzML",
            //    "nist2", 1, 0, 0);
            //SpectraFileInfo nist3 = new SpectraFileInfo(
            //    @"C:\Users\asolivai\Documents\Immuto\JD020823_TNFa_Tryp_60s_AllSamples\Task1-CalibrateTask\JD020823_TNFa_NIST_Tryp_60s_3-calib.mzML", "nist3", 1, 0, 0);
            //SpectraFileInfo nist5 = new SpectraFileInfo(
            //    @"C:\Users\asolivai\Documents\Immuto\JD020823_TNFa_Tryp_60s_AllSamples\Task1-CalibrateTask\\JD020823_TNFa_NIST_Tryp_60s_5-calib.mzML", "nist5", 1, 0, 0);

            // Alecto

            SpectraFileInfo nist2 = new SpectraFileInfo(
                @"C:\Users\Alex\Documents\Immuto\CalibratedData\JD020823_TNFa_NIST_Tryp_60s_2-calib.mzML", "nist2", 1, 0, 0);
            SpectraFileInfo nist3 = new SpectraFileInfo(
                @"C:\Users\Alex\Documents\Immuto\CalibratedData\JD020823_TNFa_NIST_Tryp_60s_3-calib.mzML", "nist3", 1, 0, 0);
            SpectraFileInfo nist5 = new SpectraFileInfo(
                @"C:\Users\Alex\Documents\Immuto\CalibratedData\JD020823_TNFa_NIST_Tryp_60s_5-calib.mzML", "nist5", 1, 0, 0);

            // create IDs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            string baseSequence = "VVSVLTVLHQDWLNGK";
            double monoisotopicMass = 1874.987;
            double peakFindingMass = 938.5006;

            // Define the modified peptide identifications
            string wModSeq = "VVSVLTVLHQDW[CF3:CF3 on W]LNGK";
            double rt2 = 30.71;
            double rt3 = 30.596;
            double rt5 = 31.50;

            // Y7 is observed in all three runs. It coelutes with y3.
            Identification vvNist2 = new Identification(nist2, baseSequence, wModSeq, monoisotopicMass,
                rt2, 3, new List<ProteinGroup> { pg });
            Identification vvNist3 = new Identification(nist3, baseSequence, wModSeq, monoisotopicMass,
                rt3, 3, new List<ProteinGroup> { pg });
            Identification vvNist5 = new Identification(nist5, baseSequence, wModSeq, monoisotopicMass,
                rt5, 2, new List<ProteinGroup> { pg });


            List<Identification> allIdentifications = new List<Identification>
            {
                vvNist2, vvNist3, vvNist5
            };

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(
                allIdentifications,
                normalize: false,
                maxThreads: 1,
                matchBetweenRuns: true,
                quantifyAmbiguousPeptides: true); // peaks are only serialized if match between runs = true

            var results = engine.Run();
            DoubleRange rtRange = new DoubleRange(30.0, 32.0);
            List<ChromatographicPeak> peaks = results.Peaks.SelectMany(kvp => kvp.Value).ToList();
            IsobarCluster testCluster = new IsobarCluster(peaks, engine, results, rtRange);

            int placeholder = 0;

            testCluster.ReassignPeakIDs();

            placeholder = 1;
        }

        [Test]
        public static void TestAmbiguousPeptideResults()
        {
            SpectraFileInfo inflix1 = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"TNFa_Inflix_Tryp_1-calib_ASQ.mzML"),
                "inflix", 1, 0, 0);
            SpectraFileInfo inflix3 = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"TNFa_Inflix_Tryp_3-calib_ASQ.mzML"),
                "inflix", 2, 0, 0);
            SpectraFileInfo inflix4 = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"TNFa_Inflix_Tryp_4-calib_ASQ.mzML"),
                "inflix", 3, 0, 0);

            // create IDs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            string baseSequence = "ASQFVGSSIHWYQQR";
            double monoisotopicMass = 1860.852;
            double peakFindingMass = 621.2912;

            // Define the modified peptide identifications
            string wModSeq = "ASQFVGSSIHW[CF3:CF3 on W]YQQR";
            double rtw1 = 23.48327;
            double rtw3 = 21.08657;
            double rtw4 = 21.22665;

            string ambigModSeq = "ASQFVGSSIHW[CF3:CF3 on W]YQQR|ASQFVGSSIHWYQQR[CF3:CF3 on R]";
            double rtambig = 21.20498;

            // Y7 is observed in all three runs. It coelutes with y3.
            Identification w1ModId = new Identification(inflix1, baseSequence, wModSeq, monoisotopicMass,
                rtw1, 3, new List<ProteinGroup> { pg });
            Identification w3ModId= new Identification(inflix3, baseSequence, wModSeq, monoisotopicMass,
                rtw3, 3, new List<ProteinGroup> { pg });
            Identification w4ModId = new Identification(inflix4, baseSequence, wModSeq, monoisotopicMass,
                rtw4, 3, new List<ProteinGroup> { pg });

            Identification ambigModId = new Identification(inflix1, baseSequence, ambigModSeq, monoisotopicMass,
                rtambig, 3, new List<ProteinGroup> { pg });

            List<Identification> allIdentifications = new List<Identification>
            {
                w1ModId, w3ModId, w4ModId, ambigModId
            };
            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(
                allIdentifications,
                normalize: false,
                maxThreads: 1,
                matchBetweenRuns: true,
                quantifyAmbiguousPeptides: false); // peaks are only serialized if match between runs = true

            var results = engine.Run();

            FlashLFQ.Peptide wModPeptide = results.PeptideModifiedSequences[wModSeq];
            List<double> wModIntensities = new();
            foreach(var file in results.SpectraFiles)
            {
                wModIntensities.Add(wModPeptide.GetIntensity(file));
            }

            // Make sure that the default behavior is picking different peaks
            double wModIntensityVariance = wModIntensities.StandardDeviation();
            Assert.Greater(wModIntensityVariance, 1e7);

            engine = new FlashLfqEngine(
                allIdentifications,
                normalize: false,
                maxThreads: 1,
                matchBetweenRuns: true,
                quantifyAmbiguousPeptides: true); // peaks are only serialized if match between runs = true

            results = engine.Run();

            wModPeptide = results.PeptideModifiedSequences[wModSeq];
            wModIntensities.Clear();
            foreach (var file in results.SpectraFiles)
            {
                wModIntensities.Add(wModPeptide.GetIntensity(file));
            }
            
            double wModIntensityVariance_Ambig = wModIntensities.StandardDeviation();
            Assert.Greater(wModIntensityVariance, wModIntensityVariance_Ambig);
            Assert.Less(wModIntensityVariance_Ambig, 1e7);

            // With QuantifyAmbiguous turned on, a new MBR peak is created by splitting the ambiguous ID
            // This peak has a higher intensity than the peak associated with the MSMS id
            // Now, I want to test the situation where the MSMS id has a higher intensity, but the retention time 
            // is shifted compared to the majority of peaks for that peptide.

            var x = results.Peaks[inflix1].First(peak => peak.Identifications.Contains(w1ModId));
            x.SetChromatographicPeakProperties(propName: "Intensity", newValue: (double)1e13);
            // IsobarClusters are used to adjust peptide quantification. We need to remove them to 
            // observed the default behaviour of CalculatePeptideResults
            var stashClusters = results.IsobarClusters;
            results.IsobarClusters = null;
            results.CalculatePeptideResults(quantifyAmbiguousPeptides: false);
            wModPeptide = results.PeptideModifiedSequences[wModSeq];
            wModIntensities.Clear();
            foreach (var file in results.SpectraFiles)
            {
                wModIntensities.Add(wModPeptide.GetIntensity(file));
            }
            double wModIntensityVarianceWithReflection = wModIntensities.StandardDeviation();
            // This tests that the MSMS id peak's intensity was changed, and as a result, 
            // that peak was used to define the peptide intensity
            Assert.Greater(wModIntensityVarianceWithReflection, wModIntensityVariance_Ambig);

            results.IsobarClusters = stashClusters;
            results.CalculatePeptideResults(quantifyAmbiguousPeptides: true);
            wModPeptide = results.PeptideModifiedSequences[wModSeq];
            wModIntensities.Clear();
            foreach (var file in results.SpectraFiles)
            {
                wModIntensities.Add(wModPeptide.GetIntensity(file));
            }
            double wModIntensityVariance_AmbigPeptideTest = wModIntensities.StandardDeviation();
            // This tests that the MSMS id peak's intensity was changed, and as a result, 
            // that peak was used to define the peptide intensity
            Assert.AreEqual(wModIntensityVariance_AmbigPeptideTest, wModIntensityVariance_Ambig, delta: 0.01);
        }

        public static void SetChromatographicPeakProperties(this ChromatographicPeak peak, string propName, Object newValue)
        {
            PropertyInfo propertyInfo = typeof(ChromatographicPeak).GetProperty(propName);
            if (propertyInfo == null || propertyInfo.PropertyType != newValue.GetType()) return;
            propertyInfo.SetValue(peak, newValue);
        }

        // In base FlashLFQ, some IDs get merged and then assigned to the wrong peaks.
        // Multi-run consensus should prevent this issue
        // This tests like 17 different things, but they all work!
        [Test]
        public static void TestReconcileExtremaRealData()
        {
            SpectraFileInfo nist = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_NIST_Tryp_60s_3-calib.mzML"),
                "nist", 1, 0, 0);
            SpectraFileInfo inflix2 = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_Inflix_Tryp_60s_2-calib.mzML"),
                "inflix", 0, 0, 0);
            SpectraFileInfo inflix3 = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_Inflix_Tryp_60s_3-calib.mzML"),
                "inflix", 0, 0, 0);
            
            // create IDs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            string baseSequence = "PWYEPIYLGGVFQLEK";
            double monoisotopicMass = 2005.980136305;
            double peakFindingMass = 670.00;

            // Define the modified peptide identifications
            string y7ModSeq = "PWYEPIY[CF3:CF3 on Y]LGGVFQLEK";
            double rty7 = 33.22371;

            // Y7 is observed in all three runs. It coelutes with y3.
            Identification y7ModId = new Identification(nist, baseSequence, y7ModSeq, monoisotopicMass,
                rty7, 3, new List<ProteinGroup> { pg });
            Identification y7ModIdInflix2 = new Identification(inflix2, baseSequence, y7ModSeq, monoisotopicMass,
                rty7, 3, new List<ProteinGroup> { pg });
            Identification y7ModIdInflix3 = new Identification(inflix3, baseSequence, y7ModSeq, monoisotopicMass,
                rty7, 3, new List<ProteinGroup> { pg });

            string y3ModSeq = "PWY[CF3:CF3 on Y]EPIYLGGVFQLEK";
            double rty3 = 33.22545;

            // The y3 modified peptide co-elutes with y7 in these runs (Actually, only in inflix 2 is that observed confidently.)
            Identification y3ModIdInflix2 = new Identification(inflix2, baseSequence, y3ModSeq, monoisotopicMass,
                rty3, 3, new List<ProteinGroup> { pg });
            Identification y3ModIdInflix3 = new Identification(inflix3, baseSequence, y3ModSeq, monoisotopicMass,
                rty3, 3, new List<ProteinGroup> { pg });

            // In NIST, the y3 mod peptide gets assigned to next peak. This peak has a FlashLFQ defined start time
            // that is after (greater than) the scan MS2 retention time. We want multi run consensus to see that the 
            // peptides coelutes in the inflix runs and adjust the peak assignment
            Identification y3ModId = new Identification(nist, baseSequence, y3ModSeq , monoisotopicMass,
                33.27658, 3, new List<ProteinGroup> { pg });


            List<Identification> allIdentifications = new List<Identification>
            {
                y3ModId, y3ModIdInflix2, y3ModIdInflix3,
                y7ModId, y7ModIdInflix2, y7ModIdInflix3
            };
            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(
                allIdentifications,
                normalize: false,
                maxThreads: 1,
                matchBetweenRuns: true,
                quantifyAmbiguousPeptides: false); // peaks are only serialized if match between runs = true

            // run the engine and grab XICs
            var results = engine.Run();
            // Check that the peptides are assigned to the same peak in the inflix conditions
            // and separate peaks in the nist condition
            Assert.AreEqual(2, results.Peaks[nist].Count);
            Assert.AreEqual(1, results.Peaks[inflix2].Count);
            Assert.AreEqual(1, results.Peaks[inflix3].Count);
           
            var indexingEngine = engine.GetIndexingEngine();
            var nistPeaks = indexingEngine.ExtractPeaks(peakFindingMass, nist);
            var inflixPeaks = indexingEngine.ExtractPeaks(peakFindingMass, inflix3);
            var inflix2Peaks = indexingEngine.ExtractPeaks(peakFindingMass, inflix2);
            double inflixAdjustment = XicProcessing.AlignPeaks(nistPeaks, inflixPeaks, 100);
            double inflix2Adjustment = XicProcessing.AlignPeaks(nistPeaks, inflix2Peaks, 100);

            Xic nistXic = new Xic(nistPeaks, peakFindingMass, nist, 0, referenceXic: true);
            Xic inflixXic = new Xic(inflixPeaks, peakFindingMass, inflix3, inflixAdjustment, referenceXic: false);
            Xic inflix2Xic = new Xic(inflix2Peaks, peakFindingMass, inflix2, inflix2Adjustment, referenceXic: false);

            var matchedExtrema = XicProcessing.ReconcileExtrema(nistXic.Extrema,
                new List<List<Extremum>> { inflix2Xic.Extrema, inflixXic.Extrema });

            Assert.AreEqual(3, matchedExtrema.GetLength(0));
            Assert.AreEqual(nistXic.Extrema.Count, matchedExtrema.GetLength(1));

            var firstRow = Enumerable.Range(0, nistXic.Extrema.Count)
                .Select(i => matchedExtrema[0, i])
                .ToArray();
            var secondRow = Enumerable.Range(0, nistXic.Extrema.Count)
                .Select(i => matchedExtrema[1, i])
                .ToArray();
            Assert.AreEqual(firstRow.Count(e => e == null), secondRow.Count(e => e == null));
            Assert.AreEqual(
                firstRow.Where(e => e != null).Count(e => e.ExtremumType == ExtremumType.Maximum),
                secondRow.Where(e => e != null).Count(e => e.ExtremumType == ExtremumType.Maximum));

            // Test Group Peaks
            //TODO: Test the ClusterPeaks function more rigorously
            var groupedPeaks = IsobarCluster.ClusterPeaks(results.Peaks.SelectMany(kvp => kvp.Value).ToList());

            Assert.AreEqual(1, groupedPeaks.Count);
            Assert.AreEqual(4, groupedPeaks.First().Count);

            // The peaks in results are all mutable, so we need to take a snapshot of the results before
            // reassignment is performed
            var nistPeaksBeforeReassignment = results.Peaks[nist]
                .OrderBy(peak => peak)
                .Select(peak => (peak.ApexRetentionTime, string.Join('|', peak.Identifications)))
                .ToList(); // Have to make sure enumeration runs here, before IDs are reassigned

            engine.IndexDict = new Dictionary<SpectraFileInfo, PeakIndexingEngine>();
            foreach (SpectraFileInfo file in new List<SpectraFileInfo>() { nist, inflix2, inflix3 })
            {
                engine.IndexDict.Add(file, new PeakIndexingEngine());
                engine.IndexDict[file].DeserializeIndex(file, deleteIndex: true);
            }
            // Test full IsobarCluster class
            var isobarCluster = IsobarCluster.FindIsobarClusters(allIdentifications, results, indexingEngine, engine).First();
            isobarCluster.ReassignPeakIDs();
            var nistPeaksAfterReassignment = results.Peaks[nist]
                .OrderBy(peak => peak)
                .Select(peak => (peak.ApexRetentionTime, string.Join('|', peak.Identifications)))
                .ToList();

            Assert.AreEqual(y7ModSeq,nistPeaksBeforeReassignment.First().Item2);
            Assert.AreEqual(y7ModSeq + "|" + y3ModSeq, nistPeaksAfterReassignment.First().Item2);

            // I'd like to sincerely apologize to my future self for the length of this test

            // We're going to remove one of the identifications from Nist and see if ReassignPeakIDs will succesfully
            // create a new chromatographic peak
            allIdentifications = new List<Identification>
            {
                y3ModId, y3ModIdInflix2, y3ModIdInflix3,
                y7ModIdInflix2, y7ModIdInflix3
            };
            // create the FlashLFQ engine
            engine = new FlashLfqEngine(
                allIdentifications,
                normalize: false,
                maxThreads: 1,
                matchBetweenRuns: true,
                quantifyAmbiguousPeptides: true);

            results = engine.Run();
            Assert.That(results.Peaks[nist].First(peak => peak.Identifications.Contains(y3ModId)).ApexRetentionTime, Is.LessThan(33.3));
            // Make sure that all chromatographic peaks with all IDs removed are filtered from the results
            Assert.AreEqual(0, results.Peaks[nist].Count(peak => peak.Identifications.Count == 0));

            // Checking for problems in reassigning the peak indexing engine and getting peak intensities
            // Make sure no peaks were accidentally duplicated/pulled from the wrong file
            Assert.AreNotEqual(
                results.Peaks[nist].First().Intensity,
                results.Peaks[inflix2].First().Intensity);
            Assert.AreNotEqual(
                results.Peaks[nist].First().Intensity,
                results.Peaks[inflix3].First().Intensity);
            Assert.AreNotEqual(
                results.Peaks[inflix2].First().Intensity,
                results.Peaks[inflix3].First().Intensity);

            var y3NistPeak = results.Peaks[nist]
                .First(peak => peak.Identifications.Contains(y3ModId));
            DoubleRange region = new DoubleRange(
                y3NistPeak.IsotopicEnvelopes.Min(e => e.IndexedPeak.RetentionTime),
                y3NistPeak.IsotopicEnvelopes.Max(e => e.IndexedPeak.RetentionTime));

            // This should be correct. This is essentially what IsobarCluster does when 
            // reassigning peaks, and this peak was reassigned
            engine.SwapPeakIndexingEngine(nist);
            var peakNist = engine.GetChromatographicPeak(y3ModId, nist, rtApex: rty3, peakRegion: region);
            Assert.AreEqual(peakNist.Intensity, results.Peaks[nist].First().Intensity, 1);

            // Assert that peaks found by the normal MS2 Quant workflow match up with those 
            // found by swapping the indexing engine and manually defining boundaries
            engine.SwapPeakIndexingEngine(inflix2);
            var peakInflix2 = engine.GetChromatographicPeak(y3ModId, nist, rtApex: rty3, peakRegion: region);
            Assert.AreEqual(peakInflix2.Intensity, results.Peaks[inflix2].First().Intensity, 1);

            engine.SwapPeakIndexingEngine(inflix3);
            var peakInflix3 = engine.GetChromatographicPeak(y3ModId, nist, rtApex: rty3, peakRegion: region);
            Assert.AreEqual(peakInflix3.Intensity, results.Peaks[inflix3].First().Intensity, 1);
        }

        [Test]
        public static void TNFaSliceMbrTest()
        {
            string psmFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", 
                "XICAlignment", @"UnambiguousPSMs_slices.psmtsv");

            SpectraFileInfo nist = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_NIST_Tryp_60s_3-calib.mzML"),
                "nist", 1, 0, 0);
            SpectraFileInfo inflix2 = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_Inflix_Tryp_60s_2-calib.mzML"),
                "inflix", 0, 0, 0);
            SpectraFileInfo inflix3 = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_Inflix_Tryp_60s_3-calib.mzML"),
                "inflix", 0, 0, 0);

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });

                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains("NIST"))
                {
                    file = nist;
                }
                else if (split[0].Contains("Inflix_Tryp_60s_2"))
                {
                    file = inflix2;
                }
                else if (split[0].Contains("Inflix_Tryp_60s_3"))
                {
                    file = inflix3;
                }

                string baseSequence = split[12];
                string fullSequence = split[13];
                double monoMass = double.Parse(split[22]);
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[25].Split(new char[] { '|' });
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

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData",
                "XICAlignment");

            // TODO: Report ambiguous but don't do multi-run consensus
            var engine = new FlashLfqEngine(ids, 
                matchBetweenRuns: true, 
                requireMsmsIdInCondition: false,
                quantifyAmbiguousPeptides: false, 
                maxThreads: 1);
            var results = engine.Run(out var exceptionList);
            Assert.IsEmpty(exceptionList);
            results.WriteResults(Path.Combine(outputFolder, "MultiRunConsensus_Peaks.tsv"), null, null, null, true);
            results.WriteResults(null, Path.Combine(outputFolder, "Default_Peptides.tsv"), null, null, true);

            var allPeaks = results.Peaks
                .SelectMany(kvp => kvp.Value)
                .ToList();

            // We only expect for the pwy3seq to be affected by the peak reassignment
            // The other two are valid isobar clusters, but they don't change
            string pwy3seq = "PWY[CF3:CF3 on Y]EPIYLGGVFQLEK";
            string pwy7seq = "PWYEPIY[CF3:CF3 on Y]LGGVFQLEK";
            string dnqyseq = "DNQLVVPSEGLY[CF3:CF3 on Y]LIYSQVLFK";
            Dictionary<string, List<double>> rtDictBefore = new Dictionary<string, List<double>>
            {
                {
                    pwy3seq,
                    allPeaks.Where(peak =>
                            peak.Identifications.Any(id => 
                                id.ModifiedSequence.Equals(pwy3seq)))
                        .Select(peak => peak.ApexRetentionTime)
                        .ToList()
                },
                {   pwy7seq,
                    allPeaks.Where(peak =>
                            peak.Identifications.Any(id => 
                                id.ModifiedSequence.Equals(pwy7seq)))
                        .Select(peak => peak.ApexRetentionTime)
                        .ToList()
                },
                {
                    dnqyseq,
                    allPeaks.Where(peak =>
                            peak.Identifications.Any(id => 
                                id.ModifiedSequence.Equals(dnqyseq)))
                        .Select(peak => peak.ApexRetentionTime)
                        .ToList()
                }
            };

            double averageRtVariance = rtDictBefore
                .Average(kvp =>
                {
                    var mean = kvp.Value.Average();
                    return kvp.Value.Sum(d => Math.Pow(d - mean, 2));
                });

            engine = new FlashLfqEngine(ids, 
                matchBetweenRuns: true, 
                requireMsmsIdInCondition: false,
                quantifyAmbiguousPeptides: true, 
                reportPeptideRetentionTimes: true,
                maxThreads: 1);
            results = engine.Run(out var exceptions);
            Assert.IsEmpty(exceptionList);

            results.WriteResults(Path.Combine(outputFolder, "MultiRunConsensus_Peaks.tsv"), null, null, null, true);
            results.WriteResults(null, Path.Combine(outputFolder, "MultiRunConsensus_Peptides.tsv"), null, null, true);

            allPeaks = results.Peaks
                .SelectMany(kvp => kvp.Value)
                .ToList();

            Dictionary<string, List<double>> rtDictAfter = new Dictionary<string, List<double>>
            {
                {
                    pwy3seq,
                    allPeaks.Where(peak =>
                            peak.Identifications.Any(id =>
                                id.ModifiedSequence.Equals(pwy3seq)))
                        .Select(peak => peak.ApexRetentionTime)
                        .ToList()
                },
                {   pwy7seq,
                    allPeaks.Where(peak =>
                            peak.Identifications
                                .Any(id => id.ModifiedSequence.Equals(pwy7seq)))
                        .Select(peak => peak.ApexRetentionTime)
                        .ToList()
                },
                {
                    dnqyseq,
                    allPeaks.Where(peak =>
                            peak.Identifications
                                .Any(id => id.ModifiedSequence.Equals(dnqyseq)))
                        .Select(peak => peak.ApexRetentionTime)
                        .ToList()
                }
            };

            // Get the average variance in retention times for all species
            double averageRtVarianceAfterAmbiguousQuant = rtDictAfter
                .Average(kvp =>
                {
                    var mean = kvp.Value.Average();
                    return kvp.Value.Sum(d => Math.Pow(d - mean, 2));
                });

            Assert.That( averageRtVarianceAfterAmbiguousQuant, Is.LessThan(averageRtVariance));
            CollectionAssert.AreEqual(rtDictBefore[dnqyseq], rtDictAfter[dnqyseq]);

            Assert.DoesNotThrow(() => results.CalculatePeptideResults(true));
        }

        [Test]
        public static void TestPeakRecognitionForCoelutingIsobars()
        {
            // In FPOP and PLIMB data, we observe cases where a modification on an aromatic residue can result in 
            // multiple chromatographic peaks belonging to the same ID ( think F modified at ortho, meta, and para positions)
            // get the raw file paths

            // This test simulates a situation where three peaks elute very close together

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
            var inflixXic = indexingEngine.ExtractPeaks(peakFindingMass, inflix);
            var nistXic = indexingEngine.ExtractPeaks(peakFindingMass, nist);

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
            SpectraFileInfo inflix2 = new SpectraFileInfo(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "XICAlignment", @"JD020823_TNFa_Inflix_Tryp_60s_2-calib.mzML"),
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
                normalize: false, maxThreads: 1, matchBetweenRuns: true); // peaks are only serialized if match between runs = true

            // run the engine and grab XICs
            var results = engine.Run();

            var indexingEngine = engine.GetIndexingEngine();
            var nistPeaks = indexingEngine.ExtractPeaks(peakFindingMass, nist);
            var inflixPeaks = indexingEngine.ExtractPeaks(peakFindingMass, inflix);
            double inflixAdjustment = XicProcessing.AlignPeaks(nistPeaks, inflixPeaks, 100);

            Xic nistXic = new Xic(nistPeaks, peakFindingMass, nist, 0, referenceXic: true);
            Xic inflixXic = new Xic(inflixPeaks, peakFindingMass, inflix, inflixAdjustment, referenceXic: false);

            var matchedExtrema = XicProcessing.ReconcileExtrema(nistXic.Extrema,
                new List<List<Extremum>> { inflixXic.Extrema });

            //var refType = matchedExtrema[0].Select(e => e.Type).ToArray();

            //Assert.AreEqual(matchedExtrema[0].Length, 16);
            //Assert.AreEqual(matchedExtrema[1].Length, 16);

            Assert.That(results.PeptideModifiedSequences[fullSequence].GetIntensity(inflix), 
                Is.EqualTo(firstPeakInflixIntensity + secondPeakInflixIntensity).Within(1));

            Assert.That(results.PeptideModifiedSequences[fullSequence].GetIntensity(nist),
                Is.EqualTo(firstPeakNistIntensity + secondPeakNistIntensity).Within(1));
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
