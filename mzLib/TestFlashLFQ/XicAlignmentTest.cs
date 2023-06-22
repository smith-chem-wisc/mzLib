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
using MassSpectrometry.MzSpectra;
using Test.FileReadingTests;
using UsefulProteomicsDatabases;
using ChromatographicPeak = FlashLFQ.ChromatographicPeak;
using Stopwatch = System.Diagnostics.Stopwatch;

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

            // run the engine
            var results = engine.Run();
            var indexingEngine = engine.GetIndexingEngine();
            var inflixXic = indexingEngine.GetXIC(peakFindingMass, inflix);
            var nistXic = indexingEngine.GetXIC(peakFindingMass, nist);

            SpectralSimilarity xicSimilarity = new SpectralSimilarity(
                P_XArray: inflixXic.Select(p => p.RetentionTime).ToArray(),
                P_YArray: inflixXic.Select(p => p.Intensity).ToArray(),
                Q_XArray: nistXic.Select(p => p.RetentionTime).ToArray(),
                Q_YArray: nistXic.Select(p => p.Intensity).ToArray(),
                SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak,
                toleranceInPpm: 250, // Not sure what the ideal ppmTolerance is here
                allPeaks: true,
                filterOutBelowThisMz: 0);
            double? xicAngle = xicSimilarity.SpectralContrastAngle();
            int placeholder = 0;
        }
    }
}
