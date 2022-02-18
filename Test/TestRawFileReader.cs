using IO.MzML;
using MassSpectrometry;
using NUnit.Framework;
using System;
using System.Diagnostics;
using System.IO;
using IO.ThermoRawFileReader;
using System.Linq;
using System.Collections.Generic;
using Easy.Common.Extensions;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestRawFileReader
    {
        [Test]
        [TestCase("testFileWMS2.raw", "a.mzML", "aa.mzML")]
        [TestCase("small.raw", "a.mzML", "aa.mzML")]
        [TestCase("05-13-16_cali_MS_60K-res_MS.raw", "a.mzML", "aa.mzML")]
        /// <summary>
        /// Tests LoadAllStaticData for ThermoRawFileReader
        /// </summary>
        public static void TestLoadAllStaticDataRawFileReader(string infile, string outfile1, string outfile2)
        {
            var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", infile);
            outfile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", outfile1);
            outfile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", outfile2);

            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
            var a = ThermoRawFileReader.LoadAllStaticData(path, maxThreads: 1);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(a, outfile1, false);
            var aa = Mzml.LoadAllStaticData(outfile1);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(aa, outfile2, true);
            Mzml.LoadAllStaticData(outfile2);
            Console.WriteLine($"Analysis time for TestLoadAllStaticDataRawFileReader({infile}): {stopwatch.Elapsed.Hours}h {stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        /// <summary>
        /// Tests the dynamic connection for ThermoRawFileReader
        /// </summary>
        public static void TestDynamicConnectionRawFileReader()
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            var path1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "small.raw");
            var dynamicConnection1 = new ThermoDynamicData(path1);

            var path2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "testFileWMS2.raw");
            var dynamicConnection2 = new ThermoDynamicData(path2);

            var msOrders = dynamicConnection1.MsOrdersByScan;
            Assert.That(msOrders != null && msOrders.Length > 0);

            var a = dynamicConnection1.GetOneBasedScanFromDynamicConnection(1);
            Assert.That(a != null);

            var b = dynamicConnection2.GetOneBasedScanFromDynamicConnection(1);
            Assert.That(b != null);

            Assert.That(a.MassSpectrum.XArray.Length != b.MassSpectrum.XArray.Length);

            a = dynamicConnection1.GetOneBasedScanFromDynamicConnection(10000);
            Assert.That(a == null);

            dynamicConnection1.CloseDynamicConnection();
            dynamicConnection2.CloseDynamicConnection();

            Console.WriteLine($"Analysis time for TestDynamicConnectionRawFileReader: {stopwatch.Elapsed.Hours}h {stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        [TestCase("testFileWMS2.raw")]
        [TestCase("small.raw")]
        [TestCase("05-13-16_cali_MS_60K-res_MS.raw")]
        /// <summary>
        /// Tests peak filtering for ThermoRawFileReader
        /// </summary>
        public static void TestPeakFilteringRawFileReader(string infile)
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();
            var filterParams = new FilteringParams(200, 0.01, 0, 1, false, true, true);

            var path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", infile);

            var a = ThermoRawFileReader.LoadAllStaticData(path, filterParams, maxThreads: 1);
            var rawScans = a.GetAllScansList();
            foreach (var scan in rawScans)
            {
                Assert.That(scan.MassSpectrum.XArray.Length <= 200);
            }

            string outfile1 = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", Path.GetFileNameWithoutExtension(infile) + ".mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(a, outfile1, false);
            var mzml = Mzml.LoadAllStaticData(outfile1, filterParams, maxThreads: 1);

            var mzmlScans = mzml.GetAllScansList();
            for (int i = 0; i < mzmlScans.Count; i++)
            {
                var mzmlScan = mzmlScans[i];
                var rawScan = rawScans[i];

                for (int j = 0; j < mzmlScan.MassSpectrum.XArray.Length; j++)
                {
                    double roundedMzmlMz = Math.Round(mzmlScan.MassSpectrum.XArray[j], 2);
                    double roundedRawMz = Math.Round(rawScan.MassSpectrum.XArray[j], 2);

                    Assert.AreEqual(roundedMzmlMz, roundedRawMz);

                    double roundedMzmlIntensity = Math.Round(mzmlScan.MassSpectrum.XArray[j], 0);
                    double roundedRawIntensity = Math.Round(rawScan.MassSpectrum.XArray[j], 0);

                    Assert.AreEqual(roundedMzmlIntensity, roundedRawIntensity);
                }
            }

            Console.WriteLine($"Analysis time for TestPeakFilteringRawFileReader: {stopwatch.Elapsed.Hours}h " +
                $"{stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        /// <summary>
        /// Just makes sure the Thermo RawFileReader licence is accessible...
        /// </summary>
        public static void TestThermoLicence()
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            var licence = ThermoRawFileReaderLicence.ThermoLicenceText;
            Assert.That(licence.Length > 100);

            Console.WriteLine($"Analysis time for TestThermoLicence: {stopwatch.Elapsed.Hours}h {stopwatch.Elapsed.Minutes}m {stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        [TestCase("small.RAW")]
        [TestCase("testFileWMS2.raw")]
        [TestCase("05-13-16_cali_MS_60K-res_MS.raw")]
        public static void TestDynamicRaw(string fileName)
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", fileName);

            ThermoRawFileReader staticRaw = ThermoRawFileReader.LoadAllStaticData(filePath);
            ThermoDynamicData dynamicRaw = new ThermoDynamicData(filePath);

            foreach (MsDataScan staticScan in staticRaw.GetAllScansList())
            {
                MsDataScan dynamicScan = dynamicRaw.GetOneBasedScanFromDynamicConnection(staticScan.OneBasedScanNumber);

                Assert.IsFalse(staticScan.MassSpectrum.YArray.Contains(0));
                Assert.IsFalse(dynamicScan.MassSpectrum.YArray.Contains(0));
                Assert.That(dynamicScan.OneBasedScanNumber == staticScan.OneBasedScanNumber);
                Assert.That(dynamicScan.MsnOrder == staticScan.MsnOrder);
                Assert.That(dynamicScan.RetentionTime == staticScan.RetentionTime);
                Assert.That(dynamicScan.Polarity == staticScan.Polarity);
                Assert.That(dynamicScan.ScanWindowRange.Minimum == staticScan.ScanWindowRange.Minimum);
                Assert.That(dynamicScan.ScanWindowRange.Maximum == staticScan.ScanWindowRange.Maximum);
                Assert.That(dynamicScan.ScanFilter == staticScan.ScanFilter);
                Assert.That(dynamicScan.NativeId == staticScan.NativeId);
                Assert.That(dynamicScan.IsCentroid == staticScan.IsCentroid);
                Assert.That(dynamicScan.IsCentroid == staticScan.IsCentroid);
                Assert.That(dynamicScan.InjectionTime == staticScan.InjectionTime);
                Assert.That(dynamicScan.NoiseData == staticScan.NoiseData);

                Assert.That(dynamicScan.IsolationMz == staticScan.IsolationMz);
                Assert.That(dynamicScan.SelectedIonChargeStateGuess == staticScan.SelectedIonChargeStateGuess);
                Assert.That(dynamicScan.SelectedIonIntensity == staticScan.SelectedIonIntensity);
                Assert.That(dynamicScan.SelectedIonMZ == staticScan.SelectedIonMZ);
                Assert.That(dynamicScan.DissociationType == staticScan.DissociationType);
                Assert.That(dynamicScan.IsolationWidth == staticScan.IsolationWidth);
                Assert.That(dynamicScan.OneBasedPrecursorScanNumber == staticScan.OneBasedPrecursorScanNumber);
                Assert.That(dynamicScan.SelectedIonMonoisotopicGuessIntensity == staticScan.SelectedIonMonoisotopicGuessIntensity);
                Assert.That(dynamicScan.SelectedIonMonoisotopicGuessMz == staticScan.SelectedIonMonoisotopicGuessMz);

                if (dynamicScan.IsolationRange != null || staticScan.IsolationRange != null)
                {
                    Assert.That(dynamicScan.IsolationRange.Minimum == staticScan.IsolationRange.Minimum);
                    Assert.That(dynamicScan.IsolationRange.Maximum == staticScan.IsolationRange.Maximum);
                }

                Assert.That(dynamicScan.MassSpectrum.XArray.Length == staticScan.MassSpectrum.XArray.Length);
                Assert.That(dynamicScan.MassSpectrum.YArray.Length == staticScan.MassSpectrum.YArray.Length);

                for (int i = 0; i < staticScan.MassSpectrum.XArray.Length; i++)
                {
                    double staticMz = staticScan.MassSpectrum.XArray[i];
                    double staticIntensity = staticScan.MassSpectrum.YArray[i];

                    double dynamicMz = dynamicScan.MassSpectrum.XArray[i];
                    double dynamicIntensity = dynamicScan.MassSpectrum.YArray[i];

                    Assert.That(dynamicMz == staticMz);
                    Assert.That(dynamicIntensity == staticIntensity);
                }
            }
        }

        [Test]
        public static void TestEthcdReading()
        {
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "sliced_ethcd.raw");
            var spectra = ThermoRawFileReader.LoadAllStaticData(filePath, null, 1);
            var hcdScan = spectra.GetOneBasedScan(5);
            Assert.That(hcdScan.DissociationType == DissociationType.HCD);
            var ethcdScan = spectra.GetOneBasedScan(6);
            Assert.That(ethcdScan.DissociationType == DissociationType.EThcD);
        }

        [Test]
        public static void FragmentationPatterns()
        {
            //string filePath = @"C:\\Users\\Nic\\OneDrive\\Research\\Source Induced Decay\\SIDTest.raw"; //Big File
            //var scans = ThermoRawFileReader.LoadAllStaticData(filePath).GetAllScansList();
            string filePath = @"C:\\Users\\Nic\\OneDrive\\Research\\Source Induced Decay\\SIDTestSmall.mzML"; //Trimmed File
            var scans = Mzml.LoadAllStaticData(filePath).GetAllScansList();

            int minAssumedChargeState = 2;
            int maxAssumedChargeState = 60;
            int deconvolutionTolerancePpm = 4;
            int intensityRatio = 3;
            var preliminaryProcessing = new List<MatchedIons>();

            //creates a list of the relavant data from the Spectrum after deconvolution and normalizing to TIC
            foreach (var scan in scans)
            {
                var tempEnvelope = scan.MassSpectrum.Deconvolute(scan.MassSpectrum.Range, minAssumedChargeState,
                    maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatio).ToList();
                var test = new MatchedIons(scan.OneBasedScanNumber, tempEnvelope, scan.TotalIonCurrent);
                preliminaryProcessing.Add(test);
            }

            var sidScans = preliminaryProcessing.Where(p => p.sidBool == true).ToList();
            var norScans = preliminaryProcessing.Where(p => p.sidBool == false).ToList();
            var matchedIons = new List<MatchedIons>();

            //Itterates through each pair of scans and finds where the monoisotopic masses are within the ppm error defined below
            int ppmError = 25;
            for (int i = 0; i < norScans.Count(); i++)
            {
                //each IsotopicEnvelope in normal scan
                foreach (var isoEnvNorm in norScans[i].isotopicEnvelope)
                {
                    //finds most intense peak in the envelope of normal scan
                    double maxIntIonNorm = 0;
                    double maxIntNorm = 0;
                    foreach (var peak in isoEnvNorm.Peaks)
                    {
                        if (peak.Item2 > maxIntNorm)
                        {
                            maxIntNorm = peak.Item2;
                            maxIntIonNorm = peak.Item1;
                        }
                    }
                    //each IsotopicEnvelope in SID scan
                    foreach (var isoEnvSid in sidScans[i].isotopicEnvelope)
                    {
                        //finds most intense peak in the envelope of SID scan
                        double maxIntIonSid = 0;
                        double maxIntSid = 0;
                        foreach (var peak in isoEnvSid.Peaks)
                        {
                            if (peak.Item2 > maxIntSid)
                            {
                                maxIntSid = peak.Item2;
                                maxIntIonSid = peak.Item1;
                            }
                        }
                        //compares the most intense peaks of every envelope and adds it to list if they match both mass and charge
                        var normTolerance = isoEnvNorm.MonoisotopicMass / 10e6 * ppmError;
                        var sidTolerance = isoEnvSid.MonoisotopicMass / 10e6 * ppmError;
                        if ((Math.Round(isoEnvNorm.MonoisotopicMass, 4) >= Math.Round(isoEnvSid.MonoisotopicMass, 4) - sidTolerance && 
                             Math.Round(isoEnvNorm.MonoisotopicMass, 4) <= Math.Round(isoEnvSid.MonoisotopicMass, 4) + sidTolerance) || 
                            (Math.Round(isoEnvNorm.MonoisotopicMass, 4) + normTolerance >= Math.Round(isoEnvSid.MonoisotopicMass, 4) &&
                             Math.Round(isoEnvNorm.MonoisotopicMass, 4) - normTolerance <= Math.Round(isoEnvSid.MonoisotopicMass, 4)))
                        {
                            double fragEfficiency = (isoEnvNorm.TotalIntensity / norScans[i].totalIonCurrent) / (isoEnvSid.TotalIntensity / sidScans[i].totalIonCurrent);
                            int sidCharge = isoEnvSid.Charge;
                            int norCharge = isoEnvNorm.Charge;
                            if (sidCharge == norCharge)
                                
                                matchedIons.Add(new MatchedIons(isoEnvSid, isoEnvNorm, Math.Round(maxIntIonNorm, 4), 
                                    norCharge, maxIntNorm, maxIntSid, norScans[i].totalIonCurrent, sidScans[i].totalIonCurrent));
                        }

                    }
                }
            }
            var sortedIons = matchedIons.OrderBy(p => p.matchedIonMZ).ToList();
            var groupedIons = sortedIons.GroupBy(p => p.chargeState).Select(p => p.ToList()).ToList(); //unused as of yet, group by other factors for future data analysis. 

            string printout;
            string filename = @"C:\\Users\\Nic\\Desktop\\OuputFolder\\aSmolScan2.txt";
            using (StreamWriter sw = new StreamWriter(filename))
            {
                sw.WriteLine("Ion MZ : Charge : Monoisotopic Mass : Fragmentation Efficiency : Total Fragmentation Efficiency");
                foreach (var ion in sortedIons)
                {
                    printout = "";
                    printout = ion.matchedIonMZ + " : " + ion.chargeState + " : " + ion.normIsotopicEnvelope.MonoisotopicMass + " : " + ion.fragmentationEfficiencyOneIon + " : " + 
                        ion.fragmentationEfficiencyTotal;
                    sw.WriteLine(printout);
                }
            }

        }
        //Discarded this as the xArray peak and what the IsotopicEnvelopes differed by a bit, so things were getting missed
        public static int FindEnvelopeIndex(double mass, List<IsotopicEnvelope> envelopes)
        {
            int i = 0;
            for (i = 0; i < envelopes.Count(); i++)
            {
                foreach (var peak in envelopes[i].Peaks)
                {
                    if (Math.Round(mass, 3) == Math.Round(peak.mz, 3))
                        return i;
                }
            }
            return -1;
        }
    }
    
}
public class MatchedIons
{
    public int chargeState;
    public IsotopicEnvelope sidIsotopicEnvelope;
    public IsotopicEnvelope normIsotopicEnvelope;
    public List<IsotopicEnvelope> isotopicEnvelope;
    public double sidIntensity;
    public double normIntensity;
    public double matchedIonMZ;
    public double fragmentationEfficiencyTotal;
    public double fragmentationEfficiencyOneIon;
    public double totalIonCurrent;
    public double sidTotalIonCurrent;
    public double normTotalIonCurrent;
    public int scanNumber;
    public bool sidBool;
    

    public MatchedIons(IsotopicEnvelope sid, IsotopicEnvelope norm, double mz, int charge, double normInt, double sidInt, double normCurrent, double sidCurrent)
    {
        sidIsotopicEnvelope = sid;
        normIsotopicEnvelope = norm;
        matchedIonMZ = mz;
        chargeState = charge;
        normIntensity = normInt;
        sidIntensity = sidInt;
        sidTotalIonCurrent = sidCurrent;
        normTotalIonCurrent = normCurrent;
        fragmentationEfficiencyOneIon = (normInt / normCurrent) / (sidInt / sidCurrent);
        fragmentationEfficiencyTotal = (norm.TotalIntensity / normCurrent) / (sid.TotalIntensity / sidCurrent);
    }
    public MatchedIons(int scanNum, List<IsotopicEnvelope> envelope, double ionCurrent)
    {
        scanNumber = scanNum;
        isotopicEnvelope = envelope;
        totalIonCurrent = ionCurrent;
        int oddEven = scanNum % 2;
        if (oddEven == 0)
            sidBool = true;
        else if (oddEven != 0)
            sidBool = false;
    }
}


