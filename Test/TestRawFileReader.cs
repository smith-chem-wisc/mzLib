using Easy.Common.Extensions;
using IO.MzML;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;

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
        public static void FragmentationPattersByXArray()
        {
            
            string filepath = @"C:\Users\Nic\OneDrive\Research\ForAustin\MS1-sidMS2_6ProtMix_R120_SID60_AGC800_20220123110514.raw"; // 6 protein standard
            //string filepath = @"c:\\users\\nic\\onedrive\\research\\source induced decay\\sidtest.raw"; //big SID file
            //ring filepath = @"C:\\Users\\Nic\\OneDrive\\Research\\Source Induced Decay\\SIDTestSmall.mzML"; //Trimmed SID File (only 2 scans)

            List<MsDataScan> scans = new();
            string scanType = filepath.Split('.')[1].Trim();
            if (scanType.Equals("mzML"))
              scans = Mzml.LoadAllStaticData(filepath).GetAllScansList();
            if (scanType.Equals("raw"))
                scans = ThermoRawFileReader.LoadAllStaticData(filepath).GetAllScansList();

            int minAssumedChargeState = 2;
            int maxAssumedChargeState = 60;
            int deconvolutionTolerancePpm = 4;
            int intensityRatio = 3;
            int roundingFactor = 2;
            int ppmError = 40;
            int relativeIntensityFilteringFactor = 100;
            bool printGrouped = false;
            var preliminaryProcessing = new List<MatchedIons>();

            //creates a list of the relavant data from the Spectrum after deconvolution and normalizing to TIC
            foreach (var scan in scans)
            {
                var tempEnvelope = scan.MassSpectrum.Deconvolute(scan.MassSpectrum.Range, minAssumedChargeState,
                    maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatio).ToList();
                scan.MassSpectrum.NormalizeToTIC(scan.TotalIonCurrent);
                scan.MassSpectrum.RoundXArray(roundingFactor);
                var test = new MatchedIons(scan.OneBasedScanNumber, tempEnvelope, scan.TotalIonCurrent, scan.MassSpectrum.XArray, scan.MassSpectrum.YArray);
                preliminaryProcessing.Add(test);
            }
            var sidScans = preliminaryProcessing.Where(p => p.sidBool == true).ToList();
            var normalScans = preliminaryProcessing.Where(p => p.sidBool == false).ToList();

            var matchedIons = new List<MatchedIons>();
            var matchedNormalIonsToEnvelopes = new List<MatchedIons>();
            var matchedSidIonsToEnvelopes = new List<MatchedIons>();
            //Goes through each pair of scans and generates a list of intersecting peaks
            for (int i = 0; i < normalScans.Count(); i++)
            {
                matchedNormalIonsToEnvelopes.Clear();
                matchedSidIonsToEnvelopes.Clear();
                List<double> intersect = null;
                List<double> intersectToRemove = new();
                intersect = normalScans[i].Xarray.Intersect(sidScans[i].Xarray).ToList();
                double normMaxIntensity = normalScans[i].Yarray.Max();
                double sidMaxIntensity = sidScans[i].Yarray.Max();
                foreach (var mass in intersect)
                {
                    if (normalScans[i].Yarray[Array.IndexOf(normalScans[i].Xarray, mass)] <= normMaxIntensity / relativeIntensityFilteringFactor)
                        intersectToRemove.Add(mass);
                    if (sidScans[i].Yarray[Array.IndexOf(sidScans[i].Xarray, mass)] <= sidMaxIntensity / relativeIntensityFilteringFactor)
                        intersectToRemove.Add(mass);
                }
                foreach (var mass in intersectToRemove)
                    if (intersect.Any(m => m.Equals(mass)))
                        intersect.Remove(mass);

                //Finds where Intersect masses can be found within normal scan isotopic envelopes
                for (int j = 0; j < normalScans[i].isotopicEnvelope.Count(); j++)
                {
                    //goes through each peak in the envelope
                    foreach (var peak in normalScans[i].isotopicEnvelope[j].Peaks)
                    {
                        var peakTolerance = peak.Item1 / 10e6 * ppmError;
                        //goes through each mass in intersect
                        foreach (var mass in intersect)
                        {
                            var massTolerance = mass / 10e6 * ppmError;
                            if ((Math.Round(peak.Item1, roundingFactor) <= mass + massTolerance &&
                                 Math.Round(peak.Item1, roundingFactor) >= mass - massTolerance) ||
                                (Math.Round(peak.Item1, roundingFactor) + peakTolerance >= mass &&
                                 Math.Round(peak.Item1, roundingFactor) - peakTolerance <= mass))
                            {
                                matchedNormalIonsToEnvelopes.Add(new MatchedIons(normalScans[i].sidBool, normalScans[i].isotopicEnvelope[j], mass, normalScans[i].totalIonCurrent, normalScans[i].Yarray[Array.IndexOf(normalScans[i].Xarray, mass)], i, normalScans[i].isotopicEnvelope[j].MonoisotopicMass));
                            }
                        }
                    }
                }
                //Finds where Intersect masses can be found within the SID scan isotopic envelopes
                for (int j = 0; j < sidScans[i].isotopicEnvelope.Count(); j++)
                {
                    //goes through each peak in the envelope
                    foreach (var peak in sidScans[i].isotopicEnvelope[j].Peaks)
                    {
                        var peakTolerance = peak.Item1 / 10e6 * ppmError;
                        //goes through each mass in intersect
                        foreach (var mass in intersect)
                        {
                            var massTolerance = mass / 10e6 * ppmError;
                            if ((Math.Round(peak.Item1, roundingFactor) <= mass + massTolerance &&
                                 Math.Round(peak.Item1, roundingFactor) >= mass - massTolerance) ||
                                (Math.Round(peak.Item1, roundingFactor) + peakTolerance >= mass &&
                                 Math.Round(peak.Item1, roundingFactor) - peakTolerance <= mass))
                            {
                                matchedSidIonsToEnvelopes.Add(new MatchedIons(sidScans[i].sidBool, sidScans[i].isotopicEnvelope[j], mass, sidScans[i].totalIonCurrent, sidScans[i].Yarray[Array.IndexOf(sidScans[i].Xarray, mass)], i, sidScans[i].isotopicEnvelope[j].MonoisotopicMass));
                                //sidTestRemaining.Remove(mass);
                            }
                        }
                    }
                }

                matchedSidIonsToEnvelopes = matchedSidIonsToEnvelopes.OrderBy(p => p.matchedIonMZ).ToList();
                matchedNormalIonsToEnvelopes = matchedNormalIonsToEnvelopes.OrderBy(p => p.matchedIonMZ).ToList();
                //Adds ions to MatchedIons List where the matched ion m/z, matched charges, and matched monoisotopic mass are equivalent between two scan types
                for (int k = 0; k < matchedNormalIonsToEnvelopes.Count(); k++)
                {
                    var normTolerance = matchedNormalIonsToEnvelopes[k].matchedIonMZ / 10e6 * ppmError;
                    foreach (var sidIon in matchedSidIonsToEnvelopes)
                    {
                        var sidTolerance = sidIon.matchedIonMZ / 10e6 * ppmError;
                        if ((matchedNormalIonsToEnvelopes[k].matchedIonMZ <= sidIon.matchedIonMZ + sidTolerance &&
                             matchedNormalIonsToEnvelopes[k].matchedIonMZ >= sidIon.matchedIonMZ - sidTolerance) ||
                            (matchedNormalIonsToEnvelopes[k].matchedIonMZ + normTolerance >= sidIon.matchedIonMZ &&
                             matchedNormalIonsToEnvelopes[k].matchedIonMZ - normTolerance <= sidIon.matchedIonMZ))
                        {
                            if (matchedNormalIonsToEnvelopes[k].normIsotopicEnvelope.Charge == sidIon.sidIsotopicEnvelope.Charge &&
                                Math.Round(matchedNormalIonsToEnvelopes[k].monoistopicMass, 0) == Math.Round(sidIon.monoistopicMass, 0))
                            {
                                matchedIons.Add(new MatchedIons(sidIon.sidIsotopicEnvelope, matchedNormalIonsToEnvelopes[k].normIsotopicEnvelope, sidIon.matchedIonMZ,
                                    sidIon.sidIsotopicEnvelope.Charge, matchedNormalIonsToEnvelopes[k].normIntensity, sidIon.sidIntensity, matchedNormalIonsToEnvelopes[k].normTotalIonCurrent,
                                    sidIon.sidTotalIonCurrent, true, sidIon.monoistopicMass));
                            }
                        }
                    }
                }
            }

            matchedIons = matchedIons.OrderBy(p => p.chargeState).ThenBy(p => p.normIsotopicEnvelope.MonoisotopicMass).ThenBy(p => p.matchedIonMZ).ToList();
            var ionsGroupedbyMonoIsotopicMass = matchedIons.GroupBy(p => p.roundedMonoisotopicMass).ToList();
            var ionsGroupedbyMonoGreaterThanSix = ionsGroupedbyMonoIsotopicMass.Where(p => p.Count() > 6).ToList();
            var ionsGroupedbyMonoGreaterThanTen = ionsGroupedbyMonoGreaterThanSix.Where(p => p.Count() > 10).ToList();
            var ionsGroupedbyMonoGreaterThanFifteen = ionsGroupedbyMonoGreaterThanSix.Where(p => p.Count() > 15).ToList();

            double[] monoIsotopicMasses = new double[6];
            string fastaLocation = @"C:\Users\Nic\OneDrive\Research\Source Induced Decay\SixProteinStandard\Six_Protein_Standard.fasta";
            string how = "mono";
            matchedIons = MatchedIons.MatchToSequence(fastaLocation, matchedIons, how);
            string printout;
            string filename;
            if (printGrouped)
            {
                filename = @"C:\\Users\\Nic\\Desktop\\OuputFolder\\AllFilesRound2Filter100GroupedSix.txt";
                using (StreamWriter sw = new StreamWriter(filename))
                {
                    sw.WriteLine("Ion MZ : Charge : Monoisotopic Mass : Rounded MonoIsotopic Mass : Fragmentation Efficiency : Total Fragmentation Efficiency");
                    foreach (var group in ionsGroupedbyMonoGreaterThanSix)
                    {
                        foreach (var ion in group)
                        {
                            printout = "";
                            printout = ion.matchedIonMZ + " : " + ion.chargeState + " : " + ion.monoistopicMass + " : " + ion.roundedMonoisotopicMass
                                + " : " + ion.fragmentationEfficiencyOneIon + " : " + ion.fragmentationEfficiencyTotal;
                            sw.WriteLine(printout);
                        }
                    }
                }
                filename = @"C:\\Users\\Nic\\Desktop\\OuputFolder\\AllFilesRound2Filter100GroupedTen.txt";
                using (StreamWriter sw = new StreamWriter(filename))
                {
                    sw.WriteLine("Ion MZ : Charge : Monoisotopic Mass : Rounded MonoIsotopic Mass : Fragmentation Efficiency : Total Fragmentation Efficiency");
                    foreach (var group in ionsGroupedbyMonoGreaterThanTen)
                    {
                        foreach (var ion in group)
                        {
                            printout = "";
                            printout = ion.matchedIonMZ + " : " + ion.chargeState + " : " + ion.monoistopicMass + " : " + ion.roundedMonoisotopicMass
                                + " : " + ion.fragmentationEfficiencyOneIon + " : " + ion.fragmentationEfficiencyTotal;
                            sw.WriteLine(printout);
                        }
                    }
                }
                filename = @"C:\\Users\\Nic\\Desktop\\OuputFolder\\AllFilesRound2Filter100GroupedFifteen.txt";
                using (StreamWriter sw = new StreamWriter(filename))
                {
                    sw.WriteLine("Ion MZ : Charge : Monoisotopic Mass : Rounded MonoIsotopic Mass : Fragmentation Efficiency : Total Fragmentation Efficiency");
                    foreach (var group in ionsGroupedbyMonoGreaterThanFifteen)
                    {
                        foreach (var ion in group)
                        {
                            printout = "";
                            printout = ion.matchedIonMZ + " : " + ion.chargeState + " : " + ion.monoistopicMass + " : " + ion.roundedMonoisotopicMass
                                + " : " + ion.fragmentationEfficiencyOneIon + " : " + ion.fragmentationEfficiencyTotal;
                            sw.WriteLine(printout);
                        }
                    }
                }
            }

            filename = @"C:\\Users\\Nic\\Desktop\\OuputFolder\\SixProteinStandardMatchToSequenceByMonoMass.txt";
            using (StreamWriter sw = new StreamWriter(filename))
            {
                sw.WriteLine("Ion MZ : Charge : Monoisotopic Mass : Rounded MonoIsotopic Mass : Fragmentation Efficiency : Total Fragmentation Efficiency");
                foreach (var ion in matchedIons)
                {
                    printout = "";
                    printout = ion.matchedIonMZ + " : " + ion.chargeState + " : " + ion.monoistopicMass + " : " + Math.Round(ion.roundedMonoisotopicMass,0)
                        + " : " + ion.fragmentationEfficiencyOneIon + " : " + ion.fragmentationEfficiencyTotal;
                    sw.WriteLine(printout);
                }
            }
        }
        [Test]
        public static void RunMatchSequence()
        {
            string filename = "file";
            List<MatchedIons> ions = new List<MatchedIons>();
            string how = "mz";
            MatchedIons.MatchToSequence(filename, ions, how);
            
        }

    }
}


public class MatchedIons
{
    public int chargeState { get; set; }
    public IsotopicEnvelope sidIsotopicEnvelope;
    public IsotopicEnvelope normIsotopicEnvelope;
    public List<IsotopicEnvelope> isotopicEnvelope;
    public double sidIntensity { get; set; }
    public double normIntensity { get; set; }
    public double matchedIonMZ;
    public double fragmentationEfficiencyTotal;
    public double fragmentationEfficiencyOneIon;
    public double totalIonCurrent;
    public double sidTotalIonCurrent;
    public double normTotalIonCurrent;
    public int scanNumber;
    public bool sidBool;
    public double[] Xarray;
    public double[] Yarray;
    public int sidIndex;
    public int normIndex;
    public double monoistopicMass;
    public double roundedMonoisotopicMass;

    public MatchedIons(IsotopicEnvelope sid, IsotopicEnvelope norm, double mz, int charge, double normInt, double sidInt, double normCurrent, double sidCurrent, bool normalized, double monoMass)
    {
        sidIsotopicEnvelope = sid;
        normIsotopicEnvelope = norm;
        matchedIonMZ = mz;
        chargeState = charge;
        normIntensity = normInt;
        sidIntensity = sidInt;
        sidTotalIonCurrent = sidCurrent;
        normTotalIonCurrent = normCurrent;
        monoistopicMass = monoMass;
        roundedMonoisotopicMass = Math.Round(monoMass, 2);
        if (!normalized)
        {
            fragmentationEfficiencyOneIon =  (sidInt / sidCurrent) / (normInt / normCurrent);
            fragmentationEfficiencyTotal =  (sid.TotalIntensity / sidCurrent) / (norm.TotalIntensity / normCurrent);
        }
        if (normalized)
        {
            fragmentationEfficiencyOneIon = sidInt / normInt;
            fragmentationEfficiencyTotal =  sid.TotalIntensity / norm.TotalIntensity;
        }
    }

    public MatchedIons(int scanNum, List<IsotopicEnvelope> envelope, double ionCurrent, double[] Xarr = null, double[] Yarr = null)
    {
        Xarray = Xarr;
        Yarray = Yarr;
        scanNumber = scanNum;
        isotopicEnvelope = envelope;
        totalIonCurrent = ionCurrent;
        int oddEven = scanNum % 2;
        if (oddEven == 0)
            sidBool = false;
        else if (oddEven != 0)
            sidBool = true;
    }

    public MatchedIons(bool sidB, IsotopicEnvelope envelope, double mz, double ionCurrent, double intensity, int index, double monomass)
    {
        sidBool = sidB;
        matchedIonMZ = mz;
        monoistopicMass = monomass;
        roundedMonoisotopicMass = Math.Round(monomass, 2);
        if (sidBool)
        {
            sidIsotopicEnvelope = envelope;
            sidTotalIonCurrent = ionCurrent;
            sidIntensity = intensity;
            sidIndex = index;
        }
        if (!sidBool)
        {
            normIsotopicEnvelope = envelope;
            normTotalIonCurrent = ionCurrent;
            normIntensity = intensity;
            normIndex = index;
        }
    }

    // Takes list of matched ions and trims it to only include monoisotopic masses from a fasta file to target proteins of interest
    public static List<MatchedIons> MatchToSequence(string filelocation, List<MatchedIons> matches, string howToProcess)
    {
        //var targets = ProteinDbLoader.LoadProteinFasta(filelocation, true, DecoyType.None, false, out var errors); This is broken with this fasta file
        // Resorting to primative methods

        //creates a deep copy of matches
        var matches2 = new List<MatchedIons>();
        foreach (var match in matches)
        {
            matches2.Add(match);
        }

        string[] sequences = new string[]
        { "MFPAMPLSSLFVNGPRTLCGAELVDALQFVCGDRGFYFNKPTGYGSSSRRAPQTGIVDECCFRSCDLRRLEMYCAPLKPAKSA",
          "TTFNIQDGPDFQDRVVNSETPVVVDFHAQWCGPCKILGPRLEKMVAKQHGKVVMAKVDIDDHTDLAIEYEVSAVPTVLAMKNGDVVDKFVGIKDEDQLEAFLKKLIG",
          "MDPYPLPKTDTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTEKPEVIDASELTPAVTTYKLVINGKTLKGETTTKAVDAETAEKAFKQYANDNGVDGVWTYDDATKTFTVTEMVTEVPGDAPTEPEKPEASIPLVPLTPATPIAKDDAKKDDTKKEDAKKPEAKKDDAKKAETAG",
          "SHHWGYGKHNGPEHWHKDFPIANGERQSPVDIDTKAVVQDPALKPLALVYGEATSRRMVNNGHSFNVEYDDSQDKAVLKDGPLTGTYRLVQFHFHWGSSDDQGSEHTVDRKKYAAELHLVHWNTKYGDFGTAAQQPDGLAVVGVFLKVGDANPALQKVLDALDSIKTKGKSTDFPNFDPGSLLPNVLNYWTYPGSLTTPPLLESVTWIVLKEPISVSSQQMLKFRTLNFNAEGEPELLMLANWRPAQPLKNRQVRGFPK",
          "AQHDEAQQNAFYQVLNMPNLNADQRNGFIQSLKDDPSQSANVLGEAQKLNDSQAPKADAQQNNFNKDQQSAFYEILNMPNLNEAQRNGFIQSLKDDPSQSTNVLGEAKKLNESQAPKADNNFNKEQQNAFYEILNMPNLNEEQRNGFIQSLKDDPSQSANLLSEAKKLNESQAPKADNKFNKEQQNAFYEILHLPNLNEEQRNGFIQSLKDDPSQSANLLAEAKKLNDAQAPKADNKFNKEQQNAFYEILHLPNLTEEQRNGFIQSLKDDPSVSKEILAEAKKLNDAQAPKEEDNNKPIEGRNSRGSVDASELTPAVTTYKLVINGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTEKPEVIDASELTPAVTTYKLVINGKTLKGETTTKAVDAETAEKAFKQYANDNGVDGVWTYDDATKTFTVTEMVTEVPLESTA",
          "MISYDNYVTILDEETLKAWIAKLEKAPVFAFATATDSLDNISANLVGLSFAIEPGVAAYIPVAHDYLDAPDQISRERALELLKPLLEDEKALKVGQNLKYDRGILANYGIELRGIAFDTMLESYILNSVAGRHDMDSLAERWLKHKTITFEEIAGKGKNQLTFNQIALEEAGRYAAEDADVTLQLHLKMWPDLQKHKGPLNVFENIEMPLVPVLSRIERNGVKIDPKVLHNHSEELTLRLAELEKKAHEIAGEEFNLSSTKQLQTILFEKQGIKPLKKTPGGAPSTSEEVLEELALDYPLPKVILEYRGLAKLKSTYTDKLPLMINPKTGRVHTSYHQAVTATGRLSSTDPNLQNIPVRNEEGRRIRQAFIAPEDYVIVSADYSQIELRIMAHLSRDKGLLTAFAEGKDIHRATAAEVFGLPLETVTSEQRRSAKAINFGLIYGMSAFGLARQLNIPRKEAQKYMDLYFERYPGVLEYMERTRAQAKEQGYVETLDGRRLYLPDIKSSNGARRAAAERAAINAPMQGTAADIIKRAMIAVDAWLQAEQPRVRMIMQVHDELVFEVHKDDVDAVAKQIHQLMENCTRLDVPLLVEVGSGENWDQAH"
        };

        // Adds each sequences to a peptide object
        List<Peptide> targets = new();
        for (int i = 0; i < sequences.Length; i++)
        {
            targets.Add(new Peptide(sequences[i]));
        }

        if (howToProcess.Equals("mz"))
        {
            // Generates list of possible m / z values from 2 - 40 charge state
            int[] charges = Enumerable.Range(2, 40).ToArray();
            List<double> massesToCheck = new();
            foreach (var peptide in targets)
            {
                foreach (var charge in charges)
                {
                    massesToCheck.Add(Math.Round(peptide.MonoisotopicMass / charge, 0));
                }
            }

            // removes any matches that are not a part of the above defined sequences
            foreach (var match in matches2)
            {
                if (!massesToCheck.Any(p => Math.Round(p, 0) == Math.Round(match.matchedIonMZ, 0)))
                    matches.Remove(match);
            }
        }

        if (howToProcess.Equals("mono"))
        {
            // removes any matches that are not a part of the above defined sequences
            foreach (var match in matches2)
            {
                bool matchFound = targets.Any(p => Math.Round(p.MonoisotopicMass, 0) >= Math.Round(match.monoistopicMass, 0) - 10 &&
                                                   Math.Round(p.MonoisotopicMass, 0) <= Math.Round(match.monoistopicMass, 0) + 10);


                if (!targets.Any(p => Math.Round(p.MonoisotopicMass, 0) == Math.Round(match.monoistopicMass, 0)))
                    matches.Remove(match);
            }
        }

        int breakpoint = 0;
        return matches.OrderBy(p => p.monoistopicMass).ThenBy(p => p.chargeState).ToList();
    }


}