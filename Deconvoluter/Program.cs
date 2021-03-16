using Chemistry;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace Deconvoluter
{
    public class Program
    {
        static void Main(string[] args)
        {
            var filepath = @"C:\Data\LVS_TD_Yeast\05-26-17_B7A_yeast_td_fract7_rep1.raw";
            var dir = Path.GetDirectoryName(filepath);
            var data = ThermoRawFileReader.LoadAllStaticData(filepath);
            double minMass = 500;

            List<string> output = new List<string> { "Scan\tMonoisotopic Mass\tPeaks m/z List\tPeaks Intensity List\t" +
                "Charge\tMS Order\tCorrelation to Averagine\tFraction Intensity Observed\tS/N\tLog2 Total Intensity\tNoise" };

            Random r = new Random(1);
            var scans = data.GetAllScansList().OrderBy(p => r.Next())
                //.Where(p => p.OneBasedScanNumber == 2466)
                .ToList();

            var envsList = new List<(int, MassSpectrometry.IsotopicEnvelope)>();

            Parallel.ForEach(Partitioner.Create(0, scans.Count),
               new ParallelOptions { MaxDegreeOfParallelism = -1 },
               (range, loopState) =>
               {
                   for (int i = range.Item1; i < range.Item2; i++)
                   {
                       var scan = scans[i];

                       var envs = scan.MassSpectrum.Deconvolute(scan.ScanWindowRange, 1, 60, 10, 3, 2, 0.4, 0.6, 1.5)
                       .Where(p => p.MonoisotopicMass > minMass).OrderBy(p => p.Peaks.First().mz);
                       //.Where(p => p.Charge >= 3)
                       //.ToList();

                       foreach (var env in envs)
                       {
                           var line = scan.OneBasedScanNumber + "\t" +
                               env.MonoisotopicMass + "\t" +
                               string.Join(";", env.Peaks.Select(p => p.mz.ToString("F3"))) + "\t" +
                               string.Join(";", env.Peaks.Select(p => p.intensity.ToString("F1"))) + "\t" +
                               env.Charge + "\t" +
                               scan.MsnOrder + "\t" +
                               env.PearsonCorrelation + "\t" +
                               env.FracIntensityObserved + "\t" +
                               env.SN + "\t" +
                               Math.Log(env.TotalIntensity, 2) + "\t" +
                               env.Noise;

                           lock (output)
                           {
                               output.Add(line);
                           }

                           lock (envsList)
                           {
                               envsList.Add((scan.OneBasedScanNumber, env));
                           }
                       }
                   }
               });

            File.WriteAllLines(Path.Combine(dir, Path.GetFileNameWithoutExtension(filepath) + "_decon.tsv"), output);

            Dictionary<int, HashSet<double>> ScanNumToClaimedMzs = new Dictionary<int, HashSet<double>>();

            output = new List<string> { "Scan\tMonoisotopic Mass\tPeaks m/z List\tPeaks Intensity List\t" +
                "Charge\tMS Order\tCorrelation to Averagine\tFraction Intensity Observed\tS/N\tLog2 Total Intensity\tNoise\tFeatureID" };

            var ms1Scans = data.GetMS1Scans().ToList();
            var ms1ScanNumbers = ms1Scans.Select(p => p.OneBasedScanNumber).ToList();
            Tolerance tol = new PpmTolerance(10);

            foreach (var scanNum in ms1ScanNumbers)
            {
                ScanNumToClaimedMzs.Add(scanNum, new HashSet<double>());
            }
            List<DeconvolutedPeak> ff = new List<DeconvolutedPeak>();

            int featureId = 0;
            foreach (var seedEnvelope in envsList.OrderByDescending(p => p.Item2.Score))
            {
                featureId++;
                int scanNumber = seedEnvelope.Item1;
                int seedScanIndex = ms1ScanNumbers.IndexOf(scanNumber);
                List<(MsDataScan, IsotopicEnvelope)> featureEnvelopes = new List<(MsDataScan, IsotopicEnvelope)>();

                if (seedScanIndex < 0)
                {
                    // not an MS1 scan
                    continue;
                }

                var theEnvelope = seedEnvelope.Item2;
                double mass = theEnvelope.Peaks.First().mz.ToMass(theEnvelope.Charge);

                int scanDirection = 1;
                int zDirection = 1;

                for (int scanIndex = seedScanIndex; scanIndex < ms1Scans.Count && scanIndex >= 0; scanIndex += scanDirection)
                {
                    var theScan = ms1Scans[scanIndex];
                    var spectrum = theScan.MassSpectrum;
                    int envsFoundInScan = 0;
                    var claimedPeaks = ScanNumToClaimedMzs[theScan.OneBasedScanNumber];

                    for (int z = theEnvelope.Charge; z <= 60 && z >= 1; z += zDirection)
                    {
                        double mainPeakMz = mass.ToMz(z);
                        int peakIndex = spectrum.GetClosestPeakIndex(mainPeakMz);
                        double experimentalMz = spectrum.XArray[peakIndex];
                        double experimentalIntensity = spectrum.YArray[peakIndex];

                        double experimentalMass = experimentalMz.ToMass(z);

                        var putativeEnv = spectrum.GetIsotopicEnvelope(experimentalMz, experimentalIntensity, z, tol, 3, ff, claimedPeaks);

                        if (tol.Within(mass, experimentalMass)
                            && putativeEnv != null
                            && putativeEnv.PearsonCorrelation > 0.4
                            && putativeEnv.FracIntensityObserved > 0.6)
                        {
                            featureEnvelopes.Add((theScan, putativeEnv));
                            envsFoundInScan++;

                            foreach (var peak in putativeEnv.Peaks)
                            {
                                claimedPeaks.Add(peak.mz);
                            }
                        }
                        else
                        {
                            if (zDirection == 1)
                            {
                                zDirection = -1;
                                z = theEnvelope.Charge;
                            }
                            else
                            {
                                break;
                            }
                        }
                    }

                    if (envsFoundInScan == 0)
                    {
                        if (scanDirection == 1)
                        {
                            scanDirection = -1;
                            zDirection = 1;
                            scanIndex = seedScanIndex;
                        }
                        else
                        {
                            break;
                        }
                    }
                }

                foreach (var envWithScan in featureEnvelopes)
                {
                    var scan = envWithScan.Item1;
                    var env = envWithScan.Item2;

                    var line = scan.OneBasedScanNumber + "\t" +
                        env.MonoisotopicMass + "\t" +
                        string.Join(";", env.Peaks.Select(p => p.mz.ToString("F3"))) + "\t" +
                        string.Join(";", env.Peaks.Select(p => p.intensity.ToString("F1"))) + "\t" +
                        env.Charge + "\t" +
                        scan.MsnOrder + "\t" +
                        env.PearsonCorrelation + "\t" +
                        env.FracIntensityObserved + "\t" +
                        env.SN + "\t" +
                        Math.Log(env.TotalIntensity, 2) + "\t" +
                        env.Noise + "\t" +
                        featureId + "\t" +
                        featureEnvelopes.Count;

                    lock (output)
                    {
                        output.Add(line);
                    }
                }
            }

            File.WriteAllLines(Path.Combine(dir, Path.GetFileNameWithoutExtension(filepath) + "_decon_features.tsv"), output);
        }
    }
}
