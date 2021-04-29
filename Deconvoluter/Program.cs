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

            // decon settings
            double minMass = 0;
            double fracIntensity = 0.4;
            double intensityRatio = 3;
            double pearsonCorr = 0.4;
            double sn = 1.5;
            double ppmTolerance = 5;



            // get "ground truth" envelopes
            HashSet<int> scansToDeconvolute = new HashSet<int>();
            Dictionary<int, List<IsotopicEnvelope>> groundTruthEnvelopes = new Dictionary<int, List<IsotopicEnvelope>>();
            var linesff = File.ReadAllLines(@"C:\Data\LVS_TD_Yeast\precursor_envelopes.txt");
            int j = 0;
            foreach (var line in linesff)
            {
                j++;

                if (j == 1)
                {
                    continue;
                }

                var split = line.Split(new char[] { '\t' });

                string filename = split[0];
                string seq = split[3];
                int charge = int.Parse(split[2]);
                int scanNum = int.Parse(split[1]);

                if (filename != "05-26-17_B7A_yeast_td_fract7_rep1" || seq.Contains("|"))
                {
                    continue;
                }

                scansToDeconvolute.Add(scanNum);

                Proteomics.AminoAcidPolymer.Peptide baseSequence = new Proteomics.AminoAcidPolymer.Peptide(seq);
                var formula = baseSequence.GetChemicalFormula();
                var isotopicDistribution = IsotopicDistribution.GetDistribution(formula, 0.125, 1e-8);
                double[] masses = isotopicDistribution.Masses.ToArray();
                double[] abundances = isotopicDistribution.Intensities.ToArray();

                double max = abundances.Max();
                int indOfMax = Array.IndexOf(abundances, max);
                double modeMass = masses[indOfMax];
                double modeMz = modeMass.ToMz(charge);

                var theScan = data.GetOneBasedScan(scanNum);
                var closestMzIndex = theScan.MassSpectrum.GetClosestPeakIndex(modeMz);
                var closestMz = theScan.MassSpectrum.XArray[closestMzIndex];
                var closestIntensity = theScan.MassSpectrum.YArray[closestMzIndex];

                var envelope = theScan.MassSpectrum.GetIsotopicEnvelope(closestMz, closestIntensity, charge, new PpmTolerance(ppmTolerance),
                    intensityRatio, new List<DeconvolutedPeak>(), new HashSet<double>(), fracIntensity, pearsonCorr);

                if (!groundTruthEnvelopes.ContainsKey(theScan.OneBasedScanNumber))
                {
                    groundTruthEnvelopes.Add(theScan.OneBasedScanNumber, new List<IsotopicEnvelope>());
                }

                if (envelope != null)
                {
                    envelope.SetMedianMonoisotopicMass(new List<double> { baseSequence.MonoisotopicMass });
                    groundTruthEnvelopes[theScan.OneBasedScanNumber].Add(envelope);
                }
            }






            List<string> output = new List<string> { IsotopicEnvelope.OutputHeader() };

            Random r = new Random(1);
            var scans = data.GetAllScansList().OrderBy(p => r.Next())
                //.Where(p => p.OneBasedScanNumber == 1437)
                .Where(p => p.MsnOrder == 1)
                .Where(p => scansToDeconvolute.Contains(p.OneBasedScanNumber))
                .ToList();

            var envsList = new List<(int, MassSpectrometry.IsotopicEnvelope)>();
            //var testt = new List<string>();

            Parallel.ForEach(Partitioner.Create(0, scans.Count),
               new ParallelOptions { MaxDegreeOfParallelism = -1 },
               (range, loopState) =>
               {
                   for (int i = range.Item1; i < range.Item2; i++)
                   {
                       var scan = scans[i];

                       var envs = scan.MassSpectrum.Deconvolute(scan.ScanWindowRange, 1, 60, ppmTolerance, intensityRatio, 2, pearsonCorr, fracIntensity, sn, scan, groundTruthEnvelopes)
                            .Where(p => p.MonoisotopicMass > minMass).OrderBy(p => p.Peaks.First().mz);

                       int identifier = 0;
                       foreach (var env in envs.OrderBy(p => p.Peaks.First().mz))
                       {
                           env.Scan = scan;
                           env.EnvelopeIdentifier = identifier;

                           var line = env.ToOutputString();

                           lock (output)
                           {
                               output.Add(line);
                           }

                           lock (envsList)
                           {
                               envsList.Add((scan.OneBasedScanNumber, env));
                           }
                           identifier++;
                       }
                   }
               });

            File.WriteAllLines(Path.Combine(dir, Path.GetFileNameWithoutExtension(filepath) + "_decon.tsv"), output);


            // test mono mass errors in decon results
            var lines = File.ReadAllLines(@"C:\Data\LVS_TD_Yeast\precursor_envelopes.txt");
            int i = 0;
            output = new List<string> { IsotopicEnvelope.OutputHeader() };
            foreach (var line in lines)
            {
                i++;

                if (i == 1)
                {
                    continue;
                }

                var split = line.Split(new char[] { '\t' });

                string filename = split[0];
                string seq = split[3];
                int charge = int.Parse(split[2]);
                int scanNum = int.Parse(split[1]);

                if (filename != "05-26-17_B7A_yeast_td_fract7_rep1" || seq.Contains("|"))
                {
                    continue;
                }

                Proteomics.AminoAcidPolymer.Peptide baseSequence = new Proteomics.AminoAcidPolymer.Peptide(seq);
                var formula = baseSequence.GetChemicalFormula();
                var isotopicDistribution = IsotopicDistribution.GetDistribution(formula, 0.125, 1e-8);
                double[] masses = isotopicDistribution.Masses.ToArray();
                double[] abundances = isotopicDistribution.Intensities.ToArray();

                double max = abundances.Max();
                int indOfMax = Array.IndexOf(abundances, max);
                double modeMass = masses[indOfMax];
                double modeMz = modeMass.ToMz(charge);

                var theScan = data.GetOneBasedScan(scanNum);
                var closestMzIndex = theScan.MassSpectrum.GetClosestPeakIndex(modeMz);
                var closestMz = theScan.MassSpectrum.XArray[closestMzIndex];
                var closestIntensity = theScan.MassSpectrum.YArray[closestMzIndex];

                // get deconvoluted envelope
                var envelopes = envsList.Where(p => p.Item1 == scanNum).ToList();
                var envelopeWithScan = envelopes.FirstOrDefault(p => p.Item2.Peaks.Select(v => v.mz).Contains(closestMz));

                // get the correct envelope according to known monoisotopic mass and charge
                var correctEnvelope = theScan.MassSpectrum.GetIsotopicEnvelope(closestMz, closestIntensity, charge,
                    new PpmTolerance(ppmTolerance), intensityRatio, new List<DeconvolutedPeak>(),
                    new HashSet<double>());

                if (correctEnvelope != null)
                {
                    correctEnvelope.SetMedianMonoisotopicMass(new List<double> { baseSequence.MonoisotopicMass });
                    correctEnvelope.Scan = theScan;
                }

                string errorType = "Correct";

                if (correctEnvelope != null && envelopeWithScan.Item2 == null)
                {
                    errorType = "Missing";
                    correctEnvelope.EnvelopeIdentifier = i;
                    output.Add(correctEnvelope.ToOutputString() + "\t" + -20 + "\t" + errorType);
                    continue;
                }
                else if (envelopeWithScan.Item2 == null)
                {
                    continue;
                }

                //if (envelopeWithScan.Item2 == null)
                //{
                //    //output.Add(seq + "\t" + baseSequence.MonoisotopicMass + "\t" + charge + "\t" + double.NaN + "\t" + double.NaN + "\t" + double.NaN + "\t" + "Not found");
                //    continue;
                //}

                var envelope = envelopeWithScan.Item2;

                var massError = envelope.MonoisotopicMass - baseSequence.MonoisotopicMass;
                double monoMassMz = baseSequence.MonoisotopicMass.ToMz(charge);


                if (charge != envelope.Charge)
                {
                    errorType = "Wrong charge (z=" + charge + ")";
                }
                else if (Math.Abs(massError) > 0.5)
                {
                    errorType = "Mono Mass Error";
                }

                envelope.EnvelopeIdentifier = i;
                output.Add(envelope.ToOutputString() + "\t" + massError + "\t" + errorType);
            }

            File.WriteAllLines(Path.Combine(dir, Path.GetFileNameWithoutExtension(filepath) + "_monoMassErrors.tsv"), output);
            //File.WriteAllLines(@"C:\Data\LVS_TD_Yeast\monoMassTest.tsv", testt);





            //Dictionary<int, HashSet<double>> ScanNumToClaimedMzs = new Dictionary<int, HashSet<double>>();

            //output = new List<string> { "Scan\tMonoisotopic Mass\tPeaks m/z List\tPeaks Intensity List\t" +
            //    "Charge\tMS Order\tCorrelation to Averagine\tFraction Intensity Observed\tS/N\tLog2 Total Intensity\tNoise\tFeatureID" };

            //var ms1Scans = data.GetMS1Scans().ToList();
            //var ms1ScanNumbers = ms1Scans.Select(p => p.OneBasedScanNumber).ToList();
            //Tolerance tol = new PpmTolerance(ppmTolerance);

            //foreach (var scanNum in ms1ScanNumbers)
            //{
            //    ScanNumToClaimedMzs.Add(scanNum, new HashSet<double>());
            //}
            //List<DeconvolutedPeak> ff = new List<DeconvolutedPeak>();

            //int featureId = 0;
            //foreach (var seedEnvelope in envsList.OrderByDescending(p => p.Item2.Peaks.Count).ThenByDescending(p => p.Item2.PearsonCorrelation))
            //{
            //    featureId++;
            //    int scanNumber = seedEnvelope.Item1;
            //    int seedScanIndex = ms1ScanNumbers.IndexOf(scanNumber);
            //    List<(MsDataScan, IsotopicEnvelope)> featureEnvelopes = new List<(MsDataScan, IsotopicEnvelope)>();

            //    if (seedScanIndex < 0)
            //    {
            //        // not an MS1 scan
            //        continue;
            //    }

            //    var theEnvelope = seedEnvelope.Item2;
            //    double mass = theEnvelope.Peaks.First().mz.ToMass(theEnvelope.Charge);

            //    int scanDirection = 1;
            //    int zDirection = 1;

            //    for (int scanIndex = seedScanIndex; scanIndex < ms1Scans.Count && scanIndex >= 0; scanIndex += scanDirection)
            //    {
            //        var theScan = ms1Scans[scanIndex];

            //        var spectrum = theScan.MassSpectrum;
            //        int envsFoundInScan = 0;
            //        var claimedPeaks = ScanNumToClaimedMzs[theScan.OneBasedScanNumber];

            //        for (int z = theEnvelope.Charge; z <= 60 && z >= 1; z += zDirection)
            //        {
            //            double mainPeakMz = mass.ToMz(z);
            //            int peakIndex = spectrum.GetClosestPeakIndex(mainPeakMz);
            //            double experimentalMz = spectrum.XArray[peakIndex];
            //            double experimentalIntensity = spectrum.YArray[peakIndex];

            //            double experimentalMass = experimentalMz.ToMass(z);

            //            var putativeEnv = spectrum.GetIsotopicEnvelope(experimentalMz, experimentalIntensity, z, tol, intensityRatio, ff, claimedPeaks);

            //            if (tol.Within(mass, experimentalMass)
            //                && putativeEnv != null
            //                && putativeEnv.PearsonCorrelation > pearsonCorr
            //                && putativeEnv.FracIntensityObserved > fracIntensity)
            //            {
            //                featureEnvelopes.Add((theScan, putativeEnv));
            //                envsFoundInScan++;
            //            }
            //            else
            //            {
            //                if (zDirection == 1)
            //                {
            //                    zDirection = -1;
            //                    z = theEnvelope.Charge;
            //                }
            //                else
            //                {
            //                    break;
            //                }
            //            }
            //        }

            //        if (envsFoundInScan == 0)
            //        {
            //            if (scanDirection == 1)
            //            {
            //                scanDirection = -1;
            //                zDirection = 1;
            //                scanIndex = seedScanIndex;
            //            }
            //            else
            //            {
            //                break;
            //            }
            //        }
            //    }

            //    if (featureEnvelopes.Count >= 3)
            //    {
            //        for (int i = 0; i < featureEnvelopes.Count; i++)
            //        {
            //            var envWithScan = featureEnvelopes[i];
            //            var scan = envWithScan.Item1;
            //            var env = envWithScan.Item2;
            //            bool seed = i == 0;

            //            var claimedPeaks = ScanNumToClaimedMzs[scan.OneBasedScanNumber];
            //            foreach (var peak in env.Peaks)
            //            {
            //                claimedPeaks.Add(peak.mz);
            //            }

            //            var line = scan.OneBasedScanNumber + "\t" +
            //                env.MonoisotopicMass + "\t" +
            //                string.Join(";", env.Peaks.Select(p => p.mz.ToString("F3"))) + "\t" +
            //                string.Join(";", env.Peaks.Select(p => p.intensity.ToString("F1"))) + "\t" +
            //                env.Charge + "\t" +
            //                scan.MsnOrder + "\t" +
            //                env.PearsonCorrelation + "\t" +
            //                env.FracIntensityObserved + "\t" +
            //                env.SN + "\t" +
            //                Math.Log(env.TotalIntensity, 2) + "\t" +
            //                env.Noise + "\t" +
            //                featureId + "\t" +
            //                featureEnvelopes.Count + "\t" +
            //                (seed ? 1 : 0);

            //            output.Add(line);
            //        }
            //    }
            //}

            //File.WriteAllLines(Path.Combine(dir, Path.GetFileNameWithoutExtension(filepath) + "_decon_features.tsv"), output);
        }
    }
}
