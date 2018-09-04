using Chemistry;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace FlashLFQ
{
    public class FlashLFQEngine
    {
        // settings
        public readonly bool Silent;
        public readonly int MaxThreads;
        public readonly double PeakfindingPpmTolerance;
        public readonly double PpmTolerance;
        public readonly double RtTol;
        public readonly double IsotopePpmTolerance;
        public readonly bool Integrate;
        public readonly int MissedScansAllowed;
        public readonly int NumIsotopesRequired;
        public readonly double MbrRtWindow;
        public readonly double InitialMbrRtWindow;
        public readonly double MbrPpmTolerance;
        public readonly bool ErrorCheckAmbiguousMatches;
        public readonly bool MatchBetweenRuns;
        public readonly bool IdSpecificChargeState;
        public readonly bool RequireMonoisotopicMass;
        public readonly bool Normalize;
        public readonly double MinDiscFactorToCutAt;
        public readonly bool AdvancedProteinQuant;

        // structures used in the FlashLFQ engine
        private List<SpectraFileInfo> spectraFileInfo;
        private Stopwatch globalStopwatch;
        private List<Identification> allIdentifications;
        private Dictionary<string, List<KeyValuePair<double, double>>> baseSequenceToIsotopicDistribution;
        private IEnumerable<int> chargeStates;
        private FlashLfqResults results;
        private int binsPerDalton = 100;
        private MsDataScan[] ms1Scans;
        private List<IndexedMassSpectralPeak>[] indexedPeaks;

        public FlashLFQEngine(List<Identification> allIdentifications, bool normalize = false, bool advancedProteinQuant = false, bool matchBetweenRuns = false, double ppmTolerance = 10.0, double isotopeTolerancePpm = 5.0, double matchBetweenRunsPpmTolerance = 5.0, bool integrate = false, int numIsotopesRequired = 2, bool idSpecificChargeState = false, bool requireMonoisotopicMass = true, bool silent = false, string optionalPeriodicTablePath = null, double maxMbrWindow = 1.5, int maxThreads = -1)
        {
            if (optionalPeriodicTablePath == null)
                optionalPeriodicTablePath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"elements.dat");
            Loaders.LoadElements(optionalPeriodicTablePath);

            globalStopwatch = new Stopwatch();
            chargeStates = new List<int>();

            this.spectraFileInfo = allIdentifications.Select(p => p.fileInfo).Distinct()
                .OrderBy(p => p.Condition)
                .ThenBy(p => p.BiologicalReplicate)
                .ThenBy(p => p.Fraction)
                .ThenBy(p => p.TechnicalReplicate).ToList();

            this.allIdentifications = allIdentifications;
            this.PpmTolerance = ppmTolerance;
            this.IsotopePpmTolerance = isotopeTolerancePpm;
            this.MatchBetweenRuns = matchBetweenRuns;
            this.MbrPpmTolerance = matchBetweenRunsPpmTolerance;
            this.Integrate = integrate;
            this.NumIsotopesRequired = numIsotopesRequired;
            this.Silent = silent;
            this.IdSpecificChargeState = idSpecificChargeState;
            this.RequireMonoisotopicMass = requireMonoisotopicMass;
            this.MbrRtWindow = maxMbrWindow;
            this.Normalize = normalize;
            this.MaxThreads = maxThreads;

            PeakfindingPpmTolerance = 20.0;
            InitialMbrRtWindow = 10.0;
            MissedScansAllowed = 1;
            RtTol = 5.0;
            ErrorCheckAmbiguousMatches = true;
            MinDiscFactorToCutAt = 0.6;
            AdvancedProteinQuant = advancedProteinQuant;
        }

        public FlashLfqResults Run()
        {
            globalStopwatch.Start();

            results = new FlashLfqResults(spectraFileInfo);

            // build m/z index keys
            CalculateTheoreticalIsotopeDistributions();

            // quantify each file
            foreach (var spectraFile in spectraFileInfo)
            {
                // fill lookup-table with peaks from the raw file
                IndexMassSpectralPeaks(spectraFile);

                if (indexedPeaks == null || indexedPeaks.Length == 0)
                {
                    // no MS1 peaks found
                    return results;
                }

                // quantify features using this file's IDs first
                QuantifyMS2IdentifiedPeptides(spectraFile);

                // find unidentified features based on other files' identification results (initial MBR peak-finding)
                if (MatchBetweenRuns)
                {
                    MatchBetweenRunsInitialPeakfinding(spectraFile);
                }

                // error checking function
                // handles features with multiple identifying scans and scans that are associated with more than one feature
                RunErrorChecking(spectraFile, false);

                if (!Silent)
                {
                    Console.WriteLine("Finished " + spectraFile.FilenameWithoutExtension);
                }

                // some memory-saving stuff
                ms1Scans = new MsDataScan[0];
                GC.Collect();
            }

            // filter initial MBR peaks with retention time calibration
            if (MatchBetweenRuns)
            {
                RetentionTimeCalibrationAndErrorCheckMatchedFeatures();
            }

            // normalize
            if (Normalize)
            {
                try
                {
                    new IntensityNormalizationEngine(results, Integrate, Silent, MaxThreads).NormalizeResults();
                }
                catch (Exception e)
                {
                    throw new MzLibException("A crash occured in FlashLFQ during the intensity normalization process:\n" + e.Message);
                }
            }

            // calculate intensities for proteins/peptides
            results.CalculatePeptideResults(true);
            results.CalculatePeptideResults(false);

            if (AdvancedProteinQuant)
            {
                new ProteinQuantificationEngine(results, MaxThreads).Run();
            }
            else
            {
                results.CalculateProteinResultsTop3();
            }

            // done
            if (!Silent)
            {
                Console.WriteLine("All done");
            }

            if (!Silent)
            {
                Console.WriteLine("Analysis time: " +
                    globalStopwatch.Elapsed.Hours + "h " +
                    globalStopwatch.Elapsed.Minutes + "m " +
                    globalStopwatch.Elapsed.Seconds + "s");
            }

            return results;
        }

        private void RetentionTimeCalibrationAndErrorCheckMatchedFeatures()
        {
            if (!Silent)
                Console.WriteLine("Running retention time calibration");

            // get all unambiguous peaks for all files
            var allFeatures = results.Peaks.SelectMany(p => p.Value).Where(p => !p.IsMbrFeature);
            var allAmbiguousFeatures = allFeatures.Where(p => p.NumIdentificationsByFullSeq > 1).ToList();
            var ambiguousFeatureSeqs = new HashSet<string>(allAmbiguousFeatures.SelectMany(p => p.Identifications.Select(v => v.ModifiedSequence)));

            foreach (var feature in allFeatures)
                if (ambiguousFeatureSeqs.Contains(feature.Identifications.First().ModifiedSequence))
                    allAmbiguousFeatures.Add(feature);

            var unambiguousPeaksGroupedByFile = allFeatures.Except(allAmbiguousFeatures).Where(v => v.Apex != null).GroupBy(p => p.RawFileInfo);

            foreach (var file in unambiguousPeaksGroupedByFile)
            {
                var allMbrFeaturesForThisFile = results.Peaks[file.Key].Where(p => p.IsMbrFeature);

                // get the best (most intense) peak for each peptide in the file
                Dictionary<string, ChromatographicPeak> pepToBestFeatureForThisFile = new Dictionary<string, ChromatographicPeak>();
                foreach (var testPeak in file)
                {
                    if (pepToBestFeatureForThisFile.TryGetValue(testPeak.Identifications.First().ModifiedSequence, out ChromatographicPeak currentBestPeak))
                    {
                        if (currentBestPeak.Intensity > testPeak.Intensity)
                            pepToBestFeatureForThisFile[testPeak.Identifications.First().ModifiedSequence] = testPeak;
                    }
                    else
                        pepToBestFeatureForThisFile.Add(testPeak.Identifications.First().ModifiedSequence, testPeak);
                }

                foreach (var otherFile in unambiguousPeaksGroupedByFile)
                {
                    // get the other files' best peak for the same peptides (to make an RT calibration curve)
                    if (otherFile.Key.Equals(file.Key))
                        continue;

                    var featuresInCommon = otherFile.Where(p => pepToBestFeatureForThisFile.ContainsKey(p.Identifications.First().ModifiedSequence));

                    Dictionary<string, ChromatographicPeak> pepToBestFeatureForOtherFile = new Dictionary<string, ChromatographicPeak>();
                    foreach (var testPeak in featuresInCommon)
                    {
                        if (pepToBestFeatureForOtherFile.TryGetValue(testPeak.Identifications.First().ModifiedSequence, out ChromatographicPeak currentBestPeak))
                        {
                            if (currentBestPeak.Intensity > testPeak.Intensity)
                                pepToBestFeatureForOtherFile[testPeak.Identifications.First().ModifiedSequence] = testPeak;
                        }
                        else
                            pepToBestFeatureForOtherFile.Add(testPeak.Identifications.First().ModifiedSequence, testPeak);
                    }

                    // create a rt-to-rt correlation for the two files' peptides
                    Dictionary<string, Tuple<double, double>> rtCalPoints = new Dictionary<string, Tuple<double, double>>();

                    foreach (var kvp in pepToBestFeatureForOtherFile)
                        rtCalPoints.Add(kvp.Key, new Tuple<double, double>(pepToBestFeatureForThisFile[kvp.Key].Apex.RetentionTime, kvp.Value.Apex.RetentionTime));

                    if (!rtCalPoints.Any())
                        continue;

                    var differencesInRt = rtCalPoints.Select(p => (p.Value.Item1 - p.Value.Item2));
                    double average = differencesInRt.Average();
                    double sd = MathNet.Numerics.Statistics.Statistics.StandardDeviation(differencesInRt);

                    // remove extreme outliers
                    if (sd > 1.0)
                    {
                        var pointsToRemove = rtCalPoints.Where(p => (p.Value.Item1 - p.Value.Item2) > (average + sd) || (p.Value.Item1 - p.Value.Item2) < (average - sd)).ToList();
                        foreach (var point in pointsToRemove)
                            rtCalPoints.Remove(point.Key);

                        differencesInRt = rtCalPoints.Select(p => (p.Value.Item1 - p.Value.Item2));
                        average = differencesInRt.Average();
                        sd = MathNet.Numerics.Statistics.Statistics.StandardDeviation(differencesInRt);
                    }

                    // find rt differences between files
                    List<Tuple<double, double>> rtCalPoints2 = rtCalPoints.Values.OrderBy(p => p.Item1).ToList();

                    int minRt = (int)Math.Floor(rtCalPoints2.First().Item1);
                    int maxRt = (int)Math.Ceiling(rtCalPoints2.Last().Item1);
                    var rtCalibrationDictionary = new Dictionary<int, List<double>>();
                    foreach (var rt in rtCalPoints2)
                    {
                        int intRt = (int)Math.Round(rt.Item1);
                        if (rtCalibrationDictionary.TryGetValue(intRt, out List<double> points))
                            points.Add(rt.Item1 - rt.Item2);
                        else
                            rtCalibrationDictionary.Add(intRt, new List<double> { rt.Item1 - rt.Item2 });
                    }

                    // create spline
                    // index is minute of source, double is rt calibration factor (in minutes) for destination file
                    double[] rtCalRunningSpline = new double[maxRt + 1];
                    double[] stdevRunningSpline = new double[maxRt + 1];
                    for (int i = 0; i < rtCalRunningSpline.Length; i++)
                    {
                        rtCalRunningSpline[i] = double.NaN;
                        stdevRunningSpline[i] = double.NaN;
                    }

                    // calculate stdev for each element in spline
                    for (int i = 1; i <= maxRt; i++)
                    {
                        if (rtCalibrationDictionary.TryGetValue(i, out List<double> rtDifferencesForThisTime))
                        {
                            if (rtDifferencesForThisTime.Count > 3)
                            {
                                rtDifferencesForThisTime.Sort();
                                rtCalRunningSpline[i] = rtDifferencesForThisTime[rtDifferencesForThisTime.Count / 2]; // median rt difference for this timepoint

                                average = rtDifferencesForThisTime.Average();
                                sd = MathNet.Numerics.Statistics.Statistics.StandardDeviation(rtDifferencesForThisTime);
                                var rtError = 3 * sd;

                                if (rtError > (MbrRtWindow / 2.0))
                                    stdevRunningSpline[i] = MbrRtWindow / 2.0;
                                else
                                    stdevRunningSpline[i] = rtError;
                            }
                        }
                    }

                    // fill gaps in spline (linear interpolation)
                    for (int i = minRt; i <= maxRt; i++)
                    {
                        if (double.IsNaN(rtCalRunningSpline[i]))
                        {
                            KeyValuePair<int, double> prevCalPoint = new KeyValuePair<int, double>(0, double.NaN);
                            KeyValuePair<int, double> nextCalPoint = new KeyValuePair<int, double>(0, double.NaN);

                            for (int j = i; j >= minRt; j--)
                            {
                                if (j < i - 3)
                                    break;
                                if (!double.IsNaN(rtCalRunningSpline[j]))
                                {
                                    prevCalPoint = new KeyValuePair<int, double>(j, rtCalRunningSpline[j]);
                                    break;
                                }
                            }
                            for (int j = i; j <= maxRt; j++)
                            {
                                if (j > i + 3)
                                    break;
                                if (!double.IsNaN(rtCalRunningSpline[j]))
                                {
                                    nextCalPoint = new KeyValuePair<int, double>(j, rtCalRunningSpline[j]);
                                    break;
                                }
                            }

                            if (!double.IsNaN(prevCalPoint.Value) && !double.IsNaN(nextCalPoint.Value))
                            {
                                var slope = (prevCalPoint.Value - nextCalPoint.Value) / (prevCalPoint.Key - nextCalPoint.Key);
                                var yint = prevCalPoint.Value - slope * prevCalPoint.Key;
                                double interpolatedRtCalPoint = slope * i + yint;

                                rtCalRunningSpline[i] = interpolatedRtCalPoint;
                                stdevRunningSpline[i] = stdevRunningSpline[i] = MbrRtWindow / 2.0;
                            }
                        }
                    }

                    // finished rt calibration for these 2 files; now use rt cal spline to find matched features
                    var allMatchedFeaturesToLookForNow = allMbrFeaturesForThisFile.Where(p => p.Identifications.First().fileInfo.Equals(otherFile.Key)).ToList();

                    // filter peak candidates with rt cal to get apex peak
                    foreach (var mbrFeature in allMatchedFeaturesToLookForNow)
                    {
                        if (mbrFeature.IsotopicEnvelopes.Any())
                        {
                            // shift = thisFileRt - otherFileRt
                            int rtSplineLookupTime = (int)Math.Round(mbrFeature.Identifications.First().ms2RetentionTimeInMinutes);

                            if (rtSplineLookupTime < rtCalRunningSpline.Length)
                            {
                                double rtShift = rtCalRunningSpline[rtSplineLookupTime];
                                double rtToleranceHere = stdevRunningSpline[rtSplineLookupTime];
                                double theoreticalRt = mbrFeature.Identifications.First().ms2RetentionTimeInMinutes + rtShift;

                                if (!double.IsNaN(rtShift))
                                    mbrFeature.IsotopicEnvelopes = mbrFeature.IsotopicEnvelopes.Where(p => Math.Abs(p.RetentionTime - theoreticalRt) < rtToleranceHere).ToList();
                                else
                                    mbrFeature.IsotopicEnvelopes = new List<IsotopicEnvelope>();
                            }
                            else
                                mbrFeature.IsotopicEnvelopes = new List<IsotopicEnvelope>();
                        }
                    }

                    foreach (var feature in allMatchedFeaturesToLookForNow)
                        if (feature.IsotopicEnvelopes.Any())
                            feature.CalculateIntensityForThisFeature(Integrate);
                }
            }

            foreach (var file in spectraFileInfo)
                RunErrorChecking(file, true);
        }

        private void CalculateTheoreticalIsotopeDistributions()
        {
            baseSequenceToIsotopicDistribution = new Dictionary<string, List<KeyValuePair<double, double>>>();

            // calculate monoisotopic masses and isotopic envelope
            foreach (var id in allIdentifications)
            {
                if (baseSequenceToIsotopicDistribution.ContainsKey(id.BaseSequence))
                    continue;

                var isotopicMassesAndNormalizedAbundances = new List<KeyValuePair<double, double>>();

                Proteomics.AminoAcidPolymer.Peptide p = new Proteomics.AminoAcidPolymer.Peptide(id.BaseSequence);
                int numCarbonsInThisPeptide = p.ElementCountWithIsotopes("C");

                var isotopicDistribution = IsotopicDistribution.GetDistribution(p.GetChemicalFormula(), 0.125, 1e-8);

                var masses = isotopicDistribution.Masses.ToArray();
                var abundances = isotopicDistribution.Intensities.ToArray();

                var monoisotopicMass = masses.Min();
                var highestAbundance = abundances.Max();

                for (int i = 0; i < masses.Length; i++)
                {
                    // expected isotopic mass shifts for peptide of this length
                    masses[i] -= monoisotopicMass;

                    // normalized abundance of each mass shift
                    abundances[i] /= highestAbundance;

                    // look for these isotopes
                    if (i < (NumIsotopesRequired - 1) || abundances[i] > 0.1)
                        isotopicMassesAndNormalizedAbundances.Add(new KeyValuePair<double, double>(masses[i], abundances[i]));
                }

                baseSequenceToIsotopicDistribution.Add(id.BaseSequence, isotopicMassesAndNormalizedAbundances);
            }

            var minChargeState = allIdentifications.Min(p => p.precursorChargeState);
            var maxChargeState = allIdentifications.Max(p => p.precursorChargeState);
            chargeStates = Enumerable.Range(minChargeState, (maxChargeState - minChargeState) + 1);

            var peptideModifiedSequences = allIdentifications.GroupBy(p => p.ModifiedSequence);
            foreach (var identifications in peptideModifiedSequences)
            {
                double lowestCommonMassShift = baseSequenceToIsotopicDistribution[identifications.First().BaseSequence].Select(p => p.Key).Min();
                var mostCommonIsotopeShift = baseSequenceToIsotopicDistribution[identifications.First().BaseSequence].Where(p => p.Value == 1).First().Key;

                var thisPeptidesLowestCommonMass = identifications.First().monoisotopicMass + lowestCommonMassShift;
                var thisPeptidesMostAbundantMass = identifications.First().monoisotopicMass + mostCommonIsotopeShift;

                foreach (var identification in identifications)
                {
                    identification.massToLookFor = RequireMonoisotopicMass ? identifications.First().monoisotopicMass : thisPeptidesMostAbundantMass;
                }
            }
        }

        private void IndexMassSpectralPeaks(SpectraFileInfo fileInfo)
        {
            if (!Silent)
            {
                Console.WriteLine("Reading spectra file");
            }

            ms1Scans = new MsDataScan[0];

            if (indexedPeaks != null)
            {
                for (int i = 0; i < indexedPeaks.Length; i++)
                {
                    if (indexedPeaks[i] == null)
                        continue;

                    indexedPeaks[i].Clear();
                }
            }

            // read spectra file
            var ext = Path.GetExtension(fileInfo.FullFilePathWithExtension).ToUpperInvariant();
            if (ext == ".MZML")
            {
                try
                {
                    ms1Scans = Mzml.LoadAllStaticData(fileInfo.FullFilePathWithExtension).GetAllScansList().OrderBy(p => p.OneBasedScanNumber).ToArray();
                }
                catch (FileNotFoundException)
                {
                    if (!Silent)
                    {
                        Console.WriteLine("\nCan't find mzml file" + fileInfo.FullFilePathWithExtension + "\n");
                    }
                    return;
                }
                catch (Exception e)
                {
                    if (!Silent)
                    {
                        Console.WriteLine("Problem opening mzml file " + fileInfo.FullFilePathWithExtension + "; " + e.Message);
                    }
                    return;
                }

                for (int i = 0; i < ms1Scans.Length; i++)
                {
                    if (ms1Scans[i].MsnOrder > 1)
                    {
                        ms1Scans[i] = null;
                    }
                }
            }
            else if (ext == ".RAW")
            {
#if NETFRAMEWORK
                using (var thermoDynamicConnection = IO.Thermo.ThermoDynamicData.InitiateDynamicConnection(fileInfo.FullFilePathWithExtension))
                {
                    var tempList = new List<MsDataScan>();

                    try
                    {
                        // use thermo dynamic connection to get the ms1 scans and then dispose of the connection
                        int[] msOrders = thermoDynamicConnection.ThermoGlobalParams.MsOrderByScan;
                        for (int i = 0; i < msOrders.Length; i++)
                        {
                            if (msOrders[i] == 1)
                            {
                                tempList.Add(thermoDynamicConnection.GetOneBasedScan(i + 1));
                            }
                            else
                            {
                                tempList.Add(null);
                            }
                        }
                    }
                    catch (FileNotFoundException)
                    {
                        thermoDynamicConnection.Dispose();

                        if (!Silent)
                        {
                            Console.WriteLine("\nCan't find raw file" + fileInfo.FullFilePathWithExtension + "\n");
                        }
                        return;
                    }
                    catch (Exception e)
                    {
                        thermoDynamicConnection.Dispose();

                        if (!Silent)
                        {
                            throw new MzLibException("FlashLFQ Error: Problem opening raw file " + fileInfo.FullFilePathWithExtension + "; " + e.Message);
                        }
                    }

                    ms1Scans = tempList.ToArray();
                }
#else
                if (!Silent)
                {
                    Console.WriteLine("Cannot open RAW with .NETStandard code - are you on Linux? " + fileInfo.FullFilePathWithExtension);
                }
                return;
#endif
            }
            else
            {
                if (!Silent)
                {
                    Console.WriteLine("Unsupported file type " + ext);
                    return;
                }
            }

            if (!Silent)
            {
                Console.WriteLine("Indexing MS1 peaks");
            }

            if (!ms1Scans.Where(p => p != null).Any())
            {
                indexedPeaks = new List<IndexedMassSpectralPeak>[0];
                return;
            }

            indexedPeaks = new List<IndexedMassSpectralPeak>[(int)Math.Ceiling(ms1Scans.Where(p => p != null && p.MassSpectrum.LastX != null).Max(p => p.MassSpectrum.LastX.Value) * binsPerDalton) + 1];

            for (int i = 0; i < ms1Scans.Length; i++)
            {
                if (ms1Scans[i] == null)
                {
                    continue;
                }

                for (int j = 0; j < ms1Scans[i].MassSpectrum.XArray.Length; j++)
                {
                    int roundedMz = (int)Math.Round(ms1Scans[i].MassSpectrum.XArray[j] * binsPerDalton, 0);
                    if (indexedPeaks[roundedMz] == null)
                    {
                        indexedPeaks[roundedMz] = new List<IndexedMassSpectralPeak>();
                    }

                    indexedPeaks[roundedMz].Add(new IndexedMassSpectralPeak(ms1Scans[i].MassSpectrum.XArray[j], ms1Scans[i].MassSpectrum.YArray[j], j, ms1Scans[i].OneBasedScanNumber));
                }
            }
        }

        private void QuantifyMS2IdentifiedPeptides(SpectraFileInfo fileInfo)
        {
            if (!Silent)
            {
                Console.WriteLine("Quantifying peptides for " + fileInfo.FilenameWithoutExtension);
            }

            var identifications = allIdentifications.Where(p => p.fileInfo.Equals(fileInfo)).ToList();

            if (!identifications.Any())
            {
                return;
            }

            Tolerance peakfindingTol = new PpmTolerance(PeakfindingPpmTolerance);
            Tolerance tol = new PpmTolerance(PpmTolerance);

            var chromatographicPeaks = new ChromatographicPeak[identifications.Count];

            Parallel.ForEach(Partitioner.Create(0, identifications.Count),
                new ParallelOptions { MaxDegreeOfParallelism = MaxThreads },
                (range, loopState) =>
            {
                List<IndexedMassSpectralPeak> binPeaks = new List<IndexedMassSpectralPeak>();
                List<IsotopicEnvelope> isotopicEnvelopes = new List<IsotopicEnvelope>();

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    //// Stop loop if canceled
                    //if (GlobalVariables.StopLoops)
                    //{
                    //    loopState.Stop();
                    //    return;
                    //}

                    binPeaks.Clear();
                    isotopicEnvelopes.Clear();
                    var identification = identifications[i];
                    ChromatographicPeak msmsFeature = new ChromatographicPeak(identification, false, fileInfo);
                    chromatographicPeaks[i] = msmsFeature;

                    foreach (var chargeState in chargeStates)
                    {
                        if (IdSpecificChargeState && chargeState != identification.precursorChargeState)
                        {
                            continue;
                        }

                        // get indexed mass spectral peaks for this ID and charge
                        double theoreticalMz = identification.massToLookFor.ToMz(chargeState);
                        int ceilingMz = (int)Math.Ceiling(peakfindingTol.GetMaximumValue(identification.massToLookFor).ToMz(chargeState) * binsPerDalton);
                        int floorMz = (int)Math.Floor(peakfindingTol.GetMinimumValue(identification.massToLookFor).ToMz(chargeState) * binsPerDalton);

                        for (int j = floorMz; j <= ceilingMz; j++)
                        {
                            if (j < indexedPeaks.Length && indexedPeaks[j] != null)
                            {
                                foreach (var peak in indexedPeaks[j])
                                {
                                    if (peakfindingTol.Within(peak.Mz.ToMass(chargeState), identification.massToLookFor) &&
                                      Math.Abs(ms1Scans[peak.OneBasedScanNumber - 1].RetentionTime - identification.ms2RetentionTimeInMinutes) < RtTol)
                                    {
                                        binPeaks.Add(peak);
                                    }
                                }
                            }
                        }

                        // do peakfinding, isotopic envelope filtering, error checking, etc
                        if (binPeaks.Any())
                        {
                            // do peakfinding
                            Peakfind(binPeaks, identification, chargeState);

                            // filter again by smaller mz tolerance
                            binPeaks.RemoveAll(p => !tol.Within(p.Mz.ToMass(chargeState), identification.massToLookFor));

                            // filter by isotopic distribution
                            isotopicEnvelopes = GetIsotopicEnvelopes(binPeaks, identification, chargeState, true);

                            // if multiple mass spectral peaks in the same scan are valid, pick the one with the smallest mass error
                            var peaksInSameScan = isotopicEnvelopes.GroupBy(p => p.IndexedPeak.OneBasedScanNumber).Where(v => v.Count() > 1);
                            if (peaksInSameScan.Any())
                            {
                                foreach (var group in peaksInSameScan)
                                {
                                    var smallestMzError = group.Min(p => Math.Abs(p.IndexedPeak.Mz - theoreticalMz));
                                    var peakToUse = group.Where(p => Math.Abs(p.IndexedPeak.Mz - theoreticalMz) == smallestMzError).First();
                                    isotopicEnvelopes.RemoveAll(p => group.Contains(p) && p != peakToUse);
                                }
                            }

                            // sort isotope envelopes by scan number (sets up for peak cutting)
                            msmsFeature.IsotopicEnvelopes.AddRange(isotopicEnvelopes);
                            msmsFeature.IsotopicEnvelopes.Sort((x, y) => x.IndexedPeak.OneBasedScanNumber.CompareTo(y.IndexedPeak.OneBasedScanNumber));
                        }
                    }

                    msmsFeature.CalculateIntensityForThisFeature(Integrate);
                    CutPeak(msmsFeature);
                }
            });

            results.Peaks.Add(fileInfo, chromatographicPeaks.ToList());
        }

        private void MatchBetweenRunsInitialPeakfinding(SpectraFileInfo fileInfo)
        {
            if (!Silent)
            {
                Console.WriteLine("Finding possible matched peptides for " + fileInfo.FilenameWithoutExtension);
            }

            if (!results.Peaks.ContainsKey(fileInfo) || results.Peaks[fileInfo].Count == 0)
            {
                return;
            }

            Tolerance mbrTol = new PpmTolerance(MbrPpmTolerance);

            // ignore IDs from proteins that have no MS/MS evidence in this biorep (or no evidence in this spectra file, if bioreps are not defined)
            var allowedFilesMatchingFrom = new HashSet<SpectraFileInfo>(spectraFileInfo.Where(
                f => f.Condition == fileInfo.Condition && f.BiologicalReplicate == fileInfo.BiologicalReplicate));

            if(spectraFileInfo.Max(p => p.BiologicalReplicate) == 0)
            {
                allowedFilesMatchingFrom = new HashSet<SpectraFileInfo> { fileInfo };
            }

            var allowedProteins = new HashSet<ProteinGroup>(allIdentifications.Where(p => allowedFilesMatchingFrom.Contains(p.fileInfo)).SelectMany(v => v.proteinGroups));

            // get IDs to look for
            var concurrentBagOfMatchedFeatures = new ConcurrentBag<ChromatographicPeak>();
            var identificationsFromOtherRunsToLookFor = new List<Identification>();
            var idsGroupedByFullSeq = allIdentifications.GroupBy(p => p.ModifiedSequence);

            foreach (var fullSequenceGroup in idsGroupedByFullSeq)
            {
                // look for peptides with no ID's in this file
                var seqsByFilename = fullSequenceGroup.GroupBy(p => p.fileInfo);

                if (!seqsByFilename.Where(p => p.Key.Equals(fileInfo)).Any())
                {
                    // ignore proteins that have no IDs in this file
                    identificationsFromOtherRunsToLookFor.AddRange(fullSequenceGroup.Where(
                        id => id.proteinGroups.Any(pg => allowedProteins.Contains(pg))));
                }
            }

            if (identificationsFromOtherRunsToLookFor.Any())
            {
                Parallel.ForEach(Partitioner.Create(0, identificationsFromOtherRunsToLookFor.Count),
                    new ParallelOptions { MaxDegreeOfParallelism = MaxThreads },
                    (range, loopState) =>
                    {
                        for (int i = range.Item1; i < range.Item2; i++)
                        {
                            //// Stop loop if canceled
                            //if (GlobalVariables.StopLoops)
                            //{
                            //    loopState.Stop();
                            //    return;
                            //}

                            var identification = identificationsFromOtherRunsToLookFor[i];

                            ChromatographicPeak mbrFeature = new ChromatographicPeak(identification, true, fileInfo);

                            foreach (var chargeState in chargeStates)
                            {
                                double theorMzHere = identification.massToLookFor.ToMz(chargeState);
                                double mzTolHere = ((MbrPpmTolerance / 1e6) * identification.monoisotopicMass) / chargeState;

                                int floorMz = (int)Math.Floor(theorMzHere * binsPerDalton);
                                int ceilingMz = (int)Math.Ceiling(theorMzHere * binsPerDalton);

                                List<IndexedMassSpectralPeak> binPeaks = new List<IndexedMassSpectralPeak>();

                                for (int j = floorMz; j <= ceilingMz && j < indexedPeaks.Length; j++)
                                {
                                    if (indexedPeaks[j] != null)
                                    {
                                        foreach (var peak in indexedPeaks[j])
                                        {
                                            if (mbrTol.Within(peak.Mz.ToMass(chargeState), identification.massToLookFor) &&
                                              Math.Abs(ms1Scans[peak.OneBasedScanNumber - 1].RetentionTime - identification.ms2RetentionTimeInMinutes) < InitialMbrRtWindow)
                                            {
                                                binPeaks.Add(peak);
                                            }
                                        }
                                    }
                                }

                                // filter by isotopic distribution
                                var validIsotopeClusters = GetIsotopicEnvelopes(binPeaks, identification, chargeState, true);

                                if (validIsotopeClusters.Any())
                                {
                                    mbrFeature.IsotopicEnvelopes.AddRange(validIsotopeClusters);
                                }
                            }

                            if (mbrFeature.IsotopicEnvelopes.Any())
                            {
                                concurrentBagOfMatchedFeatures.Add(mbrFeature);
                            }
                        }
                    }
                );
            }

            results.Peaks[fileInfo].AddRange(concurrentBagOfMatchedFeatures);
        }

        private void RunErrorChecking(SpectraFileInfo rawFile, bool postRtCal)
        {
            if (!Silent)
            {
                Console.WriteLine("Checking errors");
            }

            // remove all MBR features with intensities lower than the least-intense MS/MS-identified peak in this file
            // this is to remove MBR peaks that matched to noise
            if (MatchBetweenRuns && postRtCal)
            {
                double minMsmsIdentifiedPeakIntensity = 0;
                var msmsIdent = results.Peaks[rawFile].Where(v => !v.IsMbrFeature && v.Intensity > 0);
                if (msmsIdent.Any())
                {
                    minMsmsIdentifiedPeakIntensity = msmsIdent.Min(v => v.Intensity);
                }

                foreach (var peak in results.Peaks[rawFile])
                {
                    if (peak.IsMbrFeature && peak.Intensity < minMsmsIdentifiedPeakIntensity)
                    {
                        peak.IsotopicEnvelopes = new List<IsotopicEnvelope>();
                    }
                }

                results.Peaks[rawFile].RemoveAll(p => p.IsMbrFeature && !p.IsotopicEnvelopes.Any());
            }

            // condense duplicate features (features with same sequence and apex peak)
            var featuresWithSamePeak = results.Peaks[rawFile].Where(v => v.Apex != null).GroupBy(p => p.Apex.IndexedPeak);
            featuresWithSamePeak = featuresWithSamePeak.Where(p => p.Count() > 1);

            foreach (var duplicateFeature in featuresWithSamePeak)
            {
                duplicateFeature.First().MergeFeatureWith(duplicateFeature, Integrate);
            }
            results.Peaks[rawFile].RemoveAll(p => p.Intensity == -1);

            // condense multiple peaks for the same peptide within a time window
            var featuresToMaybeMerge = results.Peaks[rawFile].Where(p => p.NumIdentificationsByFullSeq == 1 && p.Apex != null).GroupBy(p => p.Identifications.First().ModifiedSequence).Where(p => p.Count() > 1);
            if (featuresToMaybeMerge.Any())
            {
                foreach (var group in featuresToMaybeMerge)
                {
                    if (IdSpecificChargeState)
                    {
                        var group2 = group.ToList().GroupBy(p => p.Apex.ChargeState).Where(v => v.Count() > 1);

                        foreach (var group3 in group2)
                        {
                            foreach (var feature in group3)
                            {
                                if (feature.Intensity != -1)
                                {
                                    var featuresToMerge = group3.Where(p => Math.Abs(p.Apex.RetentionTime - feature.Apex.RetentionTime) < RtTol && p.Intensity != -1);
                                    if (featuresToMerge.Any())
                                        feature.MergeFeatureWith(featuresToMerge, Integrate);
                                }
                            }
                        }
                    }
                    else
                    {
                        foreach (var feature in group)
                        {
                            if (feature.Intensity != -1)
                            {
                                var featuresToMerge = group.Where(p => Math.Abs(p.Apex.RetentionTime - feature.Apex.RetentionTime) < RtTol && p.Intensity != -1);
                                if (featuresToMerge.Any())
                                    feature.MergeFeatureWith(featuresToMerge, Integrate);
                            }
                        }
                    }
                }

                results.Peaks[rawFile].RemoveAll(p => p.Intensity == -1);
            }

            if (postRtCal || !MatchBetweenRuns)
            {
                // check for multiple peptides per feature
                var scansWithMultipleDifferentIds = results.Peaks[rawFile].Where(p => p.NumIdentificationsByFullSeq > 1);
                var ambiguousFeatures = scansWithMultipleDifferentIds.Where(p => p.NumIdentificationsByBaseSeq > 1).ToList();

                // handle ambiguous features
                foreach (var ambiguousFeature in ambiguousFeatures)
                {
                    var msmsIdentsForThisFile = ambiguousFeature.Identifications.Where(p => p.fileInfo.Equals(ambiguousFeature.RawFileInfo));

                    if (!msmsIdentsForThisFile.Any())
                    {
                        // mbr matched more than one identification to this peak - cannot resolve
                        ambiguousFeature.Intensity = -1;
                    }
                    else
                    {
                        // msms identifications take precident over mbr features
                        ambiguousFeature.Identifications.RemoveAll(p => p.fileInfo != ambiguousFeature.RawFileInfo);
                        ambiguousFeature.ResolveIdentifications();
                    }
                }

                results.Peaks[rawFile].RemoveAll(p => p.Intensity == -1);
            }
        }

        private List<IsotopicEnvelope> GetIsotopicEnvelopes(List<IndexedMassSpectralPeak> peaks, Identification identification, int chargeState, bool lookForBadIsotope)
        {
            var isotopeClusters = new List<IsotopicEnvelope>();
            var isotopeMassShifts = baseSequenceToIsotopicDistribution[identification.BaseSequence];

            if (isotopeMassShifts.Count < NumIsotopesRequired)
            {
                return isotopeClusters;
            }

            double isotopeMzTol = ((IsotopePpmTolerance / 1e6) * identification.monoisotopicMass) / chargeState;

            foreach (var peak in peaks)
            {
                // calculate theoretical isotopes relative to observed peak
                var theorIsotopeMzs = new double[isotopeMassShifts.Count];
                int isotopicPeakUnitOfPeakZeroIsMono = Convert.ToInt32(peak.Mz.ToMass(chargeState) - identification.monoisotopicMass);
                var mainpeakMz = peak.Mz;

                // left of main peak
                for (int i = 0; i < isotopicPeakUnitOfPeakZeroIsMono; i++)
                    theorIsotopeMzs[i] = mainpeakMz + (isotopeMassShifts[i].Key / chargeState);

                // main peak and right of main peak
                for (int i = isotopicPeakUnitOfPeakZeroIsMono; i < isotopeMassShifts.Count; i++)
                    theorIsotopeMzs[i] = mainpeakMz + (isotopeMassShifts[i].Key / chargeState);
                theorIsotopeMzs = theorIsotopeMzs.OrderBy(p => p).ToArray();

                var lowestMzIsotopePossible = theorIsotopeMzs.First();
                lowestMzIsotopePossible -= isotopeMzTol;
                var highestMzIsotopePossible = theorIsotopeMzs.Last();
                highestMzIsotopePossible += isotopeMzTol;

                // get possible isotope peaks from the peak's scan
                List<Tuple<double, double>> possibleIsotopePeaks = new List<Tuple<double, double>>();

                // go backwards from the peak to find the lowest-mass isotope possible
                int earliestIsotopicPeakIndexPossible = peak.ZeroBasedIndexOfPeakInScan;
                var massSpectrum = ms1Scans[peak.OneBasedScanNumber - 1].MassSpectrum;
                for (int i = earliestIsotopicPeakIndexPossible; i >= 0; i--)
                {
                    if (massSpectrum.XArray[i] < lowestMzIsotopePossible)
                        break;
                    earliestIsotopicPeakIndexPossible = i;
                }

                // find the highest-mass isotope possible
                for (int i = earliestIsotopicPeakIndexPossible; i < massSpectrum.Size; i++)
                {
                    if (massSpectrum.XArray[i] > highestMzIsotopePossible)
                        break;
                    possibleIsotopePeaks.Add(new Tuple<double, double>(massSpectrum.XArray[i], massSpectrum.YArray[i]));
                }

                if (lookForBadIsotope)
                {
                    bool badPeak = false;
                    double prevIsotopePeakMz = (theorIsotopeMzs[0] - (1.003322 / chargeState));

                    for (int i = earliestIsotopicPeakIndexPossible; i > 0; i--)
                    {
                        if (Math.Abs(massSpectrum.XArray[i] - prevIsotopePeakMz) < isotopeMzTol)
                            if (massSpectrum.YArray[i] / peak.Intensity > 1.0)
                                badPeak = true;
                        if (massSpectrum.XArray[i] < (prevIsotopePeakMz - isotopeMzTol))
                            break;
                    }

                    if (badPeak)
                        continue;
                }

                // isotopic distribution check
                bool isotopeDistributionCheck = false;
                Tuple<double, double>[] isotopePeaks = new Tuple<double, double>[isotopeMassShifts.Count];
                int isotopeIndex = 0;
                double theorIsotopeMz = theorIsotopeMzs[isotopeIndex];
                foreach (var possibleIsotopePeak in possibleIsotopePeaks)
                {
                    if (Math.Abs(possibleIsotopePeak.Item1 - theorIsotopeMz) < isotopeMzTol && possibleIsotopePeak.Item2 > 0)
                    {
                        // store the good isotope peak
                        isotopePeaks[isotopeIndex] = possibleIsotopePeak;

                        if (isotopeIndex < isotopePeaks.Length - 1)
                        {
                            // look for the next isotope
                            isotopeIndex++;
                            theorIsotopeMz = theorIsotopeMzs[isotopeIndex];
                        }
                    }
                }

                // all isotopes have been looked for - check to see if they've been observed
                if (isotopePeaks.Where(p => p != null).Count() < NumIsotopesRequired)
                    continue;

                if (isotopePeaks[0] == null && RequireMonoisotopicMass)
                    continue;

                if (RequireMonoisotopicMass)
                {
                    bool[] requiredIsotopesSeen = new bool[NumIsotopesRequired];

                    for (int i = 0; i < NumIsotopesRequired; i++)
                    {
                        if (isotopePeaks[i] != null)
                            requiredIsotopesSeen[i] = true;
                        else
                            requiredIsotopesSeen[i] = false;
                    }

                    if (requiredIsotopesSeen.Where(p => p == false).Any())
                        isotopeDistributionCheck = false;
                    else
                        isotopeDistributionCheck = true;
                }
                else
                {
                    // need to look for continuous isotope distribution when monoisotopic mass observation is not required
                    isotopeDistributionCheck = true;
                }

                if (isotopeDistributionCheck)
                {
                    double isotopeClusterIntensity = 0;

                    if (RequireMonoisotopicMass)
                    {
                        for (int i = 0; i < isotopePeaks.Length; i++)
                        {
                            if (isotopePeaks[i] != null)
                            {
                                double relIsotopeAbundance = isotopePeaks[i].Item2 / isotopePeaks[0].Item2;
                                double theorIsotopeAbundance = isotopeMassShifts[i].Value / isotopeMassShifts[0].Value;

                                // impute isotope intensity if it is very different from expected
                                if ((relIsotopeAbundance / theorIsotopeAbundance) < 2.0)
                                    isotopeClusterIntensity += isotopePeaks[i].Item2;
                                else
                                    isotopeClusterIntensity += theorIsotopeAbundance * isotopePeaks[0].Item2 * 2.0;
                            }
                            else
                                isotopeClusterIntensity += (isotopeMassShifts[i].Value / isotopeMassShifts[0].Value) * isotopePeaks[0].Item2;
                        }
                    }
                    else
                    {
                        isotopeClusterIntensity = isotopePeaks.Where(p => p != null).Sum(p => p.Item2);
                    }

                    isotopeClusters.Add(new IsotopicEnvelope(peak, chargeState, isotopeClusterIntensity, ms1Scans[peak.OneBasedScanNumber - 1].RetentionTime));
                }
            }

            return isotopeClusters;
        }

        private void Peakfind(List<IndexedMassSpectralPeak> possibleMonoisotopicPeaks, Identification identification, int chargeState)
        {
            // sort peaks by scan number
            double theorMz = identification.massToLookFor.ToMz(chargeState);
            HashSet<int> scanNumbers = new HashSet<int>(possibleMonoisotopicPeaks.Select(p => p.OneBasedScanNumber));

            // get precursor scan to start at
            int precursorScanNumber = 0;
            foreach (var ms1Scan in ms1Scans)
            {
                if (ms1Scan == null)
                {
                    continue;
                }

                if (ms1Scan.RetentionTime < identification.ms2RetentionTimeInMinutes)
                {
                    precursorScanNumber = ms1Scan.OneBasedScanNumber;
                }
                else
                {
                    break;
                }
            }
            if (precursorScanNumber == 0)
            {
                throw new MzLibException("FlashLFQ error getting precursor scan number");
            }

            // look right
            int missedScans = 0;
            for (int i = precursorScanNumber - 1; i < ms1Scans.Length; i++)
            {
                // the "ms1Scans[i].ScanWindowRange.Contains(theorMz)" part is for BoxCar compatibility
                if (ms1Scans[i] != null && (ms1Scans[i].ScanWindowRange == null || ms1Scans[i].ScanWindowRange.Contains(theorMz)))
                {
                    if (scanNumbers.Contains(ms1Scans[i].OneBasedScanNumber))
                    {
                        missedScans = 0;
                    }
                    else
                    {
                        if (i != precursorScanNumber - 1)
                        {
                            missedScans++;
                        }

                        if (missedScans > MissedScansAllowed)
                        {
                            possibleMonoisotopicPeaks.RemoveAll(p => p.OneBasedScanNumber >= ms1Scans[i].OneBasedScanNumber);
                            break;
                        }
                    }
                }
            }

            if (!possibleMonoisotopicPeaks.Any())
            {
                return;
            }

            // look left
            missedScans = 0;
            for (int i = precursorScanNumber - 1; i >= 0; i--)
            {
                if (ms1Scans[i] != null && (ms1Scans[i].ScanWindowRange == null || ms1Scans[i].ScanWindowRange.Contains(theorMz)))
                {
                    if (scanNumbers.Contains(ms1Scans[i].OneBasedScanNumber))
                    {
                        missedScans = 0;
                    }
                    else
                    {
                        if (i != precursorScanNumber - 1)
                        {
                            missedScans++;
                        }

                        if (missedScans > MissedScansAllowed)
                        {
                            possibleMonoisotopicPeaks.RemoveAll(p => p.OneBasedScanNumber <= ms1Scans[i].OneBasedScanNumber);
                            break;
                        }
                    }
                }
            }
        }

        private void CutPeak(ChromatographicPeak peak)
        {
            // find out if we need to split this peak by using the discrimination factor
            // this method assumes that the isotope envelopes in a chromatographic peak are already sorted by MS1 scan number
            bool cutThisPeak = false;

            if (peak.IsotopicEnvelopes.Count < 5)
            {
                return;
            }

            var timePointsForApexZ = peak.IsotopicEnvelopes.Where(p => p.ChargeState == peak.Apex.ChargeState).ToList();
            HashSet<int> scanNumbers = new HashSet<int>(timePointsForApexZ.Select(p => p.IndexedPeak.OneBasedScanNumber));
            int apexIndex = timePointsForApexZ.IndexOf(peak.Apex);
            IsotopicEnvelope valleyTimePoint = null;

            // -1 checks the left side, +1 checks the right side
            int[] iters = new int[] { 1, -1 };

            foreach (var iter in iters)
            {
                valleyTimePoint = null;
                int indexOfValley = 0;

                for (int i = apexIndex + iter; i < timePointsForApexZ.Count && i >= 0; i += iter)
                {
                    IsotopicEnvelope timepoint = timePointsForApexZ[i];

                    if (valleyTimePoint == null || timepoint.Intensity < valleyTimePoint.Intensity)
                    {
                        valleyTimePoint = timepoint;
                        indexOfValley = timePointsForApexZ.IndexOf(valleyTimePoint);
                    }

                    double discriminationFactor = (timepoint.Intensity - valleyTimePoint.Intensity) / timepoint.Intensity;

                    if (discriminationFactor > MinDiscFactorToCutAt && (indexOfValley + iter < timePointsForApexZ.Count && indexOfValley + iter >= 0))
                    {
                        IsotopicEnvelope secondValleyTimepoint = timePointsForApexZ[indexOfValley + iter];

                        discriminationFactor = (timepoint.Intensity - secondValleyTimepoint.Intensity) / timepoint.Intensity;

                        if (discriminationFactor > MinDiscFactorToCutAt)
                        {
                            cutThisPeak = true;
                            break;
                        }

                        int nextMs1ScanNum = -1;
                        for (int j = valleyTimePoint.IndexedPeak.OneBasedScanNumber - 1; j < ms1Scans.Length && j >= 0; j += iter)
                        {
                            if (ms1Scans[j] != null && ms1Scans[j].OneBasedScanNumber != valleyTimePoint.IndexedPeak.OneBasedScanNumber)
                            {
                                nextMs1ScanNum = j + 1;
                                break;
                            }
                        }

                        if (!scanNumbers.Contains(nextMs1ScanNum))
                        {
                            cutThisPeak = true;
                            break;
                        }
                    }
                }

                if (cutThisPeak == true)
                {
                    break;
                }
            }

            // cut
            if (cutThisPeak)
            {
                if (peak.Identifications.First().ms2RetentionTimeInMinutes > valleyTimePoint.RetentionTime)
                {
                    // MS2 identification is to the right of the valley; remove all peaks left of the valley
                    peak.IsotopicEnvelopes.RemoveAll(p => p.RetentionTime <= valleyTimePoint.RetentionTime);
                }
                else
                {
                    // MS2 identification is to the left of the valley; remove all peaks right of the valley
                    peak.IsotopicEnvelopes.RemoveAll(p => p.RetentionTime >= valleyTimePoint.RetentionTime);
                }

                // recalculate intensity for the peak
                peak.CalculateIntensityForThisFeature(Integrate);
                peak.SplitRT = valleyTimePoint.RetentionTime;

                // recursively cut
                CutPeak(peak);
            }
        }
    }
}