﻿using Chemistry;
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
        #region Public Fields

        // settings
        public readonly bool silent;
        public readonly int maxThreads;
        public readonly double peakfindingPpmTolerance;
        public readonly double ppmTolerance;
        public readonly double rtTol;
        public readonly double isotopePpmTolerance;
        public readonly bool integrate;
        public readonly int missedScansAllowed;
        public readonly int numIsotopesRequired;
        public readonly double mbrRtWindow;
        public readonly double initialMbrRtWindow;
        public readonly double mbrppmTolerance;
        public readonly bool errorCheckAmbiguousMatches;
        public readonly bool mbr;
        public readonly bool idSpecificChargeState;
        public readonly double qValueCutoff;
        public readonly bool requireMonoisotopicMass;
        public readonly bool normalize;
        public readonly double minDiscFactorToCutAt;

        #endregion Public Fields

        #region Private Fields

        private List<SpectraFileInfo> spectraFileInfo;
        private Stopwatch globalStopwatch;
        private List<Identification> allIdentifications;
        private Dictionary<string, List<KeyValuePair<double, double>>> baseSequenceToIsotopicDistribution;
        private IEnumerable<int> chargeStates;
        private FlashLFQResults results;
        private int binsPerDalton = 100;

        // these two fields will be overwritten as each file is analyzed
        private MsDataScan[] ms1Scans;
        private List<IndexedMassSpectralPeak>[] indexedPeaks;

        #endregion Private Fields

        #region Public Constructors

        public FlashLFQEngine(List<Identification> allIdentifications, bool normalize = true, double ppmTolerance = 10.0, double isotopeTolerancePpm = 5.0, bool matchBetweenRuns = false, double matchBetweenRunsPpmTolerance = 5.0, bool integrate = false, int numIsotopesRequired = 2, bool idSpecificChargeState = false, bool requireMonoisotopicMass = true, bool silent = false, string optionalPeriodicTablePath = null, double maxMbrWindow = 1.5)
        {
            if (optionalPeriodicTablePath == null)
                optionalPeriodicTablePath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"elements.dat");
            Loaders.LoadElements(optionalPeriodicTablePath);

            globalStopwatch = new Stopwatch();
            chargeStates = new List<int>();

            this.spectraFileInfo = allIdentifications.Select(p => p.fileInfo).Distinct()
                .OrderBy(p => p.condition)
                .ThenBy(p => p.biologicalReplicate)
                .ThenBy(p => p.fraction)
                .ThenBy(p => p.technicalReplicate).ToList();

            this.allIdentifications = allIdentifications;
            this.ppmTolerance = ppmTolerance;
            this.isotopePpmTolerance = isotopeTolerancePpm;
            this.mbr = matchBetweenRuns;
            this.mbrppmTolerance = matchBetweenRunsPpmTolerance;
            this.integrate = integrate;
            this.numIsotopesRequired = numIsotopesRequired;
            this.silent = silent;
            this.idSpecificChargeState = idSpecificChargeState;
            this.requireMonoisotopicMass = requireMonoisotopicMass;
            this.mbrRtWindow = maxMbrWindow;
            this.normalize = normalize;

            qValueCutoff = 0.01;
            peakfindingPpmTolerance = 20.0;
            initialMbrRtWindow = 10.0;
            missedScansAllowed = 1;
            rtTol = 5.0;
            errorCheckAmbiguousMatches = true;
            maxThreads = -1;
            minDiscFactorToCutAt = 0.6;
        }

        #endregion Public Constructors

        #region Public Methods

        public FlashLFQResults Run()
        {
            globalStopwatch.Start();

            results = new FlashLFQResults(spectraFileInfo);

            // build m/z index keys
            CalculateTheoreticalIsotopeDistributions();

            // quantify each file
            foreach (var spectraFile in spectraFileInfo)
            {
                // fill lookup-table with peaks from the raw file
                IndexMassSpectralPeaks(spectraFile);

                // quantify features using this file's IDs first
                QuantifyMS2IdentifiedPeptides(spectraFile);

                // find unidentified features based on other files' identification results (initial MBR peak-finding)
                if (mbr)
                {
                    MatchBetweenRunsInitialPeakfinding(spectraFile);
                }

                // error checking function
                // handles features with multiple identifying scans and scans that are associated with more than one feature
                RunErrorChecking(spectraFile);

                if (!silent)
                {
                    Console.WriteLine("Finished " + spectraFile.filenameWithoutExtension);
                }

                // some memory-saving stuff
                ms1Scans = new MsDataScan[0];
                GC.Collect();
            }

            // filter initial MBR peaks with retention time calibration
            if (mbr)
            {
                RetentionTimeCalibrationAndErrorCheckMatchedFeatures();
            }

            // normalize
            if (normalize)
            {
                new IntensityNormalizationEngine(results, integrate, silent).NormalizeResults();
                
                //new StatisticalAnalysisEngine(results, 0.05, 0.1).PerformStatisticalAnalysis();
            }

            // calculate intensities for proteins/peptides
            results.CalculatePeptideResults(true);
            results.CalculatePeptideResults(false);
            results.CalculateProteinResults();

            // done
            if (!silent)
                Console.WriteLine("All done");

            if (!silent)
                Console.WriteLine("Analysis time: " +
                    globalStopwatch.Elapsed.Hours + "h " +
                    globalStopwatch.Elapsed.Minutes + "m " +
                    globalStopwatch.Elapsed.Seconds + "s");

            return results;
        }

        #endregion Public Methods

        #region Private Methods

        private void RetentionTimeCalibrationAndErrorCheckMatchedFeatures()
        {
            if (!silent)
                Console.WriteLine("Running retention time calibration");

            // get all unambiguous peaks for all files
            var allFeatures = results.peaks.SelectMany(p => p.Value).Where(p => !p.isMbrFeature);
            var allAmbiguousFeatures = allFeatures.Where(p => p.NumIdentificationsByFullSeq > 1).ToList();
            var ambiguousFeatureSeqs = new HashSet<string>(allAmbiguousFeatures.SelectMany(p => p.identifications.Select(v => v.ModifiedSequence)));

            foreach (var feature in allFeatures)
                if (ambiguousFeatureSeqs.Contains(feature.identifications.First().ModifiedSequence))
                    allAmbiguousFeatures.Add(feature);

            var unambiguousPeaksGroupedByFile = allFeatures.Except(allAmbiguousFeatures).Where(v => v.apex != null).GroupBy(p => p.rawFileInfo);

            foreach (var file in unambiguousPeaksGroupedByFile)
            {
                var allMbrFeaturesForThisFile = results.peaks[file.Key].Where(p => p.isMbrFeature);

                // get the best (most intense) peak for each peptide in the file
                Dictionary<string, ChromatographicPeak> pepToBestFeatureForThisFile = new Dictionary<string, ChromatographicPeak>();
                foreach (var testPeak in file)
                {
                    if (pepToBestFeatureForThisFile.TryGetValue(testPeak.identifications.First().ModifiedSequence, out ChromatographicPeak currentBestPeak))
                    {
                        if (currentBestPeak.intensity > testPeak.intensity)
                            pepToBestFeatureForThisFile[testPeak.identifications.First().ModifiedSequence] = testPeak;
                    }
                    else
                        pepToBestFeatureForThisFile.Add(testPeak.identifications.First().ModifiedSequence, testPeak);
                }

                foreach (var otherFile in unambiguousPeaksGroupedByFile)
                {
                    // get the other files' best peak for the same peptides (to make an RT calibration curve)
                    if (otherFile.Key.Equals(file.Key))
                        continue;

                    var featuresInCommon = otherFile.Where(p => pepToBestFeatureForThisFile.ContainsKey(p.identifications.First().ModifiedSequence));

                    Dictionary<string, ChromatographicPeak> pepToBestFeatureForOtherFile = new Dictionary<string, ChromatographicPeak>();
                    foreach (var testPeak in featuresInCommon)
                    {
                        if (pepToBestFeatureForOtherFile.TryGetValue(testPeak.identifications.First().ModifiedSequence, out ChromatographicPeak currentBestPeak))
                        {
                            if (currentBestPeak.intensity > testPeak.intensity)
                                pepToBestFeatureForOtherFile[testPeak.identifications.First().ModifiedSequence] = testPeak;
                        }
                        else
                            pepToBestFeatureForOtherFile.Add(testPeak.identifications.First().ModifiedSequence, testPeak);
                    }

                    // create a rt-to-rt correlation for the two files' peptides
                    Dictionary<string, Tuple<double, double>> rtCalPoints = new Dictionary<string, Tuple<double, double>>();

                    foreach (var kvp in pepToBestFeatureForOtherFile)
                        rtCalPoints.Add(kvp.Key, new Tuple<double, double>(pepToBestFeatureForThisFile[kvp.Key].apex.retentionTime, kvp.Value.apex.retentionTime));

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

                                if (rtError > (mbrRtWindow / 2.0))
                                    stdevRunningSpline[i] = mbrRtWindow / 2.0;
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
                                stdevRunningSpline[i] = stdevRunningSpline[i] = mbrRtWindow / 2.0;
                            }
                        }
                    }

                    // finished rt calibration for these 2 files; now use rt cal spline to find matched features
                    var allMatchedFeaturesToLookForNow = allMbrFeaturesForThisFile.Where(p => p.identifications.First().fileInfo.Equals(otherFile.Key)).ToList();

                    // filter peak candidates with rt cal to get apex peak
                    foreach (var mbrFeature in allMatchedFeaturesToLookForNow)
                    {
                        if (mbrFeature.isotopicEnvelopes.Any())
                        {
                            // shift = thisFileRt - otherFileRt
                            int rtSplineLookupTime = (int)Math.Round(mbrFeature.identifications.First().ms2RetentionTimeInMinutes);

                            if (rtSplineLookupTime < rtCalRunningSpline.Length)
                            {
                                double rtShift = rtCalRunningSpline[rtSplineLookupTime];
                                double rtToleranceHere = stdevRunningSpline[rtSplineLookupTime];
                                double theoreticalRt = mbrFeature.identifications.First().ms2RetentionTimeInMinutes + rtShift;

                                if (!double.IsNaN(rtShift))
                                    mbrFeature.isotopicEnvelopes = mbrFeature.isotopicEnvelopes.Where(p => Math.Abs(p.retentionTime - theoreticalRt) < rtToleranceHere).ToList();
                                else
                                    mbrFeature.isotopicEnvelopes = new List<IsotopicEnvelope>();
                            }
                            else
                                mbrFeature.isotopicEnvelopes = new List<IsotopicEnvelope>();
                        }
                    }

                    foreach (var feature in allMatchedFeaturesToLookForNow)
                        if (feature.isotopicEnvelopes.Any())
                            feature.CalculateIntensityForThisFeature(integrate);
                }
            }

            foreach (var file in spectraFileInfo)
                RunErrorChecking(file);
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

                Proteomics.Peptide p = new Proteomics.Peptide(id.BaseSequence);
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
                    if (i < (numIsotopesRequired - 1) || abundances[i] > 0.1)
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
                    identification.massToLookFor = requireMonoisotopicMass ? identifications.First().monoisotopicMass : thisPeptidesMostAbundantMass;
                }
            }
        }

        private void IndexMassSpectralPeaks(SpectraFileInfo fileInfo)
        {
            if (!silent)
            {
                Console.WriteLine("Reading spectra file");
            }

            ms1Scans = new MsDataScan[0];

            if(indexedPeaks != null)
            {
                for(int i = 0; i < indexedPeaks.Length; i++)
                {
                    if (indexedPeaks[i] == null)
                        continue;

                    indexedPeaks[i].Clear();
                }
            }

            // read spectra file
            var ext = Path.GetExtension(fileInfo.fullFilePathWithExtension).ToUpperInvariant();
            if (ext == ".MZML")
            {
                try
                {
                    ms1Scans = Mzml.LoadAllStaticData(fileInfo.fullFilePathWithExtension).GetAllScansList().OrderBy(p => p.OneBasedScanNumber).ToArray();
                }
                catch (FileNotFoundException)
                {
                    if (!silent)
                    {
                        Console.WriteLine("\nCan't find mzml file" + fileInfo.fullFilePathWithExtension + "\n");
                    }
                    return;
                }
                catch (Exception e)
                {
                    if (!silent)
                    {
                        Console.WriteLine("Problem opening mzml file " + fileInfo.fullFilePathWithExtension + "; " + e.Message);
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
                using (var thermoDynamicConnection = IO.Thermo.ThermoDynamicData.InitiateDynamicConnection(fileInfo.fullFilePathWithExtension))
                {
                    var tempList = new List<MsDataScan>();

                    try
                    {
                        // use thermo dynamic connection to get the ms1 scans and then dispose of the connection
                        int[] msOrders = thermoDynamicConnection.ThermoGlobalParams.msOrderByScan;
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

                        if (!silent)
                        {
                            Console.WriteLine("\nCan't find raw file" + fileInfo.fullFilePathWithExtension + "\n");
                        }
                        return;
                    }
                    catch (Exception e)
                    {
                        thermoDynamicConnection.Dispose();

                        if (!silent)
                        {
                            throw new MzLibException("FlashLFQ Error: Problem opening raw file " + fileInfo.fullFilePathWithExtension + "; " + e.Message);
                        }
                    }

                    ms1Scans = tempList.ToArray();
                }
#else
                if (!silent)
                {
                    Console.WriteLine("Cannot open RAW with .NETStandard code - are you on Linux? " + fileInfo.fullFilePathWithExtension);
                }
                return;
#endif
            }
            else
            {
                if (!silent)
                {
                    Console.WriteLine("Unsupported file type " + ext);
                    return;
                }
            }

            if (!silent)
            {
                Console.WriteLine("Indexing MS1 peaks");
            }

            indexedPeaks = new List<IndexedMassSpectralPeak>[(int)Math.Ceiling(ms1Scans.Where(p => p != null).Max(p => p.MassSpectrum.LastX.Value) * binsPerDalton) + 1];

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
            if (!silent)
            {
                Console.WriteLine("Quantifying peptides for " + fileInfo.filenameWithoutExtension);
            }

            var identifications = allIdentifications.Where(p => p.fileInfo.Equals(fileInfo)).ToList();

            if (!identifications.Any())
            {
                return;
            }

            Tolerance peakfindingTol = new PpmTolerance(peakfindingPpmTolerance);
            Tolerance tol = new PpmTolerance(ppmTolerance);

            var chromatographicPeaks = new ChromatographicPeak[identifications.Count];

            Parallel.ForEach(Partitioner.Create(0, identifications.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                range =>
            {
                List<IndexedMassSpectralPeak> binPeaks = new List<IndexedMassSpectralPeak>();
                List<IsotopicEnvelope> isotopicEnvelopes = new List<IsotopicEnvelope>();

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    binPeaks.Clear();
                    isotopicEnvelopes.Clear();
                    var identification = identifications[i];
                    ChromatographicPeak msmsFeature = new ChromatographicPeak(identification, false, fileInfo);
                    chromatographicPeaks[i] = msmsFeature;

                    foreach (var chargeState in chargeStates)
                    {
                        if (idSpecificChargeState && chargeState != identification.precursorChargeState)
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
                                    if (peakfindingTol.Within(peak.mz.ToMass(chargeState), identification.massToLookFor) &&
                                      Math.Abs(ms1Scans[peak.oneBasedScanNumber - 1].RetentionTime - identification.ms2RetentionTimeInMinutes) < rtTol)
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
                            binPeaks.RemoveAll(p => !tol.Within(p.mz.ToMass(chargeState), identification.massToLookFor));

                            // filter by isotopic distribution
                            isotopicEnvelopes = GetIsotopicEnvelopes(binPeaks, identification, chargeState, true);

                            // if multiple mass spectral peaks in the same scan are valid, pick the one with the smallest mass error
                            var peaksInSameScan = isotopicEnvelopes.GroupBy(p => p.indexedPeak.oneBasedScanNumber).Where(v => v.Count() > 1);
                            if (peaksInSameScan.Any())
                            {
                                foreach (var group in peaksInSameScan)
                                {
                                    var smallestMzError = group.Min(p => Math.Abs(p.indexedPeak.mz - theoreticalMz));
                                    var peakToUse = group.Where(p => Math.Abs(p.indexedPeak.mz - theoreticalMz) == smallestMzError).First();
                                    isotopicEnvelopes.RemoveAll(p => group.Contains(p) && p != peakToUse);
                                }
                            }

                            // sort isotope envelopes by scan number (sets up for peak cutting)
                            msmsFeature.isotopicEnvelopes.AddRange(isotopicEnvelopes);
                            msmsFeature.isotopicEnvelopes.Sort((x, y) => x.indexedPeak.oneBasedScanNumber.CompareTo(y.indexedPeak.oneBasedScanNumber));
                        }
                    }

                    msmsFeature.CalculateIntensityForThisFeature(integrate);
                    CutPeak(msmsFeature);
                }
            });

            results.peaks.Add(fileInfo, chromatographicPeaks.ToList());
        }

        private void MatchBetweenRunsInitialPeakfinding(SpectraFileInfo fileInfo)
        {
            if (!silent)
                Console.WriteLine("Finding possible matched peptides for " + fileInfo.filenameWithoutExtension);

            if (!results.peaks.ContainsKey(fileInfo) || results.peaks[fileInfo].Count == 0)
                return;

            Tolerance mbrTol = new PpmTolerance(mbrppmTolerance);

            var concurrentBagOfMatchedFeatures = new ConcurrentBag<ChromatographicPeak>();
            var identificationsFromOtherRunsToLookFor = new List<Identification>();
            var idsGroupedByFullSeq = allIdentifications.GroupBy(p => p.ModifiedSequence);

            foreach (var fullSequenceGroup in idsGroupedByFullSeq)
            {
                // look for peptides with no ID's in this file
                var seqsByFilename = fullSequenceGroup.GroupBy(p => p.fileInfo);

                if (!seqsByFilename.Where(p => p.Key.Equals(fileInfo)).Any())
                    identificationsFromOtherRunsToLookFor.AddRange(fullSequenceGroup);
            }

            if (identificationsFromOtherRunsToLookFor.Any())
            {
                Parallel.ForEach(Partitioner.Create(0, identificationsFromOtherRunsToLookFor.Count),
                    new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                    (range, loopState) =>
                    {
                        for (int i = range.Item1; i < range.Item2; i++)
                        {
                            var identification = identificationsFromOtherRunsToLookFor[i];

                            ChromatographicPeak mbrFeature = new ChromatographicPeak(identification, true, fileInfo);

                            foreach (var chargeState in chargeStates)
                            {
                                double theorMzHere = identification.massToLookFor.ToMz(chargeState);
                                double mzTolHere = ((mbrppmTolerance / 1e6) * identification.monoisotopicMass) / chargeState;

                                int floorMz = (int)Math.Floor(theorMzHere * binsPerDalton);
                                int ceilingMz = (int)Math.Ceiling(theorMzHere * binsPerDalton);

                                List<IndexedMassSpectralPeak> binPeaks = new List<IndexedMassSpectralPeak>();

                                for (int j = floorMz; j <= ceilingMz && j < indexedPeaks.Length; j++)
                                {
                                    if (indexedPeaks[j] != null)
                                    {
                                        foreach (var peak in indexedPeaks[j])
                                        {
                                            if (mbrTol.Within(peak.mz.ToMass(chargeState), identification.massToLookFor) &&
                                              Math.Abs(ms1Scans[peak.oneBasedScanNumber - 1].RetentionTime - identification.ms2RetentionTimeInMinutes) < initialMbrRtWindow)
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
                                    mbrFeature.isotopicEnvelopes.AddRange(validIsotopeClusters);
                                }
                            }

                            if (mbrFeature.isotopicEnvelopes.Any())
                            {
                                concurrentBagOfMatchedFeatures.Add(mbrFeature);
                            }
                        }
                    }
                );
            }

            results.peaks[fileInfo].AddRange(concurrentBagOfMatchedFeatures);
        }

        private void RunErrorChecking(SpectraFileInfo rawFile)
        {
            results.peaks[rawFile].RemoveAll(p => p.isMbrFeature && !p.isotopicEnvelopes.Any());

            if (!silent)
                Console.WriteLine("Checking errors");
            var featuresWithSamePeak = results.peaks[rawFile].Where(v => v.apex != null).GroupBy(p => p.apex.indexedPeak);
            featuresWithSamePeak = featuresWithSamePeak.Where(p => p.Count() > 1);

            // condense duplicate features (features with same sequence and apex peak)
            foreach (var duplicateFeature in featuresWithSamePeak)
            {
                duplicateFeature.First().MergeFeatureWith(duplicateFeature, integrate);
            }
            results.peaks[rawFile].RemoveAll(p => p.intensity == -1);

            // check for multiple features per peptide within a time window
            var featuresToMaybeMerge = results.peaks[rawFile].Where(p => p.NumIdentificationsByFullSeq == 1 && p.apex != null).GroupBy(p => p.identifications.First().ModifiedSequence).Where(p => p.Count() > 1);
            if (featuresToMaybeMerge.Any())
            {
                foreach (var group in featuresToMaybeMerge)
                {
                    if (idSpecificChargeState)
                    {
                        var group2 = group.ToList().GroupBy(p => p.apex.chargeState).Where(v => v.Count() > 1);

                        foreach (var group3 in group2)
                        {
                            foreach (var feature in group3)
                            {
                                if (feature.intensity != -1)
                                {
                                    var featuresToMerge = group3.Where(p => Math.Abs(p.apex.retentionTime - feature.apex.retentionTime) < rtTol && p.intensity != -1);
                                    if (featuresToMerge.Any())
                                        feature.MergeFeatureWith(featuresToMerge, integrate);
                                }
                            }
                        }
                    }
                    else
                    {
                        foreach (var feature in group)
                        {
                            if (feature.intensity != -1)
                            {
                                var featuresToMerge = group.Where(p => Math.Abs(p.apex.retentionTime - feature.apex.retentionTime) < rtTol && p.intensity != -1);
                                if (featuresToMerge.Any())
                                    feature.MergeFeatureWith(featuresToMerge, integrate);
                            }
                        }
                    }
                }

                results.peaks[rawFile].RemoveAll(p => p.intensity == -1);
            }

            if (errorCheckAmbiguousMatches)
            {
                // check for multiple peptides per feature
                var scansWithMultipleDifferentIds = results.peaks[rawFile].Where(p => p.NumIdentificationsByFullSeq > 1);
                var ambiguousFeatures = scansWithMultipleDifferentIds.Where(p => p.NumIdentificationsByBaseSeq > 1).ToList();

                // handle ambiguous features
                foreach (var ambiguousFeature in ambiguousFeatures)
                {
                    var msmsIdentsForThisFile = ambiguousFeature.identifications.Where(p => p.fileInfo.Equals(ambiguousFeature.rawFileInfo));

                    if (!msmsIdentsForThisFile.Any())
                    {
                        // mbr matched more than one identification to this peak - cannot resolve
                        ambiguousFeature.intensity = -1;
                    }
                    else
                    {
                        // msms identifications take precident over mbr features
                        ambiguousFeature.identifications.RemoveAll(p => p.fileInfo != ambiguousFeature.rawFileInfo);
                        ambiguousFeature.ResolveIdentifications();
                    }
                }

                results.peaks[rawFile].RemoveAll(p => p.intensity == -1);
            }
        }

        private List<IsotopicEnvelope> GetIsotopicEnvelopes(List<IndexedMassSpectralPeak> peaks, Identification identification, int chargeState, bool lookForBadIsotope)
        {
            var isotopeClusters = new List<IsotopicEnvelope>();
            var isotopeMassShifts = baseSequenceToIsotopicDistribution[identification.BaseSequence];

            if (isotopeMassShifts.Count < numIsotopesRequired)
            {
                return isotopeClusters;
            }

            double isotopeMzTol = ((isotopePpmTolerance / 1e6) * identification.monoisotopicMass) / chargeState;

            foreach (var peak in peaks)
            {
                // calculate theoretical isotopes relative to observed peak
                var theorIsotopeMzs = new double[isotopeMassShifts.Count];
                int isotopicPeakUnitOfPeakZeroIsMono = Convert.ToInt32(peak.mz.ToMass(chargeState) - identification.monoisotopicMass);
                var mainpeakMz = peak.mz;

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
                int earliestIsotopicPeakIndexPossible = peak.zeroBasedIndexOfPeakInScan;
                var massSpectrum = ms1Scans[peak.oneBasedScanNumber - 1].MassSpectrum;
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
                            if (massSpectrum.YArray[i] / peak.intensity > 1.0)
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
                if (isotopePeaks.Where(p => p != null).Count() < numIsotopesRequired)
                    continue;

                if (isotopePeaks[0] == null && requireMonoisotopicMass)
                    continue;

                if (requireMonoisotopicMass)
                {
                    bool[] requiredIsotopesSeen = new bool[numIsotopesRequired];

                    for (int i = 0; i < numIsotopesRequired; i++)
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

                    if (requireMonoisotopicMass)
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

                    isotopeClusters.Add(new IsotopicEnvelope(peak, chargeState, isotopeClusterIntensity, ms1Scans[peak.oneBasedScanNumber - 1].RetentionTime));
                }
            }

            return isotopeClusters;
        }

        private void Peakfind(List<IndexedMassSpectralPeak> possibleMonoisotopicPeaks, Identification identification, int chargeState)
        {
            // sort peaks by scan number
            double theorMz = identification.massToLookFor.ToMz(chargeState);
            HashSet<int> scanNumbers = new HashSet<int>(possibleMonoisotopicPeaks.Select(p => p.oneBasedScanNumber));

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

                        if (missedScans > missedScansAllowed)
                        {
                            possibleMonoisotopicPeaks.RemoveAll(p => p.oneBasedScanNumber >= ms1Scans[i].OneBasedScanNumber);
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

                        if (missedScans > missedScansAllowed)
                        {
                            possibleMonoisotopicPeaks.RemoveAll(p => p.oneBasedScanNumber <= ms1Scans[i].OneBasedScanNumber);
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

            if (peak.isotopicEnvelopes.Count < 5)
            {
                return;
            }

            var timePointsForApexZ = peak.isotopicEnvelopes.Where(p => p.chargeState == peak.apex.chargeState).ToList();
            HashSet<int> scanNumbers = new HashSet<int>(timePointsForApexZ.Select(p => p.indexedPeak.oneBasedScanNumber));
            int apexIndex = timePointsForApexZ.IndexOf(peak.apex);
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
                    
                    if (valleyTimePoint == null || timepoint.intensity < valleyTimePoint.intensity)
                    {
                        valleyTimePoint = timepoint;
                        indexOfValley = timePointsForApexZ.IndexOf(valleyTimePoint);
                    }

                    double discriminationFactor = (timepoint.intensity - valleyTimePoint.intensity) / timepoint.intensity;

                    if (discriminationFactor > minDiscFactorToCutAt && (indexOfValley + iter < timePointsForApexZ.Count && indexOfValley + iter >= 0))
                    {
                        IsotopicEnvelope secondValleyTimepoint = timePointsForApexZ[indexOfValley + iter];

                        discriminationFactor = (timepoint.intensity - secondValleyTimepoint.intensity) / timepoint.intensity;

                        if (discriminationFactor > minDiscFactorToCutAt)
                        {
                            cutThisPeak = true;
                            break;
                        }

                        int nextMs1ScanNum = -1;
                        for (int j = valleyTimePoint.indexedPeak.oneBasedScanNumber - 1; j < ms1Scans.Length && j >= 0; j += iter)
                        {
                            if (ms1Scans[j] != null && ms1Scans[j].OneBasedScanNumber != valleyTimePoint.indexedPeak.oneBasedScanNumber)
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
                if (peak.identifications.First().ms2RetentionTimeInMinutes > valleyTimePoint.retentionTime)
                {
                    // MS2 identification is to the right of the valley; remove all peaks left of the valley
                    peak.isotopicEnvelopes.RemoveAll(p => p.retentionTime <= valleyTimePoint.retentionTime);
                }
                else
                {
                    // MS2 identification is to the left of the valley; remove all peaks right of the valley
                    peak.isotopicEnvelopes.RemoveAll(p => p.retentionTime >= valleyTimePoint.retentionTime);
                }

                // recalculate intensity for the peak
                peak.CalculateIntensityForThisFeature(integrate);
                peak.splitRT = valleyTimePoint.retentionTime;

                // recursively cut
                CutPeak(peak);
            }
        }

        #endregion Private Methods
    }
}