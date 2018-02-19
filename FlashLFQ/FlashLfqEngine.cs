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

        #endregion Public Fields

        #region Private Fields

        private List<RawFileInfo> rawFileInformation;
        private Stopwatch globalStopwatch;
        private Stopwatch fileLocalStopwatch;
        private List<Identification> allIdentifications;
        private HashSet<double> indexedMzKeys;
        private Dictionary<string, List<KeyValuePair<double, double>>> baseSequenceToIsotopicDistribution;
        private Dictionary<string, ProteinGroup> proteinGroupNameToProteinGroup;
        private IEnumerable<int> chargeStates;
        private FlashLFQResults results;

        #endregion Private Fields

        #region Public Constructors

        public FlashLFQEngine(List<Identification> allIdentifications, double ppmTolerance = 10.0, double isotopeTolerancePpm = 5.0, bool matchBetweenRuns = false, double matchBetweenRunsPpmTolerance = 5.0, bool integrate = false, int numIsotopesRequired = 2, bool idSpecificChargeState = false, bool requireMonoisotopicMass = true, bool silent = false, string optionalPeriodicTablePath = null, double maxMbrWindow = 1.5)
        {
            if (optionalPeriodicTablePath == null)
                optionalPeriodicTablePath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"elements.dat");
            PeriodicTableLoader.Load(optionalPeriodicTablePath);

            globalStopwatch = new Stopwatch();
            fileLocalStopwatch = new Stopwatch();
            chargeStates = new List<int>();
            proteinGroupNameToProteinGroup = new Dictionary<string, ProteinGroup>();

            this.rawFileInformation = allIdentifications.Select(p => p.fileInfo).Distinct().ToList();
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

            qValueCutoff = 0.01;
            peakfindingPpmTolerance = 20.0;
            initialMbrRtWindow = 10.0;
            missedScansAllowed = 1;
            rtTol = 5.0;
            errorCheckAmbiguousMatches = true;
            maxThreads = -1;

            ProteinGroup.rawFiles = rawFileInformation;
            Peptide.rawFiles = rawFileInformation;
            FlashLFQResults.rawFiles = rawFileInformation;
        }

        #endregion Public Constructors

        #region Public Methods

        public FlashLFQResults Run()
        {
            results = new FlashLFQResults();

            globalStopwatch.Start();

            // construct protein groups (intensities written later)
            foreach (var pg in allIdentifications.SelectMany(p => p.proteinGroupNames).Distinct())
                proteinGroupNameToProteinGroup.Add(pg, new ProteinGroup(pg));
            foreach (var id in allIdentifications)
                foreach (var pg in id.proteinGroupNames)
                    id.proteinGroups.Add(proteinGroupNameToProteinGroup[pg]);
            results.proteinGroups = proteinGroupNameToProteinGroup;

            // construct peptides (intensities written later)
            Dictionary<string, Peptide> baseSequences = new Dictionary<string, Peptide>();
            foreach (var baseSeq in allIdentifications.Select(p => p.BaseSequence).Distinct())
                baseSequences.Add(baseSeq, new Peptide(baseSeq));
            results.peptideBaseSequences = baseSequences;

            Dictionary<string, Peptide> modifiedSequences = new Dictionary<string, Peptide>();
            foreach (var modifiedSeq in allIdentifications.Select(p => p.ModifiedSequence).Distinct())
                modifiedSequences.Add(modifiedSeq, new Peptide(modifiedSeq));
            results.peptideModifiedSequences = modifiedSequences;

            // build m/z index keys
            ConstructIndexKeysFromIdentifications();

            // quantify each file
            foreach (var fileInfo in rawFileInformation)
            {
                GC.Collect();

                // fill lookup-table with peaks from the raw file
                var indexedMassSpectralPeaks = IndexMassSpectralPeaks(fileInfo, out Dictionary<int, IMsDataScan<IMzSpectrum<IMzPeak>>> allMs1Scans);

                // quantify features using this file's IDs first
                QuantifyMS2IdentifiedPeptides(fileInfo, indexedMassSpectralPeaks, allMs1Scans);

                // find unidentified features based on other files' identification results (initial MBR peak-finding)
                if (mbr)
                    MatchBetweenRunsInitialPeakfinding(fileInfo, indexedMassSpectralPeaks, allMs1Scans);

                // error checking function
                // handles features with multiple identifying scans and scans that are associated with more than one feature
                RunErrorChecking(fileInfo);

                if (!silent)
                    Console.WriteLine("Finished " + fileInfo.filenameWithoutExtension);

                // some memory-saving stuff
                if (fileInfo.clearAfterDone)
                    fileInfo.dataFile = null;
                allMs1Scans = null;

                fileInfo.analysisSummary = "File analysis time = " + fileLocalStopwatch.Elapsed.ToString();
            }

            // filter initial MBR peaks with retention time calibration
            if (mbr)
                RetentionTimeCalibrationAndErrorCheckMatchedFeatures();

            // sum peak intensities for peptides and proteins
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
            var ambiguousFeatureSeqs = new HashSet<string>(allAmbiguousFeatures.SelectMany(p => p.identifyingScans.Select(v => v.ModifiedSequence)));

            foreach (var feature in allFeatures)
                if (ambiguousFeatureSeqs.Contains(feature.identifyingScans.First().ModifiedSequence))
                    allAmbiguousFeatures.Add(feature);

            var unambiguousPeaksGroupedByFile = allFeatures.Except(allAmbiguousFeatures).Where(v => v.apexPeak != null).GroupBy(p => p.rawFileInfo);

            foreach (var file in unambiguousPeaksGroupedByFile)
            {
                // used later
                int fileIndex = 0;
                for (int i = 0; i < rawFileInformation.Count; i++)
                {
                    if (rawFileInformation[i].Equals(file.Key))
                    {
                        fileIndex = i;
                        break;
                    }
                }

                var allMbrFeaturesForThisFile = results.peaks[file.Key].Where(p => p.isMbrFeature);

                // get the best (most intense) peak for each peptide in the file
                Dictionary<string, ChromatographicPeak> pepToBestFeatureForThisFile = new Dictionary<string, ChromatographicPeak>();
                foreach (var testPeak in file)
                {
                    if (pepToBestFeatureForThisFile.TryGetValue(testPeak.identifyingScans.First().ModifiedSequence, out ChromatographicPeak currentBestPeak))
                    {
                        if (currentBestPeak.intensity > testPeak.intensity)
                            pepToBestFeatureForThisFile[testPeak.identifyingScans.First().ModifiedSequence] = testPeak;
                    }
                    else
                        pepToBestFeatureForThisFile.Add(testPeak.identifyingScans.First().ModifiedSequence, testPeak);
                }

                foreach (var otherFile in unambiguousPeaksGroupedByFile)
                {
                    // get the other files' best peak for the same peptides (to make an RT calibration curve)
                    if (otherFile.Key.Equals(file.Key))
                        continue;

                    var featuresInCommon = otherFile.Where(p => pepToBestFeatureForThisFile.ContainsKey(p.identifyingScans.First().ModifiedSequence));

                    Dictionary<string, ChromatographicPeak> pepToBestFeatureForOtherFile = new Dictionary<string, ChromatographicPeak>();
                    foreach (var testPeak in featuresInCommon)
                    {
                        if (pepToBestFeatureForOtherFile.TryGetValue(testPeak.identifyingScans.First().ModifiedSequence, out ChromatographicPeak currentBestPeak))
                        {
                            if (currentBestPeak.intensity > testPeak.intensity)
                                pepToBestFeatureForOtherFile[testPeak.identifyingScans.First().ModifiedSequence] = testPeak;
                        }
                        else
                            pepToBestFeatureForOtherFile.Add(testPeak.identifyingScans.First().ModifiedSequence, testPeak);
                    }

                    // create a rt-to-rt correlation for the two files' peptides
                    Dictionary<string, Tuple<double, double>> rtCalPoints = new Dictionary<string, Tuple<double, double>>();

                    foreach (var kvp in pepToBestFeatureForOtherFile)
                        rtCalPoints.Add(kvp.Key, new Tuple<double, double>(pepToBestFeatureForThisFile[kvp.Key].apexPeak.retentionTime, kvp.Value.apexPeak.retentionTime));

                    var someDoubles = rtCalPoints.Select(p => (p.Value.Item1 - p.Value.Item2));
                    double average = someDoubles.Average();
                    double sumOfSquaresOfDifferences = someDoubles.Select(val => (val - average) * (val - average)).Sum();
                    double sd = Math.Sqrt(sumOfSquaresOfDifferences / (someDoubles.Count() - 1));

                    // remove extreme outliers
                    if (sd > 1.0)
                    {
                        var pointsToRemove = rtCalPoints.Where(p => p.Value.Item1 - p.Value.Item2 > average + sd || p.Value.Item1 - p.Value.Item2 < average - sd).ToList();
                        foreach (var point in pointsToRemove)
                            rtCalPoints.Remove(point.Key);

                        someDoubles = rtCalPoints.Select(p => (p.Value.Item1 - p.Value.Item2));
                        average = someDoubles.Average();
                        sumOfSquaresOfDifferences = someDoubles.Select(val => (val - average) * (val - average)).Sum();
                        sd = Math.Sqrt(sumOfSquaresOfDifferences / (someDoubles.Count() - 1));
                    }

                    // find rt differences between files
                    List<Tuple<double, double>> rtCalPoints2 = rtCalPoints.Values.OrderBy(p => p.Item1).ToList();

                    int minRt = (int)Math.Floor(rtCalPoints2.First().Item1);
                    int maxRt = (int)Math.Ceiling(rtCalPoints2.Last().Item1);
                    var roughRts = Enumerable.Range(minRt, (maxRt - minRt) + 1);
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
                        if (rtCalibrationDictionary.TryGetValue(i, out List<double> list))
                        {
                            if (list.Count > 3)
                            {
                                list.Sort();
                                rtCalRunningSpline[i] = list[list.Count / 2]; // median rt difference for this timepoint

                                average = list.Average();
                                sumOfSquaresOfDifferences = list.Select(val => (val - average) * (val - average)).Sum();
                                sd = Math.Sqrt(sumOfSquaresOfDifferences / (list.Count - 1));

                                if (3 * sd > (mbrRtWindow / 2.0))
                                    stdevRunningSpline[i] = mbrRtWindow / 2.0;
                                else
                                    stdevRunningSpline[i] = 3 * sd;
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
                    var allMatchedFeaturesToLookForNow = allMbrFeaturesForThisFile.Where(p => p.identifyingScans.First().fileInfo.Equals(otherFile.Key)).ToList();

                    // filter peak candidates with rt cal to get apex peak
                    foreach (var mbrFeature in allMatchedFeaturesToLookForNow)
                    {
                        if (mbrFeature.isotopeClusters.Any())
                        {
                            // shift = thisFileRt - otherFileRt
                            int rtSplineLookupTime = (int)Math.Round(mbrFeature.identifyingScans.First().ms2RetentionTimeInMinutes);

                            if (rtSplineLookupTime < rtCalRunningSpline.Length)
                            {
                                double rtShift = rtCalRunningSpline[rtSplineLookupTime];
                                double rtToleranceHere = stdevRunningSpline[rtSplineLookupTime];
                                double theoreticalRt = mbrFeature.identifyingScans.First().ms2RetentionTimeInMinutes + rtShift;

                                if (!double.IsNaN(rtShift))
                                    mbrFeature.isotopeClusters = mbrFeature.isotopeClusters.Where(p => Math.Abs(p.retentionTime - theoreticalRt) < rtToleranceHere).ToList();
                                else
                                    mbrFeature.isotopeClusters = new List<IsotopeCluster>();
                            }
                            else
                                mbrFeature.isotopeClusters = new List<IsotopeCluster>();
                        }
                    }

                    foreach (var feature in allMatchedFeaturesToLookForNow)
                        if (feature.isotopeClusters.Any())
                            feature.CalculateIntensityForThisFeature(integrate);
                }
            }

            foreach (var file in rawFileInformation)
                RunErrorChecking(file);
        }

        private void ConstructIndexKeysFromIdentifications()
        {
            // start making index
            indexedMzKeys = new HashSet<double>();

            var peptideGroups = allIdentifications.GroupBy(p => p.ModifiedSequence).ToList();
            var peptideBaseSeqs = new HashSet<string>(allIdentifications.Select(p => p.BaseSequence));
            var numCarbonsToIsotopicDistribution = new Dictionary<int, IsotopicDistribution>();
            baseSequenceToIsotopicDistribution = new Dictionary<string, List<KeyValuePair<double, double>>>();

            foreach (var baseSeq in peptideBaseSeqs)
            {
                if (baseSequenceToIsotopicDistribution.ContainsKey(baseSeq))
                    continue;

                Proteomics.Peptide p = new Proteomics.Peptide(baseSeq);
                int numCarbonsInThisPeptide = p.ElementCountWithIsotopes("C");

                // get expected C13 mass shifts and abundances
                if (!numCarbonsToIsotopicDistribution.TryGetValue(numCarbonsInThisPeptide, out IsotopicDistribution isotopicDistribution))
                {
                    isotopicDistribution = IsotopicDistribution.GetDistribution(ChemicalFormula.ParseFormula("C" + numCarbonsInThisPeptide), 0.00001, 0.001);
                    numCarbonsToIsotopicDistribution.Add(numCarbonsInThisPeptide, isotopicDistribution);
                }

                var masses = isotopicDistribution.Masses.ToArray();
                var abundances = isotopicDistribution.Intensities.ToArray();
                var isotopicMassesAndNormalizedAbundances = new List<KeyValuePair<double, double>>();

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

                baseSequenceToIsotopicDistribution.Add(baseSeq, isotopicMassesAndNormalizedAbundances);
            }

            var minChargeState = allIdentifications.Select(p => p.precursorChargeState).Min();
            var maxChargeState = allIdentifications.Select(p => p.precursorChargeState).Max();
            chargeStates = Enumerable.Range(minChargeState, (maxChargeState - minChargeState) + 1);

            // build theoretical m/z bins
            foreach (var pepGroup in peptideGroups)
            {
                double lowestCommonMassShift = baseSequenceToIsotopicDistribution[pepGroup.First().BaseSequence].Select(p => p.Key).Min();
                var mostCommonIsotopeShift = baseSequenceToIsotopicDistribution[pepGroup.First().BaseSequence].Where(p => p.Value == 1).First().Key;

                var thisPeptidesLowestCommonMass = pepGroup.First().monoisotopicMass + lowestCommonMassShift;
                var thisPeptidesMostAbundantMass = pepGroup.First().monoisotopicMass + mostCommonIsotopeShift;

                foreach (var pep in pepGroup)
                {
                    pep.massToLookFor = requireMonoisotopicMass ? pepGroup.First().monoisotopicMass : thisPeptidesMostAbundantMass;
                }

                foreach (var chargeState in chargeStates)
                {
                    var t = pepGroup.First().massToLookFor.ToMz(chargeState);
                    double floorMz = Math.Floor(t * 100) / 100;
                    double ceilingMz = Math.Ceiling(t * 100) / 100;

                    if (!indexedMzKeys.Contains(floorMz))
                        indexedMzKeys.Add(floorMz);
                    if (!indexedMzKeys.Contains(ceilingMz))
                        indexedMzKeys.Add(ceilingMz);
                }
            }
        }

        private Dictionary<double, List<IndexedMassSpectralPeak>> IndexMassSpectralPeaks(RawFileInfo fileInfo, out Dictionary<int, IMsDataScan<IMzSpectrum<IMzPeak>>> allMs1Scans)
        {
            // construct bins
            var indexedMzs = indexedMzKeys.ToDictionary(v => v, v => new List<IndexedMassSpectralPeak>());
            var ms1ScanList = new List<IMsDataScan<IMzSpectrum<IMzPeak>>>();
            allMs1Scans = new Dictionary<int, IMsDataScan<IMzSpectrum<IMzPeak>>>();

            // open raw file
            var ext = Path.GetExtension(fileInfo.fullFilePathWithExtension).ToUpperInvariant();
            if (ext == ".MZML")
            {
                if (fileInfo.dataFile == null)
                {
                    try
                    {
                        ms1ScanList = Mzml.LoadAllStaticData(fileInfo.fullFilePathWithExtension).Where(p => p.MsnOrder == 1).Select(v => v as IMsDataScan<IMzSpectrum<IMzPeak>>).ToList();
                    }
                    catch (FileNotFoundException)
                    {
                        if (!silent)
                        {
                            Console.WriteLine("\nCan't find mzml file" + fileInfo.fullFilePathWithExtension + "\n");
                        }
                        return null;
                    }
                    catch (Exception e)
                    {
                        if (!silent)
                        {
                            Console.WriteLine("Problem opening mzml file " + fileInfo.fullFilePathWithExtension + "; " + e.Message);
                        }
                        return null;
                    }
                }
                else
                {
                    ms1ScanList = fileInfo.dataFile.Where(p => p.MsnOrder == 1).Select(v => v as IMsDataScan<IMzSpectrum<IMzPeak>>).ToList();
                }
            }
            else if (ext == ".RAW")
            {
#if NETFRAMEWORK
                var thermoFile = fileInfo.dataFile as IO.Thermo.ThermoFile;
                if (fileInfo.dataFile == null)
                {
                    try
                    {
                        thermoFile = IO.Thermo.ThermoDynamicData.InitiateDynamicConnection(fileInfo.fullFilePathWithExtension);
                    }
                    catch (FileNotFoundException)
                    {
                        if (!silent)
                        {
                            Console.WriteLine("\nCan't find raw file" + fileInfo.fullFilePathWithExtension + "\n");
                        }
                        return null;
                    }
                    catch (Exception e)
                    {
                        if (!silent)
                        {
                            throw new MzLibException("FlashLFQ Error: Problem opening raw file " + fileInfo.fullFilePathWithExtension + "; " + e.Message);
                        }
                    }
                }

                int[] msOrders = thermoFile.ThermoGlobalParams.msOrderByScan;
                for (int i = 0; i < msOrders.Length; i++)
                    if (msOrders[i] == 1)
                        ms1ScanList.Add(thermoFile.GetOneBasedScan(i + 1) as IMsDataScan<IMzSpectrum<IMzPeak>>);

#else
                if (!silent)
                {
                    Console.WriteLine("Cannot open RAW with .NETStandard code - are you on Linux? " + fileInfo.fullFilePathWithExtension);
                }
                return null;
#endif
            }
            else
            {
                if (!silent)
                    Console.WriteLine("Unsupported file type " + ext);
            }

            if (!silent)
                Console.WriteLine("Assigning MS1 peaks to bins");

            //multithreaded bin-filling
            var allGoodPeaks = new List<List<KeyValuePair<double, IndexedMassSpectralPeak>>>();

            Parallel.ForEach(Partitioner.Create(0, ms1ScanList.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    var threadLocalGoodPeaks = new List<KeyValuePair<double, IndexedMassSpectralPeak>>();

                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        int peakIndexInThisScan = 0;

                        for (int j = 0; j < ms1ScanList[i].MassSpectrum.XArray.Length; j++)
                        {
                            IndexedMassSpectralPeak element = null;
                            double floorMz = (Math.Floor(ms1ScanList[i].MassSpectrum.XArray[j] * 100) / 100);
                            double ceilingMz = (Math.Ceiling(ms1ScanList[i].MassSpectrum.XArray[j] * 100) / 100);

                            if (indexedMzs.ContainsKey(floorMz))
                            {
                                element = new IndexedMassSpectralPeak(ms1ScanList[i].MassSpectrum.XArray[j], ms1ScanList[i].MassSpectrum.YArray[j], peakIndexInThisScan, ms1ScanList[i].OneBasedScanNumber);
                                threadLocalGoodPeaks.Add(new KeyValuePair<double, IndexedMassSpectralPeak>(floorMz, element));
                            }
                            if (indexedMzs.ContainsKey(ceilingMz))
                            {
                                if (element == null)
                                    element = new IndexedMassSpectralPeak(ms1ScanList[i].MassSpectrum.XArray[j], ms1ScanList[i].MassSpectrum.YArray[j], peakIndexInThisScan, ms1ScanList[i].OneBasedScanNumber);
                                threadLocalGoodPeaks.Add(new KeyValuePair<double, IndexedMassSpectralPeak>(ceilingMz, element));
                            }

                            peakIndexInThisScan++;
                        }
                    }

                    lock (allGoodPeaks)
                        allGoodPeaks.Add(threadLocalGoodPeaks);
                }
            );

            Parallel.ForEach(Partitioner.Create(0, allGoodPeaks.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range) =>
                {
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        foreach (var element in allGoodPeaks[i])
                        {
                            var t = indexedMzs[element.Key];
                            lock (t)
                                t.Add(element.Value);
                        }
                    }
                });

            allMs1Scans = ms1ScanList.ToDictionary(p => p.OneBasedScanNumber, p => p);
            return indexedMzs;
        }

        private void QuantifyMS2IdentifiedPeptides(RawFileInfo fileInfo, Dictionary<double, List<IndexedMassSpectralPeak>> mzBins, Dictionary<int, IMsDataScan<IMzSpectrum<IMzPeak>>> allMs1Scans)
        {
            if (!silent)
                Console.WriteLine("Quantifying peptides for " + fileInfo.filenameWithoutExtension);

            results.peaks.Add(fileInfo, new List<ChromatographicPeak>());
            var concurrentBagOfFeatures = new ConcurrentBag<ChromatographicPeak>();

            var identifications = allIdentifications.Where(p => p.fileInfo.Equals(fileInfo)).ToList();

            if (!identifications.Any())
                return;

            // need to make this to look in RT space around a certain scan
            var listOfScans = allMs1Scans.Values.OrderBy(p => p.OneBasedScanNumber).ToList();
            var scanNumToIndex = new Dictionary<int, int>();
            for (int i = 0; i < listOfScans.Count; i++)
                scanNumToIndex.Add(listOfScans[i].OneBasedScanNumber, i);
            List<int> ms1ScanNumbers = listOfScans.Select(p => p.OneBasedScanNumber).OrderBy(p => p).ToList();

            Parallel.ForEach(Partitioner.Create(0, identifications.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        var identification = identifications[i];
                        ChromatographicPeak msmsFeature = new ChromatographicPeak();
                        msmsFeature.identifyingScans.Add(identification);
                        msmsFeature.isMbrFeature = false;
                        msmsFeature.rawFileInfo = fileInfo;

                        foreach (var chargeState in chargeStates)
                        {
                            if (idSpecificChargeState)
                                if (chargeState != identification.precursorChargeState)
                                    continue;

                            double theorMzHere = identification.massToLookFor.ToMz(chargeState);
                            double mzTolHere = ((peakfindingPpmTolerance / 1e6) * identification.monoisotopicMass) / chargeState;

                            double floorMz = Math.Floor(theorMzHere * 100) / 100;
                            double ceilingMz = Math.Ceiling(theorMzHere * 100) / 100;

                            IEnumerable<IndexedMassSpectralPeak> binPeaks = new List<IndexedMassSpectralPeak>();

                            for (double j = floorMz; j <= ceilingMz; j += 0.01)
                            {
                                if (mzBins.TryGetValue(Math.Round(j, 2), out List<IndexedMassSpectralPeak> list))
                                    binPeaks = binPeaks.Concat(list);
                            }

                            // filter by mz tolerance
                            var binPeaksHere = binPeaks.Where(p => Math.Abs(p.mz - theorMzHere) < mzTolHere);
                            // remove duplicates
                            binPeaksHere = binPeaksHere.Distinct();
                            // filter by RT
                            binPeaksHere = binPeaksHere.Where(p => Math.Abs(allMs1Scans[p.oneBasedScanNumber].RetentionTime - identification.ms2RetentionTimeInMinutes) < rtTol);

                            if (binPeaksHere.Any())
                            {
                                // get precursor scan to start at
                                int precursorScanNum = 0;
                                foreach (var ms1Scan in allMs1Scans)
                                {
                                    if (ms1Scan.Value.RetentionTime < identification.ms2RetentionTimeInMinutes)
                                        precursorScanNum = ms1Scan.Value.OneBasedScanNumber;
                                    else
                                        break;
                                }
                                if (precursorScanNum == 0)
                                    throw new MzLibException("Error getting precursor scan number");

                                // separate peaks by rt into left and right of the identification RT
                                var rightPeaks = binPeaksHere.Where(p => allMs1Scans[p.oneBasedScanNumber].RetentionTime >= identification.ms2RetentionTimeInMinutes).OrderBy(p => allMs1Scans[p.oneBasedScanNumber].RetentionTime);
                                var leftPeaks = binPeaksHere.Where(p => allMs1Scans[p.oneBasedScanNumber].RetentionTime < identification.ms2RetentionTimeInMinutes).OrderByDescending(p => allMs1Scans[p.oneBasedScanNumber].RetentionTime);

                                // store peaks on each side of the identification RT
                                var crawledRightPeaks = ScanCrawl(rightPeaks, missedScansAllowed, precursorScanNum, ms1ScanNumbers);
                                var crawledLeftPeaks = ScanCrawl(leftPeaks, missedScansAllowed, precursorScanNum, ms1ScanNumbers);

                                // filter again by smaller mz tolerance
                                mzTolHere = ((ppmTolerance / 1e6) * identification.monoisotopicMass) / chargeState;
                                var validPeaks = crawledRightPeaks.Concat(crawledLeftPeaks);
                                validPeaks = validPeaks.Where(p => Math.Abs(p.mz - theorMzHere) < mzTolHere);

                                // filter by isotopic distribution
                                var validIsotopeClusters = FilterPeaksByIsotopicDistribution(validPeaks, identification, chargeState, true, allMs1Scans);

                                // if multiple mass spectral peaks in the same scan are valid, pick the one with the smallest mass error
                                var peaksInSameScan = validIsotopeClusters.GroupBy(p => p.indexedPeak.oneBasedScanNumber).Where(v => v.Count() > 1);
                                if (peaksInSameScan.Any())
                                {
                                    foreach (var group in peaksInSameScan)
                                    {
                                        var mzToUse = group.Select(p => Math.Abs(p.indexedPeak.mz - theorMzHere)).Min();
                                        var peakToUse = group.Where(p => Math.Abs(p.indexedPeak.mz - theorMzHere) == mzToUse).First();
                                        var peaksToRemove = group.Where(p => p != peakToUse);
                                        validIsotopeClusters = validIsotopeClusters.Except(peaksToRemove);
                                    }
                                }

                                foreach (var validCluster in validIsotopeClusters)
                                    msmsFeature.isotopeClusters.Add(validCluster);
                            }
                        }

                        msmsFeature.CalculateIntensityForThisFeature(integrate);
                        CutPeak(msmsFeature, integrate, ms1ScanNumbers);
                        concurrentBagOfFeatures.Add(msmsFeature);
                    }
                }
            );

            // merge results from all threads together
            results.peaks[fileInfo] = concurrentBagOfFeatures.ToList();
        }

        private void MatchBetweenRunsInitialPeakfinding(RawFileInfo fileInfo, Dictionary<double, List<IndexedMassSpectralPeak>> mzBins, Dictionary<int, IMsDataScan<IMzSpectrum<IMzPeak>>> allMs1Scans)
        {
            if (!silent)
                Console.WriteLine("Finding possible matched peptides for " + fileInfo.filenameWithoutExtension);

            if (!results.peaks.ContainsKey(fileInfo) || results.peaks[fileInfo].Count == 0)
                return;

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

                            ChromatographicPeak mbrFeature = new ChromatographicPeak();
                            mbrFeature.identifyingScans.Add(identification);
                            mbrFeature.isMbrFeature = true;
                            mbrFeature.rawFileInfo = fileInfo;

                            foreach (var chargeState in chargeStates)
                            {
                                double theorMzHere = identification.massToLookFor.ToMz(chargeState);
                                double mzTolHere = ((mbrppmTolerance / 1e6) * identification.monoisotopicMass) / chargeState;

                                double floorMz = Math.Floor(theorMzHere * 100) / 100;
                                double ceilingMz = Math.Ceiling(theorMzHere * 100) / 100;

                                IEnumerable<IndexedMassSpectralPeak> binPeaks = new List<IndexedMassSpectralPeak>();
                                for (double j = floorMz; j <= ceilingMz; j += 0.01)
                                {
                                    if (mzBins.TryGetValue(Math.Round(j, 2), out List<IndexedMassSpectralPeak> list))
                                        binPeaks = binPeaks.Concat(list);
                                }

                                // filter by mz tolerance
                                var binPeaksHere = binPeaks.Where(p => Math.Abs(p.mz - theorMzHere) < mzTolHere);
                                // filter by rt
                                binPeaksHere = binPeaksHere.Where(p => Math.Abs(allMs1Scans[p.oneBasedScanNumber].RetentionTime - identification.ms2RetentionTimeInMinutes) < initialMbrRtWindow);
                                // remove duplicates
                                binPeaksHere = binPeaksHere.Distinct();
                                // filter by isotopic distribution
                                var validIsotopeClusters = FilterPeaksByIsotopicDistribution(binPeaksHere, identification, chargeState, true, allMs1Scans);

                                if (validIsotopeClusters.Any())
                                {
                                    mbrFeature.isotopeClusters.AddRange(validIsotopeClusters);
                                }
                            }

                            if (mbrFeature.isotopeClusters.Any())
                            {
                                concurrentBagOfMatchedFeatures.Add(mbrFeature);
                            }
                        }
                    }
                );
            }

            results.peaks[fileInfo].AddRange(concurrentBagOfMatchedFeatures);
        }

        private void RunErrorChecking(RawFileInfo rawFile)
        {
            results.peaks[rawFile].RemoveAll(p => p.isMbrFeature && !p.isotopeClusters.Any());

            if (!silent)
                Console.WriteLine("Checking errors");
            var featuresWithSamePeak = results.peaks[rawFile].Where(v => v.intensity != 0).GroupBy(p => p.apexPeak.indexedPeak);
            featuresWithSamePeak = featuresWithSamePeak.Where(p => p.Count() > 1);

            // condense duplicate features (features with same sequence and apex peak)
            foreach (var duplicateFeature in featuresWithSamePeak)
                duplicateFeature.First().MergeFeatureWith(duplicateFeature, integrate);
            results.peaks[rawFile].RemoveAll(p => p.intensity == -1);

            // check for multiple features per peptide within a time window
            var featuresToMaybeMerge = results.peaks[rawFile].Where(p => p.NumIdentificationsByFullSeq == 1 && p.apexPeak != null).GroupBy(p => p.identifyingScans.First().ModifiedSequence).Where(p => p.Count() > 1);
            if (featuresToMaybeMerge.Any())
            {
                foreach (var group in featuresToMaybeMerge)
                {
                    if (idSpecificChargeState)
                    {
                        var group2 = group.ToList().GroupBy(p => p.apexPeak.chargeState).Where(v => v.Count() > 1);

                        foreach (var group3 in group2)
                        {
                            foreach (var feature in group3)
                            {
                                if (feature.intensity != -1)
                                {
                                    var featuresToMerge = group3.Where(p => Math.Abs(p.apexPeak.retentionTime - feature.apexPeak.retentionTime) < rtTol && p.intensity != -1);
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
                                var featuresToMerge = group.Where(p => Math.Abs(p.apexPeak.retentionTime - feature.apexPeak.retentionTime) < rtTol && p.intensity != -1);
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
                    var msmsIdentsForThisFile = ambiguousFeature.identifyingScans.Where(p => p.fileInfo.Equals(ambiguousFeature.rawFileInfo));

                    if (!msmsIdentsForThisFile.Any())
                    {
                        // mbr matched more than one identification to this peak - cannot resolve
                        ambiguousFeature.intensity = -1;
                    }
                    else
                    {
                        // msms identifications take precident over mbr features
                        ambiguousFeature.identifyingScans = msmsIdentsForThisFile.ToList();
                        ambiguousFeature.ResolveIdentifications();
                    }
                }

                results.peaks[rawFile].RemoveAll(p => p.intensity == -1);
            }
        }

        private IEnumerable<IsotopeCluster> FilterPeaksByIsotopicDistribution(IEnumerable<IndexedMassSpectralPeak> peaks, Identification identification, int chargeState, bool lookForBadIsotope, Dictionary<int, IMsDataScan<IMzSpectrum<IMzPeak>>> allMs1Scans)
        {
            var isotopeClusters = new List<IsotopeCluster>();
            var isotopeMassShifts = baseSequenceToIsotopicDistribution[identification.BaseSequence];

            if (isotopeMassShifts.Count < numIsotopesRequired)
                return isotopeClusters;

            double isotopeMzTol = ((isotopePpmTolerance / 1e6) * identification.monoisotopicMass) / chargeState;

            foreach (var thisPeakWithScan in peaks)
            {
                // calculate theoretical isotopes relative to observed peak
                var theorIsotopeMzs = new double[isotopeMassShifts.Count];
                int isotopicPeakUnitOfPeakZeroIsMono = Convert.ToInt32(thisPeakWithScan.mz.ToMass(chargeState) - identification.monoisotopicMass);
                var mainpeakMz = thisPeakWithScan.mz;

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
                int earliestIsotopicPeakIndexPossible = thisPeakWithScan.zeroBasedIndexOfPeakInScan;
                for (int i = earliestIsotopicPeakIndexPossible; i >= 0; i--)
                {
                    if (allMs1Scans[thisPeakWithScan.oneBasedScanNumber].MassSpectrum.XArray[i] < lowestMzIsotopePossible)
                        break;
                    earliestIsotopicPeakIndexPossible = i;
                }

                // find the highest-mass isotope possible
                for (int i = earliestIsotopicPeakIndexPossible; i < allMs1Scans[thisPeakWithScan.oneBasedScanNumber].MassSpectrum.Size; i++)
                {
                    if (allMs1Scans[thisPeakWithScan.oneBasedScanNumber].MassSpectrum.XArray[i] > highestMzIsotopePossible)
                        break;
                    possibleIsotopePeaks.Add(new Tuple<double, double>(allMs1Scans[thisPeakWithScan.oneBasedScanNumber].MassSpectrum.XArray[i], allMs1Scans[thisPeakWithScan.oneBasedScanNumber].MassSpectrum.YArray[i]));
                }

                if (lookForBadIsotope)
                {
                    bool badPeak = false;
                    double prevIsotopePeakMz = (theorIsotopeMzs[0] - (1.003322 / chargeState));

                    for (int i = earliestIsotopicPeakIndexPossible; i > 0; i--)
                    {
                        if (Math.Abs(allMs1Scans[thisPeakWithScan.oneBasedScanNumber].MassSpectrum.XArray[i] - prevIsotopePeakMz) < isotopeMzTol)
                            if (allMs1Scans[thisPeakWithScan.oneBasedScanNumber].MassSpectrum.YArray[i] / thisPeakWithScan.intensity > 1.0)
                                badPeak = true;
                        if (allMs1Scans[thisPeakWithScan.oneBasedScanNumber].MassSpectrum.XArray[i] < (prevIsotopePeakMz - isotopeMzTol))
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

                    isotopeClusters.Add(new IsotopeCluster(thisPeakWithScan, chargeState, isotopeClusterIntensity, allMs1Scans[thisPeakWithScan.oneBasedScanNumber].RetentionTime));
                }
            }

            return isotopeClusters;
        }

        private IEnumerable<IndexedMassSpectralPeak> ScanCrawl(IOrderedEnumerable<IndexedMassSpectralPeak> peaksWithScans, int missedScansAllowed, int startingMS1ScanNumber, List<int> ms1ScanNumbers)
        {
            var validPeaksWithScans = new List<IndexedMassSpectralPeak>();

            int lastGoodIndex = ms1ScanNumbers.IndexOf(startingMS1ScanNumber);
            int ms1IndexHere = lastGoodIndex - 1;
            int missedScans = 0;

            foreach (var thisPeakWithScan in peaksWithScans)
            {
                ms1IndexHere = ms1ScanNumbers.IndexOf(thisPeakWithScan.oneBasedScanNumber);
                missedScans += Math.Abs(ms1IndexHere - lastGoodIndex) - 1;

                if (missedScans > missedScansAllowed)
                    break;

                // found a good peak; reset missed scans to 0
                missedScans = 0;
                lastGoodIndex = ms1IndexHere;

                validPeaksWithScans.Add(thisPeakWithScan);
            }

            return validPeaksWithScans;
        }

        private void CutPeak(ChromatographicPeak peak, bool integrate, List<int> ms1ScanNumbers)
        {
            bool cutThisPeak = false;
            IsotopeCluster valleyTimePoint = null;

            if (peak.isotopeClusters.Count() < 5)
                return;

            // find out if we need to split this peak by using the discrimination factor
            var timePointsForApexZ = peak.isotopeClusters.Where(p => p.chargeState == peak.apexPeak.chargeState);
            var leftTimePoints = timePointsForApexZ.Where(p => p.retentionTime <= peak.apexPeak.retentionTime).OrderByDescending(v => v.retentionTime);
            var rightTimePoints = timePointsForApexZ.Where(p => p.retentionTime >= peak.apexPeak.retentionTime).OrderBy(v => v.retentionTime);

            double mind0 = 0.6;

            foreach (var timePoint in rightTimePoints)
            {
                if (valleyTimePoint == null || timePoint.isotopeClusterIntensity < valleyTimePoint.isotopeClusterIntensity)
                    valleyTimePoint = timePoint;

                var timePointsBetweenApexAndThisTimePoint = rightTimePoints.Where(p => p.retentionTime <= timePoint.retentionTime).ToList();

                var d0 = (timePoint.isotopeClusterIntensity - valleyTimePoint.isotopeClusterIntensity) / timePoint.isotopeClusterIntensity;
                if (d0 > mind0)
                {
                    var secondValleyTimePoint = timePointsBetweenApexAndThisTimePoint[timePointsBetweenApexAndThisTimePoint.IndexOf(valleyTimePoint) + 1];

                    d0 = (timePoint.isotopeClusterIntensity - secondValleyTimePoint.isotopeClusterIntensity) / timePoint.isotopeClusterIntensity;

                    if (d0 > mind0)
                    {
                        cutThisPeak = true;
                        break;
                    }
                    else
                    {
                        // check for missed scan around valley time point
                        var tpBeforeValleyTimePoint = timePointsBetweenApexAndThisTimePoint[timePointsBetweenApexAndThisTimePoint.IndexOf(valleyTimePoint) - 1];

                        int indexOfTimepointBeforeValleyScan = ms1ScanNumbers.IndexOf(tpBeforeValleyTimePoint.indexedPeak.oneBasedScanNumber);
                        int indexOfValleyScan = ms1ScanNumbers.IndexOf(valleyTimePoint.indexedPeak.oneBasedScanNumber);
                        int indexOfSecondValleyScan = ms1ScanNumbers.IndexOf(secondValleyTimePoint.indexedPeak.oneBasedScanNumber);

                        if (Math.Abs(indexOfValleyScan - indexOfTimepointBeforeValleyScan) > 1)
                        {
                            cutThisPeak = true;
                            break;
                        }
                        else if (Math.Abs(indexOfValleyScan - indexOfSecondValleyScan) > 1)
                        {
                            cutThisPeak = true;
                            break;
                        }
                    }
                }
            }

            if (cutThisPeak == false)
            {
                valleyTimePoint = null;

                foreach (var timePoint in leftTimePoints)
                {
                    if (valleyTimePoint == null || timePoint.isotopeClusterIntensity < valleyTimePoint.isotopeClusterIntensity)
                        valleyTimePoint = timePoint;

                    var timePointsBetweenApexAndThisTimePoint = leftTimePoints.Where(p => p.retentionTime >= timePoint.retentionTime).ToList();

                    var d0 = (timePoint.isotopeClusterIntensity - valleyTimePoint.isotopeClusterIntensity) / timePoint.isotopeClusterIntensity;
                    if (d0 > mind0)
                    {
                        var secondValleyTimePoint = timePointsBetweenApexAndThisTimePoint[timePointsBetweenApexAndThisTimePoint.IndexOf(valleyTimePoint) + 1];

                        d0 = (timePoint.isotopeClusterIntensity - secondValleyTimePoint.isotopeClusterIntensity) / timePoint.isotopeClusterIntensity;

                        if (d0 > mind0)
                        {
                            cutThisPeak = true;
                            break;
                        }
                        else
                        {
                            // check for missed scan around valley time point
                            var tpBeforeValleyTimePoint = timePointsBetweenApexAndThisTimePoint[timePointsBetweenApexAndThisTimePoint.IndexOf(valleyTimePoint) - 1];

                            int indexOfTimepointBeforeValleyScan = ms1ScanNumbers.IndexOf(tpBeforeValleyTimePoint.indexedPeak.oneBasedScanNumber);
                            int indexOfValleyScan = ms1ScanNumbers.IndexOf(valleyTimePoint.indexedPeak.oneBasedScanNumber);
                            int indexOfSecondValleyScan = ms1ScanNumbers.IndexOf(secondValleyTimePoint.indexedPeak.oneBasedScanNumber);

                            if (Math.Abs(indexOfValleyScan - indexOfTimepointBeforeValleyScan) > 1)
                            {
                                cutThisPeak = true;
                                break;
                            }
                            else if (Math.Abs(indexOfValleyScan - indexOfSecondValleyScan) > 1)
                            {
                                cutThisPeak = true;
                                break;
                            }
                        }
                    }
                }
            }

            // cut
            if (cutThisPeak)
            {
                var splitLeft = peak.isotopeClusters.Where(p => p.retentionTime <= valleyTimePoint.retentionTime).ToList();
                var splitRight = peak.isotopeClusters.Where(p => p.retentionTime >= valleyTimePoint.retentionTime).ToList();

                if (peak.identifyingScans.First().ms2RetentionTimeInMinutes > splitLeft.Max(p => p.retentionTime))
                    foreach (var timePoint in splitLeft)
                        peak.isotopeClusters.Remove(timePoint);
                else
                    foreach (var timePoint in splitRight)
                        peak.isotopeClusters.Remove(timePoint);

                // recalculate intensity for the peak
                peak.CalculateIntensityForThisFeature(integrate);
                peak.splitRT = valleyTimePoint.retentionTime;

                // recursively cut
                CutPeak(peak, integrate, ms1ScanNumbers);
            }
        }

        #endregion Private Methods
    }
}