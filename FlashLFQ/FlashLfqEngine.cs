using Chemistry;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace FlashLFQ
{
    public class FlashLfqEngine
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
        public readonly bool ErrorCheckAmbiguousMatches;
        public readonly bool IdSpecificChargeState;
        public readonly bool RequireMonoisotopicMass;
        public readonly bool Normalize;
        public readonly double MinDiscFactorToCutAt;

        // MBR settings
        public readonly bool MatchBetweenRuns;
        public readonly double MbrRtWindow;
        public readonly double MbrPpmTolerance;
        public readonly double MbrQValueCutoff;

        // settings for the Bayesian protein quantification engine
        public readonly bool BayesianProteinQuant;
        public readonly string ProteinQuantBaseCondition;
        public readonly double? ProteinQuantFoldChangeCutoff;
        public readonly int McmcSteps;
        public readonly int McmcBurninSteps;
        public readonly bool UseSharedPeptidesForProteinQuant;
        public readonly int? RandomSeed;

        // structures used in the FlashLFQ engine
        private List<SpectraFileInfo> _spectraFileInfo;
        private Stopwatch _globalStopwatch;
        private List<Identification> _allIdentifications;
        private Dictionary<string, List<(double, double)>> _modifiedSequenceToIsotopicDistribution;
        private IEnumerable<int> _chargeStates;
        private FlashLfqResults _results;
        private Dictionary<SpectraFileInfo, Ms1ScanInfo[]> _ms1Scans;
        private PeakIndexingEngine _peakIndexingEngine;

        public FlashLfqEngine(
            List<Identification> allIdentifications,
            bool normalize = false,
            double ppmTolerance = 10.0,
            double isotopeTolerancePpm = 5.0,
            bool integrate = false,
            int numIsotopesRequired = 2,
            bool idSpecificChargeState = false,
            bool requireMonoisotopicMass = true,
            bool silent = false,
            int maxThreads = -1,

            // MBR settings
            bool matchBetweenRuns = false,
            double matchBetweenRunsPpmTolerance = 5.0,
            double mbrQValueCutoff = 0.10,
            double maxMbrWindow = 2.5,

            // settings for the Bayesian protein quantification engine
            bool bayesianProteinQuant = false,
            string proteinQuantBaseCondition = null,
            double? proteinQuantFoldChangeCutoff = null,
            int mcmcSteps = 3000,
            int mcmcBurninSteps = 1000,
            bool useSharedPeptidesForProteinQuant = false,
            int? randomSeed = null)
        {
            Loaders.LoadElements();

            _globalStopwatch = new Stopwatch();
            _chargeStates = new List<int>();
            _peakIndexingEngine = new PeakIndexingEngine();

            _spectraFileInfo = allIdentifications.Select(p => p.FileInfo).Distinct()
                .OrderBy(p => p.Condition)
                .ThenBy(p => p.BiologicalReplicate)
                .ThenBy(p => p.Fraction)
                .ThenBy(p => p.TechnicalReplicate).ToList();

            _allIdentifications = allIdentifications;
            PpmTolerance = ppmTolerance;
            IsotopePpmTolerance = isotopeTolerancePpm;
            MatchBetweenRuns = matchBetweenRuns;
            MbrPpmTolerance = matchBetweenRunsPpmTolerance;
            Integrate = integrate;
            NumIsotopesRequired = numIsotopesRequired;
            Silent = silent;
            IdSpecificChargeState = idSpecificChargeState;
            RequireMonoisotopicMass = requireMonoisotopicMass;
            MbrRtWindow = maxMbrWindow;
            Normalize = normalize;
            MaxThreads = maxThreads;
            BayesianProteinQuant = bayesianProteinQuant;
            ProteinQuantBaseCondition = proteinQuantBaseCondition;
            ProteinQuantFoldChangeCutoff = proteinQuantFoldChangeCutoff;
            McmcSteps = mcmcSteps;
            McmcBurninSteps = mcmcBurninSteps;
            UseSharedPeptidesForProteinQuant = useSharedPeptidesForProteinQuant;
            RandomSeed = randomSeed;
            MbrQValueCutoff = mbrQValueCutoff;

            if (MaxThreads == -1 || MaxThreads >= Environment.ProcessorCount)
            {
                MaxThreads = Environment.ProcessorCount - 1;
            }

            if (MaxThreads <= 0)
            {
                MaxThreads = 1;
            }

            PeakfindingPpmTolerance = 20.0;
            MissedScansAllowed = 1;
            RtTol = 5.0;
            ErrorCheckAmbiguousMatches = true;
            MinDiscFactorToCutAt = 0.6;
        }

        public FlashLfqResults Run()
        {
            _globalStopwatch.Start();
            _ms1Scans = new Dictionary<SpectraFileInfo, Ms1ScanInfo[]>();

            _results = new FlashLfqResults(_spectraFileInfo, _allIdentifications);

            // build m/z index keys
            CalculateTheoreticalIsotopeDistributions();

            // quantify each file
            foreach (var spectraFile in _spectraFileInfo)
            {
                // fill lookup-table with peaks from the spectra file
                if (!_peakIndexingEngine.IndexMassSpectralPeaks(spectraFile, Silent, _ms1Scans))
                {
                    // something went wrong finding/opening/indexing the file...
                    continue;
                }

                // quantify peaks using this file's IDs first
                QuantifyMs2IdentifiedPeptides(spectraFile);

                // write the indexed peaks for MBR later
                if (MatchBetweenRuns)
                {
                    _peakIndexingEngine.SerializeIndex(spectraFile);
                }

                // error checking function
                // handles features with multiple identifying scans and scans that are associated with more than one feature
                RunErrorChecking(spectraFile);

                if (!Silent)
                {
                    Console.WriteLine("Finished " + spectraFile.FilenameWithoutExtension);
                }

                // some memory-saving stuff
                _peakIndexingEngine.ClearIndex();
            }

            // do MBR
            if (MatchBetweenRuns)
            {
                foreach (var spectraFile in _spectraFileInfo)
                {
                    if (!Silent)
                    {
                        Console.WriteLine("Doing match-between-runs for " + spectraFile.FilenameWithoutExtension);
                    }

                    QuantifyMatchBetweenRunsPeaks(spectraFile);

                    _peakIndexingEngine.ClearIndex();

                    if (!Silent)
                    {
                        Console.WriteLine("Finished MBR for " + spectraFile.FilenameWithoutExtension);
                    }
                }
            }

            // normalize
            if (Normalize)
            {
                try
                {
                    new IntensityNormalizationEngine(_results, Integrate, Silent, MaxThreads).NormalizeResults();
                }
                catch (Exception e)
                {
                    throw new MzLibException("A crash occured in FlashLFQ during the intensity normalization process:\n" + e.Message);
                }
            }

            // calculate peptide intensities 
            _results.CalculatePeptideResults();

            // do top3 protein quantification
            _results.CalculateProteinResultsTop3();

            // do Bayesian protein fold-change analysis
            if (BayesianProteinQuant)
            {
                if (_spectraFileInfo.Count == 1 || _spectraFileInfo.Select(p => p.Condition).Distinct().Count() == 1)
                {
                    if (!Silent)
                    {
                        Console.WriteLine("Can't do Bayesian protein quant with only one spectra file or condition. FlashLFQ will still do a top3 protein quant");
                    }
                }
                else
                {
                    try
                    {
                        new ProteinQuantificationEngine(_results, MaxThreads, ProteinQuantBaseCondition, UseSharedPeptidesForProteinQuant,
                            ProteinQuantFoldChangeCutoff, RandomSeed, McmcBurninSteps, McmcSteps).Run();
                    }
                    catch (Exception e)
                    {
                        throw new MzLibException("A crash occured in FlashLFQ during the Bayesian protein quantification process:\n" + e.Message);
                    }
                }
            }

            // done
            if (!Silent)
            {
                Console.WriteLine("Done quantifying");
            }

            if (!Silent)
            {
                Console.WriteLine("Analysis time: " +
                                  _globalStopwatch.Elapsed.Hours + "h " +
                                  _globalStopwatch.Elapsed.Minutes + "m " +
                                  _globalStopwatch.Elapsed.Seconds + "s");
            }

            return _results;
        }

        /// <summary>
        /// Creates a theoretical isotope distribution for each of the identified sequences
        /// If the sequence is modified, this uses averagine for the modified part
        /// </summary>
        private void CalculateTheoreticalIsotopeDistributions()
        {
            _modifiedSequenceToIsotopicDistribution = new Dictionary<string, List<(double, double)>>();

            // calculate averagine (used for isotopic distributions for unknown modifications)
            double averageC = 4.9384;
            double averageH = 7.7583;
            double averageO = 1.4773;
            double averageN = 1.3577;
            double averageS = 0.0417;

            double averagineMass =
                PeriodicTable.GetElement("C").AverageMass * averageC +
                PeriodicTable.GetElement("H").AverageMass * averageH +
                PeriodicTable.GetElement("O").AverageMass * averageO +
                PeriodicTable.GetElement("N").AverageMass * averageN +
                PeriodicTable.GetElement("S").AverageMass * averageS;

            // calculate monoisotopic masses and isotopic envelope for the base sequences
            foreach (Identification id in _allIdentifications)
            {
                if (_modifiedSequenceToIsotopicDistribution.ContainsKey(id.ModifiedSequence))
                {
                    continue;
                }

                ChemicalFormula formula = id.OptionalChemicalFormula;

                var isotopicMassesAndNormalizedAbundances = new List<(double massShift, double abundance)>();

                if (formula == null)
                {
                    Proteomics.AminoAcidPolymer.Peptide baseSequence = new Proteomics.AminoAcidPolymer.Peptide(id.BaseSequence);
                    formula = baseSequence.GetChemicalFormula();

                    // add averagine for any unknown mass difference (i.e., a modification)
                    double massDiff = id.monoisotopicMass - baseSequence.MonoisotopicMass;

                    if (Math.Abs(massDiff) > 20)
                    {
                        double averagines = massDiff / averagineMass;

                        formula.Add("C", (int)Math.Round(averagines * averageC, 0));
                        formula.Add("H", (int)Math.Round(averagines * averageH, 0));
                        formula.Add("O", (int)Math.Round(averagines * averageO, 0));
                        formula.Add("N", (int)Math.Round(averagines * averageN, 0));
                        formula.Add("S", (int)Math.Round(averagines * averageS, 0));
                    }
                }

                var isotopicDistribution = IsotopicDistribution.GetDistribution(formula, 0.125, 1e-8);

                double[] masses = isotopicDistribution.Masses.ToArray();
                double[] abundances = isotopicDistribution.Intensities.ToArray();

                for (int i = 0; i < masses.Length; i++)
                {
                    masses[i] += (id.monoisotopicMass - formula.MonoisotopicMass);
                }

                double highestAbundance = abundances.Max();
                int highestAbundanceIndex = Array.IndexOf(abundances, highestAbundance);

                for (int i = 0; i < masses.Length; i++)
                {
                    // expected isotopic mass shifts for this peptide
                    masses[i] -= id.monoisotopicMass;

                    // normalized abundance of each isotope
                    abundances[i] /= highestAbundance;

                    // look for these isotopes
                    if (isotopicMassesAndNormalizedAbundances.Count < NumIsotopesRequired || abundances[i] > 0.1 || i < highestAbundanceIndex)
                    {
                        isotopicMassesAndNormalizedAbundances.Add((masses[i], abundances[i]));
                    }
                }

                _modifiedSequenceToIsotopicDistribution.Add(id.ModifiedSequence, isotopicMassesAndNormalizedAbundances);
            }

            var minChargeState = _allIdentifications.Min(p => p.PrecursorChargeState);
            var maxChargeState = _allIdentifications.Max(p => p.PrecursorChargeState);
            _chargeStates = Enumerable.Range(minChargeState, (maxChargeState - minChargeState) + 1);

            var peptideModifiedSequences = _allIdentifications.GroupBy(p => p.ModifiedSequence);
            foreach (var identifications in peptideModifiedSequences)
            {
                // isotope where normalized abundance is 1
                double mostAbundantIsotopeShift = _modifiedSequenceToIsotopicDistribution[identifications.First().ModifiedSequence].First(p => p.Item2 == 1.0).Item1;

                foreach (Identification identification in identifications)
                {
                    identification.peakfindingMass = identification.monoisotopicMass + mostAbundantIsotopeShift;
                }
            }
        }

        private void QuantifyMs2IdentifiedPeptides(SpectraFileInfo fileInfo)
        {
            if (!Silent)
            {
                Console.WriteLine("Quantifying peptides for " + fileInfo.FilenameWithoutExtension);
            }

            var ms2IdsForThisFile = _allIdentifications.Where(p => p.FileInfo.Equals(fileInfo)).ToList();

            if (!ms2IdsForThisFile.Any())
            {
                return;
            }

            Tolerance peakfindingTol = new PpmTolerance(PeakfindingPpmTolerance);
            Tolerance ppmTolerance = new PpmTolerance(PpmTolerance);

            var chromatographicPeaks = new ChromatographicPeak[ms2IdsForThisFile.Count];

            Parallel.ForEach(Partitioner.Create(0, ms2IdsForThisFile.Count),
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

                        var identification = ms2IdsForThisFile[i];
                        ChromatographicPeak msmsFeature = new ChromatographicPeak(identification, false, fileInfo);
                        chromatographicPeaks[i] = msmsFeature;

                        foreach (var chargeState in _chargeStates)
                        {
                            if (IdSpecificChargeState && chargeState != identification.PrecursorChargeState)
                            {
                                continue;
                            }

                            // get XIC (peakfinding)
                            List<IndexedMassSpectralPeak> xic = Peakfind(identification.ms2RetentionTimeInMinutes,
                                identification.peakfindingMass, chargeState, identification.fileInfo, peakfindingTol).OrderBy(p => p.RetentionTime).ToList();

                            // filter by smaller mass tolerance
                            xic.RemoveAll(p => !ppmTolerance.Within(p.Mz.ToMass(chargeState), identification.peakfindingMass));

                            // filter by isotopic distribution
                            List<IsotopicEnvelope> isotopicEnvelopes = GetIsotopicEnvelopes(xic, identification, chargeState, false);

                            // add isotopic envelopes to the chromatographic peak
                            msmsFeature.IsotopicEnvelopes.AddRange(isotopicEnvelopes);
                        }

                        msmsFeature.CalculateIntensityForThisFeature(Integrate);
                        CutPeak(msmsFeature, identification.Ms2RetentionTimeInMinutes);

                        if (!msmsFeature.IsotopicEnvelopes.Any())
                        {
                            continue;
                        }

                        var precursorXic = msmsFeature.IsotopicEnvelopes.Where(p => p.ChargeState == identification.PrecursorChargeState).ToList();

                        if (!precursorXic.Any())
                        {
                            msmsFeature.IsotopicEnvelopes.Clear();
                            continue;
                        }

                        int min = precursorXic.Min(p => p.IndexedPeak.ZeroBasedMs1ScanIndex);
                        int max = precursorXic.Max(p => p.IndexedPeak.ZeroBasedMs1ScanIndex);
                        msmsFeature.IsotopicEnvelopes.RemoveAll(p => p.IndexedPeak.ZeroBasedMs1ScanIndex < min);
                        msmsFeature.IsotopicEnvelopes.RemoveAll(p => p.IndexedPeak.ZeroBasedMs1ScanIndex > max);
                        msmsFeature.CalculateIntensityForThisFeature(Integrate);
                    }
                });

            _results.Peaks[fileInfo].AddRange(chromatographicPeaks.ToList());
        }

        private void QuantifyMatchBetweenRunsPeaks(SpectraFileInfo idAcceptorFile)
        {
            var acceptorFileIdentifiedPeaks = _results.Peaks[idAcceptorFile];

            if (!acceptorFileIdentifiedPeaks.Any())
            {
                return;
            }

            _peakIndexingEngine.DeserializeIndex(idAcceptorFile);

            var thisFilesIds = new HashSet<string>(acceptorFileIdentifiedPeaks.Where(p => p.IsotopicEnvelopes.Any())
                .SelectMany(p => p.Identifications.Select(d => d.ModifiedSequence)));

            // only match peptides from proteins that have at least one MS/MS identified peptide in the condition
            var thisFilesMsmsIdentifiedProteins = new HashSet<ProteinGroup>();
            foreach (SpectraFileInfo conditionFile in _spectraFileInfo.Where(p => p.Condition == idAcceptorFile.Condition))
            {
                foreach (ProteinGroup proteinGroup in _results.Peaks[conditionFile].Where(p => !p.IsMbrPeak).SelectMany(p => p.Identifications.SelectMany(v => v.ProteinGroups)))
                {
                    thisFilesMsmsIdentifiedProteins.Add(proteinGroup);
                }
            }

            Tolerance mbrTol = new PpmTolerance(MbrPpmTolerance);

            Dictionary<string, ChromatographicPeak> bestMbrHits = new Dictionary<string, ChromatographicPeak>();

            foreach (SpectraFileInfo idDonorFile in _spectraFileInfo)
            {
                if (idAcceptorFile.Equals(idDonorFile))
                {
                    continue;
                }

                // these peaks have no IDs in the acceptor file
                // match their IDs to this file
                var donorPeaksToMatch = _results.Peaks[idDonorFile].Where(p => p.NumIdentificationsByFullSeq == 1
                    && p.IsotopicEnvelopes.Any()
                    && !p.Identifications.Any(v => thisFilesIds.Contains(v.ModifiedSequence))
                    && p.Identifications.Any(v => v.ProteinGroups.Any(g => thisFilesMsmsIdentifiedProteins.Contains(g)))).ToList();

                if (!donorPeaksToMatch.Any())
                {
                    continue;
                }

                // generate RT calibration curve
                RetentionTimeCalibDataPoint[] rtCalibrationCurve = GetRtCalSpline(idDonorFile, idAcceptorFile);
                double[] donorRts = rtCalibrationCurve.Select(p => p.DonorFilePeak.Apex.IndexedPeak.RetentionTime).ToArray();

                // this is for calculating FDR later
                var sequenceToRealPeakInAcceptorFile = rtCalibrationCurve.ToDictionary(p => p.AcceptorFilePeak.Identifications.First().ModifiedSequence, v => v.AcceptorFilePeak);
                foreach (RetentionTimeCalibDataPoint point in rtCalibrationCurve)
                {
                    if (thisFilesIds.Contains(point.DonorFilePeak.Identifications.First().ModifiedSequence))
                    {
                        donorPeaksToMatch.Add(point.DonorFilePeak);
                    }
                }
                (ChromatographicPeak bestAcceptorPeak, bool isTruePositive)[] scoreWithKnownResults = new (ChromatographicPeak, bool)[donorPeaksToMatch.Count];

                ChromatographicPeak[] matchedPeaks = new ChromatographicPeak[donorPeaksToMatch.Count];

                Parallel.ForEach(Partitioner.Create(0, donorPeaksToMatch.Count),
                    new ParallelOptions { MaxDegreeOfParallelism = MaxThreads },
                    (range, loopState) =>
                    {
                        for (int i = range.Item1; i < range.Item2; i++)
                        {
                            ChromatographicPeak donorPeak = donorPeaksToMatch[i];
                            Identification identification = donorPeak.Identifications.First();

                            int rtHypothesisIndex = Array.BinarySearch(donorRts, donorPeak.Apex.IndexedPeak.RetentionTime);
                            if (rtHypothesisIndex < 0)
                            {
                                rtHypothesisIndex = ~rtHypothesisIndex;
                            }
                            if (rtHypothesisIndex >= donorRts.Length && rtHypothesisIndex >= 1)
                            {
                                rtHypothesisIndex = donorRts.Length - 1;
                            }

                            List<RetentionTimeCalibDataPoint> nearbyDataPoints = new List<RetentionTimeCalibDataPoint>();

                            // calculate accepted range of RTs
                            for (int r = rtHypothesisIndex; r < rtCalibrationCurve.Length; r++)
                            {
                                RetentionTimeCalibDataPoint rightDataPoint = rtCalibrationCurve[r];
                                if (Math.Abs(rightDataPoint.DonorFilePeak.Apex.IndexedPeak.RetentionTime - donorPeak.Apex.IndexedPeak.RetentionTime) < 0.5)
                                {
                                    nearbyDataPoints.Add(rtCalibrationCurve[r]);
                                }
                                else
                                {
                                    break;
                                }
                            }

                            for (int l = rtHypothesisIndex - 1; l >= 0; l--)
                            {
                                RetentionTimeCalibDataPoint leftDataPoint = rtCalibrationCurve[l];
                                if (Math.Abs(leftDataPoint.DonorFilePeak.Apex.IndexedPeak.RetentionTime - donorPeak.Apex.IndexedPeak.RetentionTime) < 0.5)
                                {
                                    nearbyDataPoints.Add(rtCalibrationCurve[l]);
                                }
                                else
                                {
                                    break;
                                }
                            }

                            double acceptorFileRtHypothesis;
                            double lowerRange;
                            double upperRange;

                            if (nearbyDataPoints.Count >= 4)
                            {
                                List<double> nearbyRts = nearbyDataPoints.Select(p => p.RtDiff).ToList();
                                double median = Statistics.Median(nearbyRts);

                                acceptorFileRtHypothesis = donorPeak.Apex.IndexedPeak.RetentionTime + median;
                                nearbyRts = nearbyRts.Select(p => p - median).ToList();

                                double firstQuartile = Statistics.LowerQuartile(nearbyRts);
                                double thirdQuartile = Statistics.UpperQuartile(nearbyRts);
                                double iqr = Statistics.InterquartileRange(nearbyRts);

                                lowerRange = firstQuartile - 1.5 * iqr;
                                upperRange = thirdQuartile + 1.5 * iqr;
                            }
                            else
                            {
                                continue;
                            }

                            double lowerBoundRt = acceptorFileRtHypothesis + lowerRange;
                            double upperBoundRt = acceptorFileRtHypothesis + upperRange;

                            if (upperBoundRt - acceptorFileRtHypothesis > MbrRtWindow)
                            {
                                upperBoundRt = acceptorFileRtHypothesis + MbrRtWindow;
                            }

                            if (acceptorFileRtHypothesis - lowerBoundRt > MbrRtWindow)
                            {
                                lowerBoundRt = acceptorFileRtHypothesis - MbrRtWindow;
                            }

                            Ms1ScanInfo[] ms1ScanInfos = _ms1Scans[idAcceptorFile];
                            Ms1ScanInfo start = ms1ScanInfos[0];
                            Ms1ScanInfo end = ms1ScanInfos[ms1ScanInfos.Length - 1];
                            for (int j = 0; j < ms1ScanInfos.Length; j++)
                            {
                                Ms1ScanInfo scan = ms1ScanInfos[j];
                                if (scan.RetentionTime <= lowerBoundRt)
                                {
                                    start = scan;
                                }

                                if (scan.RetentionTime >= upperBoundRt)
                                {
                                    end = scan;
                                    break;
                                }
                            }

                            List<ChromatographicPeak> acceptorPeakOptions = new List<ChromatographicPeak>();
                            List<IndexedMassSpectralPeak> chargeXic = new List<IndexedMassSpectralPeak>();

                            HashSet<IndexedMassSpectralPeak> claimedPeaks = new HashSet<IndexedMassSpectralPeak>();
                            foreach (int z in _chargeStates)
                            {
                                for (int j = start.ZeroBasedMs1ScanIndex; j <= end.ZeroBasedMs1ScanIndex; j++)
                                {
                                    IndexedMassSpectralPeak peak = _peakIndexingEngine.GetIndexedPeak(identification.peakfindingMass, j, mbrTol, z);
                                    IndexedMassSpectralPeak isotopePeak = _peakIndexingEngine.GetIndexedPeak(identification.peakfindingMass + Constants.C13MinusC12, j, mbrTol, z);
                                  
                                    if (peak != null && isotopePeak != null)
                                    {
                                        chargeXic.Add(peak);
                                    }
                                }

                                if (!chargeXic.Any())
                                {
                                    continue;
                                }

                                List<IsotopicEnvelope> chargeEnvelopes = GetIsotopicEnvelopes(chargeXic, identification, z, true);

                                if (!chargeEnvelopes.Any())
                                {
                                    chargeXic.Clear();
                                    continue;
                                }

                                claimedPeaks.Clear();
                                foreach (var env in chargeEnvelopes)
                                {
                                    if (claimedPeaks.Contains(env.IndexedPeak))
                                    {
                                        continue;
                                    }

                                    var acceptorPeak = new ChromatographicPeak(identification, true, idAcceptorFile);
                                    var xic = Peakfind(env.IndexedPeak.RetentionTime, identification.peakfindingMass, env.ChargeState, idAcceptorFile, mbrTol);
                                    List<IsotopicEnvelope> bestChargeEnvelopes = GetIsotopicEnvelopes(xic, identification, env.ChargeState, true);
                                    acceptorPeak.IsotopicEnvelopes.AddRange(bestChargeEnvelopes);
                                    acceptorPeak.CalculateIntensityForThisFeature(Integrate);
                                    CutPeak(acceptorPeak, env.IndexedPeak.RetentionTime);

                                    foreach (var env2 in acceptorPeak.IsotopicEnvelopes)
                                    {
                                        claimedPeaks.Add(env2.IndexedPeak);
                                    }

                                    if (!acceptorPeakOptions.Any(p => p.Apex.IndexedPeak == acceptorPeak.Apex.IndexedPeak))
                                    {
                                        acceptorPeakOptions.Add(acceptorPeak);
                                    }
                                }

                                chargeXic.Clear();
                            }

                            if (!acceptorPeakOptions.Any())
                            {
                                continue;
                            }

                            var ppmErrors = nearbyDataPoints.Select(p => p.AcceptorFilePeak.MassError).ToList();
                            var rtDiffs = nearbyDataPoints.Select(p => p.RtDiff).ToList();
                            var intensityFoldChanges = nearbyDataPoints.Select(p => Math.Log(p.AcceptorFilePeak.Intensity) - Math.Log(p.DonorFilePeak.Intensity)).ToList();

                            foreach (var peakOption in acceptorPeakOptions)
                            {
                                // matched ppm score
                                double xPpm = peakOption.MassError;
                                double muPpm = ppmErrors.Median();
                                double stddevPpm = ppmErrors.InterquartileRange() / 1.3489;
                                double ppmScore = 2 * Normal.CDF(muPpm, stddevPpm, xPpm);
                                if (ppmScore > 1)
                                {
                                    ppmScore = 2 - ppmScore;
                                }

                                // matched RT score
                                double xRt = peakOption.Apex.IndexedPeak.RetentionTime - donorPeak.Apex.IndexedPeak.RetentionTime;
                                double muRt = rtDiffs.Median();
                                double stddevRt = rtDiffs.InterquartileRange() / 1.3489;
                                double rtScore = 2 * Normal.CDF(muRt, stddevRt, xRt);
                                if (rtScore > 1)
                                {
                                    rtScore = 2 - rtScore;
                                }

                                peakOption.MbrScore = (ppmScore * rtScore) / acceptorPeakOptions.Count;

                                // matched intensity score (if the donor+acceptor files are technical replicates)
                                if (idAcceptorFile.Condition == idDonorFile.Condition
                                    && idAcceptorFile.BiologicalReplicate == idDonorFile.BiologicalReplicate
                                    && idAcceptorFile.Fraction == idDonorFile.Fraction)
                                {
                                    double xInt = Math.Log(peakOption.Intensity) - Math.Log(donorPeak.Intensity);
                                    double muInt = intensityFoldChanges.Median();
                                    double stddevInt = intensityFoldChanges.InterquartileRange() / 1.3489;
                                    double intensityFoldChangeScore = 2 * Normal.CDF(muInt, stddevInt, xInt);

                                    if (intensityFoldChangeScore > 1)
                                    {
                                        intensityFoldChangeScore = 2 - intensityFoldChangeScore;
                                    }

                                    peakOption.MbrScore = (ppmScore * rtScore * intensityFoldChangeScore) / acceptorPeakOptions.Count;
                                }
                            }

                            // accept the peak with the best (highest) score
                            acceptorPeakOptions.Sort((x, y) => y.MbrScore.CompareTo(x.MbrScore));
                            var bestAcceptorPeak = acceptorPeakOptions.First();

                            matchedPeaks[i] = bestAcceptorPeak;

                            // get alternative charge states
                            double minRt = bestAcceptorPeak.IsotopicEnvelopes.Min(p => p.IndexedPeak.RetentionTime);
                            double maxRt = bestAcceptorPeak.IsotopicEnvelopes.Max(p => p.IndexedPeak.RetentionTime);
                            var existingChromatogram = bestAcceptorPeak.IsotopicEnvelopes.Where(p => p.ChargeState == bestAcceptorPeak.Apex.ChargeState).OrderBy(p => p.IndexedPeak.RetentionTime).ToList();

                            foreach (var z in _chargeStates)
                            {
                                if (z == bestAcceptorPeak.Apex.ChargeState)
                                {
                                    continue;
                                }

                                var xic = Peakfind(bestAcceptorPeak.Apex.IndexedPeak.RetentionTime, identification.peakfindingMass, z, idAcceptorFile, mbrTol);
                                List<IsotopicEnvelope> bestChargeEnvelopes = GetIsotopicEnvelopes(xic, identification, z, true)
                                    .Where(p => p.IndexedPeak.RetentionTime <= maxRt && p.IndexedPeak.RetentionTime >= minRt).ToList();

                                if (!bestChargeEnvelopes.Any())
                                {
                                    continue;
                                }

                                double[] intensities = new double[existingChromatogram.Count];
                                for (int t = 0; t < existingChromatogram.Count; t++)
                                {
                                    var isot = bestChargeEnvelopes.FirstOrDefault(p => p.IndexedPeak.ZeroBasedMs1ScanIndex == existingChromatogram[t].IndexedPeak.ZeroBasedMs1ScanIndex);

                                    if (isot == null)
                                    {
                                        intensities[t] = 0;
                                    }
                                    else
                                    {
                                        intensities[t] = isot.Intensity;
                                    }
                                }

                                double spearman = Correlation.Spearman(existingChromatogram.Select(p => p.Intensity), intensities);

                                double max = bestChargeEnvelopes.Max(p => p.Intensity);
                                var maxEnv = bestChargeEnvelopes.First(p => p.Intensity == max);
                                bool similarTime = Math.Abs(maxEnv.IndexedPeak.ZeroBasedMs1ScanIndex - bestAcceptorPeak.Apex.IndexedPeak.ZeroBasedMs1ScanIndex) <= 1;

                                if (spearman > 0.5 || similarTime)
                                {
                                    bestAcceptorPeak.IsotopicEnvelopes.AddRange(bestChargeEnvelopes);
                                    bestAcceptorPeak.CalculateIntensityForThisFeature(Integrate);
                                    CutPeak(bestAcceptorPeak, bestAcceptorPeak.Apex.IndexedPeak.RetentionTime);
                                }
                            }

                            // used for calculating FDR
                            bool isKnownPeak = sequenceToRealPeakInAcceptorFile.ContainsKey(identification.ModifiedSequence);
                            if (isKnownPeak)
                            {
                                // get the "real" peak in this file
                                var realPeak = sequenceToRealPeakInAcceptorFile[identification.ModifiedSequence];

                                if (realPeak == null || realPeak.Intensity == 0)
                                {
                                    matchedPeaks[i] = null;
                                    continue;
                                }

                                bool isCorrectMatch = realPeak.Apex.IndexedPeak.Mz == bestAcceptorPeak.Apex.IndexedPeak.Mz
                                    && realPeak.Apex.IndexedPeak.RetentionTime == bestAcceptorPeak.Apex.IndexedPeak.RetentionTime;

                                scoreWithKnownResults[i] = (bestAcceptorPeak, isCorrectMatch);
                                matchedPeaks[i] = null;
                            }
                        }
                    });

                var res = scoreWithKnownResults.Where(p => p.bestAcceptorPeak != null).OrderBy(p => p.bestAcceptorPeak.MbrScore).ToList();

                int cumulativeTruePositives = res.Count(p => p.isTruePositive);
                int cumulativeFalsePositives = res.Count(p => !p.isTruePositive);
                double scoreCutoff = 1;

                foreach (var (bestAcceptorPeak, isTruePositive) in res)
                {
                    double qValue = (double)cumulativeFalsePositives / (cumulativeFalsePositives + cumulativeTruePositives);

                    if (qValue <= MbrQValueCutoff)
                    {
                        scoreCutoff = bestAcceptorPeak.MbrScore;
                        break;
                    }

                    if (isTruePositive)
                    {
                        cumulativeTruePositives--;
                    }
                    else
                    {
                        cumulativeFalsePositives--;
                    }
                }

                // calculate FDR
                foreach (var peak in matchedPeaks.Where(p => p != null))
                {
                    int cumulativeTrue = res.Count(p => p.isTruePositive && p.bestAcceptorPeak.MbrScore >= peak.MbrScore);
                    int cumulativeFalse = res.Count(p => !p.isTruePositive && p.bestAcceptorPeak.MbrScore >= peak.MbrScore);
                    double qValue = ((double)cumulativeFalse) / (cumulativeTrue + cumulativeFalse);

                    if (cumulativeFalse == 0)
                    {
                        qValue = 0;
                    }

                    peak.MbrQValue = qValue;
                }

                foreach (var peak in matchedPeaks.Where(p => p != null && p.MbrScore > scoreCutoff))
                {
                    if (bestMbrHits.TryGetValue(peak.Identifications.First().ModifiedSequence, out var oldPeak))
                    {
                        if (peak.MbrScore > oldPeak.MbrScore)
                        {
                            bestMbrHits[peak.Identifications.First().ModifiedSequence] = peak;
                        }
                    }
                    else
                    {
                        bestMbrHits.Add(peak.Identifications.First().ModifiedSequence, peak);
                    }
                }
            }

            // save MBR results
            _results.Peaks[idAcceptorFile].AddRange(bestMbrHits.Values);

            RunErrorChecking(idAcceptorFile);
        }

        private RetentionTimeCalibDataPoint[] GetRtCalSpline(SpectraFileInfo donor, SpectraFileInfo acceptor)
        {
            var donorFileBestMsmsPeaks = new Dictionary<string, ChromatographicPeak>();
            var acceptorFileBestMsmsPeaks = new Dictionary<string, ChromatographicPeak>();
            var rtCalibrationCurve = new List<RetentionTimeCalibDataPoint>();

            // get all peaks, not counting ambiguous peaks
            IEnumerable<ChromatographicPeak> donorPeaks = _results.Peaks[donor].Where(p => p.Apex != null && !p.IsMbrPeak && p.NumIdentificationsByFullSeq == 1);
            IEnumerable<ChromatographicPeak> acceptorPeaks = _results.Peaks[acceptor].Where(p => p.Apex != null && !p.IsMbrPeak && p.NumIdentificationsByFullSeq == 1);

            // get the best (most intense) peak for each peptide in the acceptor file
            foreach (ChromatographicPeak acceptorPeak in acceptorPeaks)
            {
                if (acceptorFileBestMsmsPeaks.TryGetValue(acceptorPeak.Identifications.First().ModifiedSequence, out ChromatographicPeak currentBestPeak))
                {
                    if (currentBestPeak.Intensity > acceptorPeak.Intensity)
                    {
                        acceptorFileBestMsmsPeaks[acceptorPeak.Identifications.First().ModifiedSequence] = acceptorPeak;
                    }
                }
                else
                {
                    acceptorFileBestMsmsPeaks.Add(acceptorPeak.Identifications.First().ModifiedSequence, acceptorPeak);
                }
            }

            // get the best (most intense) peak for each peptide in the donor file
            foreach (ChromatographicPeak donorPeak in donorPeaks)
            {
                if (donorFileBestMsmsPeaks.TryGetValue(donorPeak.Identifications.First().ModifiedSequence, out ChromatographicPeak currentBestPeak))
                {
                    if (currentBestPeak.Intensity > donorPeak.Intensity)
                    {
                        donorFileBestMsmsPeaks[donorPeak.Identifications.First().ModifiedSequence] = donorPeak;
                    }
                }
                else
                {
                    donorFileBestMsmsPeaks.Add(donorPeak.Identifications.First().ModifiedSequence, donorPeak);
                }
            }

            // create RT calibration curve
            foreach (var peak in acceptorFileBestMsmsPeaks)
            {
                ChromatographicPeak acceptorFilePeak = peak.Value;

                if (donorFileBestMsmsPeaks.TryGetValue(peak.Key, out ChromatographicPeak donorFilePeak))
                {
                    rtCalibrationCurve.Add(new RetentionTimeCalibDataPoint(donorFilePeak, acceptorFilePeak));
                }
            }

            return rtCalibrationCurve.OrderBy(p => p.DonorFilePeak.Apex.IndexedPeak.RetentionTime).ToArray();
        }

        private void RunErrorChecking(SpectraFileInfo spectraFile)
        {
            if (!Silent)
            {
                Console.WriteLine("Checking errors");
            }

            _results.Peaks[spectraFile].RemoveAll(p => p == null || p.IsMbrPeak && !p.IsotopicEnvelopes.Any());

            // merge duplicate peaks and handle MBR/MSMS peakfinding conflicts
            var peaksGroupedByApex = new Dictionary<IndexedMassSpectralPeak, ChromatographicPeak>();
            var peaks = new List<ChromatographicPeak>();
            foreach (ChromatographicPeak tryPeak in _results.Peaks[spectraFile].OrderBy(p => p.IsMbrPeak))
            {
                tryPeak.CalculateIntensityForThisFeature(Integrate);
                tryPeak.ResolveIdentifications();

                if (tryPeak.Apex == null)
                {
                    if (tryPeak.IsMbrPeak)
                    {
                        continue;
                    }

                    peaks.Add(tryPeak);
                    continue;
                }

                IndexedMassSpectralPeak apexPeak = tryPeak.Apex.IndexedPeak;
                if (peaksGroupedByApex.TryGetValue(apexPeak, out ChromatographicPeak storedPeak))
                {
                    if (tryPeak.IsMbrPeak && storedPeak == null)
                    {
                        continue;
                    }

                    if (!tryPeak.IsMbrPeak && !storedPeak.IsMbrPeak)
                    {
                        storedPeak.MergeFeatureWith(tryPeak, Integrate);
                    }
                    else if (tryPeak.IsMbrPeak && !storedPeak.IsMbrPeak)
                    {
                        continue;
                    }
                    else if (tryPeak.IsMbrPeak && storedPeak.IsMbrPeak)
                    {
                        if (tryPeak.Identifications.First().ModifiedSequence == storedPeak.Identifications.First().ModifiedSequence)
                        {
                            storedPeak.MergeFeatureWith(tryPeak, Integrate);
                        }
                        else if (storedPeak.MbrScore < tryPeak.MbrScore)
                        {
                            peaksGroupedByApex[tryPeak.Apex.IndexedPeak] = tryPeak;
                        }
                    }
                }
                else
                {
                    peaksGroupedByApex.Add(apexPeak, tryPeak);
                }
            }

            peaks.AddRange(peaksGroupedByApex.Values.Where(p => p != null));
            _results.Peaks[spectraFile] = peaks;

            // merge multiple peaks for the same peptide within a time window
            peaks = _results.Peaks[spectraFile].Where(p => p.NumIdentificationsByFullSeq > 1).ToList();
            var temp = _results.Peaks[spectraFile].Where(p => p.NumIdentificationsByFullSeq == 1).GroupBy(p => p.Identifications.First().ModifiedSequence).ToList();

            foreach (var sequence in temp)
            {
                if (sequence.Count() == 1)
                {
                    peaks.Add(sequence.First());
                    continue;
                }
                var temp2 = sequence.Where(p => p.Apex != null).ToList();
                if (temp2.Count == 0)
                {
                    peaks.Add(sequence.First());
                    continue;
                }
                var merged = new HashSet<ChromatographicPeak>();
                foreach (ChromatographicPeak peak in temp2)
                {
                    if (merged.Contains(peak))
                    {
                        continue;
                    }

                    var toMerge = temp2.Where(p => Math.Abs(p.Apex.IndexedPeak.RetentionTime - peak.Apex.IndexedPeak.RetentionTime) < RtTol && !merged.Contains(p) && p != peak);
                    foreach (ChromatographicPeak peakToMerge in toMerge)
                    {
                        peak.MergeFeatureWith(peakToMerge, Integrate);
                        merged.Add(peakToMerge);
                    }

                    peaks.Add(peak);
                }
            }

            _results.Peaks[spectraFile] = peaks;
        }

        public List<IsotopicEnvelope> GetIsotopicEnvelopes(List<IndexedMassSpectralPeak> xic, Identification identification, int chargeState, bool matchBetweenRuns)
        {
            var isotopicEnvelopes = new List<IsotopicEnvelope>();
            var isotopeMassShifts = _modifiedSequenceToIsotopicDistribution[identification.ModifiedSequence];

            if (isotopeMassShifts.Count < NumIsotopesRequired)
            {
                return isotopicEnvelopes;
            }

            PpmTolerance isotopeTolerance = new PpmTolerance(IsotopePpmTolerance);

            var experimentalIsotopeAbundances = new double[isotopeMassShifts.Count];
            var theoreticalIsotopeAbundances = isotopeMassShifts.Select(p => p.Item2).ToArray();

            foreach (IndexedMassSpectralPeak peak in xic)
            {
                // calculate theoretical isotopes relative to observed peak
                Array.Clear(experimentalIsotopeAbundances, 0, experimentalIsotopeAbundances.Length);

                double observedMass = peak.Mz.ToMass(chargeState);
                double theorMass = identification.peakfindingMass;
                double mainPeakError = observedMass - theorMass;


                if (matchBetweenRuns)
                {
                    double unexpectedIsotopeMass = observedMass - Constants.C13MinusC12;
                    var unexpectedPeak = _peakIndexingEngine.GetIndexedPeak(unexpectedIsotopeMass, peak.ZeroBasedMs1ScanIndex, isotopeTolerance, chargeState);
                    bool unexpectedIsotopePeakObserved = unexpectedPeak != null && unexpectedPeak.Intensity / peak.Intensity > 0.3;

                    if (unexpectedIsotopePeakObserved)
                    {
                        continue;
                    }
                }

                for (int t = 0; t < theoreticalIsotopeAbundances.Length; t++)
                {
                    double isotopeMass = observedMass + t * Constants.C13MinusC12;
                    var isotopePeak = _peakIndexingEngine.GetIndexedPeak(isotopeMass, peak.ZeroBasedMs1ScanIndex, isotopeTolerance, chargeState);

                    if (isotopePeak == null)
                    {
                        continue;
                    }

                    experimentalIsotopeAbundances[t] = isotopePeak.Intensity;
                }

                if (experimentalIsotopeAbundances[0] == 0 && RequireMonoisotopicMass)
                {
                    continue;
                }

                int numIsotopePeaksObserved = 0;
                int mainPeakIndex = (int)Math.Round(observedMass - identification.monoisotopicMass, 0);
                for (int i = mainPeakIndex; i < experimentalIsotopeAbundances.Length; i++)
                {
                    if (experimentalIsotopeAbundances[i] > 0)
                    {
                        numIsotopePeaksObserved++;
                    }
                    else
                    {
                        break;
                    }
                }
                for (int i = mainPeakIndex - 1; i >= 0; i--)
                {
                    if (experimentalIsotopeAbundances[i] > 0)
                    {
                        numIsotopePeaksObserved++;
                    }
                    else
                    {
                        break;
                    }
                }

                if (numIsotopePeaksObserved >= NumIsotopesRequired)
                {
                    int theoreticalMostAbundantIsotopeIndex = Array.IndexOf(theoreticalIsotopeAbundances, 1.0);
                    int replacementsAllowed = 2;

                    for (int i = 0; i < experimentalIsotopeAbundances.Length; i++)
                    {
                        double expectedIntensity = theoreticalIsotopeAbundances[i] * experimentalIsotopeAbundances[theoreticalMostAbundantIsotopeIndex];
                        double measuredIntensity = experimentalIsotopeAbundances[i];

                        if (measuredIntensity < expectedIntensity * 0.5 || measuredIntensity > expectedIntensity * 2 && replacementsAllowed >= 1)
                        {
                            experimentalIsotopeAbundances[i] = expectedIntensity;
                            replacementsAllowed--;
                        }
                    }

                    double corr = Correlation.Pearson(experimentalIsotopeAbundances, theoreticalIsotopeAbundances);

                    if (corr > 0.7)
                    {
                        isotopicEnvelopes.Add(new IsotopicEnvelope(peak, chargeState, experimentalIsotopeAbundances.Sum()));
                    }
                }
            }

            return isotopicEnvelopes;
        }

        public List<IndexedMassSpectralPeak> Peakfind(double idRetentionTime, double mass, int charge, SpectraFileInfo spectraFileInfo, Tolerance tolerance)
        {
            var xic = new List<IndexedMassSpectralPeak>();

            // get precursor scan to start at
            Ms1ScanInfo[] ms1Scans = _ms1Scans[spectraFileInfo];
            int precursorScanIndex = -1;
            foreach (Ms1ScanInfo ms1Scan in ms1Scans)
            {
                if (ms1Scan.RetentionTime < idRetentionTime)
                {
                    precursorScanIndex = ms1Scan.ZeroBasedMs1ScanIndex;
                }
                else
                {
                    break;
                }
            }

            // go right
            int missedScans = 0;
            for (int t = precursorScanIndex; t < ms1Scans.Length; t++)
            {
                var peak = _peakIndexingEngine.GetIndexedPeak(mass, t, tolerance, charge);

                if (peak == null && t != precursorScanIndex)
                {
                    missedScans++;
                }
                else if (peak != null)
                {
                    missedScans = 0;
                    xic.Add(peak);
                }

                if (missedScans > MissedScansAllowed)
                {
                    break;
                }
            }

            // go left
            missedScans = 0;
            for (int t = precursorScanIndex - 1; t >= 0; t--)
            {
                var peak = _peakIndexingEngine.GetIndexedPeak(mass, t, tolerance, charge);

                if (peak == null && t != precursorScanIndex)
                {
                    missedScans++;
                }
                else if (peak != null)
                {
                    missedScans = 0;
                    xic.Add(peak);
                }

                if (missedScans > MissedScansAllowed)
                {
                    break;
                }
            }

            // sort by RT
            xic.Sort((x, y) => x.RetentionTime.CompareTo(y.RetentionTime));

            return xic;
        }

        private void CutPeak(ChromatographicPeak peak, double identificationTime)
        {
            // find out if we need to split this peak by using the discrimination factor
            // this method assumes that the isotope envelopes in a chromatographic peak are already sorted by MS1 scan number
            bool cutThisPeak = false;

            if (peak.IsotopicEnvelopes.Count < 5)
            {
                return;
            }

            var timePointsForApexZ = peak.IsotopicEnvelopes.Where(p => p.ChargeState == peak.Apex.ChargeState).ToList();
            HashSet<int> scanNumbers = new HashSet<int>(timePointsForApexZ.Select(p => p.IndexedPeak.ZeroBasedMs1ScanIndex));
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

                    double discriminationFactor =
                        (timepoint.Intensity - valleyTimePoint.Intensity) / timepoint.Intensity;

                    if (discriminationFactor > MinDiscFactorToCutAt &&
                        (indexOfValley + iter < timePointsForApexZ.Count && indexOfValley + iter >= 0))
                    {
                        IsotopicEnvelope secondValleyTimepoint = timePointsForApexZ[indexOfValley + iter];

                        discriminationFactor =
                            (timepoint.Intensity - secondValleyTimepoint.Intensity) / timepoint.Intensity;

                        if (discriminationFactor > MinDiscFactorToCutAt)
                        {
                            cutThisPeak = true;
                            break;
                        }

                        int nextMs1ScanNum = -1;
                        for (int j = valleyTimePoint.IndexedPeak.ZeroBasedMs1ScanIndex - 1;
                            j < _ms1Scans[peak.SpectraFileInfo].Length && j >= 0;
                            j += iter)
                        {
                            if (_ms1Scans[peak.SpectraFileInfo][j].OneBasedScanNumber >= 0 &&
                                _ms1Scans[peak.SpectraFileInfo][j].OneBasedScanNumber !=
                                valleyTimePoint.IndexedPeak.ZeroBasedMs1ScanIndex)
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

                if (cutThisPeak)
                {
                    break;
                }
            }

            // cut
            if (cutThisPeak)
            {
                if (identificationTime > valleyTimePoint.IndexedPeak.RetentionTime)
                {
                    // MS2 identification is to the right of the valley; remove all peaks left of the valley
                    peak.IsotopicEnvelopes.RemoveAll(p => p.IndexedPeak.RetentionTime <= valleyTimePoint.IndexedPeak.RetentionTime);
                }
                else
                {
                    // MS2 identification is to the left of the valley; remove all peaks right of the valley
                    peak.IsotopicEnvelopes.RemoveAll(p => p.IndexedPeak.RetentionTime >= valleyTimePoint.IndexedPeak.RetentionTime);
                }

                // recalculate intensity for the peak
                peak.CalculateIntensityForThisFeature(Integrate);
                peak.SplitRT = valleyTimePoint.IndexedPeak.RetentionTime;

                // recursively cut
                CutPeak(peak, identificationTime);
            }
        }
    }
}