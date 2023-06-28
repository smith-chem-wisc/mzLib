using Chemistry;
using MassSpectrometry;
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
        public readonly double IsotopePpmTolerance;
        public readonly bool Integrate;
        public readonly int MissedScansAllowed;
        public readonly int NumIsotopesRequired;
        public readonly bool IdSpecificChargeState;
        public readonly bool Normalize;
        public readonly double DiscriminationFactorToCutPeak;
        public readonly bool QuantifyAmbiguousPeptides;

        // MBR settings
        public readonly bool MatchBetweenRuns;

        public readonly double MbrRtWindow;
        public readonly double MbrPpmTolerance;
        public readonly bool RequireMsmsIdInCondition;

        // settings for the Bayesian protein quantification engine
        public readonly bool BayesianProteinQuant;

        public readonly string ProteinQuantBaseCondition;
        public readonly double ProteinQuantFoldChangeCutoff;
        public readonly int McmcSteps;
        public readonly int McmcBurninSteps;
        public readonly bool UseSharedPeptidesForProteinQuant;
        public readonly bool PairedSamples;
        public readonly int? RandomSeed;

        // structures used in the FlashLFQ engine
        private List<SpectraFileInfo> _spectraFileInfo;

        private Stopwatch _globalStopwatch;
        private List<Identification> _allIdentifications;
        /// <summary>
        /// Dictionary linking a modified sequence to a List of tuples containing
        /// the mass shifts (isotope mass - monoisotopic mass) and normalized abundances for the
        /// isotopes for a given peptide
        /// </summary>
        private Dictionary<string, List<(double massShift, double normalizedAbundance)>> _modifiedSequenceToIsotopicDistribution;
        private List<int> _chargeStates;
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
            bool quantifyAmbiguousPeptides = false,
            bool silent = false,
            int maxThreads = -1,

            // MBR settings
            bool matchBetweenRuns = false,
            double matchBetweenRunsPpmTolerance = 10.0,
            double maxMbrWindow = 2.5,
            bool requireMsmsIdInCondition = false,

            // settings for the Bayesian protein quantification engine
            bool bayesianProteinQuant = false,
            string proteinQuantBaseCondition = null,
            double proteinQuantFoldChangeCutoff = 0.1,
            int mcmcSteps = 3000,
            int mcmcBurninSteps = 1000,
            bool useSharedPeptidesForProteinQuant = false,
            bool pairedSamples = false,
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
            QuantifyAmbiguousPeptides = quantifyAmbiguousPeptides;
            Silent = silent;
            IdSpecificChargeState = idSpecificChargeState;
            MbrRtWindow = maxMbrWindow;
            RequireMsmsIdInCondition = requireMsmsIdInCondition;
            Normalize = normalize;
            MaxThreads = maxThreads;
            BayesianProteinQuant = bayesianProteinQuant;
            PairedSamples = pairedSamples;
            ProteinQuantBaseCondition = proteinQuantBaseCondition;
            ProteinQuantFoldChangeCutoff = proteinQuantFoldChangeCutoff;
            McmcSteps = mcmcSteps;
            McmcBurninSteps = mcmcBurninSteps;
            UseSharedPeptidesForProteinQuant = useSharedPeptidesForProteinQuant;
            RandomSeed = randomSeed;

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
            DiscriminationFactorToCutPeak = 0.6;
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
                new IntensityNormalizationEngine(_results, Integrate, Silent, MaxThreads).NormalizeResults();
            }

            // calculate peptide intensities
            _results.CalculatePeptideResults(QuantifyAmbiguousPeptides);

            // do top3 protein quantification
            _results.CalculateProteinResultsMedianPolish(UseSharedPeptidesForProteinQuant);

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
                    if (!Silent)
                    {
                        Console.WriteLine("Running Bayesian protein quantification analysis");
                    }

                    new ProteinQuantificationEngine(_results, MaxThreads, ProteinQuantBaseCondition, UseSharedPeptidesForProteinQuant,
                        ProteinQuantFoldChangeCutoff, RandomSeed, McmcBurninSteps, McmcSteps, PairedSamples).Run();
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

        public PeakIndexingEngine GetIndexingEngine()
        {
            return _peakIndexingEngine;
        }
        /// <summary>
        /// Creates a theoretical isotope distribution for each of the identified sequences
        /// If the sequence is modified and the modification has an unknown chemical formula,
        /// averagine is used for the modified part
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
                    double massDiff = id.MonoisotopicMass;
                    if (!String.IsNullOrEmpty(id.BaseSequence))
                    {
                        Proteomics.AminoAcidPolymer.Peptide baseSequence = new Proteomics.AminoAcidPolymer.Peptide(id.BaseSequence);
                        formula = baseSequence.GetChemicalFormula();
                        // add averagine for any unknown mass difference (i.e., a modification)
                        massDiff -= baseSequence.MonoisotopicMass;

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
                    else
                    {
                        double averagines = massDiff / averagineMass;
                        string averagineFormulaString = 
                            "C" + (int)Math.Round(averagines * averageC, 0) +
                            "H" + (int)Math.Round(averagines * averageH, 0) +
                            "O" + (int)Math.Round(averagines * averageO, 0) +
                            "N" + (int)Math.Round(averagines * averageN, 0) +
                            "S" + (int)Math.Round(averagines * averageS, 0);
                        formula = ChemicalFormula.ParseFormula(averagineFormulaString);
                    }
                }

                var isotopicDistribution = IsotopicDistribution.GetDistribution(formula, 0.125, 1e-8);

                double[] masses = isotopicDistribution.Masses.ToArray();
                double[] abundances = isotopicDistribution.Intensities.ToArray();

                for (int i = 0; i < masses.Length; i++)
                {
                    masses[i] += (id.MonoisotopicMass - formula.MonoisotopicMass);
                }

                double highestAbundance = abundances.Max();
                int highestAbundanceIndex = Array.IndexOf(abundances, highestAbundance);

                for (int i = 0; i < masses.Length; i++)
                {
                    // expected isotopic mass shifts for this peptide
                    masses[i] -= id.MonoisotopicMass;

                    // normalized abundance of each isotope
                    abundances[i] /= highestAbundance;

                    // look for these isotopes
                    if (isotopicMassesAndNormalizedAbundances.Count < NumIsotopesRequired || abundances[i] > 0.1)
                    {
                        isotopicMassesAndNormalizedAbundances.Add((masses[i], abundances[i]));
                    }
                }

                _modifiedSequenceToIsotopicDistribution.Add(id.ModifiedSequence, isotopicMassesAndNormalizedAbundances);
            }

            var minChargeState = _allIdentifications.Min(p => p.PrecursorChargeState);
            var maxChargeState = _allIdentifications.Max(p => p.PrecursorChargeState);
            _chargeStates = Enumerable.Range(minChargeState, (maxChargeState - minChargeState) + 1).ToList();

            var peptideModifiedSequences = _allIdentifications.GroupBy(p => p.ModifiedSequence);
            foreach (var identifications in peptideModifiedSequences)
            {
                // isotope where normalized abundance is 1
                double mostAbundantIsotopeShift = _modifiedSequenceToIsotopicDistribution[identifications.First().ModifiedSequence]
                    .First(p => p.Item2 == 1.0).Item1;

                foreach (Identification identification in identifications)
                {
                    identification.PeakfindingMass = identification.MonoisotopicMass + mostAbundantIsotopeShift;
                }
            }
        }

        /// <summary>
        /// Creates an ChromatographicPeak for each MS2 ID in a given file. Works by first
        /// finding every MS1 scan that neighbors the MS2 scan the ID originated from and that
        /// contains the peak finding mass (most abundant isotope), then finds every isotope peak within that scan.
        /// Isotope peak intensities are summed and an IsotopicEnvelope object is created from the summed intensities.
        /// Multiple IsotopicEnvelopes are associated with each ChromatographicPeak (corresponding to different scans 
        /// and different charge states)
        /// </summary>
        /// <param name="fileInfo">File to be quantified</param>
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

            Tolerance peakfindingTol = new PpmTolerance(PeakfindingPpmTolerance); // Peak finding tolerance is generally higher than ppmTolerance
            Tolerance ppmTolerance = new PpmTolerance(PpmTolerance);
            ChromatographicPeak[] chromatographicPeaks = new ChromatographicPeak[ms2IdsForThisFile.Count];

            Parallel.ForEach(Partitioner.Create(0, ms2IdsForThisFile.Count),
                new ParallelOptions { MaxDegreeOfParallelism = MaxThreads },
                (range, loopState) =>
                {
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
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
                            List<IndexedMassSpectralPeak> xic = Peakfind(
                                    identification.Ms2RetentionTimeInMinutes,
                                    identification.PeakfindingMass,
                                    chargeState,
                                    identification.FileInfo,
                                    peakfindingTol)
                                .OrderBy(p => p.RetentionTime)
                                .ToList();

                            // filter by smaller mass tolerance
                            xic.RemoveAll(p => 
                                !ppmTolerance.Within(p.Mz.ToMass(chargeState), identification.PeakfindingMass));

                            // filter by isotopic distribution
                            List<IsotopicEnvelope> isotopicEnvelopes = GetIsotopicEnvelopes(xic, identification, chargeState);

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

        /// <summary>
        /// This method maps identified peaks from other chromatographic runs ("donors") onto
        /// the defined chromatographic run ("acceptor"). The goal is to reduce the number of missing
        /// intensity measurements. Missing values occur generally either because 1) the analyte is
        /// in the sample but didn't get fragmented/identified or 2) the analyte is genuinely missing
        /// from the sample.
        /// </summary>
        private void QuantifyMatchBetweenRunsPeaks(SpectraFileInfo idAcceptorFile)
        {
            bool acceptorSampleIsFractionated = _results.SpectraFiles
                .Where(p => p.Condition == idAcceptorFile.Condition && p.BiologicalReplicate == idAcceptorFile.BiologicalReplicate)
                .Select(p => p.Fraction)
                .Distinct()
                .Count() > 1;

            // acceptor file known peaks
            var acceptorFileIdentifiedPeaks = _results.Peaks[idAcceptorFile];
            var apexToAcceptorFilePeak = new Dictionary<IndexedMassSpectralPeak, ChromatographicPeak>();

            List<double> ppmErrors = new List<double>();
            foreach (var peak in acceptorFileIdentifiedPeaks.Where(p => p.Apex != null))
            {
                if (!apexToAcceptorFilePeak.ContainsKey(peak.Apex.IndexedPeak))
                {
                    apexToAcceptorFilePeak.Add(peak.Apex.IndexedPeak, peak);
                }

                ppmErrors.Add(peak.MassError);
            }

            if (ppmErrors.Count < 3)
            {
                return;
            }

            double ppmSpread = ppmErrors.Count > 30 ? ppmErrors.InterquartileRange() / 1.36 : ppmErrors.StandardDeviation();

            Normal ppmDistribution = new Normal(ppmErrors.Median(), ppmSpread);

            double filespecificMbrPpmTolerance = Math.Min(Math.Abs(ppmErrors.Median()) + ppmSpread * 4, MbrPpmTolerance);

            // match between runs PPM tolerance
            Tolerance mbrTol = new PpmTolerance(filespecificMbrPpmTolerance);

            // deserialize the file's indexed mass spectral peaks. these were stored and serialized to a file earlier
            _peakIndexingEngine.DeserializeIndex(idAcceptorFile);

            // these are the analytes already identified in this run. we don't need to try to match them from other runs
            var acceptorFileIdentifiedSequences = new HashSet<string>(acceptorFileIdentifiedPeaks.Where(p => p.IsotopicEnvelopes.Any())
                .SelectMany(p => p.Identifications.Select(d => d.ModifiedSequence)));

            var acceptorFileLogIntensities = acceptorFileIdentifiedPeaks.Where(p => p.Intensity > 0).Select(p => Math.Log(p.Intensity, 2)).ToList();
            double medianAcceptorLogIntensity = acceptorFileLogIntensities.Median();
            Normal intensityDistribution = new Normal(acceptorFileLogIntensities.Median(), acceptorFileLogIntensities.InterquartileRange() / 1.36);

            HashSet<ProteinGroup> thisFilesMsmsIdentifiedProteins = new HashSet<ProteinGroup>();
            if (RequireMsmsIdInCondition)
            {
                // only match peptides from proteins that have at least one MS/MS identified peptide in the condition
                foreach (SpectraFileInfo conditionFile in _spectraFileInfo.Where(p => p.Condition == idAcceptorFile.Condition))
                {
                    foreach (ProteinGroup proteinGroup in _results.Peaks[conditionFile].Where(p => !p.IsMbrPeak).SelectMany(p => p.Identifications.SelectMany(v => v.ProteinGroups)))
                    {
                        thisFilesMsmsIdentifiedProteins.Add(proteinGroup);
                    }
                }
            }

            // this stores the results of MBR
            var matchBetweenRunsIdentifiedPeaks = new Dictionary<string, Dictionary<IsotopicEnvelope, List<ChromatographicPeak>>>();

            // map each donor file onto this file
            foreach (SpectraFileInfo idDonorFile in _spectraFileInfo)
            {
                if (idAcceptorFile.Equals(idDonorFile))
                {
                    continue;
                }

                // this is the list of peaks identified in the other file but not in this one ("ID donor peaks")
                List<ChromatographicPeak> idDonorPeaks = _results.Peaks[idDonorFile].Where(p =>
                    !p.IsMbrPeak
                    && p.NumIdentificationsByFullSeq == 1
                    && p.IsotopicEnvelopes.Any()
                    && !acceptorFileIdentifiedSequences.Contains(p.Identifications.First().ModifiedSequence)
                    && (!RequireMsmsIdInCondition || p.Identifications.Any(v => v.ProteinGroups.Any(g => thisFilesMsmsIdentifiedProteins.Contains(g))))).ToList();

                if (!idDonorPeaks.Any())
                {
                    continue;
                }

                bool donorSampleIsFractionated = _results.SpectraFiles
                    .Where(p => p.Condition == idDonorFile.Condition && p.BiologicalReplicate == idDonorFile.BiologicalReplicate)
                    .Select(p => p.Fraction)
                    .Distinct()
                    .Count() > 1;

                var donorFileLogIntensities = idDonorPeaks.Where(p => p.Intensity > 0).Select(p => Math.Log(p.Intensity, 2)).ToList();
                double medianDonorLogIntensity = donorFileLogIntensities.Median();

                // generate RT calibration curve
                RetentionTimeCalibDataPoint[] rtCalibrationCurve = GetRtCalSpline(idDonorFile, idAcceptorFile);

                // this intensity score creates a conservative bias in MBR
                List<double> listOfFoldChangesBetweenTheFiles = new List<double>();

                if (_spectraFileInfo.Select(p => p.Condition).Distinct().Count() > 1
                    && idDonorFile.Condition != idAcceptorFile.Condition)
                {
                    var acceptorFileBestMsmsPeaks = new Dictionary<string, ChromatographicPeak>();

                    IEnumerable<ChromatographicPeak> acceptorPeaks = _results
                        .Peaks[idAcceptorFile]
                        .Where(p => p.Apex != null && !p.IsMbrPeak && p.NumIdentificationsByFullSeq == 1);

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

                    foreach (var donorPeak in idDonorPeaks)
                    {
                        double donorPeakIntensity = donorPeak.Intensity;
                        if (acceptorFileBestMsmsPeaks.TryGetValue(donorPeak.Identifications.First().ModifiedSequence, out var acceptorPeak))
                        {
                            double acceptorPeakIntensity = acceptorPeak.Intensity;

                            double intensityLogFoldChange = Math.Log(acceptorPeakIntensity, 2) - Math.Log(donorPeakIntensity, 2);

                            listOfFoldChangesBetweenTheFiles.Add(intensityLogFoldChange);
                        }
                    }
                }

                Normal foldChangeDistribution = null;
                if (listOfFoldChangesBetweenTheFiles.Count > 100)
                {
                    new Normal(listOfFoldChangesBetweenTheFiles.Median(), listOfFoldChangesBetweenTheFiles.StandardDeviation());
                }

                Parallel.ForEach(Partitioner.Create(0, idDonorPeaks.Count),
                    new ParallelOptions { MaxDegreeOfParallelism = MaxThreads },
                    (range, loopState) =>
                    {
                        var nearbyCalibrationPoints = new List<RetentionTimeCalibDataPoint>();
                        var matchBetweenRunsIdentifiedPeaksThreadSpecific = new Dictionary<string, Dictionary<IsotopicEnvelope, List<ChromatographicPeak>>>();

                        for (int i = range.Item1; i < range.Item2; i++)
                        {
                            nearbyCalibrationPoints.Clear();

                            // only compare +- 1 fraction
                            if (acceptorSampleIsFractionated && donorSampleIsFractionated)
                            {
                                int acceptorFractionNumber = idAcceptorFile.Fraction;
                                int donorFractionNumber = idDonorFile.Fraction;

                                if (Math.Abs(acceptorFractionNumber - donorFractionNumber) > 1)
                                {
                                    continue;
                                }
                            }

                            ChromatographicPeak donorPeak = idDonorPeaks[i];
                            Identification donorIdentification = donorPeak.Identifications.OrderBy(p => p.PosteriorErrorProbability).First();

                            // binary search for this donor peak in the retention time calibration spline
                            RetentionTimeCalibDataPoint testPoint = new RetentionTimeCalibDataPoint(donorPeak, null);
                            int index = Array.BinarySearch(rtCalibrationCurve, testPoint);

                            if (index < 0)
                            {
                                index = ~index;
                            }
                            if (index >= rtCalibrationCurve.Length && index >= 1)
                            {
                                index = rtCalibrationCurve.Length - 1;
                            }

                            // gather nearby data points
                            for (int r = index; r < rtCalibrationCurve.Length; r++)
                            {
                                double rtDiff = rtCalibrationCurve[r].DonorFilePeak.Apex.IndexedPeak.RetentionTime - donorPeak.Apex.IndexedPeak.RetentionTime;

                                if (Math.Abs(rtDiff) < 0.5)
                                {
                                    nearbyCalibrationPoints.Add(rtCalibrationCurve[r]);
                                }
                                else
                                {
                                    break;
                                }
                            }

                            for (int r = index - 1; r >= 0; r--)
                            {
                                double rtDiff = rtCalibrationCurve[r].DonorFilePeak.Apex.IndexedPeak.RetentionTime - donorPeak.Apex.IndexedPeak.RetentionTime;

                                if (Math.Abs(rtDiff) < 0.5)
                                {
                                    nearbyCalibrationPoints.Add(rtCalibrationCurve[r]);
                                }
                                else
                                {
                                    break;
                                }
                            }

                            if (!nearbyCalibrationPoints.Any())
                            {
                                continue;
                            }

                            // calculate difference between acceptor and donor RTs for these RT region
                            List<double> rtDiffs = nearbyCalibrationPoints
                                .Select(p => p.AcceptorFilePeak.Apex.IndexedPeak.RetentionTime - p.DonorFilePeak.Apex.IndexedPeak.RetentionTime)
                                .ToList();

                            //double medianIntensityDiff = nearbyCalibrationPoints
                            //    .Select(p => Math.Log(p.AcceptorFilePeak.Intensity, 2) - Math.Log(p.DonorFilePeak.Intensity, 2))
                            //    .Median();

                            // figure out the range of RT differences between the files that are "reasonable", centered around the median difference
                            double median = rtDiffs.Median();

                            // default range (if only 1 datapoint, or SD is 0, range is very high, etc)
                            double rtRange = MbrRtWindow;
                            double? rtStdDev = null;
                            double? rtInterquartileRange = null;

                            if (nearbyCalibrationPoints.Count < 6 && nearbyCalibrationPoints.Count > 1 && rtDiffs.StandardDeviation() > 0)
                            {
                                rtStdDev = rtDiffs.StandardDeviation();
                                rtRange = (double)rtStdDev * 6.0; // Multiplication inherited from legacy code, unsure of reason for 6
                            }
                            else if (nearbyCalibrationPoints.Count >= 6 && rtDiffs.InterquartileRange() > 0)
                            {
                                rtInterquartileRange = rtDiffs.InterquartileRange();
                                rtRange = (double)rtInterquartileRange * 4.5; // Multiplication inherited from legacy code, unsure of reason for 4.5
                            }

                            rtRange = Math.Min(rtRange, MbrRtWindow);

                            // this is the RT in the acceptor file to look around to find this analyte
                            double acceptorFileRtHypothesis = donorPeak.Apex.IndexedPeak.RetentionTime + median;
                            double lowerRangeRtHypothesis = acceptorFileRtHypothesis - (rtRange / 2.0);
                            double upperRangeRtHypothesis = acceptorFileRtHypothesis + (rtRange / 2.0);

                            Normal rtScoringDistribution = new Normal(acceptorFileRtHypothesis, rtRange / 6);

                            // get the MS1 scan info for this region so we can look up indexed peaks
                            Ms1ScanInfo[] ms1ScanInfos = _ms1Scans[idAcceptorFile];
                            Ms1ScanInfo start = ms1ScanInfos[0];
                            Ms1ScanInfo end = ms1ScanInfos[ms1ScanInfos.Length - 1];
                            for (int j = 0; j < ms1ScanInfos.Length; j++)
                            {
                                Ms1ScanInfo scan = ms1ScanInfos[j];
                                if (scan.RetentionTime <= lowerRangeRtHypothesis)
                                {
                                    start = scan;
                                }

                                if (scan.RetentionTime >= upperRangeRtHypothesis)
                                {
                                    end = scan;
                                    break;
                                }
                            }

                            // now we've identified the region in the chromatography this analyte should appear.
                            // we need to check for peaks in the region using ppm tolerance and isotope pattern matching
                            var chargesToMatch = donorPeak.Identifications.Select(p => p.PrecursorChargeState).Distinct().ToList();
                            if (!chargesToMatch.Contains(donorPeak.Apex.ChargeState))
                            {
                                chargesToMatch.Add(donorPeak.Apex.ChargeState);
                            }

                            foreach (int z in chargesToMatch)
                            {
                                List<IndexedMassSpectralPeak> chargeXic = new List<IndexedMassSpectralPeak>();

                                for (int j = start.ZeroBasedMs1ScanIndex; j <= end.ZeroBasedMs1ScanIndex; j++)
                                {
                                    IndexedMassSpectralPeak peak = _peakIndexingEngine.GetIndexedPeak(donorIdentification.PeakfindingMass, j, mbrTol, z);

                                    if (peak != null)
                                    {
                                        chargeXic.Add(peak);
                                    }
                                }

                                if (!chargeXic.Any())
                                {
                                    continue;
                                }

                                List<IsotopicEnvelope> chargeEnvelopes = GetIsotopicEnvelopes(chargeXic, donorIdentification, z);

                                // treat each isotopic envelope in the valid region as a potential seed for a chromatographic peak.
                                // remove the clustered isotopic envelopes from the list of seeds after each iteration
                                while (chargeEnvelopes.Any())
                                {
                                    var acceptorPeak = new ChromatographicPeak(donorIdentification, true, idAcceptorFile);
                                    IsotopicEnvelope seedEnv = chargeEnvelopes.First();

                                    var xic = Peakfind(seedEnv.IndexedPeak.RetentionTime, donorIdentification.PeakfindingMass, z, idAcceptorFile, mbrTol);
                                    List<IsotopicEnvelope> bestChargeEnvelopes = GetIsotopicEnvelopes(xic, donorIdentification, z);
                                    acceptorPeak.IsotopicEnvelopes.AddRange(bestChargeEnvelopes);
                                    acceptorPeak.CalculateIntensityForThisFeature(Integrate);
                                    acceptorPeak.SetRtWindow(acceptorFileRtHypothesis, rtStdDev, rtInterquartileRange);

                                    CutPeak(acceptorPeak, seedEnv.IndexedPeak.RetentionTime);

                                    var claimedPeaks = new HashSet<IndexedMassSpectralPeak>(acceptorPeak.IsotopicEnvelopes.Select(p => p.IndexedPeak));
                                    claimedPeaks.Add(seedEnv.IndexedPeak); // prevents infinite loops

                                    chargeEnvelopes.RemoveAll(p => claimedPeaks.Contains(p.IndexedPeak));

                                    // peak has already been identified by MSMS - skip it
                                    if (apexToAcceptorFilePeak.ContainsKey(seedEnv.IndexedPeak))
                                    {
                                        continue;
                                    }

                                    // score the peak hypothesis
                                    double rtScore = rtScoringDistribution.Density(acceptorPeak.Apex.IndexedPeak.RetentionTime);
                                    double ppmScore = ppmDistribution.Density(acceptorPeak.MassError);
                                    double intensityScore = 0;

                                    double logIntensity = Math.Log(acceptorPeak.Intensity, 2);

                                    if (foldChangeDistribution != null)
                                    {
                                        intensityScore = foldChangeDistribution.Density(logIntensity);
                                    }
                                    else
                                    {
                                        if (logIntensity < medianAcceptorLogIntensity)
                                        {
                                            intensityScore = intensityDistribution.Density(logIntensity);
                                        }
                                        else
                                        {
                                            intensityScore = intensityDistribution.Density(intensityDistribution.Mode);
                                        }
                                    }

                                    rtScore = Math.Log(rtScore + Math.Sqrt(Math.Pow(rtScore, 2) + 1));
                                    ppmScore = Math.Log(ppmScore + Math.Sqrt(Math.Pow(ppmScore, 2) + 1));
                                    intensityScore = Math.Log(intensityScore + Math.Sqrt(Math.Pow(intensityScore, 2) + 1));

                                    // larger scores are better
                                    acceptorPeak.MbrScore = (rtScore + ppmScore + intensityScore) * (1 - donorIdentification.PosteriorErrorProbability);

                                    // save the peak hypothesis
                                    // if this peak hypothesis already exists, sum the scores since we've mapped >1 of the same ID onto this peak
                                    if (matchBetweenRunsIdentifiedPeaksThreadSpecific.TryGetValue(donorIdentification.ModifiedSequence, out var mbrPeaks))
                                    {
                                        if (mbrPeaks.TryGetValue(acceptorPeak.Apex, out List<ChromatographicPeak> existing))
                                        {
                                            var samePeakSameSequence = existing
                                                .FirstOrDefault(p => p.Identifications.First().ModifiedSequence == acceptorPeak.Identifications.First().ModifiedSequence);

                                            if (samePeakSameSequence != null)
                                            {
                                                samePeakSameSequence.MbrScore += acceptorPeak.MbrScore;
                                                samePeakSameSequence.Identifications.Add(donorIdentification);
                                            }
                                            else
                                            {
                                                existing.Add(acceptorPeak);
                                            }
                                        }
                                        else
                                        {
                                            mbrPeaks.Add(acceptorPeak.Apex, new List<ChromatographicPeak> { acceptorPeak });
                                        }
                                    }
                                    else
                                    {
                                        matchBetweenRunsIdentifiedPeaksThreadSpecific.Add(donorIdentification.ModifiedSequence, new Dictionary<IsotopicEnvelope, List<ChromatographicPeak>>());
                                        matchBetweenRunsIdentifiedPeaksThreadSpecific[donorIdentification.ModifiedSequence].Add(acceptorPeak.Apex, new List<ChromatographicPeak> { acceptorPeak });
                                    }
                                }
                            }
                        }

                        // merge results from different threads
                        lock (matchBetweenRunsIdentifiedPeaks)
                        {
                            foreach (var kvp in matchBetweenRunsIdentifiedPeaksThreadSpecific)
                            {
                                if (matchBetweenRunsIdentifiedPeaks.TryGetValue(kvp.Key, out var list))
                                {
                                    foreach (var peak in kvp.Value)
                                    {
                                        if (list.TryGetValue(peak.Key, out List<ChromatographicPeak> existing))
                                        {
                                            foreach (var acceptorPeak in peak.Value)
                                            {
                                                var samePeakSameSequence = existing
                                                    .FirstOrDefault(p => p.Identifications.First().ModifiedSequence == acceptorPeak.Identifications.First().ModifiedSequence);

                                                if (samePeakSameSequence != null)
                                                {
                                                    samePeakSameSequence.MbrScore += acceptorPeak.MbrScore;
                                                    samePeakSameSequence.Identifications.Add(acceptorPeak.Identifications.First());
                                                }
                                                else
                                                {
                                                    existing.Add(acceptorPeak);
                                                }
                                            }
                                        }
                                        else
                                        {
                                            list.Add(peak.Key, peak.Value);
                                        }
                                    }
                                }
                                else
                                {
                                    matchBetweenRunsIdentifiedPeaks.Add(kvp.Key, kvp.Value);
                                }
                            }
                        }
                    });
            }

            // take the best result (highest scoring) for each peptide after we've matched from all the donor files
            foreach (var mbrIdentifiedPeptide in matchBetweenRunsIdentifiedPeaks.Where(p => !acceptorFileIdentifiedSequences.Contains(p.Key)))
            {
                string peptideModifiedSequence = mbrIdentifiedPeptide.Key;
                if (!mbrIdentifiedPeptide.Value.Any())
                {
                    continue;
                }

                List<ChromatographicPeak> peakHypotheses = mbrIdentifiedPeptide.Value.SelectMany(p => p.Value).OrderByDescending(p => p.MbrScore).ToList();

                ChromatographicPeak best = peakHypotheses.First();

                peakHypotheses.Remove(best);

                if (peakHypotheses.Count > 0)
                {
                    double start = best.IsotopicEnvelopes.Min(p => p.IndexedPeak.RetentionTime);
                    double end = best.IsotopicEnvelopes.Max(p => p.IndexedPeak.RetentionTime);

                    List<ChromatographicPeak> peaksToRemoveFromHypotheses = new List<ChromatographicPeak>();
                    foreach (ChromatographicPeak peak in peakHypotheses.Where(p => p.Apex.ChargeState != best.Apex.ChargeState))
                    {
                        if (peak.Apex.IndexedPeak.RetentionTime > start && peak.Apex.IndexedPeak.RetentionTime < end)
                        {
                            best.MergeFeatureWith(peak, Integrate);

                            peaksToRemoveFromHypotheses.Add(peak);
                        }
                    }
                }

                //here is where we get the cosine similarity.

                _results.Peaks[idAcceptorFile].Add(best);
            }

            RunErrorChecking(idAcceptorFile);
        }

        /// <summary>
        /// Used by the match-between-runs algorithm to determine systematic retention time drifts between
        /// chromatographic runs.
        /// </summary>
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

        /// <summary>
        /// Checks for and resolves situations where one IndexedMassSpectralPeak is defined as the apex 
        /// for multiple ChromatographicPeaks. In these situations, the two peaks are merged and the merged
        /// peak is stored in the FlashLFQ results.
        /// </summary>
        /// <param name="spectraFile"></param>
        private void RunErrorChecking(SpectraFileInfo spectraFile)
        {
            if (!Silent)
            {
                Console.WriteLine("Checking errors");
            }

            _results.Peaks[spectraFile].RemoveAll(p => p == null || p.IsMbrPeak && !p.IsotopicEnvelopes.Any());

            // merge duplicate peaks and handle MBR/MSMS peakfinding conflicts
            var errorCheckedPeaksGroupedByApex = new Dictionary<IndexedMassSpectralPeak, ChromatographicPeak>();
            var errorCheckedPeaks = new List<ChromatographicPeak>();
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

                    errorCheckedPeaks.Add(tryPeak);
                    continue;
                }

                IndexedMassSpectralPeak apexImsPeak = tryPeak.Apex.IndexedPeak;
                if (errorCheckedPeaksGroupedByApex.TryGetValue(apexImsPeak, out ChromatographicPeak storedPeak))
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
                        else if (tryPeak.MbrScore > storedPeak.MbrScore)
                        {
                            errorCheckedPeaksGroupedByApex[tryPeak.Apex.IndexedPeak] = tryPeak;
                        }
                    }
                }
                else
                {
                    errorCheckedPeaksGroupedByApex.Add(apexImsPeak, tryPeak);
                }
            }

            errorCheckedPeaks.AddRange(errorCheckedPeaksGroupedByApex.Values.Where(p => p != null));

            _results.Peaks[spectraFile] = errorCheckedPeaks;
        }

        /// <summary>
        /// Takes in a list of imsPeaks and finds all the isotopic peaks in each scan. If the experimental isotopic distribution
        /// matches the theoretical distribution, an IsotopicEnvelope object is created from the summed intensities of each isotopic peak.
        /// </summary>
        /// <param name="xic"> List of imsPeaks, where the mass of each peak is the peak finding mass (most abundant isotope) </param>
        /// <returns> A list of IsotopicEnvelopes, where each envelope contains the sum of the isotopic peak intensities from one scan </returns>
        public List<IsotopicEnvelope> GetIsotopicEnvelopes(List<IndexedMassSpectralPeak> xic, Identification identification, int chargeState)
        {
            var isotopicEnvelopes = new List<IsotopicEnvelope>();
            var isotopeMassShifts = _modifiedSequenceToIsotopicDistribution[identification.ModifiedSequence];

            if (isotopeMassShifts.Count < NumIsotopesRequired)
            {
                return isotopicEnvelopes;
            }

            PpmTolerance isotopeTolerance = new PpmTolerance(IsotopePpmTolerance);

            double[] experimentalIsotopeIntensities = new double[isotopeMassShifts.Count];
            double[] theoreticalIsotopeMassShifts = isotopeMassShifts.Select(p => p.Item1).ToArray();
            double[] theoreticalIsotopeAbundances = isotopeMassShifts.Select(p => p.Item2).ToArray();
            int peakfindingMassIndex = (int)Math.Round(identification.PeakfindingMass - identification.MonoisotopicMass, 0);

            // For each peak in the XIC, we consider the possibility that there was an off-by-one or missed monoisotopic mass
            // error in peak assignment / deconvolution. The -1 key in this dictionary corresponds to a negative off-by-one error, the 
            // +1 key corresponds to a positive off-by-one error, and the 0 key corresponds to accurate assignment/deconvolution.
            var massShiftToIsotopePeaks = new Dictionary<int, List<(double expIntensity, double theorIntensity, double theorMass)>>
            {
                { -1, new List<(double, double, double)>() },
                { 0, new List<(double, double, double)>() },
                { 1, new List<(double, double, double)>() },
            };

            List<int> directions = new List<int> { -1, 1 };

            // For each peak (most abundant mass peak), we check for the possibility that the peak was mis-assigned,
            // i.e. that the peak belongs to a species with a different mass than the identification mass
            foreach (IndexedMassSpectralPeak peak in xic)
            {
                Array.Clear(experimentalIsotopeIntensities, 0, experimentalIsotopeIntensities.Length);
                foreach (var kvp in massShiftToIsotopePeaks)
                {
                    kvp.Value.Clear();
                }

                // isotope masses are calculated relative to the observed peak
                double observedMass = peak.Mz.ToMass(chargeState);
                double observedMassError = observedMass - identification.PeakfindingMass;

                foreach (var shift in massShiftToIsotopePeaks)
                {
                    // look for each isotope peak in the data
                    // This is done by starting with the first isotope with mass less than the
                    // peak finding (most abundant) mass, then working backwards to find every isotope
                    // with mass < most abundant mass. Once an expected isotopic peak can not be found,
                    // the loop breaks and we begin working our way forward, starting with the peak finding 
                    // mass and locating every peak with mass > peak finding mass. Once an expected isotopic
                    // peak can not be found, it is assumed that we have located every isotope present and the loop breaks.
                    foreach (int direction in directions)
                    {
                        int start = direction == -1
                            ? peakfindingMassIndex - 1
                            : peakfindingMassIndex;

                        for (int i = start; i < theoreticalIsotopeAbundances.Length && i >= 0; i += direction)
                        {
                            double isotopeMass = identification.MonoisotopicMass + observedMassError +
                                                 theoreticalIsotopeMassShifts[i] + shift.Key * Constants.C13MinusC12;
                            double theoreticalIsotopeIntensity = theoreticalIsotopeAbundances[i] * peak.Intensity;

                            IndexedMassSpectralPeak isotopePeak = _peakIndexingEngine.GetIndexedPeak(isotopeMass,
                                peak.ZeroBasedMs1ScanIndex, isotopeTolerance, chargeState);

                            if (isotopePeak == null
                                || isotopePeak.Intensity < theoreticalIsotopeIntensity / 4.0
                                || isotopePeak.Intensity > theoreticalIsotopeIntensity * 4.0)
                            {
                                break;
                            }

                            shift.Value.Add((isotopePeak.Intensity, theoreticalIsotopeIntensity, isotopeMass));
                            if (shift.Key == 0)
                            {
                                experimentalIsotopeIntensities[i] = isotopePeak.Intensity;
                            }
                        }
                    }
                }

                // check number of isotope peaks observed
                if (massShiftToIsotopePeaks[0].Count < NumIsotopesRequired)
                {
                    continue;
                }

                // Check that the experimental envelope matches the theoretical
                if (CheckIsotopicEnvelopeCorrelation(massShiftToIsotopePeaks, peak, chargeState, isotopeTolerance))
                {
                    // impute unobserved isotope peak intensities
                    // TODO: Figure out why value imputation is performed. Build a toggle?
                    for (int i = 0; i < experimentalIsotopeIntensities.Length; i++)
                    {
                        if (experimentalIsotopeIntensities[i] == 0)
                        {
                            experimentalIsotopeIntensities[i] = theoreticalIsotopeAbundances[i] * experimentalIsotopeIntensities[peakfindingMassIndex];
                        }
                    }

                    isotopicEnvelopes.Add(new IsotopicEnvelope(
                        peak, chargeState, experimentalIsotopeIntensities.Sum()));
                }
            }

            return isotopicEnvelopes;
        }

        /// <summary>
        /// This function checks the correlation between experimental and actual abundances of isotopes
        /// for a given species. It returns true if the experimental data is best described by the
        /// theoretical isotope abundances, and false if there is low concordance between the theoretical
        /// and actual abundances, or if the observed data is better described by an envelope with a
        /// monoisotopic mass +/- one DA away
        /// </summary>
        /// <param name="massShiftToIsotopePeaks">Dictionary containing the experimental and theoretical abundances and expected
        /// mass for a given set of isotopic peaks, shifted -1, 0, and 1 Da (shifts = keys)</param>
        /// <returns>True if experimental data is a good match to the expected isotopic distribution </returns>
        public bool CheckIsotopicEnvelopeCorrelation(
            Dictionary<int, List<(double expIntensity, double theorIntensity, double theorMass)>> massShiftToIsotopePeaks,
            IndexedMassSpectralPeak peak,
            int chargeState,
            Tolerance isotopeTolerance)
        {
            double corr = Correlation.Pearson(
        massShiftToIsotopePeaks[0].Select(p => p.expIntensity),
        massShiftToIsotopePeaks[0].Select(p => p.theorIntensity));

            // check correlation of experimental isotope intensities to the theoretical abundances
            // check for unexpected peaks 
            foreach (var shift in massShiftToIsotopePeaks)
            {
                if (!shift.Value.Any())
                {
                    continue;
                }

                double unexpectedMass = shift.Value.Min(p => p.theorMass) - Constants.C13MinusC12;
                IndexedMassSpectralPeak unexpectedPeak = _peakIndexingEngine.GetIndexedPeak(unexpectedMass,
                            peak.ZeroBasedMs1ScanIndex, isotopeTolerance, chargeState);

                if (unexpectedPeak == null)
                {
                    shift.Value.Add((0, 0, unexpectedMass));
                }
                else
                {
                    shift.Value.Add((unexpectedPeak.Intensity, 0, unexpectedMass));
                }
            }

            double corrWithPadding = Correlation.Pearson(
                massShiftToIsotopePeaks[0].Select(p => p.expIntensity),
                massShiftToIsotopePeaks[0].Select(p => p.theorIntensity));
            double corrShiftedLeft = Correlation.Pearson(
                massShiftToIsotopePeaks[-1].Select(p => p.expIntensity),
                massShiftToIsotopePeaks[-1].Select(p => p.theorIntensity));
            double corrShiftedRight = Correlation.Pearson(
                massShiftToIsotopePeaks[1].Select(p => p.expIntensity),
                massShiftToIsotopePeaks[1].Select(p => p.theorIntensity));

            if (double.IsNaN(corrShiftedLeft))
            {
                corrShiftedLeft = -1;
            }
            if (double.IsNaN(corrShiftedRight))
            {
                corrShiftedRight = -1;
            }

            // If these conditions are true, the isotopic envelope matches the expected envelope better than 
            // either alternative (i.e., +/- missed mono-isotopic)
            return corr > 0.7 && corrShiftedLeft - corrWithPadding < 0.1 && corrShiftedRight - corrWithPadding < 0.1;
        }

        /// <summary>
        /// Finds peaks with a given mz (mass/charge + H) that occur on either side of a given
        /// retention time. Peak searching iterates backwards through MS1 scans until the peak 
        /// is no longer observed (i.e., is absent in more scans than allowed, as defined by the
        /// MissedScansAllowed property. Missed scans don't have to be sequential. The same procedure
        /// is then repeated in the forward direction.
        /// </summary>
        /// <param name="idRetentionTime"> Time where peak searching behaviour begins </param>
        /// <param name="mass"> Peakfinding mass </param>
        /// <returns></returns>
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

            // Sorts the list by RT in-place. (OrderBy does not sort in place, but creates a copy instead)
            xic.Sort((x, y) => x.RetentionTime.CompareTo(y.RetentionTime));

            return xic;
        }

        /// <summary>
        /// Recursively cuts ChromatographicPeaks, removing all IsotopicEnvelopes
        /// that occur before or after potential "valleys" surrounding the identification's
        /// MS2 retention time. Then, the peak intensity is recalculated
        /// </summary>
        /// <param name="peak"> Peak to be cut, where envelopes are sorted by MS1 scan number </param>
        /// <param name="identificationTime"> MS2 scan retention time </param>
        private void CutPeak(ChromatographicPeak peak, double identificationTime)
        {
            // find out if we need to split this peak by using the discrimination factor
            // this method assumes that the isotope envelopes in a chromatographic peak are already sorted by MS1 scan number
            bool cutThisPeak = false;

            if (peak.IsotopicEnvelopes.Count < 5)
            {
                return;
            }

            List<IsotopicEnvelope> timePointsForApexZ = peak.IsotopicEnvelopes
                .Where(p => p.ChargeState == peak.Apex.ChargeState).ToList();
            HashSet<int> scanNumbers = new HashSet<int>(timePointsForApexZ.Select(p => p.IndexedPeak.ZeroBasedMs1ScanIndex));
            int apexIndex = timePointsForApexZ.IndexOf(peak.Apex);
            IsotopicEnvelope valleyEnvelope = null;

            // -1 checks the left side, +1 checks the right side
            int[] directions = { 1, -1 };

            foreach (int direction in directions)
            {
                valleyEnvelope = null;
                int indexOfValley = 0;

                for (int i = apexIndex + direction; i < timePointsForApexZ.Count && i >= 0; i += direction)
                {
                    IsotopicEnvelope timepoint = timePointsForApexZ[i];

                    if (valleyEnvelope == null || timepoint.Intensity < valleyEnvelope.Intensity)
                    {
                        valleyEnvelope = timepoint;
                        indexOfValley = timePointsForApexZ.IndexOf(valleyEnvelope);
                    }

                    double discriminationFactor =
                        (timepoint.Intensity - valleyEnvelope.Intensity) / timepoint.Intensity;

                    if (discriminationFactor > DiscriminationFactorToCutPeak &&
                        (indexOfValley + direction < timePointsForApexZ.Count && indexOfValley + direction >= 0))
                    {
                        IsotopicEnvelope secondValleyTimepoint = timePointsForApexZ[indexOfValley + direction];

                        discriminationFactor =
                            (timepoint.Intensity - secondValleyTimepoint.Intensity) / timepoint.Intensity;

                        if (discriminationFactor > DiscriminationFactorToCutPeak)
                        {
                            cutThisPeak = true;
                            break;
                        }

                        int nextMs1ScanNum = -1;
                        for (int j = valleyEnvelope.IndexedPeak.ZeroBasedMs1ScanIndex - 1;
                            j < _ms1Scans[peak.SpectraFileInfo].Length && j >= 0;
                            j += direction)
                        {
                            if (_ms1Scans[peak.SpectraFileInfo][j].OneBasedScanNumber >= 0 &&
                                _ms1Scans[peak.SpectraFileInfo][j].OneBasedScanNumber !=
                                valleyEnvelope.IndexedPeak.ZeroBasedMs1ScanIndex)
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
                if (identificationTime > valleyEnvelope.IndexedPeak.RetentionTime)
                {
                    // MS2 identification is to the right of the valley; remove all peaks left of the valley
                    peak.IsotopicEnvelopes.RemoveAll(p => 
                        p.IndexedPeak.RetentionTime <= valleyEnvelope.IndexedPeak.RetentionTime);
                }
                else
                {
                    // MS2 identification is to the left of the valley; remove all peaks right of the valley
                    peak.IsotopicEnvelopes.RemoveAll(p => 
                        p.IndexedPeak.RetentionTime >= valleyEnvelope.IndexedPeak.RetentionTime);
                }

                // recalculate intensity for the peak
                peak.CalculateIntensityForThisFeature(Integrate);
                peak.SplitRT = valleyEnvelope.IndexedPeak.RetentionTime;

                // recursively cut
                CutPeak(peak, identificationTime);
            }
        }
    }
}