using Chemistry;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading.Tasks;
using System.Runtime.CompilerServices;
using System.IO;
using Easy.Common.Extensions;
using FlashLFQ.PEP;
using FlashLFQ.IsoTracker;
using System.Threading;
using FlashLFQ.Interfaces;
using MassSpectrometry;

[assembly: InternalsVisibleTo("Test")]

namespace FlashLFQ
{
    public class FlashLfqEngine
    {
        public FlashLfqParameters FlashParams { get; init; }
        public List<SpectraFileInfo> SpectraFileInfoList { get; init; }

        #region Constants
        public const double DiscriminationFactorToCutPeak = 0.6; // Determines behavior related to splitting bimodal peaks
        public const double PeakfindingPpmTolerance = 20.0; // Ppm tolerance used when finding XICs
        public const int MissedScansAllowed = 1; // Number of consecutive missed scans allowed when finding XICs

        // MBR Settings
        public const int NumberOfAnchorPeptidesForMbr = 3; // the number of anchor peptides used for local alignment when predicting retention times of MBR acceptor peptides
        public readonly double MbrAlignmentWindow = 2.5; // The window in which we look for anchor peptides
        public readonly double PepTrainingFraction = 0.25; // Use the top 25% of the MBR peaks as positive training examples for PEP
        #endregion

        #region Internal/Private properties and fields
        /// <summary>
        /// These peptides will be reported in the QuantifiedPeptides output and used for protein quant.
        /// Other peptides may appear in the QuantifiedPeaks output, but this list is used to enable
        /// peptide-level FDR filtering
        /// </summary>
        internal HashSet<string> PeptideModifiedSequencesToQuantify { get; init; }
        internal Dictionary<SpectraFileInfo, IFlashLfqIndexingEngine> IndexingEngineDictionary { get; private set; }
        internal Dictionary<SpectraFileInfo, List<ChromatographicPeak>> DonorFileToPeakDict { get; private set; }
        /// <summary>
        /// a flag used to indicate if the isobaric case is running, used to control the indexEngine
        /// </summary>
        internal bool IsoTrackerIsRunning { get; private set; }
        /// <summary>
        /// The list of IsoPeptide Group with the same base sequence and monoisotopic mass. Only for IsoTracker analysis
        /// </summary>
        internal List<IsobaricPeptideGroup> PeptideGroupsForIsoTracker { get; private set; } 
        /// <summary>
        /// The dictionary of isobaric peaks for each modified sequence
        /// </summary>
        internal ConcurrentDictionary<string, Dictionary<PeakRegion, List<ChromatographicPeak>>> IsobaricPeptideDict { get; private set; } // 
        /// <summary>
        /// Dictionary linking a modified sequence to a List of tuples containing
        /// the mass shifts (isotope mass - monoisotopic mass) and normalized abundances for the
        /// isotopes for a given peptide
        /// </summary>
        internal Dictionary<string, List<(double massShift, double normalizedAbundance)>> ModifiedSequenceToIsotopicDistribution { get; set; }
        private List<int> _chargeStates;
        private FlashLfqResults _results;
        private readonly List<Identification> _allIdentifications;
        private readonly Stopwatch _globalStopwatch;
        #endregion

        /// <summary>
        /// Create an instance of FlashLFQ that will quantify peptides based on their precursor intensity in MS1 spectra
        /// </summary>
        /// <param name="flashLfqParameters"> Parameters used by FlashLFQ </param>
        /// <param name="allIdentifications"></param>
        /// <param name="peptideSequencesToQuantify"></param>
        public FlashLfqEngine(FlashLfqParameters flashLfqParameters, List<Identification> allIdentifications, List<string> peptideSequencesToQuantify = null)
        {
            FlashParams = flashLfqParameters;

            _globalStopwatch = new Stopwatch();
            _chargeStates = new List<int>();
            IndexingEngineDictionary = new();

            SpectraFileInfoList = allIdentifications.Select(p => p.FileInfo).Distinct()
                .OrderBy(p => p.Condition)
                .ThenBy(p => p.BiologicalReplicate)
                .ThenBy(p => p.Fraction)
                .ThenBy(p => p.TechnicalReplicate).ToList();

            _allIdentifications = allIdentifications;
            PeptideModifiedSequencesToQuantify = peptideSequencesToQuantify.IsNotNullOrEmpty()
                ? new HashSet<string>(peptideSequencesToQuantify)
                : allIdentifications.Select(id => id.ModifiedSequence).ToHashSet();

            if (FlashParams.MaxThreads == -1 || FlashParams.MaxThreads >= Environment.ProcessorCount)
            {
                FlashParams.MaxThreads = Environment.ProcessorCount - 1;
            }

            if (FlashParams.MaxThreads <= 0)
            {
                FlashParams.MaxThreads = 1;
            }
        }

        /// <summary>
        /// Create an instance of FlashLFQ that will quantify peptides based on their precursor intensity in MS1 spectra
        /// </summary>
        /// <param name="allIdentifications">A list of identifications corresponding to MS2 peptide detections. One ID per peptide per file</param>
        /// <param name="peptideSequencesToQuantify">Optional. A list of strings corresponding to the modified sequences of peptides that should be quantified/used for
        /// protein level quant. Reccommended use is to pass in the full sequence of every peptide at 1% peptide-level FDR</param>
        ///  /// <param name="integrate">Optional. Bool indicating whether peaks should be integrated before quantification. It is HIGHLY recommended this is set to FALSE</param>
        public FlashLfqEngine(
            List<Identification> allIdentifications,
            List<string> peptideSequencesToQuantify = null,
            bool normalize = false,
            double ppmTolerance = 10.0,
            double isotopeTolerancePpm = 5.0,
            bool integrate = false,
            int numIsotopesRequired = 2,
            bool idSpecificChargeState = false,
            bool quantifyAmbiguousPeptides = false,
            bool silent = false,
            int maxThreads = -1,

            // IsoTracker settings
            bool isoTracker = false,
            List<char> motifsList = null, // The target motifs for the IsoTracker
            bool requireMultipleIdsInOneFiles = true,

            // MBR settings
            bool matchBetweenRuns = false,
            double matchBetweenRunsPpmTolerance = 10.0,
            double maxMbrWindow = 1.0,
            bool requireMsmsIdInCondition = false,
            double matchBetweenRunsFdrThreshold = 0.05,
            double donorQValueThreshold = 0.01,
            DonorCriterion donorCriterion = DonorCriterion.Score,

            // settings for the Bayesian protein quantification engine
            bool bayesianProteinQuant = false,
            string proteinQuantBaseCondition = null,
            double proteinQuantFoldChangeCutoff = 0.1,
            int mcmcSteps = 3000,
            int mcmcBurninSteps = 1000,
            bool useSharedPeptidesForProteinQuant = false,
            bool pairedSamples = false,
            int? randomSeed = null) :
            this(
                new FlashLfqParameters()
                {
                    PpmTolerance = ppmTolerance,
                    IsotopePpmTolerance = isotopeTolerancePpm,
                    Integrate = integrate,
                    NumIsotopesRequired = numIsotopesRequired,
                    IdSpecificChargeState = idSpecificChargeState,
                    QuantifyAmbiguousPeptides = quantifyAmbiguousPeptides,
                    Silent = silent,
                    MaxThreads = maxThreads,
                    Normalize = normalize,
                    IsoTracker = isoTracker,
                    IsoTrackerIdFilter = new IsoTrackerIdFilter(motifsList),
                    RequireMultipleIdsInOneFiles = requireMultipleIdsInOneFiles,
                    MatchBetweenRuns = matchBetweenRuns,
                    MaxMbrRtWindow = maxMbrWindow,
                    MbrPpmTolerance = matchBetweenRunsPpmTolerance,
                    MbrQValueThreshold = matchBetweenRunsFdrThreshold,
                    DonorQValueThreshold = donorQValueThreshold,
                    DonorCriterion = donorCriterion,
                    RequireMsmsIdInCondition = requireMsmsIdInCondition,
                    BayesianProteinQuant = bayesianProteinQuant,
                    ProteinQuantBaseCondition = proteinQuantBaseCondition,
                    ProteinQuantFoldChangeCutoff = proteinQuantFoldChangeCutoff,
                    McmcSteps = mcmcSteps,
                    McmcBurninSteps = mcmcBurninSteps,
                    UseSharedPeptidesForProteinQuant = useSharedPeptidesForProteinQuant,
                    PairedSamples = pairedSamples,
                    RandomSeed = randomSeed
                },
                allIdentifications,
                peptideSequencesToQuantify
            )
        { }

        public FlashLfqResults Run()
        {
            _globalStopwatch.Start();
            _results = new FlashLfqResults(SpectraFileInfoList, _allIdentifications, FlashParams.MbrQValueThreshold, PeptideModifiedSequencesToQuantify, FlashParams.IsoTracker);
            // build m/z index keys
            CalculateTheoreticalIsotopeDistributions();
            List<float> targetMzs = new List<float>();
            if (FlashParams.IsoTracker) // If IsoTracker is on, we need to set the target mass for pruning
            {
                PeptideGroupsForIsoTracker = GroupPeptideForIsoTracker(); // Generate the isoPeptideGroups for IsoTracker
                targetMzs = GetTargetMz();
            }

            // quantify each file
            foreach (var spectraFile in SpectraFileInfoList)
            {
                if (!FlashParams.Silent) Console.WriteLine("Reading spectra file");
                var indexingEngine = PeakIndexingEngine.InitializeIndexingEngine(spectraFile);
                if(indexingEngine == null)
                {
                    // something went wrong finding/opening/indexing the file...
                    if( !FlashParams.Silent ) Console.WriteLine("FlashLFQ Error: The file " + spectraFile.FilenameWithoutExtension + " contained no MS1 peaks!");
                    continue;
                }
                IndexingEngineDictionary[spectraFile] = indexingEngine;

                // quantify peaks using this file's IDs first
                QuantifyMs2IdentifiedPeptides(spectraFile);

                // Retain, Serialize and Clear, or Clear the indexing engine
                // If IsoTracker is on, we don't need to serialize the index for each file. The indexed peaks are retained.
                if (FlashParams.IsoTracker) 
                {
                    if (FlashParams.MatchBetweenRuns) // If MBR is on then we need to serialize the index.
                    {
                        IndexingEngineDictionary[spectraFile].SerializeIndex(); 
                    }
                    IndexingEngineDictionary[spectraFile].PruneIndex(targetMzs); // Remove some unused data from the index to save memory for IsoTracker
                }
                // If IsoTracker is off, we  need to clear the indexing engine array to save memory.
                else
                {
                    if (FlashParams.MatchBetweenRuns) // If MBR is turned on then we need to serialize the index before clearing the array.
                    {
                        IndexingEngineDictionary[spectraFile].SerializeIndex(); 
                    }
                    IndexingEngineDictionary[spectraFile].ClearIndex();
                }
                // error checking function
                // handles features with multiple identifying scans and scans that are associated with more than one feature
                RunErrorChecking(spectraFile);

                if (!FlashParams.Silent)
                {
                    Console.WriteLine("Finished " + spectraFile.FilenameWithoutExtension);
                }
            }

            //IsoTracker
            if (FlashParams.IsoTracker)
            {
                IsoTrackerIsRunning = true; // Turn on the flag, then we will use the separate indexEngine for each files
                IsobaricPeptideDict = new ConcurrentDictionary<string, Dictionary<PeakRegion, List<ChromatographicPeak>>>();
                QuantifyIsobaricPeaks();
                _results.IsobaricPeptideDict = IsobaricPeptideDict;
                AddIsoPeaks();
                foreach (var file in _results.SpectraFiles)
                {
                    RunErrorChecking(file);
                }
            }
            IndexingEngineDictionary.ForEach(p=>p.Value.ClearIndex()); // Clear the indexEngine for each file after the quantification is done
            IsoTrackerIsRunning = false;

            // do MBR
            if (FlashParams.MatchBetweenRuns)
            {
                Console.WriteLine("Find the best donors for match-between-runs");
                FindPeptideDonorFiles();
                foreach (var spectraFile in SpectraFileInfoList)
                {            
                    if ( !IndexingEngineDictionary.ContainsKey(spectraFile) ) continue;
                    if (!FlashParams.Silent)
                    {
                        Console.WriteLine("Doing match-between-runs for " + spectraFile.FilenameWithoutExtension);
                    }

                    //Deserialize the relevant index prior to MBR
                    IndexingEngineDictionary[spectraFile].DeserializeIndex();
                    QuantifyMatchBetweenRunsPeaks(spectraFile);
                    IndexingEngineDictionary[spectraFile].ClearIndex();

                    if (!FlashParams.Silent)
                    {
                        Console.WriteLine("Finished MBR for " + spectraFile.FilenameWithoutExtension);
                    }
                }

                Console.WriteLine("Computing PEP for MBR Transfers");
                bool pepSuccesful = RunPEPAnalysis();

                foreach (var spectraFile in SpectraFileInfoList)
                {
                    CalculateFdrForMbrPeaks(spectraFile, pepSuccesful);
                }
            }

            // normalize
            if (FlashParams.Normalize)
            {
                new IntensityNormalizationEngine(_results, FlashParams.Integrate, FlashParams.Silent, FlashParams.MaxThreads).NormalizeResults();
            }

            // calculate peptide intensities
            _results.CalculatePeptideResults(FlashParams.QuantifyAmbiguousPeptides);

            // do top3 protein quantification
            _results.CalculateProteinResultsMedianPolish(FlashParams.UseSharedPeptidesForProteinQuant);

            // do Bayesian protein fold-change analysis
            if (FlashParams.BayesianProteinQuant)
            {
                if (SpectraFileInfoList.Count == 1 || SpectraFileInfoList.Select(p => p.Condition).Distinct().Count() == 1)
                {
                    if (!FlashParams.Silent)
                    {
                        Console.WriteLine("Can't do Bayesian protein quant with only one spectra file or condition. FlashLFQ will still do a top3 protein quant");
                    }
                }
                else
                {
                    if (!FlashParams.Silent)
                    {
                        Console.WriteLine("Running Bayesian protein quantification analysis");
                    }

                    new ProteinQuantificationEngine(_results, FlashParams.MaxThreads, FlashParams.ProteinQuantBaseCondition, FlashParams.UseSharedPeptidesForProteinQuant,
                        FlashParams.ProteinQuantFoldChangeCutoff, FlashParams.RandomSeed, FlashParams.McmcBurninSteps, FlashParams.McmcSteps, FlashParams.PairedSamples).Run();
                }
            }

            // done
            if (!FlashParams.Silent)
            {
                Console.WriteLine("Done quantifying");
            }

            if (!FlashParams.Silent)
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
        /// If the sequence is modified and the modification has an unknown chemical formula,
        /// averagine is used for the modified part
        /// </summary>
        internal void CalculateTheoreticalIsotopeDistributions()
        {
            ModifiedSequenceToIsotopicDistribution = new Dictionary<string, List<(double, double)>>();

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
                if (ModifiedSequenceToIsotopicDistribution.ContainsKey(id.ModifiedSequence))
                {
                    continue;
                }

                ChemicalFormula formula = id.OptionalChemicalFormula;

                var isotopicMassesAndNormalizedAbundances = new List<(double massShift, double abundancence)>();

                if(formula is null)
                {
                    formula = new ChemicalFormula();
                    if (id.BaseSequence.AllSequenceResiduesAreValid())
                    {
                        // there are sometimes non-parsable sequences in the base sequence input
                        formula = new Proteomics.AminoAcidPolymer.Peptide(id.BaseSequence).GetChemicalFormula();
                        double massDiff = id.MonoisotopicMass;
                        massDiff -= formula.MonoisotopicMass;

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
                        double averagines = id.MonoisotopicMass / averagineMass;

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
                    if (isotopicMassesAndNormalizedAbundances.Count < FlashParams.NumIsotopesRequired || abundances[i] > 0.1)
                    {
                        isotopicMassesAndNormalizedAbundances.Add((masses[i], abundances[i]));
                    }
                }

                ModifiedSequenceToIsotopicDistribution.Add(id.ModifiedSequence, isotopicMassesAndNormalizedAbundances);
            }

            var minChargeState = _allIdentifications.Min(p => p.PrecursorChargeState);
            var maxChargeState = _allIdentifications.Max(p => p.PrecursorChargeState);
            _chargeStates = Enumerable.Range(minChargeState, (maxChargeState - minChargeState) + 1).ToList();

            var peptideModifiedSequences = _allIdentifications.GroupBy(p => p.ModifiedSequence);
            foreach (var identifications in peptideModifiedSequences)
            {
                // isotope where normalized abundance is 1
                double mostAbundantIsotopeShift = ModifiedSequenceToIsotopicDistribution[identifications.First().ModifiedSequence]
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
            if (!FlashParams.Silent)
            {
                Console.WriteLine("Quantifying peptides for " + fileInfo.FilenameWithoutExtension);
            }

            var ms2IdsForThisFile = _allIdentifications.Where(p => p.FileInfo.Equals(fileInfo)).ToList();

            if (!ms2IdsForThisFile.Any())
            {
                return;
            }

            PpmTolerance peakfindingTol = new PpmTolerance(PeakfindingPpmTolerance); // Peak finding tolerance is generally higher than ppmTolerance
            PpmTolerance ppmTolerance = new PpmTolerance(FlashParams.PpmTolerance);
            ChromatographicPeak[] chromatographicPeaks = new ChromatographicPeak[ms2IdsForThisFile.Count];

            Parallel.ForEach(Partitioner.Create(0, ms2IdsForThisFile.Count),
                new ParallelOptions { MaxDegreeOfParallelism = FlashParams.MaxThreads },
                (range, loopState) =>
                {
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        var identification = ms2IdsForThisFile[i];
                        ChromatographicPeak msmsFeature = new ChromatographicPeak(identification, fileInfo);
                        chromatographicPeaks[i] = msmsFeature;

                        foreach (var chargeState in _chargeStates)
                        {
                            if (FlashParams.IdSpecificChargeState && chargeState != identification.PrecursorChargeState)
                            {
                                continue;
                            }

                            // get XIC (peakfinding)
                            List<IIndexedPeak> xic = IndexingEngineDictionary[fileInfo].GetXic(
                                    identification.PeakfindingMass.ToMz(chargeState),
                                    identification.Ms2RetentionTimeInMinutes,
                                    peakfindingTol,
                                    MissedScansAllowed)
                                .OrderBy(p => p.RetentionTime)
                                .ToList();
                            
                            // filter by smaller mass tolerance
                            xic.RemoveAll(p => 
                                !ppmTolerance.Within(p.M.ToMass(chargeState), identification.PeakfindingMass));

                            // filter by isotopic distribution
                            List<IsotopicEnvelope> isotopicEnvelopes = GetIsotopicEnvelopes(xic, identification, chargeState, fileInfo);

                            // add isotopic envelopes to the chromatographic peak
                            msmsFeature.IsotopicEnvelopes.AddRange(isotopicEnvelopes);
                        }

                        msmsFeature.CalculateIntensityForThisFeature(FlashParams.Integrate);
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

                        int min = precursorXic.Min(p => p.IndexedPeak.ZeroBasedScanIndex);
                        int max = precursorXic.Max(p => p.IndexedPeak.ZeroBasedScanIndex);
                        msmsFeature.IsotopicEnvelopes.RemoveAll(p => p.IndexedPeak.ZeroBasedScanIndex < min);
                        msmsFeature.IsotopicEnvelopes.RemoveAll(p => p.IndexedPeak.ZeroBasedScanIndex > max);
                        msmsFeature.CalculateIntensityForThisFeature(FlashParams.Integrate);
                    }
                });

            _results.Peaks[fileInfo].AddRange(chromatographicPeaks.ToList());
        }

        #region MatchBetweenRuns
        /// <summary>
        /// Used by the match-between-runs algorithm to determine systematic retention time drifts between
        /// chromatographic runs.
        /// </summary>
        private RetentionTimeCalibDataPoint[] GetRtCalSpline(SpectraFileInfo donor, SpectraFileInfo acceptor, MbrScorer scorer,
            out List<ChromatographicPeak> donorFileBestMsmsPeaksOrderedByMass)
        {
            Dictionary<string, ChromatographicPeak> donorFileBestMsmsPeaks = new();
            Dictionary<string, ChromatographicPeak> acceptorFileBestMsmsPeaks = new();
            List<RetentionTimeCalibDataPoint> rtCalibrationCurve = new();
            List<double> anchorPeptideRtDiffs = new(); // anchor peptides are peptides that were MS2 detected in both the donor and acceptor runs

            Dictionary<string, List<ChromatographicPeak>> donorFileAllMsmsPeaks = _results.Peaks[donor]
                .Where(peak => peak.NumIdentificationsByFullSeq == 1
                    && peak.DetectionType == DetectionType.MSMS
                    && peak.IsotopicEnvelopes.Any()
                    && peak.Identifications.Min(id => id.QValue) < FlashParams.DonorQValueThreshold)
                .GroupBy(peak => peak.Identifications.First().ModifiedSequence)
                .ToDictionary(group => group.Key, group => group.ToList());

            // iterate through each unique donor sequence
            foreach (var sequencePeakListKvp in donorFileAllMsmsPeaks)
            {
                List<ChromatographicPeak> peaksForPeptide = sequencePeakListKvp.Value;
                if (!peaksForPeptide.Any())
                    continue;

                ChromatographicPeak bestPeak = ChooseBestPeak(peaksForPeptide);

                if (bestPeak == null) continue;
                donorFileBestMsmsPeaks.Add(sequencePeakListKvp.Key, bestPeak);
            }

            Dictionary<string, List<ChromatographicPeak>> acceptorFileAllMsmsPeaks = _results.Peaks[acceptor]
                .Where(peak => peak.NumIdentificationsByFullSeq == 1
                    && peak.DetectionType == DetectionType.MSMS
                    && peak.IsotopicEnvelopes.Any()
                    && peak.Identifications.Min(id => id.QValue) < FlashParams.DonorQValueThreshold)
                .GroupBy(peak => peak.Identifications.First().ModifiedSequence)
                .ToDictionary(group => group.Key, group => group.ToList());

            // iterate through each acceptor sequence
            foreach (var sequencePeakListKvp in acceptorFileAllMsmsPeaks)
            {
                List<ChromatographicPeak> peaksForPeptide = sequencePeakListKvp.Value;
                if (!peaksForPeptide.Any())
                    continue;

                ChromatographicPeak bestPeak = ChooseBestPeak(peaksForPeptide);

                if (bestPeak == null) continue;
                acceptorFileBestMsmsPeaks.Add(sequencePeakListKvp.Key, bestPeak);
            }

            // create RT calibration curve
            foreach (var peak in acceptorFileBestMsmsPeaks)
            {
                ChromatographicPeak acceptorFilePeak = peak.Value;

                if (donorFileBestMsmsPeaks.TryGetValue(peak.Key, out ChromatographicPeak donorFilePeak))
                {
                    rtCalibrationCurve.Add(new RetentionTimeCalibDataPoint(donorFilePeak, acceptorFilePeak));
                    if (donorFilePeak.ApexRetentionTime > 0 && acceptorFilePeak.ApexRetentionTime > 0)
                    {
                        anchorPeptideRtDiffs.Add(donorFilePeak.ApexRetentionTime - acceptorFilePeak.ApexRetentionTime);
                    }
                }
            }

            scorer.AddRtPredErrorDistribution(donor, anchorPeptideRtDiffs, NumberOfAnchorPeptidesForMbr);
            donorFileBestMsmsPeaksOrderedByMass = donorFileBestMsmsPeaks.Select(kvp => kvp.Value).OrderBy(p => p.Identifications.First().PeakfindingMass).ToList();

            return rtCalibrationCurve.OrderBy(p => p.DonorFilePeak.Apex.IndexedPeak.RetentionTime).ToArray();
        }

        /// <summary>
        /// For every MSMS identified peptide, selects one file that will be used as the donor
        /// by finding files that contain the most peaks in the local neighborhood,
        /// then writes the restults to the DonorFileToIdsDict.
        /// WARNING! Strong assumption that this is called BEFORE MBR peaks are identified/assigned to the results
        /// </summary>
        private void FindPeptideDonorFiles()
        {
            DonorFileToPeakDict = new Dictionary<SpectraFileInfo, List<ChromatographicPeak>>();

            Dictionary<string, List<ChromatographicPeak>> seqPeakDict = _results.Peaks
                    .SelectMany(kvp => kvp.Value)
                    .Where(peak => peak.NumIdentificationsByFullSeq == 1
                        && peak.IsotopicEnvelopes.Any()
                        && peak.Identifications.Min(id => id.QValue) < FlashParams.DonorQValueThreshold)
                    .GroupBy(peak => peak.Identifications.First().ModifiedSequence)
                    .Where(group => PeptideModifiedSequencesToQuantify.Contains(group.Key))
                    .ToDictionary(group => group.Key, group => group.ToList());

            // iterate through each unique sequence
            foreach (var sequencePeakListKvp in seqPeakDict)
            {
                List<ChromatographicPeak> peaksForPeptide = sequencePeakListKvp.Value;
                if (!peaksForPeptide.Any())
                    continue;

                ChromatographicPeak bestPeak = ChooseBestPeak(peaksForPeptide);

                if (bestPeak == null) continue;
                if (DonorFileToPeakDict.ContainsKey(bestPeak.SpectraFileInfo))
                {
                    DonorFileToPeakDict[bestPeak.SpectraFileInfo].Add(bestPeak);
                }
                else
                {
                    DonorFileToPeakDict.Add(bestPeak.SpectraFileInfo, new List<ChromatographicPeak> { bestPeak });
                }
            }
        }

        internal ChromatographicPeak ChooseBestPeak(List<ChromatographicPeak> peaks)
        {
            ChromatographicPeak bestPeak = null;
            switch (FlashParams.DonorCriterion)
            {
                case DonorCriterion.Score: // Select best peak by the PSM score
                    bestPeak = peaks.MaxBy(peak => peak.Identifications.Max(id => id.PsmScore));
                    if (bestPeak.Identifications.First().PsmScore > 0)
                        break;
                    else // if every ID has a score of zero, let it fall through to the default case
                        goto default;
                case DonorCriterion.Neighbors: // Select peak with the most neighboring peaks
                    int maxPeaks = 0;
                    foreach (var donorPeak in peaks)
                    {
                        // Count the number of neighboring peaks with unique peptides
                        int neighboringPeaksCount = _results.Peaks[donorPeak.SpectraFileInfo]
                            .Where(peak => Math.Abs(peak.ApexRetentionTime - donorPeak.ApexRetentionTime) < MbrAlignmentWindow)
                            .Select(peak => peak.Identifications.First().ModifiedSequence)
                            .Distinct()
                            .Count();

                        if (neighboringPeaksCount > maxPeaks)
                        {
                            maxPeaks = neighboringPeaksCount;
                            bestPeak = donorPeak;
                        }
                    }
                    break;
                case DonorCriterion.Intensity: // Select the peak with the highest intensity
                default:
                    bestPeak = peaks.MaxBy(peak => peak.Intensity);
                    break;
            }

            return bestPeak;
        }

        /// <summary>
        /// Used by MBR. Predicts the retention time of a peak in an acceptor file based on the 
        /// retention time of the peak in the donor file. This is done with a local alignment
        /// where all peaks within 30 seconds of the donor peak are matched to peaks with the same associated peptide in the acceptor file,
        /// if such a peak exists.
        /// </summary>
        /// <param name="rtCalibrationCurve">Array of all shared peaks between the donor and the acceptor file</param>
        /// <returns> RtInfo object containing the predicted retention time of the acceptor peak and the width of the predicted retention time window </returns>
        internal RtInfo PredictRetentionTime(
            RetentionTimeCalibDataPoint[] rtCalibrationCurve,
            ChromatographicPeak donorPeak,
            SpectraFileInfo acceptorFile,
            bool acceptorSampleIsFractionated,
            bool donorSampleIsFractionated)
        {
            var nearbyCalibrationPoints = new List<RetentionTimeCalibDataPoint>(); // The number of anchor peptides to be used for local alignment (on either side of the donor peptide)

            // only compare +- 1 fraction
            if (acceptorSampleIsFractionated && donorSampleIsFractionated)
            {
                int acceptorFractionNumber = acceptorFile.Fraction;
                int donorFractionNumber = donorPeak.SpectraFileInfo.Fraction;

                if (Math.Abs(acceptorFractionNumber - donorFractionNumber) > 1)
                {
                    return null;
                }
            }

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

            int numberOfForwardAnchors = 0;
            // gather nearby data points
            for (int r = index + 1; r < rtCalibrationCurve.Length; r++)
            {
                double rtDiff = rtCalibrationCurve[r].DonorFilePeak.Apex.IndexedPeak.RetentionTime - donorPeak.Apex.IndexedPeak.RetentionTime;
                if (rtCalibrationCurve[r].AcceptorFilePeak != null
                    && rtCalibrationCurve[r].AcceptorFilePeak.ApexRetentionTime > 0)
                {
                    if (Math.Abs(rtDiff) > 0.5) // If the rtDiff is too large, it's no longer local alignment
                    {
                        break;
                    }
                    nearbyCalibrationPoints.Add(rtCalibrationCurve[r]);
                    numberOfForwardAnchors++;
                    if (numberOfForwardAnchors >= NumberOfAnchorPeptidesForMbr) // We only want a handful of anchor points
                    {
                        break;
                    }
                }
            }

            int numberOfBackwardsAnchors = 0;
            for (int r = index - 1; r >= 0; r--)
            {
                double rtDiff = rtCalibrationCurve[r].DonorFilePeak.Apex.IndexedPeak.RetentionTime - donorPeak.Apex.IndexedPeak.RetentionTime;
                if (rtCalibrationCurve[r].AcceptorFilePeak != null
                    && rtCalibrationCurve[r].AcceptorFilePeak.ApexRetentionTime > 0)
                {
                    if (Math.Abs(rtDiff) > 0.5) // If the rtDiff is too large, it's no longer local alignment
                    {
                        break;
                    }
                    nearbyCalibrationPoints.Add(rtCalibrationCurve[r]);
                    numberOfBackwardsAnchors++;
                    if (numberOfBackwardsAnchors >= NumberOfAnchorPeptidesForMbr) // We only want a handful of anchor points
                    {
                        break;
                    }
                }
            }

            if (!nearbyCalibrationPoints.Any())
            {
                // If there are no nearby calibration points, return the donor peak's RT and a width of 15 seconds
                return new RtInfo(predictedRt: donorPeak.Apex.IndexedPeak.RetentionTime, width: 0.25);
            }

            // calculate difference between acceptor and donor RTs for these RT region
            List<double> rtDiffs = nearbyCalibrationPoints
                .Select(p => p.DonorFilePeak.ApexRetentionTime - p.AcceptorFilePeak.ApexRetentionTime)
                .ToList();

            double medianRtDiff = rtDiffs.Median();
            if(rtDiffs.Count == 1)
            {
                // If there are no nearby calibration points, return the donor peak's RT and a width of 15 seconds
                return new RtInfo(predictedRt: donorPeak.Apex.IndexedPeak.RetentionTime - medianRtDiff, width: 0.25);
            }

            double rtRange = rtDiffs.StandardDeviation() * 6;

            rtRange = Math.Min(rtRange, FlashParams.MaxMbrRtWindow);

            return new RtInfo(predictedRt: donorPeak.Apex.IndexedPeak.RetentionTime - medianRtDiff, width: rtRange);
        }

        /// <summary>
        /// Returns a pseudo-randomly selected peak that does not have the same mass as the donor
        /// </summary>
        /// <param name="peaksOrderedByMass"></param>
        /// <param name="donorPeakPeakfindingMass"> Will search for a peak at least 5 Da away from the peakfinding mass </param>
        /// <returns></returns>
        internal ChromatographicPeak GetRandomPeak(
            List<ChromatographicPeak> peaksOrderedByMass,
            double donorPeakRetentionTime,
            double retentionTimeMinDiff,
            Identification donorIdentification)
        {
            double minDiff = 5 * PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
            double maxDiff = 11 * PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
            double donorPeakPeakfindingMass = donorIdentification.PeakfindingMass;

            // Theoretically we could do a binary search but we're just going to iterate through the whole list of donor peaks
            List<ChromatographicPeak> randomPeakCandidates = peaksOrderedByMass
                .Where(p => 
                    p.ApexRetentionTime > 0
                    && Math.Abs(p.ApexRetentionTime - donorPeakRetentionTime) > retentionTimeMinDiff
                    && p.Identifications.First().BaseSequence != donorIdentification.BaseSequence
                    && Math.Abs(p.Identifications.First().PeakfindingMass - donorPeakPeakfindingMass) > minDiff
                    && Math.Abs(p.Identifications.First().PeakfindingMass - donorPeakPeakfindingMass) < maxDiff)
                .ToList();

            while (!randomPeakCandidates.Any() & maxDiff < 1e5)
            {
                // Increase the search space by a factor of 10 and try again
                maxDiff *= 10;
                randomPeakCandidates = peaksOrderedByMass
                .Where(p =>
                    p.ApexRetentionTime > 0
                    && Math.Abs(p.ApexRetentionTime - donorPeakRetentionTime) > retentionTimeMinDiff
                    && p.Identifications.First().BaseSequence != donorIdentification.BaseSequence
                    && Math.Abs(p.Identifications.First().PeakfindingMass - donorPeakPeakfindingMass) > minDiff
                    && Math.Abs(p.Identifications.First().PeakfindingMass - donorPeakPeakfindingMass) < maxDiff)
                .ToList();
            }

            if (!randomPeakCandidates.Any())
            {
                return null;
            }

            // Generates a pseudo-random number based on the donor peak finding mass + retention time
            int pseudoRandomNumber = (int)(1e5 * (donorIdentification.PeakfindingMass % 1.0) * (donorIdentification.Ms2RetentionTimeInMinutes % 1.0)) % randomPeakCandidates.Count;
            return randomPeakCandidates[pseudoRandomNumber];
        }

        /// <summary>
        /// This method maps identified peaks from other chromatographic runs ("donors") onto
        /// the defined chromatographic run ("acceptor"). The goal is to reduce the number of missing
        /// intensity measurements. Missing values occur generally either because 1) the analyte is
        /// in the sample but didn't get fragmented/identified or 2) the analyte is genuinely missing
        /// from the sample.
        /// </summary>
        private void QuantifyMatchBetweenRunsPeaks(SpectraFileInfo acceptorFile)
        {
            bool acceptorSampleIsFractionated = _results.SpectraFiles
                .Where(p => p.Condition == acceptorFile.Condition && p.BiologicalReplicate == acceptorFile.BiologicalReplicate)
                .Select(p => p.Fraction)
                .Distinct()
                .Count() > 1;

            // acceptor file known peaks
            var acceptorFileIdentifiedPeaks = _results.Peaks[acceptorFile]
                .Where(p => p.Identifications.Any(id => PeptideModifiedSequencesToQuantify.Contains(id.ModifiedSequence)))
                .ToList();

            // these are the analytes already identified in this run. we don't need to try to match them from other runs
            var acceptorFileIdentifiedSequences = new HashSet<string>(acceptorFileIdentifiedPeaks
                .Where(peak => peak.IsotopicEnvelopes.Any() && peak.Identifications.Min(id => id.QValue) < 0.01)
                .SelectMany(p => p.Identifications.Select(d => d.ModifiedSequence)));

            MbrScorer scorer = MbrScorerFactory.BuildMbrScorer(acceptorFileIdentifiedPeaks, FlashParams, out var mbrTol);
            if (scorer == null)
                return;

            mbrTol = new PpmTolerance(FlashParams.MbrPpmTolerance);
            HashSet<ProteinGroup> thisFilesMsmsIdentifiedProteins = new HashSet<ProteinGroup>();
            if (FlashParams.RequireMsmsIdInCondition)
            {
                // only match peptides from proteins that have at least one MS/MS identified peptide in the condition
                foreach (SpectraFileInfo conditionFile in SpectraFileInfoList.Where(p => p.Condition == acceptorFile.Condition))
                {
                    foreach (ProteinGroup proteinGroup in _results.Peaks[conditionFile].Where(p => p.DetectionType == DetectionType.MSMS).SelectMany(p => p.Identifications.SelectMany(v => v.ProteinGroups)))
                    {
                        thisFilesMsmsIdentifiedProteins.Add(proteinGroup);
                    }
                }
            }

            // this stores the results of MBR
            ConcurrentDictionary<string, ConcurrentDictionary<IsotopicEnvelope, List<MbrChromatographicPeak>>> matchBetweenRunsIdentifiedPeaks = new();

            // map each donor file onto this file
            foreach (var donorFilePeakListKvp in DonorFileToPeakDict)
            {
                if (acceptorFile.Equals(donorFilePeakListKvp.Key))
                {
                    continue;
                }

                // this is the list of peaks identified in the other file but not in this one ("ID donor peaks")
                List<ChromatographicPeak> idDonorPeaks = donorFilePeakListKvp.Value
                    .Where(p => 
                        !acceptorFileIdentifiedSequences.Contains(p.Identifications.First().ModifiedSequence)
                        && (!FlashParams.RequireMsmsIdInCondition 
                            || p.Identifications.Any(v => v.ProteinGroups.Any(g => thisFilesMsmsIdentifiedProteins.Contains(g))))
                        && this.PeptideModifiedSequencesToQuantify.Contains(p.Identifications.First().ModifiedSequence))
                    .ToList();

                if (!idDonorPeaks.Any())
                {
                    continue;
                }

                bool donorSampleIsFractionated = _results.SpectraFiles
                    .Where(p => p.Condition == donorFilePeakListKvp.Key.Condition && p.BiologicalReplicate == donorFilePeakListKvp.Key.BiologicalReplicate)
                    .Select(p => p.Fraction)
                    .Distinct()
                    .Count() > 1;

                // We're only interested in the fold change if the conditions are different. Otherwise, we score based off of the intensities
                // of the acceptor file
                if (SpectraFileInfoList.Select(p => p.Condition).Distinct().Count() > 1
                    && donorFilePeakListKvp.Key.Condition != acceptorFile.Condition)
                {
                    scorer.CalculateFoldChangeBetweenFiles(idDonorPeaks);
                }

                // generate RT calibration curve
                RetentionTimeCalibDataPoint[] rtCalibrationCurve = GetRtCalSpline(donorFilePeakListKvp.Key, acceptorFile, scorer, out var donorPeaksMassOrdered);

                // break if MBR transfers can't be scored
                if (!scorer.IsValid(donorFilePeakListKvp.Key)) continue;

                // Loop through every MSMS id in the donor file
                Parallel.ForEach(Partitioner.Create(0, idDonorPeaks.Count),
                    new ParallelOptions { MaxDegreeOfParallelism = FlashParams.MaxThreads },
                    (range, loopState) =>
                    {
                        for (int i = range.Item1; i < range.Item2; i++)
                        {
                            ChromatographicPeak donorPeak = idDonorPeaks[i];
                            // TODO: Add a toggle that set rtRange to be maximum width
                            RtInfo rtInfo = PredictRetentionTime(rtCalibrationCurve, donorPeak, acceptorFile, acceptorSampleIsFractionated, donorSampleIsFractionated);
                            if (rtInfo == null) continue;

                            // Look for MBR target (predicted-RT peak)
                            FindAllAcceptorPeaks(acceptorFile, scorer, rtInfo, mbrTol, donorPeak, out var bestAcceptor);
                            AddPeakToConcurrentDict(matchBetweenRunsIdentifiedPeaks, bestAcceptor, donorPeak.Identifications.First());

                            //Draw a random donor that has an rt sufficiently far enough away
                            double minimumRtDifference = rtInfo.Width*2;
                            ChromatographicPeak randomDonor = GetRandomPeak(donorPeaksMassOrdered,
                                donorPeak.ApexRetentionTime,
                                minimumRtDifference,
                                donorPeak.Identifications.First());

                            // Look for MBR decoy (random-RT peak) 
                            MbrChromatographicPeak bestDecoy = null;
                            RtInfo decoyRtInfo = null;
                            if (randomDonor != null)
                            {
                                decoyRtInfo = PredictRetentionTime(rtCalibrationCurve, randomDonor, acceptorFile, acceptorSampleIsFractionated, donorSampleIsFractionated);
                                if (decoyRtInfo != null)
                                {
                                    //Find a decoy peak using the randomly drawn retention time
                                    FindAllAcceptorPeaks(acceptorFile, scorer, rtInfo, mbrTol, donorPeak, out bestDecoy,
                                        randomRt: decoyRtInfo.PredictedRt);
                                    AddPeakToConcurrentDict(matchBetweenRunsIdentifiedPeaks, bestDecoy, donorPeak.Identifications.First());
                                }
                            }

                            double windowWidth = Math.Max(0.5, rtInfo.Width);
                            // If the search turned up empty, try again with a wider search window
                            while (bestAcceptor == null && bestDecoy == null)
                            {
                                windowWidth = Math.Min(windowWidth, FlashParams.MaxMbrRtWindow);
                                rtInfo.Width = windowWidth;
                                FindAllAcceptorPeaks(acceptorFile, scorer, rtInfo, mbrTol, donorPeak, out bestAcceptor);
                                AddPeakToConcurrentDict(matchBetweenRunsIdentifiedPeaks, bestAcceptor, donorPeak.Identifications.First());

                                if(decoyRtInfo != null)
                                {
                                    decoyRtInfo.Width = windowWidth;
                                    FindAllAcceptorPeaks(acceptorFile, scorer, rtInfo, mbrTol, donorPeak, out bestDecoy,
                                    randomRt: decoyRtInfo.PredictedRt);
                                    AddPeakToConcurrentDict(matchBetweenRunsIdentifiedPeaks, bestDecoy, donorPeak.Identifications.First());
                                }
                                if (windowWidth >= FlashParams.MaxMbrRtWindow)
                                {
                                    break;
                                }
                                else
                                {
                                    windowWidth *= 2;
                                }
                            }

                        }
                    });
            }

            // Eliminate duplicate peaks (not sure where they come from)
            foreach (var seqDictionaryKvp in matchBetweenRunsIdentifiedPeaks)
            {
                // Each isotopic envelope is linked to a list of ChromatographicPeaks
                // Here, we remove instances where the same envelope is associated with multiple chromatographic peaks but the peaks correspond to the same donor peptide
                // I don't know why this happens lol
                // If multiple peaks are associated with the same envelope, and they have different associated peptide identifications, then they're kept separate.
                foreach (var envelopePeakListKvp in seqDictionaryKvp.Value)
                {
                    List<MbrChromatographicPeak> bestPeaks = new();
                    foreach (var peakGroup in envelopePeakListKvp.Value.GroupBy(peak => peak.Identifications.First().ModifiedSequence))
                    {
                        bestPeaks.Add(peakGroup.MaxBy(peak => peak.MbrScore));
                    }
                    envelopePeakListKvp.Value.Clear();
                    envelopePeakListKvp.Value.AddRange(bestPeaks);
                }
            }

            // Create a dictionary that stores imsPeak associated with an ms/ms identified peptide
            Dictionary<int, List<IIndexedPeak>> msmsImsPeaks = _results.Peaks[acceptorFile]
                .Where(peak => 
                        !peak.DecoyPeptide 
                        && peak.Apex?.IndexedPeak != null 
                        && PeptideModifiedSequencesToQuantify.Contains(peak.Identifications.First().ModifiedSequence))
                .Select(peak => peak.Apex.IndexedPeak)
                .GroupBy(imsPeak => imsPeak.ZeroBasedScanIndex)
                .ToDictionary(g => g.Key, g => g.ToList());

            // take the best result (highest scoring) for each peptide after we've matched from all the donor files
            foreach (var mbrIdentifiedPeptide in matchBetweenRunsIdentifiedPeaks.Where(p => !acceptorFileIdentifiedSequences.Contains(p.Key)))
            {
                string peptideModifiedSequence = mbrIdentifiedPeptide.Key;
                if (!mbrIdentifiedPeptide.Value.Any())
                {
                    continue;
                }

                foreach (var peakHypothesisGroup in mbrIdentifiedPeptide.Value.SelectMany(kvp => kvp.Value).OrderByDescending(p => p.MbrScore).GroupBy(p => p.RandomRt))
                {
                    var peakHypotheses = peakHypothesisGroup.ToList();
                    MbrChromatographicPeak best = peakHypotheses.First();
                    peakHypotheses.Remove(best);

                    // Discard any peaks that are already associated with an ms/ms identified peptide
                    while (best?.Apex?.IndexedPeak != null && msmsImsPeaks.TryGetValue(best.Apex.IndexedPeak.ZeroBasedScanIndex, out var peakList))
                    {
                        if (peakList.Contains(best.Apex.IndexedPeak))
                        {
                            if (!peakHypotheses.Any())
                            {
                                best = null;
                                break;
                            }
                            best = peakHypotheses.First();
                            peakHypotheses.Remove(best);
                        }
                        else
                        {
                            break;
                        }
                    }
                    if (best == null) continue;

                    // merge peaks with different charge states
                    if (peakHypotheses.Count > 0)
                    {
                        double start = best.IsotopicEnvelopes.Min(p => p.IndexedPeak.RetentionTime);
                        double end = best.IsotopicEnvelopes.Max(p => p.IndexedPeak.RetentionTime);

                        _results.Peaks[acceptorFile].Add(best);
                        foreach (ChromatographicPeak peak in peakHypotheses.Where(p => p.Apex.ChargeState != best.Apex.ChargeState))
                        {
                            if (peak.Apex.IndexedPeak.RetentionTime >= start
                                && peak.Apex.IndexedPeak.RetentionTime <= end)
                                //&& Math.Abs(peak.MbrScore - best.MbrScore) / best.MbrScore < 0.25)// 25% difference is a rough heuristic, but I don't want super shitty peaks being able to supercede the intensity of a good peak!
                            {
                                if (msmsImsPeaks.TryGetValue(peak.Apex.IndexedPeak.ZeroBasedScanIndex, out var peakList) && peakList.Contains(peak.Apex.IndexedPeak))
                                {
                                    continue; // If the peak is already accounted for, skip it.
                                }
                                else
                                {
                                    best.MergeFeatureWith(peak, FlashParams.Integrate);
                                }
                            }
                        }
                    }
                    _results.Peaks[acceptorFile].Add(best);
                }
            }

            RunErrorChecking(acceptorFile);
        }

        /// <summary>
        /// A concurrent dictionary is used to keep track of MBR peaks that have been identified in the acceptor file. This function updates that dictionary
        /// </summary>
        /// <param name="matchBetweenRunsIdentifiedPeaks"> concurrent dictionary. Key = Peptide sequence. Value = ConcurrentDictionary mapping where keys are isotopic envelopes and values are list of associated peaks</param>
        /// <param name="peakToSave">Peak to add to the dictionary</param>
        /// <param name="donorIdentification">The donor ID associated with the MBR peaks</param>
        private void AddPeakToConcurrentDict(ConcurrentDictionary<string, ConcurrentDictionary<IsotopicEnvelope, List<MbrChromatographicPeak>>> matchBetweenRunsIdentifiedPeaks,
            MbrChromatographicPeak peakToSave,
            Identification donorIdentification)
        {
            if(peakToSave == null)
            {
                return;
            }
            // save the peak hypothesis
            matchBetweenRunsIdentifiedPeaks.AddOrUpdate
            (
                // new key
                key: donorIdentification.ModifiedSequence,
                // if we are adding a value for the first time, we simply create a new dictionatry with one entry
                addValueFactory: (sequenceKey) =>
                new ConcurrentDictionary<IsotopicEnvelope, List<MbrChromatographicPeak>>(
                    new Dictionary<IsotopicEnvelope, List<MbrChromatographicPeak>>
                    {
                        { peakToSave.Apex, new List<MbrChromatographicPeak> { peakToSave } }
                    }),
                // if the key (sequence) already exists, we have to add the new peak to the existing dictionary
                updateValueFactory: (sequenceKey, envelopePeakListDict) =>
                {
                    envelopePeakListDict.AddOrUpdate(
                        key: peakToSave.Apex,
                        addValueFactory: (envelopeKey) => new List<MbrChromatographicPeak> { peakToSave }, // if the key (envelope) doesnt exist, just create a new list
                        updateValueFactory: (envelopeKey, peakList) => { peakList.Add(peakToSave); return peakList; }); // if the key (envelope) already exists, add the peak to the associated list
                    return envelopePeakListDict;
                }
            );
        }

        /// <summary>
        /// Finds MBR acceptor peaks by looping  through every possible peak for every possible charge state
        /// in a given retention time range. Identified peaks are added to the matchBetweenRunsIdentifiedPeaks dictionary.
        /// </summary>
        /// <param name="scorer"> The MbrScorer object used to score acceptor peaks</param>
        /// <param name="rtInfo"> RtInfo object containing the predicted retention time for the acceptor peak and the width of the expected RT window</param>
        /// <param name="fileSpecificTol"> Ppm Tolerance specific to the acceptor file</param>
        /// <param name="donorPeak"> The donor peak. Acceptor peaks are presumed to represent the same peptide ast he donor peak</param>
        /// <param name="matchBetweenRunsIdentifiedPeaksThreadSpecific"> A dictionary containing peptide sequences and their associated mbr peaks </param>
        internal void FindAllAcceptorPeaks(
            SpectraFileInfo acceptorFile, 
            MbrScorer scorer,
            RtInfo rtInfo,
            PpmTolerance fileSpecificTol,
            ChromatographicPeak donorPeak,
            out MbrChromatographicPeak bestAcceptor,
            double? randomRt = null)
        {
            // get the MS1 scan info for this region so we can look up indexed peaks
            ScanInfo[] ms1ScanInfos = IndexingEngineDictionary[acceptorFile].ScanInfoArray;
            ScanInfo start = ms1ScanInfos[0];
            ScanInfo end = ms1ScanInfos[ms1ScanInfos.Length - 1];
            double rtStartHypothesis = randomRt == null ? rtInfo.RtStartHypothesis : (double)randomRt - (rtInfo.Width / 2.0);
            double rtEndHypothesis = randomRt == null ? rtInfo.RtEndHypothesis : (double)randomRt + (rtInfo.Width / 2.0);

            // Try to snip the MS1 scans to the region where the analyte should appear
            for (int j = 0; j < ms1ScanInfos.Length; j++)
            {
                ScanInfo scan = ms1ScanInfos[j];
                if (scan.RetentionTime <= rtStartHypothesis)
                {
                    start = scan;
                }
                if (scan.RetentionTime >= rtEndHypothesis)
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

            Identification donorIdentification = donorPeak.Identifications.First();
            bestAcceptor = null;

            foreach (int z in chargesToMatch)
            {
                List<IIndexedPeak> chargeXic = new List<IIndexedPeak>();

                for (int j = start.ZeroBasedScanIndex; j <= end.ZeroBasedScanIndex; j++)
                {
                    IIndexedPeak peak = IndexingEngineDictionary[acceptorFile]
                        .GetIndexedPeak(donorIdentification.PeakfindingMass.ToMz(z), j, fileSpecificTol);
                    if (peak != null)
                        chargeXic.Add(peak);
                }
                if (!chargeXic.Any())
                    continue;

                List<IsotopicEnvelope> chargeEnvelopes = GetIsotopicEnvelopes(chargeXic, donorIdentification, z, acceptorFile).OrderBy(env => env.Intensity).ToList();

                // treat each isotopic envelope in the valid region as a potential seed for a chromatographic peak.
                // remove the clustered isotopic envelopes from the list of seeds after each iteration
                while (chargeEnvelopes.Any())
                {
                    MbrChromatographicPeak acceptorPeak = FindIndividualAcceptorPeak(acceptorFile, scorer, donorPeak,
                        fileSpecificTol, rtInfo, z, chargeEnvelopes, randomRt);
                    if (acceptorPeak == null)
                        continue;
                    if (bestAcceptor == null || bestAcceptor.MbrScore < acceptorPeak.MbrScore)
                    {
                        acceptorPeak.ChargeList = chargesToMatch;
                        bestAcceptor = acceptorPeak;
                    }
                }
            }
        }

        /// <summary>
        /// Grabs the first isotopic envelope in the list of chargeEnvelopes as a potential seed for a chromatographic peak.
        /// remove the isotopic envelope from chargeEnvelopes afterward.
        /// </summary>
        /// <param name="acceptorFile"></param>
        /// <param name="mbrTol"></param>
        /// <param name="rtInfo"></param>
        /// <param name="rtScoringDistribution"></param>
        /// <param name="z"></param>
        /// <param name="chargeEnvelopes"></param>
        /// <returns> An acceptor chromatographic peak, unless the peak found was already linked to an MS/MS id, in which case it return null. </returns>
        internal MbrChromatographicPeak FindIndividualAcceptorPeak(
            SpectraFileInfo acceptorFile,
            MbrScorer scorer,
            ChromatographicPeak donorPeak,
            PpmTolerance mbrTol,
            RtInfo rtInfo,
            int z,
            List<IsotopicEnvelope> chargeEnvelopes,
            double? randomRt = null)
        {
            var donorId = donorPeak.Identifications.OrderBy(p => p.QValue).First();
            var acceptorPeak = new MbrChromatographicPeak(donorId, acceptorFile, randomRt ?? rtInfo.PredictedRt, randomRt != null);

            // Grab the first scan/envelope from charge envelopes. This should be the most intense envelope in the list
            IsotopicEnvelope seedEnv = chargeEnvelopes.First();
            var xic = IndexingEngineDictionary[acceptorFile].GetXic(donorId.PeakfindingMass.ToMz(z), seedEnv.IndexedPeak.RetentionTime, mbrTol, MissedScansAllowed);
            List<IsotopicEnvelope> bestChargeEnvelopes = GetIsotopicEnvelopes(xic, donorId, z, acceptorFile);
            acceptorPeak.IsotopicEnvelopes.AddRange(bestChargeEnvelopes);
            acceptorPeak.CalculateIntensityForThisFeature(FlashParams.Integrate);

            CutPeak(acceptorPeak, seedEnv.IndexedPeak.RetentionTime);

            var claimedPeaks = new HashSet<IIndexedPeak>(acceptorPeak.IsotopicEnvelopes.Select(p => p.IndexedPeak))
            {
                seedEnv.IndexedPeak // prevents infinite loops
            };
            chargeEnvelopes.RemoveAll(p => claimedPeaks.Contains(p.IndexedPeak));

            // peak has already been identified by MSMS - skip it
            if (scorer.ApexToAcceptorFilePeakDict.ContainsKey(seedEnv.IndexedPeak))
            {
                return null;
            }

            acceptorPeak.MbrScore = scorer.ScoreMbr(acceptorPeak, donorPeak, randomRt ?? rtInfo.PredictedRt);

            return acceptorPeak;
        }


        /// <summary>
        /// Checks for and resolves situations where one IndexedMassSpectralPeak is defined as the apex 
        /// for multiple ChromatographicPeaks. In these situations, the two peaks are merged and the merged
        /// peak is stored in the FlashLFQ results.
        /// </summary>
        /// <param name="spectraFile"></param>
        private void RunErrorChecking(SpectraFileInfo spectraFile)
        {
            if (!FlashParams.Silent)
            {
                Console.WriteLine("Checking errors");
            }

            _results.Peaks[spectraFile].RemoveAll(p => p == null || p.DetectionType == DetectionType.MBR && !p.IsotopicEnvelopes.Any());

            // merge duplicate peaks and handle MBR/MSMS peakfinding conflicts
            var errorCheckedPeaksGroupedByApex = new Dictionary<IIndexedPeak, ChromatographicPeak>();
            var errorCheckedPeaks = new List<ChromatographicPeak>();
            
            foreach (ChromatographicPeak tryPeak in _results.Peaks[spectraFile].OrderBy(p => p.DetectionType == DetectionType.MBR))
            {
                tryPeak.CalculateIntensityForThisFeature(FlashParams.Integrate);
                tryPeak.ResolveIdentifications();

                if (tryPeak.Apex == null)
                {
                    if (tryPeak.DetectionType == DetectionType.MBR)
                    {
                        continue;
                    }

                    errorCheckedPeaks.Add(tryPeak);
                    continue;
                }

                IIndexedPeak apexImsPeak = tryPeak.Apex.IndexedPeak;
                if (errorCheckedPeaksGroupedByApex.TryGetValue(apexImsPeak, out ChromatographicPeak storedPeak) && storedPeak != null)
                {
                    // At here, the peaks detected from the IsoTracker will be confident, then we don't want to eliminate.
                    // Logically, we view the IsoTracker_MBR and IsoTracker_Ambiguity as the MSMS.
                    // Therefore, we only need to check the MBR peaks.
                    if (tryPeak.DetectionType != DetectionType.MBR && storedPeak.DetectionType != DetectionType.MBR)
                    {
                        if (PeptideModifiedSequencesToQuantify.Contains(tryPeak.Identifications.First().ModifiedSequence))
                        {
                            if (PeptideModifiedSequencesToQuantify.Contains(storedPeak.Identifications.First().ModifiedSequence))
                            {
                                // if the try peak is merge into stored peaks, the try peak should be labeled as MSMSAmbiguousPeakfinding and then the intensity should be set as 0 in the peptide result
                                if (tryPeak.DetectionType == DetectionType.IsoTrack_MSMS ||
                                    tryPeak.DetectionType == DetectionType.IsoTrack_MBR)
                                {
                                    tryPeak.DetectionType = DetectionType.MSMSAmbiguousPeakfinding;
                                }
                                storedPeak.MergeFeatureWith(tryPeak, FlashParams.Integrate);
                            }
                            else
                            {
                                // If the stored peak id isn't in the list of peptides to quantify, overwrite it
                                errorCheckedPeaksGroupedByApex[tryPeak.Apex.IndexedPeak] = tryPeak;
                            }
                        }
                    }
                    else if (tryPeak.DetectionType == DetectionType.MBR && storedPeak.DetectionType != DetectionType.MBR)
                    {
                        // Default to MSMS peaks over MBR Peaks.
                        // Most of these have already been eliminated
                        // However, sometimes merging MBR peaks with different charge states reveals that
                        // The MBR peak conflicts with an MSMS peak
                        // Removing the peak when this happens is a conservative step.
                        // Sometimes the MSMS peak is a decoy, or has a peptides level Q-value < 0.01 (i.e., the modified sequence isn't in PeptideModifiedSequencesToQuantify).
                        // In this case, we keep the MBR peak.
                        if (storedPeak.DecoyPeptide || !PeptideModifiedSequencesToQuantify.Contains(storedPeak.Identifications.First().ModifiedSequence))
                        {
                            errorCheckedPeaksGroupedByApex[tryPeak.Apex.IndexedPeak] = tryPeak;
                        }
                        else
                        {
                            continue;
                        }
                    }
                    else if (tryPeak.DetectionType == DetectionType.MBR && storedPeak.DetectionType == DetectionType.MBR)
                    {
                        if (tryPeak.Identifications.First().ModifiedSequence == storedPeak.Identifications.First().ModifiedSequence)
                        {
                            storedPeak.MergeFeatureWith(tryPeak, FlashParams.Integrate);
                        }
                        else if (((MbrChromatographicPeak)tryPeak).MbrScore > ((MbrChromatographicPeak)storedPeak).MbrScore)
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

        private bool RunPEPAnalysis()
        {
            List<MbrChromatographicPeak> mbrPeaks = _results.Peaks.SelectMany(kvp => kvp.Value)
                .Where(peak => peak.DetectionType == DetectionType.MBR)
                .Cast<MbrChromatographicPeak>()
                .OrderByDescending(peak => peak.MbrScore)
                .ToList();

            if (!mbrPeaks.IsNotNullOrEmpty()) return false;
            int decoyPeakTotal = mbrPeaks.Count(peak => peak.RandomRt);

            List<double> tempPepQs = new();
            List<double> tempQs = new();
            if (mbrPeaks.Count > 100 && decoyPeakTotal > 20)
            {
                PepAnalysisEngine pepAnalysisEngine = new PepAnalysisEngine(
                    mbrPeaks,
                    outputFolder: Path.GetDirectoryName(SpectraFileInfoList.First().FullFilePathWithExtension),
                    maxThreads: FlashParams.MaxThreads,
                    pepTrainingFraction: PepTrainingFraction);
                var pepOutput = pepAnalysisEngine.ComputePEPValuesForAllPeaks();

                _results.PepResultString = pepOutput;

                return true;
            }
            return false;
        }

        /// <summary>
        /// Calculates the FDR for each MBR-detected peak using decoy peaks and decoy peptides,
        /// Then filters out all peaks below a given FDR threshold
        /// </summary>
        private void CalculateFdrForMbrPeaks(SpectraFileInfo acceptorFile, bool usePep)
        {
            List<MbrChromatographicPeak> mbrPeaks;
            if (usePep)
            {
                // Take only the top scoring acceptor for each donor (acceptor can be target or decoy!)
                // Maybe we're sorting twice when we don't have to but idk if order is preserved using group by
                mbrPeaks = _results.Peaks[acceptorFile]
                    .Where(peak => peak.DetectionType == DetectionType.MBR)
                    .Cast<MbrChromatographicPeak>()
                    .GroupBy(peak => peak.Identifications.First())
                    .Select(group => group.OrderBy(peak => peak.MbrPep)
                    .ThenByDescending(peak => peak.MbrScore).First())
                    .OrderBy(peak => peak.MbrPep)
                    .ThenByDescending(peak => peak.MbrScore)
                    .ToList();

                _results.Peaks[acceptorFile] = mbrPeaks.Concat(_results.Peaks[acceptorFile].Where(peak => peak.DetectionType != DetectionType.MBR)).ToList();
            }
            else
            {
                // If PEP wasn't performed, things probably aren't calibrated very well, and so it's better
                // To err on the safe side and not remove the decoys
                mbrPeaks = _results.Peaks[acceptorFile]
                    .Where(peak => peak.DetectionType == DetectionType.MBR)
                    .Cast<MbrChromatographicPeak>()
                    .OrderByDescending(peak => peak.MbrScore)
                    .ToList();
            }

            if (!mbrPeaks.IsNotNullOrEmpty()) return;

            List<double> tempQs = new();
            int totalPeaks = 0;
            int decoyPeptides = 0;
            int decoyPeaks = 0;
            int doubleDecoys = 0;
            for (int i = 0; i < mbrPeaks.Count; i++)
            {
                totalPeaks++;
                switch((mbrPeaks[i].DecoyPeptide, mbrPeaks[i].RandomRt))
                {
                    case (false, false):
                        break;
                    case (true, false):
                        decoyPeptides++;
                        break;
                    case (false, true):
                        decoyPeaks++;
                        break;
                    case (true, true):
                        doubleDecoys++;
                        break;
                }

                // There are two parts to this score. We're summing the PEPs of peaks derived from target peptides. For peaks derived from decoy peptides,
                // We do the double decoy things where we count decoyPeptidePeaks - doubleDecoypeaks
                tempQs.Add(Math.Round(EstimateFdr(doubleDecoys, decoyPeptides, decoyPeaks, totalPeaks), 6));
            }

            // Set the q-value for each peak
            double[] correctedQs = CorrectQValues(tempQs);
            for (int i = 0; i < correctedQs.Length; i++)
            {
                mbrPeaks[i].MbrQValue = correctedQs[i];
            }
        }

        private int EstimateDecoyPeptideErrors(int decoyPeptideCount, int doubleDecoyCount)
        {
            return Math.Max(0, decoyPeptideCount - doubleDecoyCount);
        }

        private double EstimateFdr(int doubleDecoyCount, int decoyPeptideCount, int decoyPeakCount, int totalPeakCount)
        {
            return (double)(1 + decoyPeakCount + EstimateDecoyPeptideErrors(decoyPeptideCount, doubleDecoyCount)) / totalPeakCount;
        }

        /// <summary>
        /// Standard q-value correction, ensures that in a list of temporary q-values, a q-value is equal to
        /// Min(q-values, every q-value below in the list). As you work your way down a list of q-values, the value should only increase or stay the same.
        /// </summary>
        /// <param name="tempQs"></param>
        /// <returns></returns>
        private double[] CorrectQValues(List<double> tempQs)
        {
            if (!tempQs.IsNotNullOrEmpty()) return null;
            double[] correctedQValues = new double[tempQs.Count];
            correctedQValues[tempQs.Count - 1] = tempQs.Last();
            for(int i = tempQs.Count-2; i >=0; i--)
            {
                if (tempQs[i] > correctedQValues[i+1])
                {
                    correctedQValues[i] = correctedQValues[i + 1];
                }
                else
                {
                    correctedQValues[i] = tempQs[i];
                }
            }

            return correctedQValues;
        }

        #endregion

        /// <summary>
        /// Takes in a list of imsPeaks and finds all the isotopic peaks in each scan. If the experimental isotopic distribution
        /// matches the theoretical distribution, an IsotopicEnvelope object is created from the summed intensities of each isotopic peak.
        /// </summary>
        /// <param name="xic"> List of imsPeaks, where the mass of each peak is the peak finding mass (most abundant isotope) </param>
        /// <returns> A list of IsotopicEnvelopes, where each envelope contains the sum of the isotopic peak intensities from one scan </returns>
        public List<IsotopicEnvelope> GetIsotopicEnvelopes(
            List<IIndexedPeak> xic,
            Identification identification,
            int chargeState,
            SpectraFileInfo spectraFile)
        {
            var isotopicEnvelopes = new List<IsotopicEnvelope>();
            var isotopeMassShifts = ModifiedSequenceToIsotopicDistribution[identification.ModifiedSequence];

            if (isotopeMassShifts.Count < FlashParams.NumIsotopesRequired)
            {
                return isotopicEnvelopes;
            }

            PpmTolerance isotopeTolerance = new PpmTolerance(FlashParams.IsotopePpmTolerance);

            double[] experimentalIsotopeIntensities = new double[isotopeMassShifts.Count];
            double[] theoreticalIsotopeMassShifts = isotopeMassShifts.Select(p => p.Item1).ToArray();
            double[] theoreticalIsotopeAbundances = isotopeMassShifts.Select(p => p.Item2).ToArray();
            int peakfindingMassIndex = (int)Math.Round(identification.PeakfindingMass - identification.MonoisotopicMass, 0);
            List<IIndexedPeak> isotopologuePeaks = new List<IIndexedPeak>();

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
            foreach (IIndexedPeak peak in xic)
            {
                Array.Clear(experimentalIsotopeIntensities, 0, experimentalIsotopeIntensities.Length);
                foreach (var kvp in massShiftToIsotopePeaks)
                {
                    kvp.Value.Clear();
                }

                // isotope masses are calculated relative to the observed peak
                double observedMass = peak.M.ToMass(chargeState);
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

                            IIndexedPeak isotopePeak = IndexingEngineDictionary[spectraFile]
                                .GetIndexedPeak(isotopeMass.ToMz(chargeState), peak.ZeroBasedScanIndex, isotopeTolerance);

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
                                isotopologuePeaks.Add(isotopePeak);
                            }
                        }
                    }
                }

                // check number of isotope peaks observed
                if (massShiftToIsotopePeaks[0].Count < FlashParams.NumIsotopesRequired)
                {
                    continue;
                }

                // Check that the experimental envelope matches the theoretical
                if (CheckIsotopicEnvelopeCorrelation(massShiftToIsotopePeaks, peak, chargeState, isotopeTolerance, spectraFile, out var pearsonCorr))
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

                    isotopicEnvelopes.Add(new IsotopicEnvelope(peak, chargeState, experimentalIsotopeIntensities.Sum(), pearsonCorr));
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
            IIndexedPeak peak,
            int chargeState,
            PpmTolerance isotopeTolerance,
            SpectraFileInfo spectraFile,
            out double pearsonCorrelation)
        {
            pearsonCorrelation = Correlation.Pearson(
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
                IIndexedPeak unexpectedPeak = IndexingEngineDictionary[spectraFile]
                    .GetIndexedPeak(unexpectedMass.ToMz(chargeState), peak.ZeroBasedScanIndex, isotopeTolerance);

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
            return pearsonCorrelation > 0.7 && corrShiftedLeft - corrWithPadding < 0.1 && corrShiftedRight - corrWithPadding < 0.1;
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

            // Ordered list of all time points where the apex charge state had a valid isotopic envelope
            List<IsotopicEnvelope> timePointsForApexZ = peak.IsotopicEnvelopes
                .Where(p => p.ChargeState == peak.Apex.ChargeState).ToList();
            HashSet<int> scanNumbers = new HashSet<int>(timePointsForApexZ.Select(p => p.IndexedPeak.ZeroBasedScanIndex));
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

                    // Valley envelope is the lowest intensity point that has been encountered thus far
                    if (valleyEnvelope == null || timepoint.Intensity < valleyEnvelope.Intensity)
                    {
                        valleyEnvelope = timepoint;
                        indexOfValley = timePointsForApexZ.IndexOf(valleyEnvelope);
                    }

                    double discriminationFactor =
                        (timepoint.Intensity - valleyEnvelope.Intensity) / timepoint.Intensity;

                    // If the time point is at least discriminationFactor times more intense than the valley
                    // We perform an additional check to see if the time point is more intense than the point next to the valley
                    if (discriminationFactor > DiscriminationFactorToCutPeak &&
                        (indexOfValley + direction < timePointsForApexZ.Count && indexOfValley + direction >= 0))
                    {

                        IsotopicEnvelope secondValleyTimepoint = timePointsForApexZ[indexOfValley + direction];

                        discriminationFactor =
                            (timepoint.Intensity - secondValleyTimepoint.Intensity) / timepoint.Intensity;

                        // If the current timepoint is more intense than the second valley, we cut the peak
                        // If the scan following the valley isn't in the timePointsForApexZ list (i.e., no isotopic envelope is observed in the scan immediately after the valley), we also cut the peak
                        if (discriminationFactor > DiscriminationFactorToCutPeak || !scanNumbers.Contains(valleyEnvelope.IndexedPeak.ZeroBasedScanIndex + direction))
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
                peak.CalculateIntensityForThisFeature(FlashParams.Integrate);
                peak.SplitRT = valleyEnvelope.IndexedPeak.RetentionTime;

                // recursively cut
                CutPeak(peak, identificationTime);
            }
        }

        private void QuantifyIsobaricPeaks()
        {
            if(!FlashParams.Silent)
                Console.WriteLine("Quantifying isobaric species...");
            int isoGroupsSearched = 0;
            double lastReportedProgress = 0;
            double currentProgress = 0;


            Parallel.ForEach(Partitioner.Create(0, PeptideGroupsForIsoTracker.Count),
                new ParallelOptions { MaxDegreeOfParallelism = FlashParams.MaxThreads },
                (range, loopState) =>
                {
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        var idGroup = PeptideGroupsForIsoTracker[i];
                        List<XIC> xicGroup = new List<XIC>();
                        var mostCommonChargeIdGroup = idGroup.Identifications.GroupBy(p => p.PrecursorChargeState).OrderBy(p => p.Count()).Last();
                        var id = mostCommonChargeIdGroup.First();

                        // Try to get the primitive window for the XIC , from the (firstOne - 2) ->  (lastOne + 2)
                        double rightWindow = mostCommonChargeIdGroup.Select(p => p.Ms2RetentionTimeInMinutes).Max() + 2;
                        double leftWindow = mostCommonChargeIdGroup.Select(p => p.Ms2RetentionTimeInMinutes).Min() - 2;

                        //generate XIC from each file
                        foreach (var spectraFile in SpectraFileInfoList)
                        {
                            var xicIds = mostCommonChargeIdGroup.Where(p => p.FileInfo.Equals(spectraFile)).ToList();

                            //If in this run, there is no any ID, we will borrow the ID from other run to build the XIC.
                            if (xicIds.Count == 0)
                            {
                                xicIds = mostCommonChargeIdGroup.Where(p => !p.IsDecoy).ToList();
                            }

                            XIC xic = BuildXIC(xicIds, spectraFile, false, leftWindow, rightWindow);
                            if (xic != null && xic.Ms1Peaks.Count() > 5) // each XIC should have at least 5 points
                            {
                                xicGroup.Add(xic);
                            }
                        }

                        // In order to eliminate the bad XIC case that only one id for each file.
                        // We need to check if the XICGroup has more than one ID in one file.
                        if ( (!FlashParams.RequireMultipleIdsInOneFiles || MoreThanOneID(xicGroup)) && xicGroup.Count > 1)
                        {
                            // Step 1: Find the XIC with most IDs then, set as reference XIC
                            xicGroup.OrderByDescending(p => p.Ids.Count()).First().Reference = true;

                            //Step 2: Build the XICGroups
                            XICGroups xICGroups = new XICGroups(xicGroup);

                            //Step 3: Build the peakPairChrom dictionary
                            //The shared peak dictionary, each sharedPeak included serval ChromatographicPeak for each run [SharedPeak <-> ChromatographicPeaks]
                            Dictionary<PeakRegion, List<ChromatographicPeak>> sharedPeaksDict = new Dictionary<PeakRegion, List<ChromatographicPeak>>();

                            //Step 4: Iterate the shared peaks in the XICGroups
                            foreach (var peak in xICGroups.SharedPeaks)
                            {
                                double peakWindow = peak.Width;
                                if (peakWindow > 0.001) //make sure we have enough length of the window (0.001 min) for peak searching
                                {
                                    List<ChromatographicPeak> chromPeaksInSharedPeak = new List<ChromatographicPeak>();
                                    CollectChromPeakInRuns(peak, chromPeaksInSharedPeak, xICGroups);

                                    var allChromPeakInDict = sharedPeaksDict.SelectMany(pair => pair.Value).ToList();

                                    // ensuring no duplicates are added and only non-null peaks are considered.
                                    if (chromPeaksInSharedPeak.Where(p => p != null).Any() &&
                                        !chromPeaksInSharedPeak.Any(chromPeak => chromPeak != null && allChromPeakInDict.Any(peakInDict => peakInDict != null && peakInDict.Equals(chromPeak))))
                                    {
                                        sharedPeaksDict[peak] = chromPeaksInSharedPeak;
                                    }
                                }
                            }
                            IsobaricPeptideDict.TryAdd(id.ModifiedSequence, sharedPeaksDict);
                        }

                        // report search progress (proteins searched so far out of total proteins in database)
                        if (!FlashParams.Silent)
                        {
                            Interlocked.Increment(ref isoGroupsSearched);

                            double percentProgress = ((double)isoGroupsSearched / PeptideGroupsForIsoTracker.Count * 100);
                            currentProgress = Math.Max(percentProgress, currentProgress);

                            if (currentProgress > lastReportedProgress + 10)
                            {
                                Console.WriteLine("{0:0.}% of isobaric species quantified", currentProgress);
                                lastReportedProgress = currentProgress;
                            }
                        }
                    }
                });
            if (!FlashParams.Silent)
                Console.WriteLine("Finished quantifying isobaric species!");
        }

        /// <summary>
        /// Check if any XIC has more than one ID
        /// </summary>
        /// <param name="xics"></param>
        /// <returns></returns>
        internal bool MoreThanOneID(List<XIC> xics)
        {
            return xics.Any(p=> p.Ids.Count > 1);
        }

        /// <summary>
        /// Build the XIC from the given IDs.
        /// </summary>
        /// <param name="ids"> The IDs that produce the XIC</param>
        /// <param name="spectraFile">The file name</param>
        /// <param name="isReference"></param>
        /// <param name="start">The starting time</param>
        /// <param name="end">The end time</param>
        /// <returns></returns>
        internal XIC BuildXIC(List<Identification> ids, SpectraFileInfo spectraFile, bool isReference, double start, double end)
        {
            Identification id = ids.FirstOrDefault(); 
            var peakIndexingEngine = IndexingEngineDictionary[spectraFile];
            PpmTolerance isotopeTolerance = new PpmTolerance(FlashParams.PpmTolerance);
            ScanInfo[] ms1ScanInfos = peakIndexingEngine.ScanInfoArray;

            ScanInfo startScan = ms1ScanInfos
                .Where(p => p.RetentionTime < start)
                .OrderBy(p => p.RetentionTime)
                .LastOrDefault()
                ?? ms1ScanInfos.OrderBy(p => p.RetentionTime).First(); // If the start time is before the first scan, use the first scan

            ScanInfo endScan = ms1ScanInfos
                .Where(p => p.RetentionTime > end)
                .OrderBy(p => p.RetentionTime)
                .FirstOrDefault()
                ?? ms1ScanInfos.OrderBy(p => p.RetentionTime).Last(); // If the end time is after the last scan, use the last scan

            // Collect all peaks from the Ms1 scans in the given time window, then build the XIC
            List<IIndexedPeak> peaks = new List<IIndexedPeak>();
            for (int j = startScan.ZeroBasedScanIndex; j <= endScan.ZeroBasedScanIndex; j++)
            {
                double mz = id.PeakfindingMass.ToMz(id.PrecursorChargeState);
                IIndexedPeak peak = peakIndexingEngine.GetIndexedPeak(mz , j, isotopeTolerance);
                if (peak != null)
                {
                    peaks.Add(peak);
                }
            }

            XIC xIC = null;
            if (peaks.Count > 0)
            {
                xIC = new XIC(peaks, id.PeakfindingMass, spectraFile, isReference, ids);
            }
            return xIC;
        }

        /// <summary>
        /// Using the time window to track the sharedPeak across each run, and try to find the corresponding ChromatographicPeak.
        /// The ChromPeaks will be stored in chromPeaksInThisSequence
        /// For example, we have four invaild peak in the XIC, then we look at the first peak, generate chormPeak from each file. The dictionary will be like: peak1: [chromPeak (run 1), chromPeak (run 2), chromPeak (run 3)...]
        /// </summary>
        /// <param name="peaksInOneXIC"> The time window for the valid peak, format is <time < start, end>> </start></param>
        /// <param name="chromPeaksInThisSequence"> The resilt will store the inforamtion</param>
        internal void CollectChromPeakInRuns(PeakRegion sharedPeak, List<ChromatographicPeak> chromPeaksInSharedPeak, XICGroups xICGroups)
        {
            foreach (var xic in xICGroups)
            {
                double peakStart = sharedPeak.StartRT ;
                double PeakEnd = sharedPeak.EndRT ;
                bool isMBR = false;
                DetectionType detectionType = DetectionType.IsoTrack_MSMS;

                // Check is there any Id in the XIC within the time window.
                // Yes: use the Id from the same file. No: use the Id from other file, then set the detection type property as MBR.
                var idsForThisPeak = xICGroups.IdList
                    .Where(p=>Within(p.Ms2RetentionTimeInMinutes, peakStart, PeakEnd) && p.FileInfo.Equals(xic.SpectraFile))
                    .DistinctBy(p=>p.ModifiedSequence)
                    .ToList();
                switch (idsForThisPeak.Count)
                {
                    case 0: // If there is no any Id from the same file in the time window, then we borrow one from others files.
                        idsForThisPeak = xICGroups.IdList
                            .Where(p => Within(p.Ms2RetentionTimeInMinutes, peakStart, PeakEnd))
                            .DistinctBy(p=>p.ModifiedSequence)
                            .ToList();
                        detectionType = DetectionType.IsoTrack_MBR;
                        // If there are more than one Id from other files in the time window, then detectionType should be IsoTrack_Ambiguous.
                        if (idsForThisPeak.Count > 1) 
                        {
                            detectionType = DetectionType.IsoTrack_Ambiguous;
                        }

                        break;
                    case 1: // If there is one Id from the same file in the time window, then detectionType should be MSMS.
                        detectionType = DetectionType.IsoTrack_MSMS;
                        break;
                    case > 1: // If there are more than one Id from the same file in the time window, then detectionType should be IsoTrack_Ambiguous.
                        detectionType = DetectionType.IsoTrack_Ambiguous;
                        break;
                }

                // If there is no any existed Id in the time window, then skip this XIC.
                if (idsForThisPeak.Count == 0) { continue; }

                // Use the RT shift to calculate the realRT in this XIC, then searching the ChromatographicPeak.
                double timeShift = xic.RtShift;
                double rt = sharedPeak.ApexRT - timeShift;
                double start = sharedPeak.StartRT - timeShift;
                double end = sharedPeak.EndRT - timeShift;

                // Generate the practical time for searching. Time info: predicted RT, RtStartHypothesis, RtEndHypothesis
                Tuple<double, double, double> rtInfo = new Tuple<double, double, double>(rt, start, end); 

                ChromatographicPeak peak = FindChromPeak(rtInfo, xic, idsForThisPeak, detectionType);
                chromPeaksInSharedPeak.Add(peak);
            }
        }

        public bool Within(double time, double start, double end)
        {
            if (time >= start && time <= end)
            {
                return true;
            } 
            return false;
        }

        /// <summary>
        /// Using the RtInfo to find the ChromatographicPeak in the XIC.
        /// </summary>
        /// <param name="rtInfo">Retention information</param>
        /// <param name="xic">The searched XIC</param>
        /// <param name="idForChrom"></param>
        /// <param name="isMBR"></param>
        /// <returns></returns>
        internal ChromatographicPeak FindChromPeak(Tuple<double, double, double> rtInfo, XIC xic, List<Identification> idsForChrom, DetectionType detectionType) 
        {
            // Get the snippedPeaks from the window, then used for finding the isotopic envelope.
            List<IIndexedPeak> snippedPeaks = new ();
            Identification id = idsForChrom.FirstOrDefault();
            SpectraFileInfo spectraFile = xic.SpectraFile;

            foreach (var peak in xic.Ms1Peaks)
            {
                if (peak.RetentionTime >= rtInfo.Item2 && peak.RetentionTime <= rtInfo.Item3)
                {
                    snippedPeaks.Add(peak);
                }
            }

            List<IsotopicEnvelope> chargeEnvelopes = GetIsotopicEnvelopes(snippedPeaks, id, id.PrecursorChargeState, spectraFile).OrderBy(env => env.Intensity).ToList();
            if (chargeEnvelopes.Count == 0)
            {
                return null; // If there is no envelope, then return null.
            }

            // Find envelope with the highest intensity
            // envelopes = Find all envelopes that are adjacent to the highest intensity envelope (similar to PeakFind function)
            ChromatographicPeak acceptorPeak = null;
            if (idsForChrom.Count == 1)
            {
                acceptorPeak = new ChromatographicPeak(id, xic.SpectraFile, detectionType);
            }
            else
            {
                acceptorPeak = new ChromatographicPeak(idsForChrom, xic.SpectraFile, detectionType);
            }

            IsotopicEnvelope bestEnvelopes = chargeEnvelopes.OrderByDescending(p => p.Intensity).First();
            acceptorPeak.IsotopicEnvelopes.Add(bestEnvelopes);
            acceptorPeak.CalculateIntensityForThisFeature(FlashParams.Integrate);

            return acceptorPeak;
        }

        /// <summary>
        /// Add the isobaric peaks into the result dictionary (_peaks)
        /// </summary>
        internal void AddIsoPeaks()
        {
            foreach (var fileInfo in SpectraFileInfoList)
            {
                var allIsoTrackerPeaksInFile = IsobaricPeptideDict
                    .SelectMany(p => p.Value)
                    .SelectMany(p => p.Value)
                    .Where(peak => peak != null && peak.SpectraFileInfo.Equals(fileInfo))
                    .ToList();

                //remove the repeated peaks from FlashLFQ with the same identification list, the priority is IsoTrack > MSMS
                foreach (var peak in allIsoTrackerPeaksInFile)
                {
                    _results.Peaks[fileInfo].RemoveAll(p => IDsEqual(p.Identifications,peak.Identifications));
                }

                // Add the peaks into the result dictionary, and remove the duplicated peaks.
                // And we choice the IsoTrack_MSMS as the priority.
                _results.Peaks[fileInfo].AddRange(allIsoTrackerPeaksInFile);
                _results.Peaks[fileInfo] = _results.Peaks[fileInfo]
                   .GroupBy(peak => new { peak.ApexRetentionTime, peak.SpectraFileInfo, peak.Identifications.First().BaseSequence })
                   .Select(group =>
                   {
                       // Prioritize IsoTrack_MSMS over other detection types
                       var prioritizedPeak = group
                           .OrderByDescending(p => p.DetectionType == DetectionType.IsoTrack_MSMS)
                           .ThenByDescending(p => p.DetectionType == DetectionType.IsoTrack_MBR)
                            .ThenByDescending(p => p.DetectionType == DetectionType.IsoTrack_Ambiguous)
                           .ThenByDescending(p => p.DetectionType == DetectionType.MSMS)
                           .First();
                       return prioritizedPeak;
                   })
                   .ToList();
            }
        }

        /// <summary>
        /// Check if two lists of identifications are equal, to be compatible for MovId as well, we only consider modified sequence, and fileInfo.
        /// Only for chromatographicPeaks comparison.
        /// </summary>
        /// <param name="idList1"></param>
        /// <param name="idList2"></param>
        /// <returns></returns>
        private bool IDsEqual(List<Identification> idList1, List<Identification> idList2)
        {
            if (idList1.Count != idList2.Count)
            {
                return false;
            }

            var sortedIdList1 = idList1.OrderBy(id => id.ModifiedSequence).ToList();
            var sortedIdList2 = idList2.OrderBy(id => id.ModifiedSequence).ToList();

            for (int i = 0; i < sortedIdList1.Count; i++)
            {
                if (!sortedIdList1[i].ModifiedSequence.Equals(sortedIdList2[i].ModifiedSequence))
                {
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Calculates the target m/z (mass-to-charge) values based on isotopic distributions and modifications for the
        /// peptides in the current dataset.
        /// </summary>
        /// <returns>A list of doubles representing the calculated m/z values of interest.</returns>
        internal List<float> GetTargetMz()
        {
            // we only keep the masses from the PeptideListForIsoTracker.
            // Then we will calculate the masses of interest from the isotopic distributions of the unique peptides.
            IEnumerable<string> uniqueFullSequences = PeptideGroupsForIsoTracker
                .SelectMany(g => g.Identifications)
                .Select(p=>p.ModifiedSequence)
                .Distinct();
            HashSet<float> interestedMass = new HashSet<float>();

            // Collect all the masses of interest from the isotopic distributions of the peptides
            foreach (var fullSeq in uniqueFullSequences)
            {
                var monoMass = _allIdentifications.FirstOrDefault(id => id.ModifiedSequence == fullSeq)?.MonoisotopicMass;
                var peptideIsotopicDistribution= ModifiedSequenceToIsotopicDistribution[fullSeq];
                List<float> massesOfInterest = peptideIsotopicDistribution.Select(p => (float)(p.massShift + monoMass)).ToList();
                var minMass = massesOfInterest.Min();
                var maxMass = massesOfInterest.Max();
                massesOfInterest.Add(minMass - (float)Constants.C13MinusC12); // Except the isotopic distribution, we also add three mass (one before the lowest, two over the biggest)
                massesOfInterest.Add(maxMass + (float)Constants.C13MinusC12);
                massesOfInterest.Add(maxMass + 2f * (float)Constants.C13MinusC12);
                interestedMass.AddRange(massesOfInterest);
            }
            return GetAllPossibleMzs(interestedMass);
        }

        // Helper method for float version
        private List<float> GetAllPossibleMzs(HashSet<float> targetMasses)
        {
            HashSet<float> targetMzs = new();
            foreach (var mass in targetMasses)
            {
                _chargeStates.ForEach(chargeState =>
                {
                    targetMzs.Add(mass.ToMz(chargeState));
                });
            }
            var sortedTargetMzs = targetMzs.ToList();
            sortedTargetMzs.Sort();
            return sortedTargetMzs;
        }

        /// <summary>
        /// Groups the peptide identifications for IsoTracker based on their base sequence and monoisotopic mass.
        /// Each peptide is filtered by motif and modification, excluding decoys.
        /// </summary>
        /// <returns> the isobaricPeptideGroup</returns>
        private List<IsobaricPeptideGroup> GroupPeptideForIsoTracker()
        {
            // Filter IDs by motif and modification, exclude decoys
            var filteredIds = _allIdentifications
                .Where(p => FlashParams.IsoTrackerIdFilter.ContainsAcceptableModifiedResidue(p.ModifiedSequence))
                .Where(p => p.BaseSequence != p.ModifiedSequence && !p.IsDecoy);

            // Group by base sequence and rounded monoisotopic mass
            var grouped = filteredIds
                .GroupBy(p => (p.BaseSequence, MonoisotopicMassGroup: Math.Round(p.MonoisotopicMass / 0.001)));

            // Create IsobaricPeptideGroup objects for each group
            var isoGroups = grouped
                .Select(g => new IsobaricPeptideGroup(g.Key.BaseSequence, g.Key.MonoisotopicMassGroup,g.ToList()))
                .ToList();

            return isoGroups;
        }

    }

}