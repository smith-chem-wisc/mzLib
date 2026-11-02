using Easy.Common.EasyComparer;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    /// <summary>
    /// This class takes in an intensity distribution, a log foldchange distribution, and a ppm distribution 
    /// unique to each donor file - acceptor file pair
    /// </summary>
    internal class MbrScorer
    {
        // Intensity and ppm distributions are specific to each acceptor file
        private Normal _logIntensityDistribution;
        private Normal _ppmDistribution;
        private Normal _scanCountDistribution;
        private Gamma _isotopicCorrelationDistribution;

        // Keep raw ppm stats to support exact tolerance calculation even when the fitted distribution must be made valid
        private double _ppmMedianRaw;
        private double _ppmStdDevRaw;

        // Set a floor on the minimum StdDev for RtErrorDistributions and add a default StdDev for when a donor file has no anchor peptides. This prevents degenerate distributions from being created and allows MBR scoring to proceed.
        public readonly double RtStandardDeviationMin = 0.0167; // 0.0167 minutes = 1 second
        public readonly double RtStandardDeviationDefault = 1;

        // The logFcDistributions and rtDifference distributions are unique to each donor file - acceptor file pair
        private readonly Dictionary<SpectraFileInfo, Normal> _logFcDistributionDictionary;
        private readonly Dictionary<SpectraFileInfo, Normal> _rtPredictionErrorDistributionDictionary;

        internal Dictionary<IIndexedPeak, ChromatographicPeak> ApexToAcceptorFilePeakDict { get; }
        internal List<ChromatographicPeak> UnambiguousMsMsAcceptorPeaks { get; }
        internal double MaxNumberOfScansObserved { get; }

        /// <summary>
        /// Takes in an intensity distribution, a log foldchange distribution, and a ppm distribution 
        /// unique to each donor file - acceptor file pair. These are used to score MBR matches
        /// </summary>
        internal MbrScorer(
            Dictionary<IIndexedPeak, ChromatographicPeak> apexToAcceptorFilePeakDict,
            List<ChromatographicPeak> unambiguousAcceptorFilePeaks)
        {
            ApexToAcceptorFilePeakDict = apexToAcceptorFilePeakDict ?? new Dictionary<IIndexedPeak, ChromatographicPeak>();
            UnambiguousMsMsAcceptorPeaks = unambiguousAcceptorFilePeaks ?? new List<ChromatographicPeak>();

            if (UnambiguousMsMsAcceptorPeaks.Count > 0)
            {
                MaxNumberOfScansObserved = UnambiguousMsMsAcceptorPeaks.Max(peak => peak.ScanCount);
            }
            else
            {
                MaxNumberOfScansObserved = 0;
            }

            _logFcDistributionDictionary = new();
            _rtPredictionErrorDistributionDictionary = new();

            _ppmMedianRaw = 0;
            _ppmStdDevRaw = 0;
        }

        /// <summary>
        /// Constructs the distributions that are used to score MBR matches
        /// </summary>
        /// <returns>Returns true if the scorer was initialized successfully, false otherwise</returns>
        internal bool InitializeScorer()
        {
            if (UnambiguousMsMsAcceptorPeaks == null || UnambiguousMsMsAcceptorPeaks.Count < 3)
                return false;

            // Populate distributions for scoring MBR matches
            _logIntensityDistribution = GetLogIntensityDistribution();
            _ppmDistribution = GetPpmErrorDistribution();
            _isotopicCorrelationDistribution = GetIsotopicEnvelopeCorrDistribution();
            _scanCountDistribution = GetScanCountDistribution();

            return IsValid();
        }

        #region BuildDistributions
        private Normal GetPpmErrorDistribution()
        {
            // Construct a distribution of ppm errors for all MSMS peaks in the acceptor file
            List<double> ppmErrors = UnambiguousMsMsAcceptorPeaks
                .Select(p => p.MassError)
                .Where(e => !double.IsNaN(e) && !double.IsInfinity(e))
                .ToList();

            if (ppmErrors.Count < 2)
                return null;

            _ppmMedianRaw = ppmErrors.Median();
            double spread = ppmErrors.Count > 30 ? ppmErrors.InterquartileRange() / 1.36 : ppmErrors.StandardDeviation();
            _ppmStdDevRaw = double.IsNaN(spread) || spread < 0 ? 0 : spread;

            // MathNet requires StdDev > 0. If the spread is zero (identical errors), use a safe stddev for the Normal instance,
            // but preserve the raw (zero) stddev for GetPpmErrorTolerance to return exactly median + 4*0 as per tests.
            double safeStdDev = _ppmStdDevRaw > 0 && Normal.IsValidParameterSet(_ppmMedianRaw, _ppmStdDevRaw) ? _ppmStdDevRaw : 1.0;

            return new Normal(_ppmMedianRaw, safeStdDev);
        }

        private Normal GetLogIntensityDistribution()
        {
            var logIntensities = UnambiguousMsMsAcceptorPeaks
                .Where(p => p.Intensity > 0 && !double.IsNaN(p.Intensity) && !double.IsInfinity(p.Intensity))
                .Select(p => Math.Log(p.Intensity, 2))
                .ToList();

            if (logIntensities.Count < 2)
                return null;

            double mean = logIntensities.Median();
            double stdDev = logIntensities.InterquartileRange() / 1.36;
            if (double.IsNaN(stdDev) || stdDev <= 0 || !Normal.IsValidParameterSet(mean, stdDev))
            {
                // Fall back to a safe, non-zero std-dev to avoid invalid Normal parameterization
                stdDev = 1.0;
            }

            return new Normal(mean, stdDev);
        }

        /// <summary>
        /// Find the difference in peptide intensities between donor and acceptor files
        /// this intensity score creates a conservative bias in MBR
        /// </summary>
        /// <param name="idDonorPeaks"> List of peaks in the donoro file. </param>
        internal void CalculateFoldChangeBetweenFiles(List<ChromatographicPeak> idDonorPeaks)
        {
            if (idDonorPeaks == null || idDonorPeaks.Count == 0)
                return;

            // Find the difference in peptide intensities between donor and acceptor files
            // this intensity score creates a conservative bias in MBR
            List<double> listOfFoldChangesBetweenTheFiles = new List<double>();
            var acceptorFileBestMsmsPeaks = new Dictionary<string, ChromatographicPeak>();

            // get the best (most intense) peak for each peptide in the acceptor file
            foreach (ChromatographicPeak acceptorPeak in UnambiguousMsMsAcceptorPeaks)
            {
                string key = acceptorPeak.Identifications.First().ModifiedSequence;
                if (acceptorFileBestMsmsPeaks.TryGetValue(key, out ChromatographicPeak currentBestPeak))
                {
                    if (currentBestPeak.Intensity > acceptorPeak.Intensity)
                    {
                        acceptorFileBestMsmsPeaks[key] = acceptorPeak;
                    }
                }
                else
                {
                    acceptorFileBestMsmsPeaks.Add(key, acceptorPeak);
                }
            }

            foreach (var donorPeak in idDonorPeaks)
            {
                double donorPeakIntensity = donorPeak.Intensity;
                if (acceptorFileBestMsmsPeaks.TryGetValue(donorPeak.Identifications.First().ModifiedSequence, out var acceptorPeak))
                {
                    double acceptorPeakIntensity = acceptorPeak.Intensity;

                    double intensityLogFoldChange = Math.Log(acceptorPeakIntensity, 2) - Math.Log(donorPeakIntensity, 2);

                    if (!double.IsNaN(intensityLogFoldChange) && !double.IsInfinity(intensityLogFoldChange))
                    {
                        listOfFoldChangesBetweenTheFiles.Add(intensityLogFoldChange);
                    }
                }
            }

            listOfFoldChangesBetweenTheFiles = listOfFoldChangesBetweenTheFiles
                .Where(d => !(double.IsNaN(d) | double.IsInfinity(d)))
                .ToList();

            if (listOfFoldChangesBetweenTheFiles.Count < 100)
                return;

            double medianFC = listOfFoldChangesBetweenTheFiles.Median();
            double stdDevFC = listOfFoldChangesBetweenTheFiles.StandardDeviation();

            if (double.IsNaN(medianFC) || double.IsNaN(stdDevFC) || !Normal.IsValidParameterSet(medianFC, stdDevFC))
                return;

            _logFcDistributionDictionary[idDonorPeaks.First().SpectraFileInfo] = new Normal(medianFC, stdDevFC);
        }


        // This is kludgey, because scan counts are discrete
        private Normal GetScanCountDistribution()
        {
            List<double> scanList = UnambiguousMsMsAcceptorPeaks.Select(peak => (double)peak.ScanCount).ToList();

            // Should be at least 3 peaks to be here; still handle degenerate cases robustly
            double mean = scanList.Count > 0 ? scanList.Average() : 0.0;
            double stdDev = scanList.Count > 30 ? scanList.StandardDeviation() : scanList.InterquartileRange() / 1.36;

            if (double.IsNaN(stdDev) || stdDev <= 0 || !Normal.IsValidParameterSet(mean, stdDev))
            {
                stdDev = 1.0;
            }

            return new Normal(mean, stdDev);
        }

        /// <summary>
        /// This distribution represents (1 - Pearson Correlation) for isotopic envelopes of MS/MS acceptor peaks
        /// </summary>
        /// <returns></returns>
        private Gamma GetIsotopicEnvelopeCorrDistribution()
        {
            var pearsonCorrs = UnambiguousMsMsAcceptorPeaks
                .Select(p => 1 - p.IsotopicPearsonCorrelation)
                .Where(p => p > 0 && !double.IsNaN(p) && !double.IsInfinity(p))
                .ToList();

            if (pearsonCorrs.Count < 2) return null;

            double mean = pearsonCorrs.Mean();
            double variance = pearsonCorrs.Variance();

            if (variance <= 0 || double.IsNaN(mean) || double.IsNaN(variance))
                return null;

            var alpha = Math.Pow(mean, 2) / variance;
            var beta = mean / variance;

            if (!Gamma.IsValidParameterSet(alpha, beta))
                return null;

            return new Gamma(alpha, beta);
        }

        /// <summary>
        /// Uses the anchor peptides shared between the donor and acceptor files to calculate the distribution
        /// of prediction errors of the local RT alignment strategy employed by match-between-runs for the
        /// specified donor file. The curve carries its own donor-RT ordering, which the local alignment below
        /// relies on, so callers do not have to sort before passing it in.
        /// </summary>
        /// <param name="calibrationCurve">The anchor peptides shared between the donor and acceptor files, ordered by donor apex retention time</param>
        /// <param name="numberOfAnchorPeptides">The number of neighboring anchor peptides used on either side to predict each anchor's RT shift</param>
        internal void AddRtPredErrorDistribution(SpectraFileInfo donorFile, RetentionTimeCalibrationCurve calibrationCurve, int numberOfAnchorPeptides)
        {
            // Default distribution: safe, non-degenerate. Also the distribution used when a donor file has too
            // few anchor peptides to estimate a prediction-error spread.
            Normal rtPredictionErrorDist = new Normal(0, RtStandardDeviationDefault);
            RetentionTimeCalibDataPoint[] validCalibrationDataPoints = calibrationCurve.DataPoints.Where(x => !double.IsNaN(x.RtDiff)).ToArray();

            // in MBR, we use anchor peptides on either side of the donor to predict the retention time
            // here, we're going to repeat the same process, using neighboring anchor peptides to predict the Rt shift for each
            // individual anchor peptide 
            // then, we'll check how close our predicted rt shift was to the observed rt shift
            // and build a distribution based on the predicted v actual rt diffs
            if (validCalibrationDataPoints != null && numberOfAnchorPeptides >= 0 && validCalibrationDataPoints.Length >= (2 * numberOfAnchorPeptides + 1))
            {
                List<double> rtPredictionErrors = new();
                List<double> rtDiffs = new(2*numberOfAnchorPeptides);

                for (int i = numberOfAnchorPeptides; i < (validCalibrationDataPoints.Length - numberOfAnchorPeptides); i++)
                {
                    rtDiffs.Clear();
                    for (int j = 1; j <= numberOfAnchorPeptides; j++)
                    {
                        rtDiffs.Add(validCalibrationDataPoints[i - j].RtDiff);
                        rtDiffs.Add(validCalibrationDataPoints[i + j].RtDiff);
                    }
                    double avgDiff = rtDiffs.Average();
                    double err = avgDiff - validCalibrationDataPoints[i].RtDiff;
                    if (!double.IsNaN(err) && !double.IsInfinity(err))
                    {
                        rtPredictionErrors.Add(err);
                    }
                }

                if (rtPredictionErrors.Count >= 2)
                {
                    double medianRtError = rtPredictionErrors.Median();
                    double stdDevRtError = rtPredictionErrors.InterquartileRange() / 1.36; // Use IQR to estimate stddev for robustness. IQR/1.36 is a robust estimator of stddev for normal distributions, and is less sensitive to outliers than the standard deviation.

                    if (!double.IsNaN(medianRtError))
                    {
                        if(stdDevRtError < RtStandardDeviationMin)
                            stdDevRtError = RtStandardDeviationMin;
                        if(Normal.IsValidParameterSet(medianRtError, stdDevRtError))
                           rtPredictionErrorDist = new Normal(medianRtError, stdDevRtError);
                    }
                }
            }

            _rtPredictionErrorDistributionDictionary[donorFile] = rtPredictionErrorDist;
        }

        #endregion

        /// <summary>
        /// Scores an MBR candidate peak.
        /// </summary>
        /// <returns> An MBR Score ranging between 0 and 100. Higher scores are better. </returns>
        internal double ScoreMbr(MbrChromatographicPeak acceptorPeak, ChromatographicPeak donorPeak, double predictedRt)
        {
            acceptorPeak.IntensityScore = CalculateIntensityScore(acceptorPeak.Intensity, donorPeak);
            acceptorPeak.RtPredictionError = predictedRt - acceptorPeak.ApexRetentionTime;

            // Use a safe default RT distribution if none was computed for this donor file
            if (!_rtPredictionErrorDistributionDictionary.TryGetValue(donorPeak.SpectraFileInfo, out var rtDist) || rtDist == null)
            {
                rtDist = new Normal(0, RtStandardDeviationDefault);
            }

            acceptorPeak.RtScore = CalculateScore(rtDist, acceptorPeak.RtPredictionError);
            acceptorPeak.PpmScore = CalculateScore(_ppmDistribution, acceptorPeak.MassError);
            acceptorPeak.ScanCountScore = CalculateScore(_scanCountDistribution, acceptorPeak.ScanCount);
            acceptorPeak.IsotopicDistributionScore = CalculateScore(_isotopicCorrelationDistribution, 1 - acceptorPeak.IsotopicPearsonCorrelation);

            // Returns 100 times the geometric mean of the five scores
            return 100 * Math.Pow(
                acceptorPeak.IntensityScore
                * acceptorPeak.RtScore
                * acceptorPeak.PpmScore
                * acceptorPeak.ScanCountScore
                * acceptorPeak.IsotopicDistributionScore, 0.20);
        }

        // Setting a minimum score prevents the MBR score from going to zero if one component of that score is 0
        // 3e-7 is the fraction of a normal distribution that lies at least 5 stdDev away from the mean
        private readonly double _minScore = 3e-7;

        internal double CalculateScore(Normal distribution, double value)
        {
            if (distribution == null || double.IsNaN(value) || double.IsInfinity(value))
            {
                return _minScore;
            }

            double absoluteDiffFromMean = Math.Abs(distribution.Mean - value);
            // Returns a value between (0, 1] where 1 means the value was equal to the distribution mean
            // The score represents the fraction of the distribution that lies absoluteDiffFromMean away from the mean or further
            // i.e., what fraction of the distribution is more extreme than value
            double score = 2 * distribution.CumulativeDistribution(distribution.Mean - absoluteDiffFromMean);
            return (double.IsNaN(score) || score == 0) ? _minScore : score;
        }

        internal double CalculateScore(Gamma distribution, double value)
        {
            if (value < 0 || distribution == null || double.IsNaN(value) || double.IsInfinity(value))
            {
                return _minScore;
            }

            // For the gamma distribtuion, the CDF is 0 when the pearson correlation is equal to 1 (value = 0)
            // The CDF then rapidly rises, reaching ~1 at a value of 0.3 (corresponding to a pearson correlation of 0.7)
            double cdf = distribution.CumulativeDistribution(value);
            if (double.IsNaN(cdf))
                return _minScore;

            return 1 - cdf;
        }

        internal double CalculateIntensityScore(double acceptorIntensity, ChromatographicPeak donorPeak)
        {
            if (donorPeak != null && acceptorIntensity != 0 && donorPeak.Intensity != 0 &&
                _logFcDistributionDictionary.TryGetValue(donorPeak.SpectraFileInfo, out var logFcDistribution))
            {
                var logFoldChange = Math.Log(acceptorIntensity, 2) - Math.Log(donorPeak.Intensity, 2);
                return CalculateScore(logFcDistribution, logFoldChange);
            }
            else
            {
                var logIntensity = Math.Log(acceptorIntensity, 2);
                return CalculateScore(_logIntensityDistribution, logIntensity);
            }

        }

        /// <summary>
        /// Determines whether or not the scorer is validly paramaterized and capable 
        /// of scoring MBR transfers originating from the given donorFile
        /// </summary>
        internal bool IsValid(SpectraFileInfo donorFile)
        {
            return _rtPredictionErrorDistributionDictionary.TryGetValue(donorFile, out var rtDist)
                && rtDist != null
                && IsValid();
        }

        /// <summary>
        /// This method checks whether the scorer is validly parameterized and capable of scoring MBR transfers
        /// Notably, it is indifferent to the isotopic correlation distribution being null, as a null isotopic distribution correlation
        /// results in all MBR transfers receiving the minimum score for the isotopic distribution component.
        /// This could be changed in the future, but currently multiple tests results in null isotopic distributions, and will break if they can't do MBR
        /// </summary>
        /// <returns></returns>
        internal bool IsValid()
        {
            return _ppmDistribution != null
                && _scanCountDistribution != null
                && _logIntensityDistribution != null;
        }

        /// <summary>
        /// Returns the standard deviation of the Ppm error distribution + the median of the Ppm error distribution.
        /// Uses raw ppm stats to preserve exact behavior for degenerate distributions (e.g., stddev = 0).
        /// </summary>
        internal double GetPpmErrorTolerance()
        {
            return Math.Abs(_ppmMedianRaw) + 4 * _ppmStdDevRaw;
        }
    }
}
