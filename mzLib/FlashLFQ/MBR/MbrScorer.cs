using Easy.Common.EasyComparer;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Data;
using System.Data.Entity.ModelConfiguration.Conventions;
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
        private readonly Normal _logIntensityDistribution;
        private readonly Normal _ppmDistribution;
        private readonly Normal _scanCountDistribution;
        private readonly Gamma _isotopicCorrelationDistribution;
        // The logFcDistributions and rtDifference distributions are unique to each donor file - acceptor file pair
        private Dictionary<SpectraFileInfo, Normal> _logFcDistributionDictionary;
        private Dictionary<SpectraFileInfo, Normal> _rtPredictionErrorDistributionDictionary;

        internal Dictionary<IIndexedMzPeak, ChromatographicPeak> ApexToAcceptorFilePeakDict { get; }
        internal List<ChromatographicPeak> UnambiguousMsMsAcceptorPeaks { get; }
        internal double MaxNumberOfScansObserved { get; }

        /// <summary>
        /// Takes in an intensity distribution, a log foldchange distribution, and a ppm distribution 
        /// unique to each donor file - acceptor file pair. These are used to score MBR matches
        /// </summary>
        internal MbrScorer(
            Dictionary<IIndexedMzPeak, ChromatographicPeak> apexToAcceptorFilePeakDict,
            List<ChromatographicPeak> acceptorFileMsmsPeaks,
            Normal ppmDistribution, 
            Normal logIntensityDistribution)
        {
            ApexToAcceptorFilePeakDict = apexToAcceptorFilePeakDict;
            UnambiguousMsMsAcceptorPeaks = acceptorFileMsmsPeaks.Where(p => p.Apex != null && p.DetectionType != DetectionType.MBR && p.NumIdentificationsByFullSeq == 1).ToList();
            MaxNumberOfScansObserved = acceptorFileMsmsPeaks.Max(peak => peak.ScanCount);
            _logIntensityDistribution = logIntensityDistribution;
            _ppmDistribution = ppmDistribution;
            _isotopicCorrelationDistribution = GetIsotopicEnvelopeCorrDistribution();
            _logFcDistributionDictionary = new();
            _rtPredictionErrorDistributionDictionary = new();

            // This is kludgey, because scan counts are discrete
            List<double> scanList = UnambiguousMsMsAcceptorPeaks.Select(peak => (double)peak.ScanCount).ToList();
            // build a normal distribution for the scan list of the acceptor peaks
            _scanCountDistribution = new Normal(scanList.Average(), scanList.Count > 30 ? scanList.StandardDeviation() : scanList.InterquartileRange() / 1.36);
        }

        /// <summary>
        /// This distribution represents (1 - Pearson Correlation) for isotopic envelopes of MS/MS acceptor peaks
        /// </summary>
        /// <returns></returns>
        private Gamma GetIsotopicEnvelopeCorrDistribution()
        {
            var pearsonCorrs = UnambiguousMsMsAcceptorPeaks.Select(p => 1 - p.IsotopicPearsonCorrelation).Where(p => p > 0).ToList();
            if (pearsonCorrs.Count <= 1) return null;
            double mean = pearsonCorrs.Mean();
            double variance = pearsonCorrs.Variance();
            var alpha = Math.Pow(mean, 2) / variance;
            var beta = mean / variance;
            return new Gamma(alpha, beta);
        }

        /// <summary>
        /// Takes in a list of retention time differences for anchor peptides (donor RT - acceptor RT) and uses
        /// this list to calculate the distribution of prediction errors of the local RT alignment strategy employed by
        /// match-between-runs for the specified donor file
        /// </summary>
        /// <param name="anchorPeptideRtDiffs">List of retention time differences (doubles) calculated as donor file RT - acceptor file RT</param>
        internal void AddRtPredErrorDistribution(SpectraFileInfo donorFile, List<double> anchorPeptideRtDiffs, int numberOfAnchorPeptides)
        {
            // in MBR, we use anchor peptides on either side of the donor to predict the retention time
            // here, we're going to repeat the same process, using neighboring anchor peptides to predicte the Rt shift for each
            // individual anchor peptide 
            // then, we'll check how close our predicted rt shift was to the observed rt shift
            // and build a distribution based on the predicted v actual rt diffs

            double cumSumRtDiffs;
            List<double> rtPredictionErrors = new();

            for (int i = numberOfAnchorPeptides; i < (anchorPeptideRtDiffs.Count - numberOfAnchorPeptides); i++)
            {
                cumSumRtDiffs = 0;
                for(int j = 1; j <= numberOfAnchorPeptides; j++)
                {
                    cumSumRtDiffs += anchorPeptideRtDiffs[i - j];
                    cumSumRtDiffs += anchorPeptideRtDiffs[i + j];
                }
                double avgDiff = cumSumRtDiffs / (2 * numberOfAnchorPeptides);
                rtPredictionErrors.Add(avgDiff - anchorPeptideRtDiffs[i]);
            }

            Normal rtPredictionErrorDist = new Normal(0, 0); 
            // Default distribution. Effectively assigns a RT Score of zero if no alignment can be performed
            // between the donor and acceptor based on shared MS/MS IDs

            if(rtPredictionErrors.Any())
            {
                double medianRtError = rtPredictionErrors.Median();
                double stdDevRtError = rtPredictionErrors.StandardDeviation();
                if(stdDevRtError >= 0.0 && !double.IsNaN(medianRtError))
                {
                    rtPredictionErrorDist = new Normal(medianRtError, 1);
                }
            }
            
            _rtPredictionErrorDistributionDictionary.Add(donorFile, rtPredictionErrorDist);
        }

        /// <summary>
        /// Takes in a list of retention time differences for anchor peptides (donor RT - acceptor RT) and uses
        /// this list to calculate the distribution of prediction errors of the local RT alignment strategy employed by
        /// match-between-runs for the specified donor file
        /// </summary>
        /// <returns> An MBR Score ranging between 0 and 100. Higher scores are better. </returns>
        internal double ScoreMbr(MbrChromatographicPeak acceptorPeak, ChromatographicPeak donorPeak, double predictedRt)
        {
            acceptorPeak.IntensityScore = CalculateIntensityScore(acceptorPeak.Intensity, donorPeak);
            acceptorPeak.RtPredictionError = predictedRt - acceptorPeak.ApexRetentionTime;
            acceptorPeak.RtScore = CalculateScore(_rtPredictionErrorDistributionDictionary[donorPeak.SpectraFileInfo],
                acceptorPeak.RtPredictionError);
            acceptorPeak.PpmScore = CalculateScore(_ppmDistribution, acceptorPeak.MassError);
            acceptorPeak.ScanCountScore = CalculateScore(_scanCountDistribution, acceptorPeak.ScanCount);
            acceptorPeak.IsotopicDistributionScore = CalculateScore(_isotopicCorrelationDistribution, 1 - acceptorPeak.IsotopicPearsonCorrelation);

            // Returns 100 times the geometric mean of the four scores (scan count, intensity score, rt score, ppm score)
            return 100 * Math.Pow(acceptorPeak.IntensityScore 
                * acceptorPeak.RtScore
                * acceptorPeak.PpmScore 
                * acceptorPeak.ScanCountScore
                * acceptorPeak.IsotopicDistributionScore, 0.20);
        }

        // Setting a minimum score prevents the MBR score from going to zero if one component of that score is 0
        // 3e-7 is the fraction of a normal distribution that lies at least 5 stdDev away from the mean
        private double _minScore = 3e-7;

        internal double CalculateScore(Normal distribution, double value)
        {
            // new method
            double absoluteDiffFromMean = Math.Abs(distribution.Mean - value);
            // Returns a value between (0, 1] where 1 means the value was equal to the distribution mean
            // The score represents the fraction of the distribution that lies absoluteDiffFromMean away from the mean or further
            // i.e., what fraction of the distribution is more extreme than value
            double score = 2 * distribution.CumulativeDistribution(distribution.Mean - absoluteDiffFromMean);
            return (double.IsNaN(score) || score == 0) ? _minScore : score;
        }

        internal double CalculateScore(Gamma distribution, double value)
        {
            if (value < 0 || distribution == null)
            {
                return _minScore;
            }

            // For the gamma distribtuion, the CDF is 0 when the pearson correlation is equal to 1 (value = 0)
            // The CDF then rapidly rises, reaching ~1 at a value of 0.3 (corresponding to a pearson correlation of 0.7)
            return 1 - distribution.CumulativeDistribution(value);
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
        /// Find the difference in peptide intensities between donor and acceptor files
        /// this intensity score creates a conservative bias in MBR
        /// </summary>
        /// <param name="idDonorPeaks"> List of peaks in the donoro file. </param>
        internal void CalculateFoldChangeBetweenFiles(List<ChromatographicPeak> idDonorPeaks)
        {

            var donorFileLogIntensities = idDonorPeaks.Where(p => p.Intensity > 0).Select(p => Math.Log(p.Intensity, 2)).ToList();
            double medianDonorLogIntensity = donorFileLogIntensities.Median();

            // Find the difference in peptide intensities between donor and acceptor files
            // this intensity score creates a conservative bias in MBR
            List<double> listOfFoldChangesBetweenTheFiles = new List<double>();
            var acceptorFileBestMsmsPeaks = new Dictionary<string, ChromatographicPeak>();

            // get the best (most intense) peak for each peptide in the acceptor file
            foreach (ChromatographicPeak acceptorPeak in UnambiguousMsMsAcceptorPeaks)
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
            Normal foldChangeDistribution = listOfFoldChangesBetweenTheFiles.Count > 100
                ? new Normal(listOfFoldChangesBetweenTheFiles.Median(), listOfFoldChangesBetweenTheFiles.StandardDeviation())
                : null;

            if (foldChangeDistribution != null)
            {
                _logFcDistributionDictionary.Add(idDonorPeaks.First().SpectraFileInfo, foldChangeDistribution);
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
                && _ppmDistribution != null
                && _scanCountDistribution != null
                && _logIntensityDistribution != null;
        }

    }
}
