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
        // The logFcDistributions and rtDifference distributions are unique to each donor file - acceptor file pair
        private Dictionary<SpectraFileInfo, Normal> _logFcDistributionDictionary;
        private Dictionary<SpectraFileInfo, Normal> _rtPredictionErrorDistributionDictionary;

        internal Dictionary<IndexedMassSpectralPeak, ChromatographicPeak> ApexToAcceptorFilePeakDict { get; }
        internal List<ChromatographicPeak> UnambiguousMsMsAcceptorPeaks { get; }
        internal double MaxNumberOfScansObserved { get; }

        /// <summary>
        /// Takes in an intensity distribution, a log foldchange distribution, and a ppm distribution 
        /// unique to each donor file - acceptor file pair. These are used to score MBR matches
        /// </summary>
        internal MbrScorer(
            Dictionary<IndexedMassSpectralPeak, ChromatographicPeak> apexToAcceptorFilePeakDict,
            List<ChromatographicPeak> acceptorPeaks,
            Normal ppmDistribution, 
            Normal logIntensityDistribution)
        {
            ApexToAcceptorFilePeakDict = apexToAcceptorFilePeakDict;
            UnambiguousMsMsAcceptorPeaks = acceptorPeaks.Where(p => p.Apex != null && !p.IsMbrPeak && p.NumIdentificationsByFullSeq == 1).ToList();
            MaxNumberOfScansObserved = acceptorPeaks.Max(peak => peak.ScanCount);
            _logIntensityDistribution = logIntensityDistribution;
            _ppmDistribution = ppmDistribution;
            _logFcDistributionDictionary = new();
            _rtPredictionErrorDistributionDictionary = new();

            // This is kludgey, because scan counts are discrete
            List<double> scanList = acceptorPeaks.Select(peak => (double)peak.ScanCount).ToList();
            // build a normal distribution for the scan list of the acceptor peaks
            _scanCountDistribution = new Normal(scanList.Average(), scanList.Count > 30 ? scanList.StandardDeviation() : scanList.InterquartileRange() / 1.36);
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

            for (int i = numberOfAnchorPeptides; i < anchorPeptideRtDiffs.Count - (numberOfAnchorPeptides); i++)
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

            if(!rtPredictionErrors.Any())
            {
                _rtPredictionErrorDistributionDictionary.Add(donorFile, new Normal(0, 1));
                return;
            }

            double medianRtError = rtPredictionErrors.Median();
            double stdDevRtError = rtPredictionErrors.StandardDeviation();

            _rtPredictionErrorDistributionDictionary.Add(donorFile, new Normal(medianRtError, stdDevRtError));
        }

        /// <summary>
        /// Get the RT window width for a given donor file,
        /// where RT window width is equal to 4*stdDev of the rtDiffs for all anchor peptides
        /// </summary>
        /// <returns>The width of the retention time window in minutes</returns>
        internal double GetRTWindowWidth(SpectraFileInfo donorFile)
        {
            // 95% of all peaks are expected to fall within six standard deviations
            return _rtPredictionErrorDistributionDictionary[donorFile].StdDev * 4;
        }

        internal double GetMedianRtDiff(SpectraFileInfo donorFile)
        {
            return _rtPredictionErrorDistributionDictionary[donorFile].Median;
        }

        /// <summary>
        /// Scores a MBR peak based on it's retention time, ppm error, and intensity
        /// </summary>
        /// <returns> An MBR Score ranging between 0 and 100. Higher scores are better. </returns>
        internal double ScoreMbr(ChromatographicPeak acceptorPeak, ChromatographicPeak donorPeak, double predictedRt)
        {
            acceptorPeak.IntensityScore = CalculateIntensityScore(acceptorPeak.Intensity, donorPeak);
            acceptorPeak.RtScore = CalculateScore(_rtPredictionErrorDistributionDictionary[donorPeak.SpectraFileInfo],
                predictedRt - acceptorPeak.ApexRetentionTime);
            acceptorPeak.PpmScore = CalculateScore(_ppmDistribution, acceptorPeak.MassError);
            acceptorPeak.ScanCountScore = CalculateScore(_scanCountDistribution, acceptorPeak.ScanCount);
            
            // Returns 100 times the geometric mean of the three scores (scan count, intensity score, rt score)
            return 100 * Math.Pow(acceptorPeak.IntensityScore * acceptorPeak.IntensityScore 
                * acceptorPeak.RtScore * acceptorPeak.RtScore 
                * acceptorPeak.PpmScore * acceptorPeak.ScanCountScore, 1.0/6.0);
        }

        internal double CalculateScore(Normal distribution, double value)
        {
            // new method
            double absoluteDiffFromMean = Math.Abs(distribution.Mean - value);
            // Returns a value between (0, 1] where 1 means the value was equal to the distribution mean
            return 2 * distribution.CumulativeDistribution(distribution.Mean - absoluteDiffFromMean);
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
      
    }
}
