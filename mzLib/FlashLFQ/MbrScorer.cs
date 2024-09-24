using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;

namespace FlashLFQ
{
    /// <summary>
    /// This class takes in an intensity distribution, a log foldchange distribution, and a ppm distribution 
    /// unique to each donor file - acceptor file pair
    /// </summary>
    internal class MbrScorer
    {
        internal SpectraFileInfo AcceptorFile { get; init; }
        // Intensity and ppm distributions are specific to each acceptor file
        private readonly Normal _logIntensityDistribution;
        private readonly Normal _ppmDistribution;
        private readonly Normal _scanCountDistribution;
        // The logFcDistributions and rtDifference distributions are unique to each donor file - acceptor file pair
        private Dictionary<SpectraFileInfo, Normal> _logFcDistributionDictionary;
        private Dictionary<SpectraFileInfo, Normal> _rtPredictionErrorDistributionDictionary;
        internal Dictionary<IndexedMassSpectralPeak, ChromatographicPeak> ApexToAcceptorFilePeakDict { get; }
        internal List<ChromatographicPeak> UnambiguousMsMsPeaks { get; }

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
            AcceptorFile = acceptorPeaks.First().SpectraFileInfo;
            ApexToAcceptorFilePeakDict = apexToAcceptorFilePeakDict;
            UnambiguousMsMsPeaks = acceptorPeaks.Where(p => p.Apex != null && !p.IsMbrPeak && p.NumIdentificationsByFullSeq == 1).ToList();
            _logIntensityDistribution = logIntensityDistribution;
            _ppmDistribution = ppmDistribution;
            _logFcDistributionDictionary = new();
            _rtPredictionErrorDistributionDictionary = new();

            // This is kludgey, because scan counts are discrete
            List<double> scanList = acceptorPeaks.Select(peak => (double)peak.ScanCount).ToList();
            // build a normal distribution for the scan list of the acceptor peaks
            // InterQuartileRange / 1.35 = StandardDeviation for a normal distribution
            _scanCountDistribution = new Normal(scanList.Average(), scanList.Count > 30 ? scanList.StandardDeviation() : scanList.InterquartileRange() / 1.36);
        }

        /// <summary>
        /// Takes in a list of retention time differences for anchor peptides (donor RT - acceptor RT) and uses
        /// this list to calculate the distribution of prediction errors of the local RT alignment strategy employed by
        /// match-between-runs for the specified donor file
        /// </summary>
        /// <param name="anchorPeptideRtDiffs">List of retention time differences (doubles) calculated as donor file RT - acceptor file RT</param>
        internal void AddRtPredErrorDistribution(SpectraFileInfo donorFile, List<double> anchorPeptideRtDiffs, int numberOfAnchorPeptidesPerSide)
        {
            // in MBR, we use anchor peptides on either side of the donor to predict the retention time
            // here, we're going to repeat the same process, using neighboring anchor peptides to predicte the Rt shift for each
            // individual anchor peptide 
            // then, we'll check how close our predicted rt shift was to the observed rt shift
            // and build a distribution based on the predicted v actual rt diffs

            double cumSumRtDiffs;
            List<double> rtPredictionErrors = new();

            for (int i = numberOfAnchorPeptidesPerSide; i < anchorPeptideRtDiffs.Count - (numberOfAnchorPeptidesPerSide) ; i++)
            {
                cumSumRtDiffs = 0;
                for(int j = 1; j <= numberOfAnchorPeptidesPerSide; j++)
                {
                    cumSumRtDiffs += anchorPeptideRtDiffs[i - j];
                    cumSumRtDiffs += anchorPeptideRtDiffs[i + j];
                }
                double avgDiff = cumSumRtDiffs / (2 * numberOfAnchorPeptidesPerSide);
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

        private double CalculateScore(Normal distribution, double value)
        {
            double absoluteDiffFromMean = Math.Abs(distribution.Mean - value);
            // Returns a value between (0, 1] where 1 means the value was equal to the distribution mean
            return 2 * distribution.CumulativeDistribution(distribution.Mean - absoluteDiffFromMean);
        }

        internal double CalculateIntensityScore(ChromatographicPeak acceptorPeak, ChromatographicPeak donorPeak)
        {
            double acceptorIntensity = acceptorPeak.Intensity;
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
        /// Calculates the retention time score for a given MbrAcceptor by comparing to the 
        /// distribution of all retention time prediction errors for all anchor peptides shared between 
        /// the donor and acceptor files
        /// </summary>
        /// <returns>Score bounded by 0 and 1, where higher scores are better</returns>
        internal double CalculateRetentionTimeScore(ChromatographicPeak acceptorPeak, ChromatographicPeak donorPeak)
        {
            double rtPredictionError = acceptorPeak.PredictedRetentionTime - acceptorPeak.ApexRetentionTime;
            return CalculateScore(_rtPredictionErrorDistributionDictionary[donorPeak.SpectraFileInfo], rtPredictionError);
        }

        /// <summary>
        /// Calculates the Ppm error score for a given acceptor by comparing the ppm error for the given peak
        /// to the ppm error of all non-MBR peaks in the acceptor file
        /// </summary>
        /// <returns>Score bounded by 0 and 1, where higher scores are better</returns>
        internal double CalculatePpmErrorScore(ChromatographicPeak acceptorPeak)
        {
            return CalculateScore(_ppmDistribution, acceptorPeak.MassError);
        }

        /// <summary>
        /// Calculates the scan count score for a given acceptor by comparing the number of scans observed for the given peak
        /// to the ppm error of all non-MBR peaks in the acceptor file
        /// </summary>
        /// <returns>Score bounded by 0 and 1, where higher scores are better</returns>
        internal double CalculateScanCountScore(ChromatographicPeak acceptorPeak)
        {
            return CalculateScore(_scanCountDistribution, acceptorPeak.ScanCount);
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
            foreach (ChromatographicPeak acceptorPeak in UnambiguousMsMsPeaks)
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
