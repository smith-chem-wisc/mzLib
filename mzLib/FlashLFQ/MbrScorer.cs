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
        // Intensity and ppm distribution are specific to each acceptor file
        private readonly Normal _logIntensityDistribution;
        private readonly Normal _ppmDistribution;
        private readonly Normal _scanCountDistribution;
        // The logFcDistributions are unique to each donor file - acceptor file pair
        private Dictionary<SpectraFileInfo, Normal> _logFcDistributionDictionary;

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
            // This is kludgey, because scan counts are discrete
            List<double> scanList = acceptorPeaks.Select(peak => (double)peak.ScanCount).ToList();
            // build a normal distribution for the scan list of the acceptor peaks
            _scanCountDistribution = new Normal(scanList.Average(), scanList.Count > 30 ? scanList.StandardDeviation() : scanList.InterquartileRange() / 1.36);
        }

        /// <summary>
        /// Scores a MBR peak based on it's retention time, ppm error, and intensity
        /// </summary>
        /// <returns> The MBR score as a double. Higher scores are better. </returns>
        internal double ScoreMbr(ChromatographicPeak acceptorPeak, Normal rtDistribution, ChromatographicPeak? donorPeak = null)
        {
            acceptorPeak.IntensityScore = CalculateIntensityScore(acceptorPeak.Intensity, donorPeak);
            acceptorPeak.RtScore = CalculateScore(rtDistribution, acceptorPeak.ApexRetentionTime);
            acceptorPeak.PpmScore = CalculateScore(_ppmDistribution, acceptorPeak.MassError);
            acceptorPeak.ScanCountScore = CalculateScore(_scanCountDistribution, acceptorPeak.ScanCount);
            //acceptorPeak.ScanCountScore = (double)acceptorPeak.ScanCount / _scanCountDistribution.Median;

            double donorIdPEP = donorPeak.Identifications.OrderBy(p => p.PosteriorErrorProbability).First().PosteriorErrorProbability;

            return (acceptorPeak.IntensityScore + acceptorPeak.RtScore + acceptorPeak.PpmScore + acceptorPeak.ScanCountScore) * (1 - donorIdPEP);
        }

        internal double CalculateScore(Normal distribution, double value)
        {
            // new method
            double absoluteDiffFromMean = Math.Abs(distribution.Mean - value);
            // Returns a value between (0, 1] where 1 means the value was equal to the distribution mean
            return 2 * distribution.CumulativeDistribution(distribution.Mean - absoluteDiffFromMean);

            // old method
            //return DensityScoreConversion(distribution.Density(value));
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


                // I don't know what the if/else statement accomplishes. It feels like we should take the density regardless
                // As it is, the score is artifically inflated for very intense peaks
                //if (logIntensity < _logIntensityDistribution.Median)
                //    intensityDensity = _logIntensityDistribution.Density(logIntensity);
                //else
                //    intensityDensity = _logIntensityDistribution.Density(_logIntensityDistribution.Mode);

                //alternate, more straightforward approach
            }

        }

        internal double CalculatePpmScore(double ppmError)
        {
            return DensityScoreConversion(_ppmDistribution.Density(ppmError));
        }

        internal double CalculateRtScore(double retentionTime, Normal rtDistribution)
        {
            return DensityScoreConversion(rtDistribution.Density(retentionTime));
        }

        internal double CalculateScanCountScore(int scanCount)
        {
            return (double)scanCount / (double)MaxNumberOfScansObserved;
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
        /// Takes in the density of a normal distribution at a given point, and transforms it
        /// by taking the log of the density plus the square root of the squared density plus one
        /// This transformation was implemented in the original code, and we're unsure of the rationale
        /// </summary>
        /// <param name="density"> A Normal distribution</param>
        /// <returns> The transformed score</returns>
        private double DensityScoreConversion(double density)
        {
            return Math.Log(density + Math.Sqrt(Math.Pow(density, 2) + 1));
        }
      
    }
}
