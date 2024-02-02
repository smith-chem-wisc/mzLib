using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
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
        // Intensity and ppm distribution are specific to each acceptor file
        private readonly Normal _logIntensityDistribution;
        private readonly Normal _ppmDistribution;
        // The logFcDistributions are unique to each donor file - acceptor file pair
        private Dictionary<SpectraFileInfo, Normal> _logFcDistributionDictionary;

        internal Dictionary<IndexedMassSpectralPeak, ChromatographicPeak> ApexToAcceptorFilePeakDict { get; }
        internal List<ChromatographicPeak> UnambiguousMsMsAcceptorPeaks { get; }


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
            _logIntensityDistribution = logIntensityDistribution;
            _ppmDistribution = ppmDistribution;
            _logFcDistributionDictionary = new();
        }

        /// <summary>
        /// Scores a MBR peak based on it's retention time, ppm error, and intensity
        /// </summary>
        /// <returns> The MBR score as a double. Higher scores are better. </returns>
        internal double ScoreMbr(Normal rtDistribution, double retentionTime, double ppmError, double acceptorIntensity, ChromatographicPeak? donorPeak = null)
        {
            double intensityDensity;
            if (donorPeak != null && acceptorIntensity != 0 && donorPeak.Intensity != 0 &&
                _logFcDistributionDictionary.TryGetValue(donorPeak.SpectraFileInfo, out var logFcDistribution))
            {
                intensityDensity = logFcDistribution.Density(
                    Math.Log(acceptorIntensity, 2) - Math.Log(donorPeak.Intensity, 2)
                    );
            }
            else
            {
                var logIntensity = Math.Log(acceptorIntensity, 2);

                // I don't know what the if/else statement accomplishes. It feels like we should take the density regardless
                // As it is, the score is artifically inflated for very intense peaks
                //if (logIntensity < _logIntensityDistribution.Median)
                //    intensityDensity = _logIntensityDistribution.Density(logIntensity);
                //else
                //    intensityDensity = _logIntensityDistribution.Density(_logIntensityDistribution.Mode);

                //alternate, more straightforward approach
                intensityDensity = _logIntensityDistribution.Density(logIntensity);
            }

            double intensityScore = DensityScoreConversion(intensityDensity);
            double ppmScore = DensityScoreConversion(_ppmDistribution.Density(ppmError));
            double rtScore = DensityScoreConversion(rtDistribution.Density(retentionTime));

            double donorIdPEP = donorPeak.Identifications.OrderBy(p => p.PosteriorErrorProbability).First().PosteriorErrorProbability;

            return (ppmScore + rtScore + intensityScore) * (1 - donorIdPEP);
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
