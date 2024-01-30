using MathNet.Numerics.Distributions;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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

        /// <summary>
        /// Takes in an intensity distribution, a log foldchange distribution, and a ppm distribution 
        /// unique to each donor file - acceptor file pair. These are used to score MBR matches
        /// </summary>
        internal MbrScorer(Dictionary<IndexedMassSpectralPeak, ChromatographicPeak> apexToAcceptorFilePeakDict,
            Normal ppmDistribution, Normal logIntensityDistribution)
        {
            _logIntensityDistribution = intensityDistribution;
            _ppmDistribution = ppmDistribution;
            _logFcDistributionDictionary = new();
        }

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
                if (logIntensity < _logIntensityDistribution.Median)
                    intensityDensity = _logIntensityDistribution.Density(logIntensity);
                else
                    intensityDensity = _logIntensityDistribution.Density(_logIntensityDistribution.Mode);
            }

            double intensityScore = DensityScoreConversion(intensityDensity);
            double ppmScore = DensityScoreConversion(_ppmDistribution.Density(ppmError));
            double rtScore = DensityScoreConversion(rtDistribution.Density(retentionTime));

            return ppmScore + rtScore + intensityScore;
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
