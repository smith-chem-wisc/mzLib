using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    internal static class MbrScorerFactory
    {
        /// <summary>
        /// Constructs a MbrScorer object that is used to score all MBR peaks for a given acceptor file
        /// </summary>
        /// <param name="acceptorFileMsmsPeaks"> All MSMS identified peaks in the acceptor file that contain quantifiable peptides </param>
        /// <param name="fileSpecificMbrTolerance">A ppm tolerance specific to the given file</param>
        /// <returns> A MbrScorer object </returns>
        public static MbrScorer BuildMbrScorer(List<ChromatographicPeak> acceptorFileMsmsPeaks,
            FlashLfqParameters flashParams, out PpmTolerance fileSpecificMbrTolerance)
        {
            // Construct a dictionary linking each MSMS peaks to the indexed peak of its apex.
            // This is to ensure MBR doesn't assign a peptide to a peak that is already claimed by one or more other peptides.
            var apexToAcceptorFilePeakDict = acceptorFileMsmsPeaks
                .Where(p => p.Apex != null)
                .DistinctBy(p => p.Apex.IndexedPeak)
                .ToDictionary(p => p.Apex.IndexedPeak, p => p);

            // Then, filter down so that only unambiguous peaks are used to construct scoring distributions
            List<ChromatographicPeak> unambiguousAcceptorFilePeaks = apexToAcceptorFilePeakDict
                .Select(kvp => kvp.Value)
                .Where(p => p.NumIdentificationsByFullSeq == 1 && p.DetectionType != DetectionType.MBR && p.DetectionType != DetectionType.IsoTrack_MBR)
                .ToList();

            // If we don't have enough unambiguous peaks, return null
            if (unambiguousAcceptorFilePeaks.Count < 3)
            {
                fileSpecificMbrTolerance = null;
                return null;
            }

            MbrScorer scorer = new MbrScorer(apexToAcceptorFilePeakDict, acceptorFileMsmsPeaks);
            double mbrPpmTolerance = Math.Min(scorer.GetPpmErrorTolerance(), flashParams.MbrPpmTolerance);
            fileSpecificMbrTolerance = new PpmTolerance(mbrPpmTolerance); // match between runs PPM tolerance

            return scorer;
        }

    }
}
