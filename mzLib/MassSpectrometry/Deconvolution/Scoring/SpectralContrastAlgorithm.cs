using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry.MzSpectra;
using MzLibUtil;

namespace MassSpectrometry.Deconvolution.Scoring
{
    public class SpectralContrastAlgorithm : ScoringAlgorithm
    {
        public SpectralContrastAlgorithm(PpmTolerance tolerance) : base(tolerance)
        {

        }

        public override double GetScore(IScoreArgs args)
        {
            switch (args)
            {
                case MinimalSpectraArgs spectraArgs:
                    SpectralSimilarity spectralSimilarity =
                        new(spectraArgs.ExperimentalSpectrum.MzArray, spectraArgs.ExperimentalSpectrum.IntensityArray,
                            spectraArgs.TheoreticalSpectrum.MzArray, spectraArgs.TheoreticalSpectrum.IntensityArray,
                            SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum, PpmTolerance.Value,
                            allPeaks: true, filterOutBelowThisMz: 1);
                    return spectralSimilarity.SpectralContrastAngle() ?? 0;
                default:
                    throw new ArgumentException();
            }
        }

    }
}
