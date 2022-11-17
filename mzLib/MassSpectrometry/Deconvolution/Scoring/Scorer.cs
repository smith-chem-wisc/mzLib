using System;
using MzLibUtil;

namespace MassSpectrometry.Deconvolution.Scoring;

// Context class for scoring deconvolution hypotheses
public class Scorer
{
    public enum ScoringMethods
    {
        KullbackLeibler,
        SpectralContrastAngle
    }
    public ScoringAlgorithm ScoringAlgorithm { get; private set; }
    public ScoringMethods ScoringMethod { get; }

    public Scorer(ScoringMethods scoringMethod, PpmTolerance tolerance)
    {
        ScoringMethod = scoringMethod;
        ConstructScoringAlgorithm(scoringMethod, tolerance);
    }

    public double Score(IScoreArgs args)
    {
        return ScoringAlgorithm.GetScore(args);
    }

    public double Score(MinimalSpectrum experimentalSpectrum, MinimalSpectrum theoreticalSpectrum)
    {
        IScoreArgs args = new MinimalSpectraArgs(experimentalSpectrum, theoreticalSpectrum);
        return ScoringAlgorithm.GetScore(args);
    }

    private void ConstructScoringAlgorithm(ScoringMethods method, PpmTolerance tolerance)
    {
        switch (method)
        {
            case ScoringMethods.KullbackLeibler:
                throw new NotImplementedException();
            case ScoringMethods.SpectralContrastAngle:
                ScoringAlgorithm = new SpectralContrastAlgorithm(tolerance);
                break;
            default:
                throw new NotImplementedException();
        }
    }

}