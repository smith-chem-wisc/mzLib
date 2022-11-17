using System;

namespace MassSpectrometry.Deconvolution.Scoring;

// Context class for scoring deconvolution hypotheses
public class Scorer
{
    public enum ScoringMethods
    {
        KullbackLeibler,
        SpectralContrastAngle
    }
    public ScoringAlgorithm ScoringAlgorithm { get; }
    public ScoringMethods ScoringMethod { get; }

    public Scorer(ScoringMethods scoringMethod)
    {
        ScoringMethod = scoringMethod;
        ScoringAlgorithm = ConstructScoringAlgorithm(scoringMethod);
    }

    public static ScoringAlgorithm ConstructScoringAlgorithm(ScoringMethods method)
    {
        switch (method)
        {
            case ScoringMethods.KullbackLeibler:
                throw new NotImplementedException();
            case ScoringMethods.SpectralContrastAngle:
                throw new NotImplementedException();
            default:
                throw new NotImplementedException();
        }

        return null;
    }

}