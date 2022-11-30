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
    private double? _poorScore;

    public double PoorScore
    {
        get
        {
            if (_poorScore.HasValue) return (double)_poorScore;
            switch (ScoringMethod)
            {
                case ScoringMethods.KullbackLeibler:
                    _poorScore = Double.MaxValue;
                    return (double)_poorScore;
                case ScoringMethods.SpectralContrastAngle:
                    _poorScore = 0;
                    return (double)_poorScore;
                default:
                    _poorScore = Double.MinValue;
                    return (double)_poorScore;
            }
        }
    }

    public Scorer(ScoringMethods scoringMethod, PpmTolerance tolerance)
    {
        ScoringMethod = scoringMethod;
        ConstructScoringAlgorithm(tolerance);
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

    /// <summary>
    /// Compares two scores in a method specific fashion. Returns true if the instanceScore (first)
    /// is better than the argumentScore (second). Outputs the better of the two. This method is necessary
    /// because there are some metrics where lower scores are better.
    /// </summary>
    /// <param name="instanceScore"></param>
    /// <param name="argumentScore"></param>
    /// <param name="betterScore"></param>
    /// <returns></returns>
    /// <exception cref="NotImplementedException"></exception>
    public bool TestForScoreImprovement(double instanceScore, double argumentScore, out double betterScore)
    {
        switch (ScoringMethod)
        {
            case ScoringMethods.KullbackLeibler:
                if (instanceScore < argumentScore)
                {
                    betterScore = instanceScore;
                    return true;
                }
                else
                {
                    betterScore = argumentScore;
                    return false;
                }
            case ScoringMethods.SpectralContrastAngle:
                return DefaultCompare(instanceScore, argumentScore, out betterScore);
            default:
                return DefaultCompare(instanceScore, argumentScore, out betterScore);
        }
    }

    /// <summary>
    /// The default score comparison, where higher scores are better. Compares two scores, returns true
    /// if the instance score is higher than the argument score, returns false if instance score is lower.
    /// </summary>
    /// <param name="instanceScore"></param>
    /// <param name="argumentScore"></param>
    /// <param name="betterScore"> The higher of the two scores</param>
    /// <returns></returns>
    private bool DefaultCompare(double instanceScore, double argumentScore, out double betterScore)
    {
        if (instanceScore > argumentScore)
        {
            betterScore = instanceScore;
            return true;
        }
        else
        {
            betterScore = argumentScore;
            return false;
        }
    }

    private void ConstructScoringAlgorithm(PpmTolerance tolerance)
    {
        switch (ScoringMethod)
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