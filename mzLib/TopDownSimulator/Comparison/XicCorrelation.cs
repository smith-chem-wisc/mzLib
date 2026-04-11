using System;
using System.Collections.Generic;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;

namespace TopDownSimulator.Comparison;

public sealed record ChargeCorrelation(int Charge, double Correlation);

/// <summary>
/// Per-charge Pearson correlation between observed and simulated XICs.
/// </summary>
public static class XicCorrelation
{
    public static ChargeCorrelation[] ComputePerCharge(ProteoformGroundTruth truth, IReadOnlyList<ProteoformModel> proteoforms, double sigmaMz)
    {
        var projected = ModelProjection.ProjectIsotopologueIntensities(truth, proteoforms, sigmaMz);
        var projectedXics = ModelProjection.SumIsotopologuesByCharge(projected);
        var correlations = new ChargeCorrelation[truth.ChargeCount];

        for (int c = 0; c < truth.ChargeCount; c++)
        {
            correlations[c] = new ChargeCorrelation(
                Charge: truth.MinCharge + c,
                Correlation: ComputePearson(truth.ChargeXics[c], projectedXics[c]));
        }

        return correlations;
    }

    public static double ComputePearson(IReadOnlyList<double> observed, IReadOnlyList<double> predicted)
    {
        if (observed.Count != predicted.Count)
            throw new ArgumentException("Observed and predicted vectors must have the same length.");

        int n = observed.Count;
        if (n == 0)
            return 0;

        double meanObs = 0;
        double meanPred = 0;
        for (int i = 0; i < n; i++)
        {
            meanObs += observed[i];
            meanPred += predicted[i];
        }
        meanObs /= n;
        meanPred /= n;

        double num = 0;
        double denObs = 0;
        double denPred = 0;
        for (int i = 0; i < n; i++)
        {
            double obs = observed[i] - meanObs;
            double pred = predicted[i] - meanPred;
            num += obs * pred;
            denObs += obs * obs;
            denPred += pred * pred;
        }

        if (denObs <= 0 || denPred <= 0)
            return 0;

        return num / Math.Sqrt(denObs * denPred);
    }
}
