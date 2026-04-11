using System.Collections.Generic;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;

namespace TopDownSimulator.Comparison;

public sealed record ResidualSummary(double ResidualEnergyFraction, double ObservedEnergy, double ResidualEnergy);

/// <summary>
/// Quantifies the model-unexplained energy in the extracted tensor.
/// </summary>
public static class ResidualAnalyzer
{
    public static ResidualSummary Analyze(ProteoformGroundTruth truth, IReadOnlyList<ProteoformModel> proteoforms, double sigmaMz)
    {
        var projected = ModelProjection.ProjectIsotopologueIntensities(truth, proteoforms, sigmaMz);
        double observedEnergy = 0;
        double residualEnergy = 0;

        for (int c = 0; c < truth.ChargeCount; c++)
        for (int i = 0; i < truth.CentroidMzs[c].Length; i++)
        for (int s = 0; s < truth.ScanCount; s++)
        {
            double observed = truth.IsotopologueIntensities[c][i][s];
            double residual = observed - projected[c][i][s];
            observedEnergy += observed * observed;
            residualEnergy += residual * residual;
        }

        double fraction = observedEnergy > 0 ? residualEnergy / observedEnergy : 0;
        return new ResidualSummary(fraction, observedEnergy, residualEnergy);
    }
}
