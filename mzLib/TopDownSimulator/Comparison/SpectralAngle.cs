using System;
using System.Collections.Generic;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;

namespace TopDownSimulator.Comparison;

/// <summary>
/// Per-scan spectral-angle and cosine-style similarity utilities.
/// </summary>
public static class SpectralAngle
{
    public static double ComputeCosineSimilarity(IReadOnlyList<double> observed, IReadOnlyList<double> predicted)
    {
        if (observed.Count != predicted.Count)
            throw new ArgumentException("Observed and predicted vectors must have the same length.");

        double dot = 0;
        double obsNorm = 0;
        double predNorm = 0;
        for (int i = 0; i < observed.Count; i++)
        {
            dot += observed[i] * predicted[i];
            obsNorm += observed[i] * observed[i];
            predNorm += predicted[i] * predicted[i];
        }

        if (obsNorm <= 0 && predNorm <= 0)
            return 1;
        if (obsNorm <= 0 || predNorm <= 0)
            return 0;

        return dot / Math.Sqrt(obsNorm * predNorm);
    }

    public static double Compute(IReadOnlyList<double> observed, IReadOnlyList<double> predicted)
    {
        double cosine = Math.Clamp(ComputeCosineSimilarity(observed, predicted), -1.0, 1.0);
        return 1.0 - 2.0 * Math.Acos(cosine) / Math.PI;
    }

    public static double[] ComputePerScan(ProteoformGroundTruth truth, IReadOnlyList<ProteoformModel> proteoforms, double sigmaMz)
    {
        var projected = ModelProjection.ProjectIsotopologueIntensities(truth, proteoforms, sigmaMz);
        int nCharges = truth.ChargeCount;
        int nIso = truth.CentroidMzs[0].Length;
        int nScans = truth.ScanCount;
        int vectorLength = nCharges * nIso;
        var result = new double[nScans];

        for (int s = 0; s < nScans; s++)
        {
            var observed = new double[vectorLength];
            var predicted = new double[vectorLength];
            int k = 0;
            for (int c = 0; c < nCharges; c++)
            for (int i = 0; i < nIso; i++)
            {
                observed[k] = truth.IsotopologueIntensities[c][i][s];
                predicted[k] = projected[c][i][s];
                k++;
            }

            result[s] = Compute(observed, predicted);
        }

        return result;
    }
}
