using System.Collections.Generic;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;

namespace TopDownSimulator.Comparison;

internal static class ModelProjection
{
    public static double[][][] ProjectIsotopologueIntensities(
        ProteoformGroundTruth truth,
        IReadOnlyList<ProteoformModel> proteoforms,
        double sigmaMz)
    {
        var forwardModel = new ForwardModel(proteoforms, truth.MinCharge, truth.MaxCharge, sigmaMz);
        int nCharges = truth.ChargeCount;
        int nIso = truth.CentroidMzs[0].Length;
        int nScans = truth.ScanCount;

        var projected = new double[nCharges][][];
        for (int c = 0; c < nCharges; c++)
        {
            projected[c] = new double[nIso][];
            for (int i = 0; i < nIso; i++)
            {
                projected[c][i] = new double[nScans];
                double mz = truth.CentroidMzs[c][i];
                for (int s = 0; s < nScans; s++)
                    projected[c][i][s] = forwardModel.Evaluate(truth.ScanTimes[s], mz);
            }
        }

        return projected;
    }

    public static double[][] SumIsotopologuesByCharge(double[][][] isotopologueTensor)
    {
        int nCharges = isotopologueTensor.Length;
        int nIso = isotopologueTensor[0].Length;
        int nScans = isotopologueTensor[0][0].Length;
        var xics = new double[nCharges][];

        for (int c = 0; c < nCharges; c++)
        {
            xics[c] = new double[nScans];
            for (int i = 0; i < nIso; i++)
            for (int s = 0; s < nScans; s++)
                xics[c][s] += isotopologueTensor[c][i][s];
        }

        return xics;
    }
}
