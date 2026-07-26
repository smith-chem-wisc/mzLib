using System;
using System.Collections.Generic;
using System.Linq;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;

namespace TopDownSimulator.Comparison;

public sealed record ComparisonReport(
    double[] PerScanSpectralAngles,
    ChargeCorrelation[] PerChargeXicCorrelations,
    ResidualSummary Residuals)
{
    public double MeanSpectralAngle => PerScanSpectralAngles.Length > 0 ? PerScanSpectralAngles.Average() : 0;
    public double MeanChargeCorrelation => PerChargeXicCorrelations.Length > 0 ? PerChargeXicCorrelations.Average(p => p.Correlation) : 0;
}

public static class ComparisonReportBuilder
{
    public static ComparisonReport Create(ProteoformGroundTruth truth, IReadOnlyList<ProteoformModel> proteoforms, double sigmaMz) =>
        Create(truth, proteoforms, new ConstantPeakWidth(sigmaMz));

    public static ComparisonReport Create(ProteoformGroundTruth truth, IReadOnlyList<ProteoformModel> proteoforms, IPeakWidthModel widthModel)
    {
        return new ComparisonReport(
            PerScanSpectralAngles: SpectralAngle.ComputePerScan(truth, proteoforms, widthModel),
            PerChargeXicCorrelations: XicCorrelation.ComputePerCharge(truth, proteoforms, widthModel),
            Residuals: ResidualAnalyzer.Analyze(truth, proteoforms, widthModel));
    }
}
