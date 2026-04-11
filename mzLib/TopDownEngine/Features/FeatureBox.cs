using MzLibUtil;

namespace TopDownEngine.Features;

public sealed record FeatureBox(
    MzRange MzRange,
    DoubleRange RtRange,
    double SeedIntensity,
    double TotalIntensity,
    string SourceFile,
    int PeakCount);
