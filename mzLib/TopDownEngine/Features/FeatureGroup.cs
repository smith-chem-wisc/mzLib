using MassSpectrometry;
using System.Collections.Generic;

namespace TopDownEngine.Features;

public sealed record FeatureGroup(
    FeatureBox DonorBox,
    Dictionary<SpectraFileInfo, (int KTrue, double PValue, double Intensity)> Matches,
    int TotalSupport);
