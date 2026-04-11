namespace TopDownEngine.Alignment;

public readonly record struct AnchorBin(int BinIndex, double MzCenter, int FileCount, double MeanMaximumIntensity);
