using System;
using System.Collections.Generic;
using System.Text;

namespace SpectralAveraging;

public class ModWtOutput
{
    public ModWtOutput(int maxScale, WaveletType waveletType, BoundaryType boundaryType)
    {
        Levels = new List<Level>();
        MaxScale = maxScale;
        WaveletType = waveletType;
        BoundaryType = boundaryType;
    }

    public List<Level> Levels { get; private set; }
    public int MaxScale { get; private set; }
    public WaveletType WaveletType { get; }
    public BoundaryType BoundaryType { get; }

    public void AddLevel(double[] waveletCoeff, double[] scalingCoeff, int scale,
        BoundaryType boundaryType, int originalSignalLength, int filterLength)
    {
        if (boundaryType == BoundaryType.Reflection)
        {
            int startIndex = ((int)Math.Pow(2, scale) - 1) * (filterLength - 1);
            int stopIndex = Math.Min(startIndex + originalSignalLength, waveletCoeff.Length - 1);
            Levels.Add(new Level(scale,
                waveletCoeff[startIndex..stopIndex],
                scalingCoeff[startIndex..stopIndex]));
        }
    }
}