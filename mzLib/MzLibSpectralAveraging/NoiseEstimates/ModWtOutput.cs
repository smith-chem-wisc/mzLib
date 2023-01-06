using System;
using System.Collections.Generic;

namespace MzLibSpectralAveraging;

/// <summary>
/// Class used to organize and unify access to relevant information resulting from the
/// maximum overlap discrete wavelet transform of a signal. Exposes a a list of levels and the
/// maximum scale calculated. 
/// </summary>
/// <see cref="WaveletMath"/>
/// <seealso cref="Level"/>
internal class ModWtOutput
{
    /// <summary>
    /// Creates a ModWtOutput object given a maximum number of scales. 
    /// </summary>
    /// <param name="maxScale">The highest scale to be calculated.</param>
    internal ModWtOutput(int maxScale)
    {
        Levels = new List<Level>();
        MaxScale = maxScale;
    }

    internal List<Level> Levels { get; private set; }
    internal int MaxScale { get; private set; }
    /// <summary>
    /// Creates a new Level object and stores it in the Levels property. If BoundaryType == Reflection,
    /// it calculates the true start and stop indices and subsets the results arrays so that no extra coefficients are stored. 
    /// </summary>
    /// <remarks>When using a BoundaryType.Reflection MODWT, which is the default and no others are currently implemented,
    /// the resulting wavelet and scaling coefficients are 2x longer than the original signal length. This method subsets the
    /// results arrays so that the coefficients corresponding to only the original signal length is stored.</remarks>
    /// <param name="waveletCoeff">Wavelet coefficients resulting from a ModWt</param>
    /// <param name="scalingCoeff">Scaling coefficients resulting from a ModWt</param>
    /// <param name="scale">The scale of the coefficients being stored.</param>
    /// <param name="boundaryType">The boundary type used in the ModWt calculation.</param>
    /// <param name="originalSignalLength">The original signal length (before it was multiplied by 2 in the case of reflected boundaries.</param>
    /// <param name="filterLength">The length of the wavelet or scaling filter (both are equivalent).</param>
    internal void AddLevel(double[] waveletCoeff, double[] scalingCoeff, int scale,
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