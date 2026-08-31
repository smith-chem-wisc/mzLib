using System;

namespace MzLibUtil.NoiseEstimation;
internal static class WaveletMathUtils
{
    /// <summary>
    /// Calculate the quadrature mirror filter for x, an array of filter coefficients.
    /// The alternating sign runs from the end of the array towards the front.
    /// </summary>
    /// <param name="x">Filter coefficients to mirror.</param>
    /// <returns>The mirrored coefficients.</returns>
    internal static double[] QMF(double[] x)
    {
        double[] y = new double[x.Length];
        Buffer.BlockCopy(x, 0, y, 0, x.Length * sizeof(double));
        Array.Reverse(y);

        int firstIndex = 1;
        for (int i = x.Length - 1; i >= 0; i--)
        {
            y[i] *= Math.Pow(-1d, firstIndex);
            firstIndex++;
        }

        return y;
    }
}