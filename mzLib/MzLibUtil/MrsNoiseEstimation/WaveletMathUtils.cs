using System;

namespace MzLibUtil.NoiseEstimation;
internal static class WaveletMathUtils
{
    /// <summary>
    /// Calcualte the quadruture mirror filter for x, an array of filter coefficients
    /// </summary>
    /// <param name="x"></param>
    /// <param name="inverse"></param>
    /// <returns></returns>
    internal static double[] QMF(double[] x, bool inverse = false)
    {
        double[] y = new double[x.Length];
        Buffer.BlockCopy(x, 0, y, 0, x.Length * sizeof(double));
        Array.Reverse(y);
        if (inverse)
        {
            // start the for loop from the back
            int firstIndex = 1;
            for (int i = x.Length - 1; i >= 0; i--)
            {
                y[i] *= Math.Pow(-1d, firstIndex);
                firstIndex++;
            }
        }
        else
        {
            for (int i = 0; i < x.Length; i++)
            {
                y[i] *= Math.Pow(-1d, i + 1);
            }
        }
        return y;
    }
}