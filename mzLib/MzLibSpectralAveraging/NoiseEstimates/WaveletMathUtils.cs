using System;

namespace MzLibSpectralAveraging;
internal static class WaveletMathUtils
{
    /// <summary>
    /// Calcualte the quadruture mirror filter for x, an array of filter coefficients
    /// </summary>
    /// <param name="x">An array of filter coefficients, either generated via a program like MATLAB or found in a reference.</param>
    /// <param name="inverse">Determines whether the inverse filter should be returned. Set to true if the filter is being generated
    /// for a maximal overlap discrete wavelet transform. Set to false if generating a scaling filter for the discrete wavelet transform.</param>
    /// <returns>A quadrature mirror filter transformed set of wavelet coefficients.</returns>
    internal static double[] QMF(double[] x, bool inverse = false)
    {

        // QMF only uses the inverse to create the MODWT filter inputs. If you want to 
        // expand into any other type of wavelet or wavelet transform, you may need to 
        // use the inverse. Specifically, you may need to use the inverse for the discrete wavelet
        // transform. Uncomment the else and then it should be good to good. -AVC 

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
        //else
        //{
        //    for (int i = 0; i < x.Length; i++)
        //    {
        //        y[i] *= Math.Pow(-1d, i + 1);
        //    }
        //}
        return y;
    }
}