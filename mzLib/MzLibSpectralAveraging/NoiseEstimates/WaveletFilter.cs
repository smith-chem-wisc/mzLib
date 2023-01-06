using System;

namespace MzLibSpectralAveraging;
/// <summary>
/// Object to facilitate creating and storing wavelets. Calculates and stores the Wavelet and Scaling filter coefficients
/// for subsequent use in the ModWt methods. 
/// </summary>
/// <remarks>Wavelet filter coefficients are stored as private readonly fields in this object.
/// If you want more wavelets, you need to add the coefficients to this class.</remarks>
public class WaveletFilter
{
    public double[] WaveletCoefficients { get; private set; }
    public double[] ScalingCoefficients { get; private set; }
    public WaveletType WaveletType { get; private set; }
    /// <summary>
    /// Private method that facilitates creation of the wavelet coefficients. Calculates the
    /// scaling coefficients from the wavelet coefficients using a quadrature mirror filter. 
    /// </summary>
    /// <see cref="WaveletMathUtils.QMF"/>
    /// <param name="filterCoeffs">filter coefficients that are stored as readonly fields in this class.</param>
    private void CreateFiltersFromCoeffs(double[] filterCoeffs)
    {
        WaveletCoefficients = new double[filterCoeffs.Length];
        ScalingCoefficients = new double[filterCoeffs.Length];

        // calculate wavelet coefficients
        for (int i = 0; i < ScalingCoefficients.Length; i++)
        {
            ScalingCoefficients[i] = filterCoeffs[i] / Math.Sqrt(2d);
        }
        WaveletCoefficients = WaveletMathUtils.QMF(ScalingCoefficients, inverse: true);
    }
    /// <summary>
    /// Calculates the wavelet and scaling filter coefficients given the type of wavelet and assigns to the properties.
    /// </summary>
    /// <remarks>Need to add a case to the switch for each type of supported wavelet. If you want to add more wavelets,
    /// you need to add a case to the switch in this method as well as adding the wavelet filter coefficients as
    /// private readonly fields. </remarks>
    /// <param name="waveletType">Wavelet option</param>
    public void CreateFiltersFromCoeffs(WaveletType waveletType)
    {
        switch (waveletType)
        {
            case WaveletType.Haar:
            {
                WaveletType = WaveletType.Haar;
                CreateFiltersFromCoeffs(_haarCoefficients);
                return;
            }
            case WaveletType.Db4:
            {
                WaveletType = WaveletType.Db4; 
                CreateFiltersFromCoeffs(_db4Coefficients);
                return; 
            }
        }
    }
    #region Wavelet Coefficients
    private readonly double[] _haarCoefficients =
    {
        0.7071067811865475,
        0.7071067811865475
    };

    private readonly double[] _db4Coefficients =
    {
        0.2304, 
        0.7148, 
        0.6309, 
        -0.0280, 
        -0.1870, 
        0.0308, 
        0.0329, 
        -0.0106
    };
#endregion
}