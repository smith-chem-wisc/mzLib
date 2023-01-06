namespace MzLibSpectralAveraging;
/// <summary>
/// Object used to unify wavelet and scaling coefficients calculated at a particular scale. 
/// </summary>
public class Level
{
    /// <summary>
    /// Creates a Level object that contains the scale and the wavelet
    /// and scaling coefficients calculated at that scale. 
    /// </summary>
    /// <param name="scale">Integer scale value for this level.</param>
    /// <param name="waveletCoeff">The wavelet coefficients at this level.</param>
    /// <param name="scalingCoeff">The scaling coefficients at this level.</param>
    public Level(int scale, double[] waveletCoeff, double[] scalingCoeff)
    {
        Scale = scale;
        WaveletCoeff = waveletCoeff;
        ScalingCoeff = scalingCoeff;
    }
    public int Scale { get; private set; }
    public double[] WaveletCoeff { get; private set; }
    public double[] ScalingCoeff { get; private set; }
}