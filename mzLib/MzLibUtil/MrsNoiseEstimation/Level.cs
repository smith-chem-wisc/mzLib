namespace MzLibUtil.NoiseEstimation;

public class Level
{
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