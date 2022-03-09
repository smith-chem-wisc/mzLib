namespace MassSpectrometry
{
    internal class GeneratedMzSpectrum : MzSpectrum
    {
        public GeneratedMzSpectrum(double[] mz, double[] intensities, bool shouldCopy) : base(mz, intensities, shouldCopy)
        {
        }

        protected MzPeak GeneratePeak(int index)
        {
            return new MzPeak(XArray[index], YArray[index]);
        }
    }
}