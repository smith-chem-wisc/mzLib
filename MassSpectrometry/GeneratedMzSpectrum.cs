namespace MassSpectrometry
{
    internal class GeneratedMzSpectrum : MzSpectrumZR
    {
        #region Public Constructors

        public GeneratedMzSpectrum(double[] mz, double[] intensities, bool shouldCopy) : base(mz, intensities, shouldCopy)
        {
        }

        #endregion Public Constructors

        #region Protected Methods

        protected MzPeakZR GeneratePeak(int index)
        {
            return new MzPeakZR(XArray[index], YArray[index]);
        }

        #endregion Protected Methods
    }
}