namespace MassSpectrometry
{
    internal class GeneratedMzSpectrum : MzSpectrum<IMzPeak>
    {
        #region Public Constructors

        public GeneratedMzSpectrum(double[] mz, double[] intensities, bool shouldCopy) : base(mz, intensities, shouldCopy)
        {
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override IMzPeak GeneratePeak(int index)
        {
            return new MzPeak(XArray[index], YArray[index]);
        }

        #endregion Protected Methods
    }
}