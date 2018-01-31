namespace FlashLFQ
{
    public class MassSpectralPeak
    {
        #region Public Fields

        public readonly double Mz;
        public readonly double Intensity;

        #endregion Public Fields

        #region Public Constructors

        public MassSpectralPeak(double Mz, double Intensity)
        {
            this.Mz = Mz;
            this.Intensity = Intensity;
        }

        #endregion Public Constructors
    }
}