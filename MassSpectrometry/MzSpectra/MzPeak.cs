using Spectra;

namespace MassSpectrometry
{
    public class MzPeak : Peak, IMzPeak
    {
        #region Public Constructors

        public MzPeak(double mz, double intensity)
            : base(mz, intensity)
        {
        }

        #endregion Public Constructors

        #region Public Properties

        public double Intensity
        {
            get
            {
                return Y;
            }
        }

        public double Mz
        {
            get
            {
                return X;
            }
        }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            return string.Format("({0:G7},{1:G7})", X, Y);
        }

        #endregion Public Methods
    }
}