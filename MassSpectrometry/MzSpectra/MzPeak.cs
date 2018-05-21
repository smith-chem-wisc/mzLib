namespace MassSpectrometry
{
    public class MzPeak
    {

        #region Public Properties

        public double Mz { get; protected set; }
        public double Intensity { get; protected set; }

        #endregion Public Properties

        #region Public Constructors

        public MzPeak(double mz, double intensity)
        {
            Mz = mz;
            Intensity = intensity;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            return string.Format("({0:G7},{1:G7})", Mz, Intensity);
        }

        #endregion Public Methods
    }
}