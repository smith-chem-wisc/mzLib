namespace Spectra
{
    public abstract class Peak : IPeak
    {
        #region Public Constructors

        public Peak(double x, double y)
        {
            X = x;
            Y = y;
        }

        #endregion Public Constructors

        #region Public Properties

        public double X { get; private set; }
        public double Y { get; private set; }

        #endregion Public Properties
    }
}