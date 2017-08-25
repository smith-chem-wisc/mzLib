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

        public double X { get; protected set; }
        public double Y { get; protected set; }

        #endregion Public Properties
    }
}