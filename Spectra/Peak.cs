namespace Spectra
{
    public abstract class Peak : IPeak
    {
        public Peak(double x, double y)
        {
            X = x;
            Y = y;
        }

        public double X { get; protected set; }
        public double Y { get; protected set; }
    }
}