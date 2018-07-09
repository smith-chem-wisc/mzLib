namespace MassSpectrometry
{
    public class MzPeak
    {
        public double Mz { get; protected set; }
        public double Intensity { get; protected set; }

        public MzPeak(double mz, double intensity)
        {
            Mz = mz;
            Intensity = intensity;
        }

        public override string ToString()
        {
            return string.Format("({0:G7},{1:G7})", Mz, Intensity);
        }
    }
}