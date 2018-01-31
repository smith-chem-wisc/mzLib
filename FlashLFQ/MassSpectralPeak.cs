namespace FlashLFQ
{
    public class MassSpectralPeak
    {
        public readonly double Mz;
        public readonly double Intensity;

        public MassSpectralPeak(double Mz, double Intensity)
        {
            this.Mz = Mz;
            this.Intensity = Intensity;
        }
    }
}
