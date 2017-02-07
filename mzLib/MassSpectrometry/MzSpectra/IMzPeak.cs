using Spectra;

namespace MassSpectrometry
{
    public interface IMzPeak : IPeak
    {
        double Intensity { get; }

        double Mz { get; }
    }
}