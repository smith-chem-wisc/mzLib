using Spectra;

namespace MassSpectrometry
{
    public interface IMzPeak : IPeak
    {
        #region Public Properties

        double Intensity { get; }

        double Mz { get; }

        #endregion Public Properties
    }
}