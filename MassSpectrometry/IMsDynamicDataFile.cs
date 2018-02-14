using System;

namespace MassSpectrometry
{
    public interface IMsDynamicDataFile<TScan> : IDisposable
        where TScan : IMsDataScan<IMzSpectrum<IMzPeak>>
    {
        #region Public Methods

        void ClearCachedScans();

        #endregion Public Methods
    }
}