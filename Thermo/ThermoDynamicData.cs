using MassSpectrometry;
using MSFileReaderLib;
using System;

namespace IO.Thermo
{
    public class ThermoDynamicData : ThermoFile, IMsDynamicDataFile<IThermoScan>
    {

        #region Private Fields

        private IXRawfile5 _rawConnection;

        #endregion Private Fields

        #region Private Constructors

        private ThermoDynamicData(IXRawfile5 _rawConnection, int numSpectra, ManagedThermoHelperLayer.PrecursorInfo[] couldBePrecursor) : base(_rawConnection, numSpectra, couldBePrecursor)
        {
            this._rawConnection = _rawConnection;
        }

        #endregion Private Constructors

        #region Public Methods

        public static ThermoDynamicData InitiateDynamicConnection(string fileName)
        {
            var ok = new ManagedThermoHelperLayer.HelperClass();
            var nice = ok.GetAllPrecursorInfos(fileName);
            IXRawfile5 _rawConnection = (IXRawfile5)new MSFileReader_XRawfile();
            _rawConnection.Open(fileName);
            _rawConnection.SetCurrentController(0, 1);

            int lastspectrumNumber = -1;
            _rawConnection.GetLastSpectrumNumber(ref lastspectrumNumber);
            int firstspectrumNumber = -1;
            _rawConnection.GetFirstSpectrumNumber(ref firstspectrumNumber);

            return new ThermoDynamicData(_rawConnection, lastspectrumNumber - firstspectrumNumber + 1, nice);
        }

        public override IThermoScan GetOneBasedScan(int oneBasedScanNumber)
        {
            if (Scans[oneBasedScanNumber - 1] == null)
                Scans[oneBasedScanNumber - 1] = GetMsDataOneBasedScanFromThermoFile(oneBasedScanNumber, _rawConnection, ThermoGlobalParams);
            return Scans[oneBasedScanNumber - 1];
        }

        public virtual void ClearCachedScans()
        {
            Array.Clear(Scans, 0, Scans.Length);
        }

        public override int GetClosestOneBasedSpectrumNumber(double retentionTime)
        {
            int spectrumNumber = 0;
            _rawConnection.ScanNumFromRT(retentionTime, ref spectrumNumber);
            return spectrumNumber;
        }

        public void Dispose()
        {
            Dispose(true);
            GC.SuppressFinalize(this);
        }

        #endregion Public Methods

        #region Protected Methods

        protected virtual void Dispose(bool disposing)
        {
            if (disposing)
                if (_rawConnection != null)
                    _rawConnection.Close();
        }

        #endregion Protected Methods

    }
}