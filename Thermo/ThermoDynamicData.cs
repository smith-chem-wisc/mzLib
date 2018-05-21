using MassSpectrometry;
using MSFileReaderLib;
using MzLibUtil;
using System;
using System.IO;
using System.Security.Cryptography;

namespace IO.Thermo
{
    public class ThermoDynamicData : ThermoDataFile, IDisposable
    {
        #region Private Fields

        private readonly IXRawfile5 _rawConnection;
        private readonly IFilteringParams filterParams;

        #endregion Private Fields

        #region Private Constructors

        private ThermoDynamicData(IXRawfile5 _rawConnection, IFilteringParams filterParams, int numSpectra, SourceFile sourceFile, ThermoGlobalParams thermoGlobalParams) : base(numSpectra, sourceFile)
        {
            this._rawConnection = _rawConnection;
            this.filterParams = filterParams;
            ThermoGlobalParams = thermoGlobalParams;
        }

        #endregion Private Constructors

        #region Public Methods

        public static ThermoDynamicData InitiateDynamicConnection(string filePath, IFilteringParams filterParams = null)
        {
            if (!File.Exists(filePath))
            {
                throw new FileNotFoundException();
            }

            if (!ThermoStaticData.CheckForMsFileReader())
            {
                throw new MzLibException("MsFileReader Not Installed");
            }

            IXRawfile5 _rawConnection = (IXRawfile5)new MSFileReader_XRawfile();
            _rawConnection.Open(filePath);
            _rawConnection.SetCurrentController(0, 1);

            int lastspectrumNumber = -1;
            _rawConnection.GetLastSpectrumNumber(ref lastspectrumNumber);
            int firstspectrumNumber = -1;
            _rawConnection.GetFirstSpectrumNumber(ref firstspectrumNumber);

            var precursorInfoArray = new ManagedThermoHelperLayer.PrecursorInfo[lastspectrumNumber - firstspectrumNumber + 1];

            string sendCheckSum;
            using (FileStream stream = File.OpenRead(filePath))
            {
                using (SHA1Managed sha = new SHA1Managed())
                {
                    byte[] checksum = sha.ComputeHash(stream);
                    sendCheckSum = BitConverter.ToString(checksum)
                        .Replace("-", string.Empty);
                }
            }
            SourceFile sourceFile = new SourceFile(
                @"Thermo nativeID format",
                @"Thermo RAW format",
                sendCheckSum,
                @"SHA-1",
                filePath,
                Path.GetFileNameWithoutExtension(filePath));

            var thermoGlobalParams = ThermoStaticData.GetAllGlobalStuff(_rawConnection, precursorInfoArray, filePath);

            // if the spectra file only contains 1 scan and its MS order is 0, this indicates an errored read result
            if (thermoGlobalParams.msOrderByScan.Length == 1 && thermoGlobalParams.msOrderByScan[0] == 0)
            {
                throw new MzLibException("Could not read data from file " + Path.GetFileName(filePath));
            }

            return new ThermoDynamicData(_rawConnection, filterParams, lastspectrumNumber - firstspectrumNumber + 1, sourceFile, thermoGlobalParams);
        }

        public new MsDataScan GetOneBasedScan(int oneBasedScanNumber)
        {
            if (Scans[oneBasedScanNumber - 1] == null)
                Scans[oneBasedScanNumber - 1] = ThermoStaticData.GetMsDataOneBasedScanFromThermoFile(_rawConnection, oneBasedScanNumber, ThermoGlobalParams, filterParams);
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
            {
                if (_rawConnection != null)
                {
                    _rawConnection.Close();
                }
            }
        }

        #endregion Protected Methods
    }
}