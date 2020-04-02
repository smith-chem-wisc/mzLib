using MassSpectrometry;
using MzLibUtil;
using System.IO;
using System.Linq;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.Interfaces;
using ThermoFisher.CommonCore.RawFileReader;
using UsefulProteomicsDatabases;

namespace ThermoRawFileReader
{
    public class ThermoDynamicData
    {
        private IRawDataPlus dynamicConnection;
        public readonly int[] MsOrdersByScan;

        public ThermoDynamicData(string filePath)
        {
            InitiateDynamicConnection(filePath);
            MsOrdersByScan = GetMsOrdersByScanInDynamicConnection();

            if (MsOrdersByScan.Length == 0)
            {
                throw new MzLibException("The file contained zero scans and could not be loaded: " + filePath);
            }
        }

        /// <summary>
        /// Gets all the MS orders of all scans in a dynamic connection. This is useful if you want to open all MS1 scans
        /// without loading all of the other MSn scans.
        /// </summary>
        private int[] GetMsOrdersByScanInDynamicConnection()
        {
            if (dynamicConnection != null)
            {
                int lastSpectrum = dynamicConnection.RunHeaderEx.LastSpectrum;
                var scanEvents = dynamicConnection.GetScanEvents(1, lastSpectrum);

                int[] msorders = scanEvents.Select(p => (int)p.MSOrder).ToArray();

                return msorders;
            }

            return null;
        }

        /// <summary>
        /// Initiates a dynamic connection with a Thermo .raw file. Data can be "streamed" instead of loaded all at once. Use 
        /// GetOneBasedScanFromDynamicConnection to get data from a particular scan. Use CloseDynamicConnection to close the 
        /// dynamic connection after all desired data has been retrieved from the dynamic connection.
        /// </summary>
        private void InitiateDynamicConnection(string filePath)
        {
            Loaders.LoadElements();

            if (dynamicConnection != null)
            {
                dynamicConnection.Dispose();
            }

            dynamicConnection = RawFileReaderAdapter.FileFactory(filePath);

            if (!File.Exists(filePath))
            {
                throw new FileNotFoundException();
            }

            if (!dynamicConnection.IsOpen)
            {
                throw new MzLibException("Unable to access RAW file!");
            }

            if (dynamicConnection.IsError)
            {
                throw new MzLibException("Error opening RAW file!");
            }

            if (dynamicConnection.InAcquisition)
            {
                throw new MzLibException("RAW file still being acquired!");
            }

            dynamicConnection.SelectInstrument(Device.MS, 1);
        }

        /// <summary>
        /// Allows access to a .raw file one scan at a time via an open dynamic connection. Returns null if the raw file does not contain the 
        /// scan number specified. Use InitiateDynamicConnection to open a dynamic connection before using this method.
        /// </summary>
        public MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
        {
            if (dynamicConnection == null)
            {
                throw new MzLibException("The dynamic connection has not been created yet!");
            }

            if (oneBasedScanNumber > dynamicConnection.RunHeaderEx.LastSpectrum || oneBasedScanNumber < dynamicConnection.RunHeaderEx.FirstSpectrum)
            {
                return null;
            }

            return ThermoRawFileReaderData.GetOneBasedScan(dynamicConnection, filterParams, oneBasedScanNumber);
        }

        /// <summary>
        /// Disposes of the dynamic connection, if one is open.
        /// </summary>
        public void CloseDynamicConnection()
        {
            if (dynamicConnection != null)
            {
                dynamicConnection.Dispose();
            }
        }
    }
}
