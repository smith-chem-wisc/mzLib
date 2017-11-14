// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (ThermoRawFile.cs) is part of MassSpecFiles.
//
// MassSpecFiles is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MassSpecFiles is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with MassSpecFiles. If not, see <http://www.gnu.org/licenses/>.

using MassSpectrometry;
using MSFileReaderLib;
using MzLibUtil;
using System;
using System.IO;
using System.Runtime.InteropServices;
using System.Security.Cryptography;

namespace IO.Thermo
{
    public class ThermoStaticData : ThermoFile, IMsStaticDataFile<IThermoScan>
    {
        #region Private Constructors

        private ThermoStaticData(IThermoScan[] scans, ThermoGlobalParams p, SourceFile sourceFile) : base(scans, p, sourceFile)
        {
        }

        #endregion Private Constructors

        #region Public Methods
        private const string THERMO_READER_CLSID = "{1d23188d-53fe-4c25-b032-dc70acdbdc02}";

        public static ThermoStaticData LoadAllStaticData(string filePath, int? topNpeaks = null, double? minRatio = null, bool trimMs1Peaks = true, bool trimMsMsPeaks = true)
        {

            try
            {
                var thermoReader = Type.GetTypeFromCLSID(Guid.Parse(THERMO_READER_CLSID));
                Console.WriteLine("Instantiated Type object from CLSID {0}",
                               THERMO_READER_CLSID);
                Object wordObj = Activator.CreateInstance(thermoReader);
                Console.WriteLine("\n\n\n\n\n\n\n\n\n\n\n Instantiated {0}",
                                  wordObj.GetType().FullName, THERMO_READER_CLSID);
            }
            catch(COMException ex)
            {
                if(ex.ErrorCode == -2147287036)
                {
                    throw new MzLibException("MS File Reader Not Installed");
                }
            }


            var ok = new ManagedThermoHelperLayer.HelperClass();
            IXRawfile5 theConnection = (IXRawfile5)new MSFileReader_XRawfile();
            theConnection.Open(filePath);
            int pbSMData = 0;
            theConnection.IsThereMSData(ref pbSMData);
            if (pbSMData == 0)
                throw new MzLibException("File not found");

            theConnection.SetCurrentController(0, 1);

            var precursorInfosArray = ok.GetAllPrecursorInfos(filePath);
            for (int i = 0; i < precursorInfosArray.Length; i++)
            {
                if (precursorInfosArray[i].nScanNumber == 0)
                    precursorInfosArray[i].nScanNumber = -1;
            }

            int pnFirstSpectrum = 0;
            theConnection.GetFirstSpectrumNumber(ref pnFirstSpectrum);
            int pnLastSpectrum = 0;
            theConnection.GetLastSpectrumNumber(ref pnLastSpectrum);

            ThermoGlobalParams p = GetAllGlobalStuff(theConnection, precursorInfosArray, filePath);

            IThermoScan[] scans = new IThermoScan[pnLastSpectrum - pnFirstSpectrum + 1];
            for (int nScanNumber = pnFirstSpectrum; nScanNumber <= pnLastSpectrum; nScanNumber++)
                scans[nScanNumber - pnFirstSpectrum] = GetMsDataOneBasedScanFromThermoFile(nScanNumber, theConnection, p, topNpeaks, minRatio, trimMs1Peaks, trimMsMsPeaks);

            theConnection.Close();

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

            return new ThermoStaticData(scans, p, sourceFile);
        }

        public override IThermoScan GetOneBasedScan(int oneBasedScanNumber)
        {
            return Scans[oneBasedScanNumber - 1];
        }

        #endregion Public Methods
    }
}