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

namespace IO.Thermo
{
    public class ThermoStaticData : ThermoFile, IMsStaticDataFile<IThermoScan>
    {
        #region Private Constructors

        private ThermoStaticData(IThermoScan[] scans, ThermoGlobalParams p) : base(scans, p)
        {
        }

        #endregion Private Constructors

        #region Public Methods

        public static ThermoStaticData LoadAllStaticData(string filePath)
        {
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
                scans[nScanNumber - pnFirstSpectrum] = GetMsDataOneBasedScanFromThermoFile(nScanNumber, theConnection, p);

            theConnection.Close();

            return new ThermoStaticData(scans, p);
        }

        public override IThermoScan GetOneBasedScan(int oneBasedScanNumber)
        {
            return Scans[oneBasedScanNumber - 1];
        }

        #endregion Public Methods
    }
}