// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (MsDataScan.cs) is part of MassSpectrometry.
//
// MassSpectrometry is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MassSpectrometry is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with MassSpectrometry. If not, see <http://www.gnu.org/licenses/>.

using MzLibUtil;
using System.Collections.Generic;

namespace MassSpectrometry
{
    public class MsDataScan<TSpectrum> : IMsDataScan<TSpectrum>
        where TSpectrum : IMzSpectrum<IMzPeak>
    {
        #region Public Constructors

        public MsDataScan(TSpectrum massSpectrum, int oneBasedScanNumber, int msnOrder, bool isCentroid, Polarity polarity, double retentionTime, MzRange scanWindowRange, string scanFilter, MZAnalyzerType mzAnalyzer, double totalIonCurrent, double? injectionTime, double[,] noiseData, string nativeId)
        {
            OneBasedScanNumber = oneBasedScanNumber;
            MsnOrder = msnOrder;
            IsCentroid = isCentroid;
            Polarity = polarity;
            RetentionTime = retentionTime;
            ScanWindowRange = scanWindowRange;
            ScanFilter = scanFilter;
            MzAnalyzer = mzAnalyzer;
            TotalIonCurrent = totalIonCurrent;
            InjectionTime = injectionTime;
            NoiseData = noiseData;
            MassSpectrum = massSpectrum;
            NativeId = nativeId;
        }

        #endregion Public Constructors

        #region Public Properties

        /// <summary>
        /// The mass spectrum associated with the scan
        /// </summary>
        public TSpectrum MassSpectrum { get; protected set; }

        public int OneBasedScanNumber { get; }

        public int MsnOrder { get; }

        public double RetentionTime { get; }

        public Polarity Polarity { get; }

        public MZAnalyzerType MzAnalyzer { get; }

        public MzRange ScanWindowRange { get; }

        public string ScanFilter { get; }

        public string NativeId { get; }

        public bool IsCentroid { get; }

        public double TotalIonCurrent { get; }

        public double? InjectionTime { get; }
        public double[,] NoiseData { get; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            return string.Format("Scan #{0}", OneBasedScanNumber);
        }

        public byte[] Get64BitNoiseDataMass()
        {
            return MzSpectrum<IMzPeak>.Get64Bitarray(GetNoiseDataMass(NoiseData));
        }

        public byte[] Get64BitNoiseDataNoise()
        {
            return MzSpectrum<IMzPeak>.Get64Bitarray(GetNoiseDataNoise(NoiseData));
        }

        public byte[] Get64BitNoiseDataBaseline()
        {
            return MzSpectrum<IMzPeak>.Get64Bitarray(GetNoiseDataBaseline(NoiseData));
        }

        #endregion Public Methods

        #region Private Methods

        private IEnumerable<double> GetNoiseDataMass(double[,] noiseData)
        {
            for (int i = 0; i < noiseData.Length / 3; i++)
                yield return noiseData[0, i];
        }

        private IEnumerable<double> GetNoiseDataNoise(double[,] noiseData)
        {
            for (int i = 0; i < noiseData.Length / 3; i++)
                yield return noiseData[1, i];
        }

        private IEnumerable<double> GetNoiseDataBaseline(double[,] noiseData)
        {
            for (int i = 0; i < noiseData.Length / 3; i++)
                yield return noiseData[2, i];
        }

        #endregion Private Methods
    }
}