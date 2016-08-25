// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
// 
// This file (MsDataFile.cs) is part of MassSpectrometry.
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

using Spectra;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;

namespace MassSpectrometry
{
    /// <summary>
    /// A data file for storing data collected from a Mass Spectrometer
    /// </summary>
    public abstract class MsDataFile<TSpectrum> : IMsDataFile<TSpectrum>
        where TSpectrum : IMzSpectrum<MzPeak>
    {
        /// <summary>
        /// Defines if MS scans should be cached for quicker retrieval. Cached scans are held in an internal
        /// array and don't get cleared until the file is disposed or the ClearCacheScans() method is called.
        /// Of course, if you store the scans somewhere else, they will persist. The default value is True.
        /// </summary>
        public bool CacheScans { get; private set; }

        internal MsDataScan<TSpectrum>[] Scans = null;

        private string _filePath;

        private string _name;

        protected MsDataFile(string filePath, bool cacheScans, MsDataFileType filetype = MsDataFileType.UnKnown)
        {
            FilePath = filePath;
            FileType = filetype;
            CacheScans = cacheScans;
        }

        public string FilePath
        {
            get { return _filePath; }
            private set
            {
                _filePath = value;
                _name = Path.GetFileNameWithoutExtension(value);
            }
        }

        public MsDataFileType FileType { get; private set; }

        bool _firstSpectrumNumberSet = false;
        int _firstSpectrumNumber;
        public virtual int FirstSpectrumNumber
        {
            get
            {
                if (_firstSpectrumNumberSet)
                    return _firstSpectrumNumber;
                _firstSpectrumNumberSet = true;
                _firstSpectrumNumber = GetFirstSpectrumNumber();
                return _firstSpectrumNumber;
            }
        }

        bool _lastSpectrumNumberSet = false;
        int _lastSpectrumNumber;
        public virtual int LastSpectrumNumber
        {
            get
            {
                if (_lastSpectrumNumberSet)
                    return _lastSpectrumNumber;
                _lastSpectrumNumberSet = true;
                _lastSpectrumNumber = GetLastSpectrumNumber();
                return _lastSpectrumNumber;
            }
        }

        public string Name
        {
            get { return _name; }
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }

        public IEnumerator<IMsDataScan<TSpectrum>> GetEnumerator()
        {
            return GetMsScans().GetEnumerator();
        }

        /// <summary>
        /// Get the MS Scan at the specific spectrum number.
        /// </summary>
        /// <param name="scanNumber">The spectrum number to get the MS Scan at</param>
        /// <returns></returns>
        public virtual IMsDataScan<TSpectrum> GetScan(int scanNumber)
        {
            if (!CacheScans)
                return GetMsDataScanFromFile(scanNumber);
            if (Scans == null)
                Scans = new MsDataScan<TSpectrum>[LastSpectrumNumber - FirstSpectrumNumber + 1];
            if (Scans[scanNumber - FirstSpectrumNumber] == null)
            {
                Scans[scanNumber - FirstSpectrumNumber] = GetMsDataScanFromFile(scanNumber);
            }

            return Scans[scanNumber - FirstSpectrumNumber];
        }

        public virtual void LoadAllScansInMemory()
        {
            if (Scans == null)
            {
                Scans = new MsDataScan<TSpectrum>[LastSpectrumNumber - FirstSpectrumNumber + 1];
            }

            for (int scanNumber = FirstSpectrumNumber; scanNumber <= LastSpectrumNumber; scanNumber++)
            {
                if (Scans[scanNumber - FirstSpectrumNumber] == null)
                {
                    Scans[scanNumber - FirstSpectrumNumber] = GetMsDataScanFromFile(scanNumber);
                }
            }
        }

        public virtual void ClearCachedScans()
        {
            if (Scans == null)
                return;
            Array.Clear(Scans, 0, Scans.Length);
        }

        protected abstract MsDataScan<TSpectrum> GetMsDataScanFromFile(int spectrumNumber);

        public virtual IEnumerable<IMsDataScan<TSpectrum>> GetMsScans()
        {
            return GetMsScansInIndexRange(FirstSpectrumNumber, LastSpectrumNumber);
        }

        public virtual IEnumerable<IMsDataScan<TSpectrum>> GetMsScansInIndexRange(int FirstSpectrumNumber, int LastSpectrumNumber)
        {
            for (int spectrumNumber = FirstSpectrumNumber; spectrumNumber <= LastSpectrumNumber; spectrumNumber++)
            {
                yield return GetScan(spectrumNumber);
            }
        }

        public virtual IEnumerable<IMsDataScan<TSpectrum>> GetMsScansInTimeRange(double firstRT, double lastRT)
        {
            int spectrumNumber = GetSpectrumNumber(firstRT - 0.0000001);
            while (spectrumNumber <= LastSpectrumNumber)
            {
                IMsDataScan<TSpectrum> scan = GetScan(spectrumNumber++);
                double rt = scan.RetentionTime;
                if (rt < firstRT)
                    continue;
                if (rt > lastRT)
                    yield break;
                yield return scan;
            }
        }

        public abstract int GetSpectrumNumber(double retentionTime);

        public override string ToString()
        {
            return string.Format("{0} ({1})", Name, Enum.GetName(typeof(MsDataFileType), FileType));
        }

        protected abstract int GetFirstSpectrumNumber();

        protected abstract int GetLastSpectrumNumber();
        public abstract void Open();
    }
}