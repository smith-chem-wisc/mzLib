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

        #region Internal Fields

        internal MsDataScan<TSpectrum>[] Scans = null;

        #endregion Internal Fields

        #region Private Fields

        private string _filePath;

        private string _name;

        private bool _numSpectraSet = false;

        private int _numSpectra;

        #endregion Private Fields

        #region Protected Constructors

        protected MsDataFile(string filePath, MsDataFileType filetype = MsDataFileType.UnKnown)
        {
            FilePath = filePath;
            FileType = filetype;
        }

        #endregion Protected Constructors

        #region Public Properties

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

        public virtual int NumSpectra
        {
            get
            {
                if (_numSpectraSet)
                    return _numSpectra;
                _numSpectraSet = true;
                _numSpectra = GetNumSpectra();
                return _numSpectra;
            }
        }

        public string Name
        {
            get { return _name; }
        }

        #endregion Public Properties

        #region Public Methods

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
        /// <param name="oneBasedScanNumber">The spectrum number to get the MS Scan at</param>
        /// <returns></returns>
        public virtual IMsDataScan<TSpectrum> GetOneBasedScan(int oneBasedScanNumber)
        {
            if (Scans == null)
                Scans = new MsDataScan<TSpectrum>[NumSpectra];
            if (Scans[oneBasedScanNumber - 1] == null)
            {
                Scans[oneBasedScanNumber - 1] = GetMsDataOneBasedScanFromFile(oneBasedScanNumber);
            }

            return Scans[oneBasedScanNumber - 1];
        }

        public void LoadAllScansInMemory()
        {
            if (Scans == null)
            {
                Scans = new MsDataScan<TSpectrum>[NumSpectra];
            }

            for (int scanNumber = 1; scanNumber <= NumSpectra; scanNumber++)
            {
                if (Scans[scanNumber - 1] == null)
                {
                    Scans[scanNumber - 1] = GetMsDataOneBasedScanFromFile(scanNumber);
                }
            }
        }

        public virtual void ClearCachedScans()
        {
            Array.Clear(Scans, 0, Scans.Length);
        }

        public virtual IEnumerable<IMsDataScan<TSpectrum>> GetMsScans()
        {
            return GetMsScansInIndexRange(1, NumSpectra);
        }

        public virtual IEnumerable<IMsDataScan<TSpectrum>> GetMsScansInIndexRange(int FirstSpectrumNumber, int LastSpectrumNumber)
        {
            for (int oneBasedSpectrumNumber = FirstSpectrumNumber; oneBasedSpectrumNumber <= LastSpectrumNumber; oneBasedSpectrumNumber++)
            {
                yield return GetOneBasedScan(oneBasedSpectrumNumber);
            }
        }

        public virtual IEnumerable<IMsDataScan<TSpectrum>> GetMsScansInTimeRange(double firstRT, double lastRT)
        {
            int oneBasedSpectrumNumber = GetClosestOneBasedSpectrumNumber(firstRT);
            while (oneBasedSpectrumNumber <= NumSpectra)
            {
                IMsDataScan<TSpectrum> scan = GetOneBasedScan(oneBasedSpectrumNumber);
                double rt = scan.RetentionTime;
                oneBasedSpectrumNumber++;
                if (rt < firstRT)
                    continue;
                if (rt > lastRT)
                    yield break;
                yield return scan;
            }
        }

        public abstract int GetClosestOneBasedSpectrumNumber(double retentionTime);

        public override string ToString()
        {
            return string.Format("{0} ({1})", Name, Enum.GetName(typeof(MsDataFileType), FileType));
        }

        public abstract void Open();

        public abstract void Close();

        #endregion Public Methods

        #region Protected Methods

        protected abstract MsDataScan<TSpectrum> GetMsDataOneBasedScanFromFile(int oneBasedSpectrumNumber);

        protected abstract int GetNumSpectra();

        #endregion Protected Methods

    }
}