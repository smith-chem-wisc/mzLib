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

using System;
using System.Collections;
using System.Collections.Generic;

namespace MassSpectrometry
{
    /// <summary>
    /// A class for interacting with data collected from a Mass Spectrometer, and stored in a file
    /// </summary>
    public abstract class MsDataFile<TScan> : IMsDataFile<TScan>
        where TScan : IMsDataScan<IMzSpectrum<IMzPeak>>
    {
        #region Protected Fields

        protected TScan[] Scans;

        #endregion Protected Fields

        #region Public Constructors

        public MsDataFile(int numSpectra)
        {
            Scans = new TScan[numSpectra];
        }

        public MsDataFile(TScan[] scans)
        {
            Scans = scans;
        }

        #endregion Public Constructors

        #region Public Properties

        public int NumSpectra
        {
            get
            {
                return Scans.Length;
            }
        }

        #endregion Public Properties

        #region Public Methods

        public abstract TScan GetOneBasedScan(int scanNumber);

        public IEnumerable<TScan> GetMsScansInIndexRange(int FirstSpectrumNumber, int LastSpectrumNumber)
        {
            for (int oneBasedSpectrumNumber = FirstSpectrumNumber; oneBasedSpectrumNumber <= LastSpectrumNumber; oneBasedSpectrumNumber++)
            {
                yield return GetOneBasedScan(oneBasedSpectrumNumber);
            }
        }

        public IEnumerable<TScan> GetMsScansInTimeRange(double firstRT, double lastRT)
        {
            int oneBasedSpectrumNumber = GetClosestOneBasedSpectrumNumber(firstRT);
            while (oneBasedSpectrumNumber <= NumSpectra)
            {
                TScan scan = GetOneBasedScan(oneBasedSpectrumNumber);
                double rt = scan.RetentionTime;
                if (rt < firstRT)
                {
                    oneBasedSpectrumNumber++;
                    continue;
                }
                if (rt > lastRT)
                    yield break;
                yield return scan;
                oneBasedSpectrumNumber++;
            }
        }

        public virtual int GetClosestOneBasedSpectrumNumber(double retentionTime)
        {
            // TODO need to convert this to a binary search of some sort. Or if the data is indexedMZML see if the indices work better.
            double bestDiff = double.MaxValue;
            for (int i = 0; i < NumSpectra; i++)
            {
                double diff = Math.Abs(GetOneBasedScan(i + 1).RetentionTime - retentionTime);
                if (diff > bestDiff)
                    return i;
                bestDiff = diff;
            }
            return NumSpectra;
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetMsScansInIndexRange(1, NumSpectra).GetEnumerator();
        }

        public IEnumerator<TScan> GetEnumerator()
        {
            return GetMsScansInIndexRange(1, NumSpectra).GetEnumerator();
        }

        #endregion Public Methods
    }
}