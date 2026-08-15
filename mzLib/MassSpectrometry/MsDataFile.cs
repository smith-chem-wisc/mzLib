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
using System.Linq;

namespace MassSpectrometry
{
    /// <summary>
    /// A class for interacting with data collected from a Mass Spectrometer, and stored in a file
    /// </summary>
    public abstract class MsDataFile : IEnumerable<MsDataScan>
    {
        protected readonly object DynamicReadingLock = new();
        public MsDataScan[] Scans { get; protected set; }
        public SourceFile SourceFile { get; set; }
        public int NumSpectra => Scans?.Length ?? 0;
        public string FilePath { get; }
        protected MsDataFile(int numSpectra, SourceFile sourceFile)
        {
            Scans = new MsDataScan[numSpectra];
            SourceFile = sourceFile;
        }

        protected MsDataFile(MsDataScan[] scans, SourceFile sourceFile)
        {
            Scans = scans;
            SourceFile = sourceFile;
        }

        protected MsDataFile(string filePath)
        {
            FilePath = filePath;
        }

        public MsDataScan this[int index] => Scans[index];

        #region Abstract members

        // static connection
        public abstract MsDataFile LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1);

        public abstract SourceFile GetSourceFile();

        // Dynamic Connection
        public abstract MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber,
            IFilteringParams filterParams = null);

        public abstract void CloseDynamicConnection();
        public abstract void InitiateDynamicConnection();

        #endregion

        #region Utilities

        public virtual MsDataScan[] GetMsDataScans()
        {
            if (!CheckIfScansLoaded())
                LoadAllStaticData();
            return Scans;
        }

        public virtual List<MsDataScan> GetAllScansList()
        {
            if (!CheckIfScansLoaded())
                LoadAllStaticData();

            return Scans.ToList();
        }

        public virtual IEnumerable<MsDataScan> GetMS1Scans()
        {
            if (!CheckIfScansLoaded())
                LoadAllStaticData();
            for (int i = 1; i <= NumSpectra; i++)
            {
                var scan = GetOneBasedScan(i);
                if (scan.MsnOrder == 1)
                {
                    yield return scan;
                }
            }
        }

        public virtual MsDataScan GetOneBasedScan(int scanNumber)
        {
            if (!CheckIfScansLoaded())
                LoadAllStaticData();

            return Scans[scanNumber - 1];
        }

        public virtual IEnumerable<MsDataScan> GetMsScansInIndexRange(int firstSpectrumNumber, int lastSpectrumNumber)
        {
            if (!CheckIfScansLoaded())
                LoadAllStaticData();

            for (int oneBasedSpectrumNumber = firstSpectrumNumber;
                 oneBasedSpectrumNumber <= lastSpectrumNumber;
                 oneBasedSpectrumNumber++)
            {
                yield return GetOneBasedScan(oneBasedSpectrumNumber);
            }
        }

        public virtual IEnumerable<MsDataScan> GetMsScansInTimeRange(double firstRT, double lastRT)
        {
            if (!CheckIfScansLoaded())
                LoadAllStaticData();

            int oneBasedSpectrumNumber = GetClosestOneBasedSpectrumNumber(firstRT);
            while (oneBasedSpectrumNumber <= NumSpectra)
            {
                MsDataScan scan = GetOneBasedScan(oneBasedSpectrumNumber);
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
            if (!CheckIfScansLoaded())
                LoadAllStaticData();

            if (NumSpectra == 0)
                return 0;

            // Binary search rather than a walk. Retention times ascending is not a new assumption: the
            // linear version this replaces stopped as soon as a scan's distance to the target grew, so it
            // never looked past the first local minimum either. Making it explicit turns O(n) into
            // O(log n), which matters because callers such as GetMsScansInTimeRange start here.
            int low = 1;
            int high = NumSpectra;
            while (low < high)
            {
                int mid = low + (high - low) / 2;
                if (GetOneBasedScan(mid).RetentionTime < retentionTime)
                    low = mid + 1;
                else
                    high = mid;
            }

            // low is now the first scan at or after retentionTime, so the closest is either it or its
            // predecessor. Strict less-than keeps the previous behaviour of preferring the later scan
            // when a target falls exactly between two of them.
            if (low > 1)
            {
                double distanceAfter = Math.Abs(GetOneBasedScan(low).RetentionTime - retentionTime);
                double distanceBefore = Math.Abs(GetOneBasedScan(low - 1).RetentionTime - retentionTime);
                if (distanceBefore < distanceAfter)
                    return low - 1;
            }

            // Scans can share a retention time. The linear version kept walking while the distance stayed
            // equal, so it came to rest on the last scan of such a run; binary search lands on the first.
            // Advancing to the end of the run keeps the answer identical to the previous implementation.
            double closestRetentionTime = GetOneBasedScan(low).RetentionTime;
            while (low < NumSpectra && GetOneBasedScan(low + 1).RetentionTime == closestRetentionTime)
                low++;

            return low;
        }

        public virtual int[] GetMsOrderByScanInDynamicConnection()
        {
            throw new NotImplementedException();
        }

        #endregion

        public virtual bool CheckIfScansLoaded()
        {
            return (Scans != null && Scans.Length > 0);
        }

        public IEnumerator<MsDataScan> GetEnumerator()
        {
            return Scans.Where(scan => scan is not null).GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
        }
    }

    
}