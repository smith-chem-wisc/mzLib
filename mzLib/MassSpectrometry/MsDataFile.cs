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
            //
            // Indexes Scans rather than calling GetOneBasedScan. Mgf and MsAlign override that accessor to
            // read a separate null-padded array sized to the largest scan number, so probing an arbitrary
            // slot in [1, NumSpectra] throws on a file whose scan numbers have gaps. Scans is the array
            // NumSpectra is measured from, so searching it directly cannot land on a hole those two readers
            // introduced.
            //
            // The result is converted back through OneBasedScanNumber rather than returned as a position,
            // because for those same two readers the two differ. Everywhere else they coincide: the base
            // GetOneBasedScan is Scans[n - 1], and Mzml keeps that true even for an out-of-order file by
            // rebuilding Scans indexed by scan number. So this changes no answer that was previously
            // correct, and starts giving a usable one where the position was meaningless.
            //
            // Chemistry.ClassExtensions.GetClosestIndex is the existing helper for this shape of search and
            // was considered. It does not fit: it takes a double[], so it would need a materialised array of
            // retention times, which is O(n) per call unless cached on the file and invalidated on load; and
            // Array.BinarySearch is documented as unspecified for duplicate keys, so it cannot reproduce the
            // last-scan-of-an-equal-run answer below. Mslindex.BinarySearchLowerBound hand-rolls its loop for
            // the same reason.
            int low = 0;
            int high = NumSpectra - 1;
            while (low < high)
            {
                int mid = low + (high - low) / 2;
                if (Scans[mid].RetentionTime < retentionTime)
                    low = mid + 1;
                else
                    high = mid;
            }

            // low is now the first scan at or after retentionTime, so the closest is either it or its
            // predecessor. Strict less-than keeps the previous behaviour of preferring the later scan
            // when a target falls exactly between two of them. No run-advance is needed on this branch:
            // the search has already landed past the run, so its predecessor is the last member.
            if (low > 0)
            {
                double distanceAfter = Math.Abs(Scans[low].RetentionTime - retentionTime);
                double distanceBefore = Math.Abs(Scans[low - 1].RetentionTime - retentionTime);
                if (distanceBefore < distanceAfter)
                    return Scans[low - 1].OneBasedScanNumber;
            }

            // Scans can share a retention time. The linear version kept walking while the distance stayed
            // equal, so it came to rest on the last scan of such a run; binary search lands on the first.
            // A second search for the end of the run keeps the answer identical to the previous
            // implementation without reintroducing a linear step: the run is unbounded, because MsAlign
            // assigns DefaultErrorValue to every scan whose retention time field fails to parse, which makes
            // the whole file one run.
            //
            // The exact == is deliberate. The walk being reproduced continued only while the distance was
            // not strictly greater than the best so far, which is itself an exact-equality test, so a
            // tolerance here would merge scans a few ULPs apart that the old code kept separate.
            double closestRetentionTime = Scans[low].RetentionTime;
            int runEnd = low;
            int runHigh = NumSpectra - 1;
            while (runEnd < runHigh)
            {
                int mid = runEnd + (runHigh - runEnd + 1) / 2;
                if (Scans[mid].RetentionTime == closestRetentionTime)
                    runEnd = mid;
                else
                    runHigh = mid - 1;
            }

            return Scans[runEnd].OneBasedScanNumber;
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