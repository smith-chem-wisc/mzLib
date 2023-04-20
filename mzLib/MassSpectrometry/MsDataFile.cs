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

using Chemistry;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    // TODO: Define scope of class 
    // Class scope is to provide to the data loaded from the DataFile. 

    /// <summary>
    /// A class for interacting with data collected from a Mass Spectrometer, and stored in a file
    /// </summary>
    public abstract class MsDataFile
    {
        public MsDataScan[] Scans { get; set; }
        public SourceFile SourceFile { get; set; }
        public int NumSpectra => Scans.Length;
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
                if (scan != null && scan.MsnOrder == 1)
                {
                    yield return scan;
                }
            }
        }

        public virtual MsDataScan GetOneBasedScan(int scanNumber)
        {
            if (!CheckIfScansLoaded())
                LoadAllStaticData();
            return Scans.SingleOrDefault(i => i.OneBasedScanNumber == scanNumber);
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

        public virtual IEnumerator<MsDataScan> GetEnumerator()
        {
            return GetMsScansInIndexRange(1, NumSpectra).GetEnumerator();
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
    }
}