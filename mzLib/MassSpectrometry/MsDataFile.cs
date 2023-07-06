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
using System.Collections.Generic;
using System.Collections.Immutable;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

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
            {
                LoadAllStaticData();
            }
            return Scans;
        }

        public virtual List<MsDataScan> GetAllScansList()
        {
            if (!CheckIfScansLoaded())
            {
                LoadAllStaticData();
            }

            return Scans.ToList();
        }

        public virtual IEnumerable<MsDataScan> GetMS1Scans()
        {
            if (!CheckIfScansLoaded())
            {
                LoadAllStaticData();
            }

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
            {
                LoadAllStaticData();
            }

            return Scans.SingleOrDefault(i => i.OneBasedScanNumber == scanNumber);
        }

        public virtual IEnumerable<MsDataScan> GetMsScansInIndexRange(int firstSpectrumNumber, int lastSpectrumNumber)
        {
            if (!CheckIfScansLoaded())
            {
                LoadAllStaticData();
            }

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
            {
                LoadAllStaticData();
            }

            int oneBasedSpectrumNumber = GetClosestOneBasedSpectrumNumber(firstRT);
            
            while (oneBasedSpectrumNumber <= NumSpectra)
            {
                MsDataScan scan = GetOneBasedScan(oneBasedSpectrumNumber);

                if (scan == null)
                {
                    break;
                }

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

        /// <summary>
        /// 
        /// Performs a Binary Search and returns the closest Spectrum Number to the given Retention Time (RT).
        /// If RT is negative or bigger than last RT recorded in file, it will return 0.
        /// 
        /// </summary>
        /// <param name="retentionTime"></param>
        /// <returns></returns>
        public virtual int GetClosestOneBasedSpectrumNumber(double retentionTime)
        {
            if (!CheckIfScansLoaded())
            {
                LoadAllStaticData();
            }

            ImmutableSortedDictionary<int, double> scansAndRetentionTime =
                Scans.ToImmutableSortedDictionary(x => x.OneBasedScanNumber, x => x.RetentionTime);


            if (retentionTime <= scansAndRetentionTime.ElementAt(scansAndRetentionTime.Keys.Count() - 1).Key && retentionTime >= 0)
            {
                int search = Array.BinarySearch(scansAndRetentionTime.Values.ToArray(), index: 0, length: scansAndRetentionTime.Values.Count(), retentionTime);
                int indexFromSearch = ~search;

                if (search < 0)
                {
                    if (indexFromSearch < scansAndRetentionTime.Keys.Count())
                    {
                        return scansAndRetentionTime.ElementAt(indexFromSearch).Key;
                    }

                    if (indexFromSearch == scansAndRetentionTime.Keys.Count())
                    {
                        return scansAndRetentionTime.ElementAt(indexFromSearch - 1).Key;
                    }
                }
                else if (search < scansAndRetentionTime.Keys.Count())
                {
                    return scansAndRetentionTime.ElementAt(search).Key;
                }
            }

            return 0;
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