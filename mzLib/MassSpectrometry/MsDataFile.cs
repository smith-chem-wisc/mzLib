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
using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    // TODO: Define scope of class 
    // Class scope is to provide to the data loaded from the DataFile. 

    /// <summary>
    /// A class for interacting with data collected from a Mass Spectrometer, and stored in a file
    /// </summary>
    public abstract class MsDataFile : IEnumerable<MsDataScan>
    {
        public MsDataScan[] Scans { get; protected set; }
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
           
            return Scans.SingleOrDefault(i => i.OneBasedScanNumber == scanNumber);
        }

        public virtual IEnumerable<MsDataScan> GetMsScansInIndexRange(int firstSpectrumNumber, int lastSpectrumNumber)
        {
            if (!CheckIfScansLoaded())
                LoadAllStaticData();
                
            for (int oneBasedSpectrumNumber = firstSpectrumNumber;
                 oneBasedSpectrumNumber <= lastSpectrumNumber; oneBasedSpectrumNumber++)
                yield return GetOneBasedScan(oneBasedSpectrumNumber);
        }

        public virtual IEnumerable<MsDataScan> GetMsScansInTimeRange(double firstRT, double lastRT)
        {
            if (!CheckIfScansLoaded())
                LoadAllStaticData();
                
            int oneBasedSpectrumNumber = GetClosestOneBasedSpectrumNumber(firstRT) + 1;

            while (oneBasedSpectrumNumber < NumSpectra + 1)
            {
                MsDataScan scan = GetOneBasedScan(oneBasedSpectrumNumber);

                double rt = scan.RetentionTime;

                if (rt < firstRT)
                {
                    oneBasedSpectrumNumber++;
                    continue;
                }
                else if (rt > lastRT)
                    yield break;

                yield return scan;
                oneBasedSpectrumNumber++;
            }
        }

        /// <summary>
        /// Returns index of the scan with the closest retention time to the target retention time.
        /// If retention time target is out of range, it returns the
        /// scan with the closest retention time to target. 
        ///
        /// For example: 
        /// if fileRetentionTimeRange = [2,3,4] and targetRetentionTime = 1, returns 0 and
        /// if fileRetentionTimeRange = [2,3,4] and targetRetentionTime = 5, returns 2.
        /// </summary>
        /// <param name="retentionTime"></param>
        /// <returns></returns>
        public virtual int GetClosestOneBasedSpectrumNumber(double retentionTime)
        {
            if (!CheckIfScansLoaded())
                LoadAllStaticData();

            return ClassExtensions.GetClosestIndex(Scans
                .Select(scan => scan.RetentionTime)
                .ToArray(), retentionTime);
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