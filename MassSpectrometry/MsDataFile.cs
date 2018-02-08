﻿// Copyright 2012, 2013, 2014 Derek J. Bailey
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

using MzLibUtil;
using System;
using System.Collections;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

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

        #region Protected Constructors

        protected MsDataFile(int numSpectra, SourceFile sourceFile) : this(sourceFile)
        {
            Scans = new TScan[numSpectra];
        }

        protected MsDataFile(TScan[] scans, SourceFile sourceFile) : this(sourceFile)
        {
            Scans = scans;
        }

        #endregion Protected Constructors

        #region Private Constructors

        private MsDataFile(SourceFile sourceFile)
        {
            this.SourceFile = sourceFile;
        }

        #endregion Private Constructors

        #region Public Properties

        public SourceFile SourceFile { get; }

        public int NumSpectra
        {
            get
            {
                return Scans.Length;
            }
        }

        #endregion Public Properties

        #region Public Methods

        public static int TopNpeakHelper(double[] intensities, double[] mArray, IFilteringParams filteringParams)
        {
            IComparer<double> c = new ReverseComparer();
            Array.Sort(intensities, mArray, c);

            int numPeaks = intensities.Length;
            if (filteringParams.MinimumAllowedIntensityRatioToBasePeakM.HasValue)
            {
                double minIntensity = filteringParams.MinimumAllowedIntensityRatioToBasePeakM.Value * intensities[0];
                numPeaks = Math.Min(intensities.Count(b => b >= minIntensity), numPeaks);
            }

            if (filteringParams.NumberOfPeaksToKeepPerWindow.HasValue)
                numPeaks = Math.Min(filteringParams.NumberOfPeaksToKeepPerWindow.Value, numPeaks);
            return numPeaks;
        }

        public static void WindowModeHelper(ref double[] intensities, ref double[] mArray, IFilteringParams filteringParams)
        {
            List<double> mzResults = new List<double>();
            List<double> intensityResults = new List<double>();
            //windowNum is the number of windows
            for (int i = 0; i < filteringParams.NumberOfWindows; i++)
            {
                int temp = (i == filteringParams.NumberOfWindows) ? intensities.Length - i * (intensities.Length / filteringParams.NumberOfWindows.Value) : intensities.Length / filteringParams.NumberOfWindows.Value;
                var mzTemp = new double[temp];
                var intensityTemp = new double[temp];

                Buffer.BlockCopy(mArray, sizeof(double) * temp * i, mzTemp, 0, sizeof(double) * temp);
                Buffer.BlockCopy(intensities, sizeof(double) * temp * i, intensityTemp, 0, sizeof(double) * temp);

                int numPeaks = TopNpeakHelper(intensities, mArray, filteringParams);
                Array.Resize(ref intensityTemp, numPeaks);
                Array.Resize(ref mzTemp, numPeaks);
                mzResults.AddRange(mzTemp);
                intensityResults.AddRange(intensityTemp);
            }
            mArray = mzResults.ToArray();
            intensities = intensityResults.ToArray();
        }

        public abstract IEnumerable<TScan> GetMS1Scans();

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

        public IEnumerable<DeconvolutionFeatureWithMassesAndScans> Deconvolute(int? minScan, int? maxScan, int minAssumedChargeState, int maxAssumedChargeState, double deconvolutionTolerancePpm, double intensityRatioLimit, double aggregationTolerancePpm, Func<TScan, bool> scanFilterFunc)
        {
            minScan = minScan ?? 1;
            maxScan = maxScan ?? NumSpectra;

            var allAggregateGroups = new List<IsotopicEnvelope>[maxScan.Value - minScan.Value + 1];
            Parallel.ForEach(Partitioner.Create(minScan.Value, maxScan.Value + 1), fff =>
            {
                for (int scanIndex = fff.Item1; scanIndex < fff.Item2; scanIndex++)
                {
                    var theScan = GetOneBasedScan(scanIndex);
                    if (scanFilterFunc(theScan))
                        allAggregateGroups[scanIndex - minScan.Value] = theScan.MassSpectrum.Deconvolute(new MzRange(0, double.PositiveInfinity), minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatioLimit).ToList();
                }
            });

            List<DeconvolutionFeatureWithMassesAndScans> currentListOfGroups = new List<DeconvolutionFeatureWithMassesAndScans>();
            for (int scanIndex = minScan.Value; scanIndex <= maxScan.Value; scanIndex++)
            {
                if (allAggregateGroups[scanIndex - minScan.Value] == null)
                    continue;
                foreach (var isotopicEnvelope in allAggregateGroups[scanIndex - minScan.Value])
                {
                    DeconvolutionFeatureWithMassesAndScans matchingGroup = null;
                    var mass = isotopicEnvelope.monoisotopicMass;
                    foreach (var possibleGroup in currentListOfGroups)
                    {
                        var possibleGroupMass = possibleGroup.Mass;
                        if (Math.Abs(mass - possibleGroupMass) / possibleGroupMass * 1e6 <= aggregationTolerancePpm ||
                            Math.Abs(mass + 1.002868314 - possibleGroupMass) / possibleGroupMass * 1e6 <= aggregationTolerancePpm ||
                            Math.Abs(mass + 2.005408917 - possibleGroupMass) / possibleGroupMass * 1e6 <= aggregationTolerancePpm ||
                            Math.Abs(mass + 3.007841294 - possibleGroupMass) / possibleGroupMass * 1e6 <= aggregationTolerancePpm ||
                            Math.Abs(mass - 1.002868314 - possibleGroupMass) / possibleGroupMass * 1e6 <= aggregationTolerancePpm ||
                            Math.Abs(mass - 2.005408917 - possibleGroupMass) / possibleGroupMass * 1e6 <= aggregationTolerancePpm ||
                            Math.Abs(mass - 3.007841294 - possibleGroupMass) / possibleGroupMass * 1e6 <= aggregationTolerancePpm)
                        {
                            matchingGroup = possibleGroup;
                            matchingGroup.AddEnvelope(isotopicEnvelope, scanIndex, GetOneBasedScan(scanIndex).RetentionTime);
                            break;
                        }
                    }

                    if (matchingGroup == null)
                    {
                        var newGroupScans = new DeconvolutionFeatureWithMassesAndScans();
                        newGroupScans.AddEnvelope(isotopicEnvelope, scanIndex, GetOneBasedScan(scanIndex).RetentionTime);
                        currentListOfGroups.Add(newGroupScans);
                    }
                }
                foreach (var ok in currentListOfGroups.Where(b => b.MaxScanIndex < scanIndex))
                    yield return ok;
                currentListOfGroups.RemoveAll(b => b.MaxScanIndex < scanIndex);
            }
            foreach (var ok in currentListOfGroups)
                yield return ok;
        }

        #endregion Public Methods

        #region Protected Classes

        protected class ReverseComparer : IComparer<double>
        {
            #region Public Methods

            public int Compare(double x, double y)
            {
                return y.CompareTo(x);
            }

            #endregion Public Methods
        }

        #endregion Protected Classes
    }
}