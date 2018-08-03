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

using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    /// <summary>
    /// A class for interacting with data collected from a Mass Spectrometer, and stored in a file
    /// </summary>
    public class MsDataFile
    {
        protected MsDataScan[] Scans;

        public MsDataFile(int numSpectra, SourceFile sourceFile)
        {
            Scans = new MsDataScan[numSpectra];
            SourceFile = sourceFile;
        }

        public MsDataFile(MsDataScan[] scans, SourceFile sourceFile)
        {
            Scans = scans;
            SourceFile = sourceFile;
        }

        public SourceFile SourceFile { get; }

        public int NumSpectra
        {
            get
            {
                return Scans.Length;
            }
        }

        public virtual List<MsDataScan> GetAllScansList()
        {
            return Scans.ToList();
        }

        public static int TopNpeakHelper(ref double[] intensities, ref double[] mArray, IFilteringParams filteringParams)
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

        /// <summary>
        /// This method is designed to break a scan up into windows and take the top N peaks (by intensity)
        /// from each window, then merge the results as the scan's new mass spectrum
        /// </summary>
        /// <param name="intensities"></param>
        /// <param name="mArray"></param>
        /// <param name="filteringParams"></param>
        public static void WindowModeHelper(ref double[] intensities, ref double[] mArray, IFilteringParams filteringParams)
        {
            List<double> mzResults = new List<double>();
            List<double> intensityResults = new List<double>();

            int windowSize = intensities.Length / filteringParams.NumberOfWindows.Value;

            // window must always have at least one peak in it
            if (windowSize < 1)
            {
                windowSize = 1;
            }

            int windowPeakIndexMinimum = 0;
            int windowPeakIndexMaximum = windowSize - 1;

            // this loop breaks a scan up into "windows" (e.g., a scan with 100 peaks
            // divided into 10 windows would have 10 peaks per window) and takes the top N peaks per window.
            // the results of each trimmed window are concatenated into mzResults and intensityResults
            for (int i = 0; i < filteringParams.NumberOfWindows; i++)
            {
                // make the last window end at the end of the spectrum
                // this is to prevent rounding errors in windowSize from cutting off the end of the spectrum
                if (i == filteringParams.NumberOfWindows - 1)
                {
                    windowPeakIndexMaximum = intensities.Length - 1;
                }

                // avoid index out of range problems
                if (windowPeakIndexMinimum >= intensities.Length)
                {
                    break;
                }

                // determine the valid peaks given filtering conditions for this window
                var windowIntensities = new double[windowPeakIndexMaximum - windowPeakIndexMinimum + 1];
                Array.Copy(intensities, windowPeakIndexMinimum, windowIntensities, 0, windowIntensities.Length);

                var windowMzs = new double[windowPeakIndexMaximum - windowPeakIndexMinimum + 1];
                Array.Copy(mArray, windowPeakIndexMinimum, windowMzs, 0, windowMzs.Length);

                int numPeaksToTakeInThisWindow = TopNpeakHelper(ref windowIntensities, ref windowMzs, filteringParams);

                // merge results of this window into the global results
                intensityResults.AddRange(windowIntensities.Take(numPeaksToTakeInThisWindow));
                mzResults.AddRange(windowMzs.Take(numPeaksToTakeInThisWindow));

                // set up for the next window
                windowPeakIndexMinimum = windowPeakIndexMaximum + 1;
                windowPeakIndexMaximum = windowPeakIndexMinimum + windowSize;
            }

            // convert merged results to array and sort by m/z
            mArray = mzResults.ToArray();
            intensities = intensityResults.ToArray();

            Array.Sort(mArray, intensities);
        }

        public virtual IEnumerable<MsDataScan> GetMS1Scans()
        {
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
            return Scans[scanNumber - 1];
        }

        public IEnumerable<MsDataScan> GetMsScansInIndexRange(int FirstSpectrumNumber, int LastSpectrumNumber)
        {
            for (int oneBasedSpectrumNumber = FirstSpectrumNumber; oneBasedSpectrumNumber <= LastSpectrumNumber; oneBasedSpectrumNumber++)
            {
                yield return GetOneBasedScan(oneBasedSpectrumNumber);
            }
        }

        public IEnumerable<MsDataScan> GetMsScansInTimeRange(double firstRT, double lastRT)
        {
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

        public IEnumerator<MsDataScan> GetEnumerator()
        {
            return GetMsScansInIndexRange(1, NumSpectra).GetEnumerator();
        }

        public IEnumerable<DeconvolutionFeatureWithMassesAndScans> Deconvolute(int? minScan, int? maxScan, int minAssumedChargeState, int maxAssumedChargeState, double deconvolutionTolerancePpm, double intensityRatioLimit, double aggregationTolerancePpm, Func<MsDataScan, bool> scanFilterFunc, int maxThreads = -1)
        {
            minScan = minScan ?? 1;
            maxScan = maxScan ?? NumSpectra;

            var allAggregateGroups = new List<IsotopicEnvelope>[maxScan.Value - minScan.Value + 1];
            Parallel.ForEach(Partitioner.Create(minScan.Value, maxScan.Value + 1), new ParallelOptions { MaxDegreeOfParallelism = maxThreads }, fff =>
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

        protected class ReverseComparer : IComparer<double>
        {
            public int Compare(double x, double y)
            {
                return y.CompareTo(x);
            }
        }
    }
}