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

        /// <summary>
        /// This method is designed to break a scan up into windows and take the top N peaks (by intensity)
        /// from each window, then merge the results as the scan's new mass spectrum
        /// </summary>
        /// <param name="intensities"></param>
        /// <param name="mArray"></param>
        /// <param name="filteringParams"></param>
        public static void WindowModeHelper(ref double[] intensities, ref double[] mArray, IFilteringParams filteringParams, double scanRangeMinMz, double scanRangeMaxMz, bool keepZeroPeaks = false)
        {
            Array.Sort(intensities, mArray);
            double TIC = intensities.Sum();

            //filter low intensites based on a percent for the whole spectrum.
            if (filteringParams.MinimumAllowedIntensityRatioToBasePeakM.HasValue)
            {
                double maxIntensity = intensities.Max();
                double cutOff = maxIntensity * filteringParams.MinimumAllowedIntensityRatioToBasePeakM.Value;
                for (int i = 0; i < intensities.Length; i++)
                {
                    if (intensities[i] <= cutOff && intensities[i] > 0)
                    {
                        intensities[i] = 0;
                    }
                }
            }

            const double shiftToMakeRangeInclusive = 0.000000001;

            Chemistry.ClassExtensions.TupleList<double, double> ranges = new Chemistry.ClassExtensions.TupleList<double, double>();

            if (filteringParams.WindowWidthThomsons != null && filteringParams.WindowWidthThomsons > 0)
            {
                double scanRangeToUse = Math.Min(scanRangeMaxMz - scanRangeMinMz, filteringParams.WindowWidthThomsons.Value);

                List<double> ends = new List<double>();
                double end = 0;
                bool first = true;
                while (end < scanRangeMaxMz)
                {
                    if ((end + scanRangeToUse) > scanRangeMinMz)
                    {
                        if (first)
                        {
                            ends.Add(scanRangeMinMz - shiftToMakeRangeInclusive);
                            first = false;
                        }
                        else
                        {
                            ends.Add(end);
                        }
                    }
                    end += scanRangeToUse;
                }

                for (int i = 0; i < ends.Count; i++)
                {
                    if (i == 0)
                    {
                        ranges.Add(ends[i], ends[i + 1]);
                    }
                    else if (i != (ends.Count - 1))
                    {
                        ranges.Add(ends[i] + shiftToMakeRangeInclusive, ends[i + 1]);
                    }
                    else
                    {
                        ranges.Add(ends[i] + shiftToMakeRangeInclusive, scanRangeMaxMz + shiftToMakeRangeInclusive);
                    }
                }
            }
            else if (filteringParams.NumberOfWindows != null && filteringParams.NumberOfWindows > 1)
            {
                double mzRangeInOneWindow = (scanRangeMaxMz - scanRangeMinMz) / filteringParams.NumberOfWindows.Value;

                ranges.Add(scanRangeMinMz - shiftToMakeRangeInclusive, (scanRangeMinMz + mzRangeInOneWindow));
                scanRangeMinMz += mzRangeInOneWindow;

                for (int i = 2; i < filteringParams.NumberOfWindows; i++)
                {
                    ranges.Add(scanRangeMinMz, (scanRangeMinMz + mzRangeInOneWindow));
                    scanRangeMinMz += mzRangeInOneWindow;
                }
                ranges.Add(scanRangeMinMz, (scanRangeMinMz + mzRangeInOneWindow) + shiftToMakeRangeInclusive);
                scanRangeMinMz += mzRangeInOneWindow;
            }
            else
            {
                ranges.Add(scanRangeMinMz - shiftToMakeRangeInclusive, scanRangeMaxMz + shiftToMakeRangeInclusive);
            }

            Dictionary<int, List<int>> mzInRange = new Dictionary<int, List<int>>(); //index of range and  list of index values in mArray
            Dictionary<int, List<double>> mzRangeIntensities = new Dictionary<int, List<double>>(); //index of range and  list of index values in mArray
            for (int i = 0; i < ranges.Count; i++)
            {
                mzInRange.Add(i, new List<int>());
                mzRangeIntensities.Add(i, new List<double>());
            }

            //we're going to keep track of the array indicies b/c that's easier than the m/zs b/c rounding issues
            //we're only keeping peaks with intensity greater than 1 (assuming those <= 1 has "zero" intensity)
            for (int j = mArray.Length - 1; j >= 0; j--)
            {
                foreach (int rangeIndex in Enumerable.Range(0, ranges.Count))
                {
                    if (mArray[j] > ranges[rangeIndex].Item1 && mArray[j] <= ranges[rangeIndex].Item2 && (intensities[j] > 0.000000001 || keepZeroPeaks))
                    {
                        mzInRange[rangeIndex].Add(j);
                        break;
                    }
                }
            }

            int countOfPeaksToKeepPerWindow = filteringParams.NumberOfPeaksToKeepPerWindow ?? int.MaxValue;

            foreach (int rangeIndex in mzInRange.Keys)
            {
                List<double> tempIntList = new List<double>();
                foreach (int arrayIndex in mzInRange[rangeIndex])
                {
                    tempIntList.Add(intensities[arrayIndex]);
                }
                mzRangeIntensities[rangeIndex] = tempIntList;
            }

            int countOfRangesWithIntensities = 0;
            foreach (int range in mzRangeIntensities.Keys)
            {
                if (mzRangeIntensities[range].Sum() > 0)
                {
                    countOfRangesWithIntensities++;
                }
            }

            List<double> reducedMzList = new List<double>();
            List<double> reducedIntensityList = new List<double>();

            foreach (int rangeIndex in mzInRange.Keys)
            {
                List<double> tempMzList = new List<double>();
                foreach (int arrayIndex in mzInRange[rangeIndex])
                {
                    tempMzList.Add(mArray[arrayIndex]);
                }
                //There is no need to do any normalization unless there are multiple windows
                if (filteringParams.NormalizePeaksAcrossAllWindows)
                {

                    double max = mzRangeIntensities[rangeIndex].Max();
                    if (max == 0)
                    {
                        max = 1;
                    }
                    mzRangeIntensities[rangeIndex] = mzRangeIntensities[rangeIndex].Select(x => x / max * 50.0000).ToList();

                    //Saving b/c I might want to put it back.

                    //double sum = mzRangeIntensities[rangeIndex].Sum();
                    //if (sum > 0)
                    //{
                    //    double normalizationFactor = TIC / sum / countOfRangesWithIntensities;
                    //    mzRangeIntensities[rangeIndex] = mzRangeIntensities[rangeIndex].Select(x => x * normalizationFactor).ToList();
                    //}
                }

                if (tempMzList.Count > 0 && mzRangeIntensities[rangeIndex].Count > 0)
                {
                    reducedMzList.AddRange(tempMzList.GetRange(0, Math.Min(tempMzList.Count, countOfPeaksToKeepPerWindow)));
                    reducedIntensityList.AddRange(mzRangeIntensities[rangeIndex].GetRange(0, Math.Min(mzRangeIntensities[rangeIndex].Count, countOfPeaksToKeepPerWindow)));
                }
            }

            intensities = reducedIntensityList.ToArray();
            mArray = reducedMzList.ToArray();
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
                    var mass = isotopicEnvelope.MonoisotopicMass;
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

        /// <summary>
        /// Extracts an ion chromatogram from the spectra file, given a mass, charge, retention time, and mass tolerance.
        /// </summary>
        public ExtractedIonChromatogram ExtractIonChromatogram(double neutralMass, int charge, Tolerance massTolerance,
            double retentionTimeInMinutes, int msOrder = 1, double retentionTimeWindowWidthInMinutes = 5)
        {
            double theorMz = neutralMass.ToMz(charge);
            double startRt = retentionTimeInMinutes - retentionTimeWindowWidthInMinutes / 2;
            double endRt = retentionTimeInMinutes + retentionTimeWindowWidthInMinutes / 2;
            List<Datum> xicData = new List<Datum>();

            IEnumerable<MsDataScan> scansInRtWindow = GetMsScansInTimeRange(startRt, endRt);

            foreach (MsDataScan scan in scansInRtWindow.Where(p => p.MsnOrder == msOrder))
            {
                int ind = scan.MassSpectrum.GetClosestPeakIndex(theorMz);

                double expMz = scan.MassSpectrum.XArray[ind];

                if (massTolerance.Within(expMz.ToMass(charge), neutralMass))
                {
                    xicData.Add(new Datum(scan.RetentionTime, scan.MassSpectrum.YArray[ind]));
                }
                else
                {
                    xicData.Add(new Datum(scan.RetentionTime, 0));
                }
            }

            return new ExtractedIonChromatogram(xicData);
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