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
        public static void WindowModeHelper(ref double[] intensities, ref double[] mArray, IFilteringParams filteringParams, double scanRangeMinMz, double scanRangeMaxMz, double? WindowMaxNormalizationToValue = null)
        {
            Array.Sort(mArray, intensities);
            const double shiftToMakeRangeInclusive = 0.000000001;
            List<MzPeak> mzIntensites = new List<MzPeak>();
            List<MzPeak> mzIntensites_reduced = new List<MzPeak>();
            int numberOfWindows = 1;
            if (filteringParams.NumberOfWindows != null)
            {
                numberOfWindows = filteringParams.NumberOfWindows.Value;
            }

            //make MzPeaks of each mz/intensity pair and create a list of pairs for everything with intensity above the minimum cutoff
            double maximumIntensityInArray = intensities.Max();
            for (int i = 0; i < mArray.Length; i++)
            {
                if ((!filteringParams.MinimumAllowedIntensityRatioToBasePeakM.HasValue) || (intensities[i] > (maximumIntensityInArray * filteringParams.MinimumAllowedIntensityRatioToBasePeakM.Value)))
                {
                    mzIntensites.Add(new MzPeak(mArray[i], intensities[i]));
                }
            }

            double mzRangeInOneWindow = scanRangeMaxMz - scanRangeMinMz;
            Chemistry.ClassExtensions.TupleList<double, double> ranges = new Chemistry.ClassExtensions.TupleList<double, double>();

            double scanMin = scanRangeMinMz;
            if (filteringParams.NumberOfWindows.HasValue && filteringParams.NumberOfWindows.Value > 0)
            {
                mzRangeInOneWindow = mzRangeInOneWindow / filteringParams.NumberOfWindows.Value;

                for (int i = 1; i <= numberOfWindows; i++)
                {
                    if (i == 1) // first
                    {
                        ranges.Add(scanMin - shiftToMakeRangeInclusive, (scanMin + mzRangeInOneWindow));
                    }
                    else if (i == (numberOfWindows))//last
                    {
                        ranges.Add(scanMin, (scanMin + mzRangeInOneWindow) + shiftToMakeRangeInclusive);
                    }
                    else//middle
                    {
                        ranges.Add(scanMin, (scanMin + mzRangeInOneWindow));
                    }
                    scanMin += mzRangeInOneWindow;
                }
            }
            else
            {
                ranges.Add(scanRangeMinMz - shiftToMakeRangeInclusive, scanRangeMaxMz + shiftToMakeRangeInclusive);
            }

            foreach (Tuple<double, double> range in ranges)
            {
                List<MzPeak> handyLittleListofMzPeaksThatIsOnlyNeededTemporarily = new List<MzPeak>();
                foreach (MzPeak mzIntensityPair in mzIntensites)
                {
                    if (mzIntensityPair.Mz > range.Item1 && mzIntensityPair.Mz <= range.Item2)
                    {
                        handyLittleListofMzPeaksThatIsOnlyNeededTemporarily.Add(mzIntensityPair);
                    }
                }
                if (handyLittleListofMzPeaksThatIsOnlyNeededTemporarily.Count > 0)
                {
                    handyLittleListofMzPeaksThatIsOnlyNeededTemporarily.Sort((x, y) => y.Intensity.CompareTo(x.Intensity)); //sort tuple in place decending by item 2, reverse by changing x and y
                    double maxIntensity = handyLittleListofMzPeaksThatIsOnlyNeededTemporarily.Max(x => x.Intensity);
                    int countOfPeaksToKeep = handyLittleListofMzPeaksThatIsOnlyNeededTemporarily.Count;
                    if (filteringParams.NumberOfPeaksToKeepPerWindow != null)
                    {
                        countOfPeaksToKeep = Math.Min(handyLittleListofMzPeaksThatIsOnlyNeededTemporarily.Count, filteringParams.NumberOfPeaksToKeepPerWindow.Value);
                    }
                    for (int i = 0; i < countOfPeaksToKeep; i++)
                    {
                        if (WindowMaxNormalizationToValue.HasValue)
                        {
                            double normalizedIntensity = handyLittleListofMzPeaksThatIsOnlyNeededTemporarily[i].Intensity / maxIntensity * WindowMaxNormalizationToValue.Value;
                            mzIntensites_reduced.Add(new MzPeak(handyLittleListofMzPeaksThatIsOnlyNeededTemporarily[i].Mz, normalizedIntensity));
                        }
                        else
                        {
                            mzIntensites_reduced.Add(handyLittleListofMzPeaksThatIsOnlyNeededTemporarily[i]);
                        }
                    }
                }
            }

            // convert merged results to array and sort by m/z
            mzIntensites_reduced.Sort((x, y) => x.Mz.CompareTo(y.Mz));
            mArray = mzIntensites_reduced.Select(i => i.Mz).ToArray();
            intensities = mzIntensites_reduced.Select(i => i.Intensity).ToArray();
        }

        public static void XCorrPrePreprocessing(ref double[] intensities, ref double[] mArray, double scanRangeMinMz, double scanRangeMaxMz, double precursorMz, double precursorDiscardRange = 1.5, double discreteMassBin = 1.0005079, double percentMaxThreshold = 5)
        {
            //The discrete bin value 1.0005079 was from J. Proteome Res., 2018, 17 (11), pp 3644–3656

            Array.Sort(mArray, intensities);
            int numberOfWindows = (int)Math.Round((scanRangeMaxMz - scanRangeMinMz + discreteMassBin) / discreteMassBin, 0);

            FilteringParams firstFilter = new FilteringParams(1, percentMaxThreshold / 100, numberOfWindows, false, false);

            for (int i = 0; i < mArray.Length; i++)
            {
                if (Math.Abs(intensities[i]) > 0.000001)//some small number
                {
                    intensities[i] = Math.Sqrt(intensities[i]);
                }
            }

            //set intensity of precursor ions to zero
            for (int i = 0; i < mArray.Length; i++)
            {
                if (mArray[i] > (precursorMz - precursorDiscardRange) && mArray[i] < (precursorMz + precursorDiscardRange))
                {
                    intensities[i] = 0;
                }
            }

            //We filter twice below: the second filter, which is 10 windows with 10 peaks, you could end up with 10 peaks that are not all separated minimally by 1.0005079 daltons.T
            //The first pass gets one peak(maximally) at every dalton and the second pass keeps up to 10 of those in each of the 10 windows.

            //filter peaks into discrete mass bins. makes one bin for every single integer m/z.
            double minMz = scanRangeMinMz - discreteMassBin / 2;
            double maxMz = scanRangeMaxMz + discreteMassBin / 2;

            MsDataFile.WindowModeHelper(ref intensities, ref mArray, firstFilter, minMz, maxMz);

            //normalize intensity in segments to max 50 intensity keeping only 10 peaks per window.
            FilteringParams secondFilter = new FilteringParams((int)Math.Ceiling((decimal)numberOfWindows / 10), null, 10, false, false);
            MsDataFile.WindowModeHelper(ref intensities, ref mArray, secondFilter, minMz, maxMz, 50);

            for (int i = 0; i < mArray.Length; i++)
            {
                mArray[i] = Math.Round(mArray[i] / discreteMassBin, 0) * discreteMassBin;
            }

            //fill in missing mz values in array with zero intensity
            List<double> allPossibleMzValues = new List<double>();
            double start = Math.Round(scanRangeMinMz / discreteMassBin, 0) * discreteMassBin;//don't use minMz here, that would be WRONG
            while (start < (maxMz))
            {
                allPossibleMzValues.Add(start);
                start += discreteMassBin;
            }

            //remove any that are already in the mArray
            foreach (double remainingMzValue in mArray)
            {
                allPossibleMzValues.RemoveAll(x => x > (remainingMzValue - discreteMassBin / 2) && x < (remainingMzValue + discreteMassBin / 2));
            }

            double[] zeroArray = new double[allPossibleMzValues.Count];
            for (int i = 0; i < zeroArray.Length; i++)
            {
                zeroArray[i] = 0;
            }

            int original_mArray_count = mArray.Length;
            Array.Resize<double>(ref mArray, mArray.Length + allPossibleMzValues.Count);
            Array.Copy(allPossibleMzValues.ToArray(), 0, mArray, original_mArray_count, allPossibleMzValues.Count);

            Array.Resize<double>(ref intensities, intensities.Length + zeroArray.Length);
            Array.Copy(zeroArray.ToArray(), 0, intensities, original_mArray_count, zeroArray.Length);

            Array.Sort(mArray, intensities);

            //Scale the intensities

            const int rangeEnd = 75; //from J. Proteome Res., 2018, 17 (11), pp 3644–3656

            double[] scaledIntensities = new double[intensities.Length];
            for (int i = 0; i < intensities.Length; i++)
            {
                double scaleValue = 0;

                int low = Math.Max(0, i - rangeEnd);
                int high = Math.Min(intensities.Length - 1, i + rangeEnd);
                int denominator = high - low + 1;

                for (int j = low; j <= high; j++)
                {
                    scaleValue += intensities[j];
                }
                scaledIntensities[i] = Math.Max(0, intensities[i] - 1d / (denominator) * scaleValue);
            }
            intensities = scaledIntensities;
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