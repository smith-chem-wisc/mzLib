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

        /// <summary>
        /// This method is designed to break a scan up into windows and take the top N peaks (by intensity)
        /// from each window, then merge the results as the scan's new mass spectrum
        /// </summary>
        /// <param name="intensities"></param>
        /// <param name="mArray"></param>
        /// <param name="filteringParams"></param>
        public static void WindowModeHelper(ref double[] intensities, ref double[] mArray, IFilteringParams filteringParams, double scanRangeMinMz, double scanRangeMaxMz, double? WindowMaxNormalizationToValue = null, bool keepZeroPeaks = false)
        {
            Array.Sort(intensities, mArray);

            //filter low intensites based on a percent for the whole spectrum.
            if (filteringParams.MinimumAllowedIntensityRatioToBasePeakM.HasValue)
            {
                double maxIntensity = intensities.Max();
                for (int i = 0; i < intensities.Length; i++)
                {
                    if (intensities[i] <= maxIntensity * filteringParams.MinimumAllowedIntensityRatioToBasePeakM.Value)
                    {
                        intensities[i] = 0;
                    }
                }
            }

            const double shiftToMakeRangeInclusive = 0.000000001;

            int numberOfWindows = filteringParams.NumberOfWindows ?? 1;

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

            Dictionary<int, List<int>> mzInRange = new Dictionary<int, List<int>>(); //index of range and  list of index values in mArray
            for (int i = 0; i < ranges.Count; i++)
            {
                mzInRange.Add(i, new List<int>());
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

            List<double> reducedMzList = new List<double>();
            List<double> reducedIntensityList = new List<double>();

            foreach (int rangeIndex in mzInRange.Keys)
            {
                List<double> tempIntList = new List<double>();
                List<double> tempMzList = new List<double>();
                foreach (int arrayIndex in mzInRange[rangeIndex])
                {
                    tempIntList.Add(intensities[arrayIndex]);
                    tempMzList.Add(mArray[arrayIndex]);
                }
                if (WindowMaxNormalizationToValue.HasValue)
                {
                    double max = tempIntList.Max();
                    tempIntList = tempIntList.Select(x => x / max * WindowMaxNormalizationToValue.Value).ToList();
                }
                if (tempMzList.Count > 0 && tempIntList.Count > 0)
                {
                    reducedMzList.AddRange(tempMzList.GetRange(0, Math.Min(tempMzList.Count, countOfPeaksToKeepPerWindow)));
                    reducedIntensityList.AddRange(tempIntList.GetRange(0, Math.Min(tempIntList.Count, countOfPeaksToKeepPerWindow)));
                }
            }
            intensities = reducedIntensityList.ToArray();
            mArray = reducedMzList.ToArray();
            Array.Sort(mArray, intensities);
        }

        public static void XCorrPrePreprocessing(ref double[] intensities, ref double[] mArray, double scanRangeMinMz, double scanRangeMaxMz, double precursorMz, double precursorDiscardRange = 1.5, double discreteMassBin = 1.0005079, double minimumAllowedIntensityRatioToBasePeak = 0.05)
        {
            //The discrete bin value 1.0005079 was from J. Proteome Res., 2018, 17 (11), pp 3644–3656

            Array.Sort(mArray, intensities);
            int numberOfNominalMasses = (int)Math.Round((scanRangeMaxMz - scanRangeMinMz) / 1.0005079, 0) + 1;

            double[] genericIntensityArray = new double[numberOfNominalMasses];
            double[] genericMzArray = new double[numberOfNominalMasses];

            double lowestMz = Math.Round(scanRangeMinMz / discreteMassBin, 0) * discreteMassBin;

            for (int i = 0; i < numberOfNominalMasses; i++)
            {
                genericMzArray[i] = i * discreteMassBin + lowestMz;
            }

            int lowTheortetical = (int)Math.Round((precursorMz - precursorDiscardRange) / 1.0005079, 0);
            int hiTheortetical = (int)Math.Round((precursorMz + precursorDiscardRange) / 1.0005079, 0);

            //this leaves you with one intensity per nominal mass, even if they come in as more than one. Intensity is Square-rooted
            for (int i = 0; i < mArray.Length; i++)
            {
                //the nominalMass is really just the integer index
                int nominalMass = (int)Math.Round((mArray[i] - scanRangeMinMz) / 1.0005079, 0);

                //this might have be be exclusive (i.e. get rid of the = sign) should eliminate unfragmented precursors
                if (nominalMass < numberOfNominalMasses && (mArray[i] <= lowTheortetical || mArray[i] >= hiTheortetical))
                {
                    genericIntensityArray[nominalMass] = Math.Max(genericIntensityArray[nominalMass], Math.Sqrt(intensities[i]));
                }
            }

            //we've already filtered for when multiple mzs appear in a single nominal mass bin
            FilteringParams secondFilter = new FilteringParams(null, minimumAllowedIntensityRatioToBasePeak, 10, false, false);

            MsDataFile.WindowModeHelper(ref genericIntensityArray, ref genericMzArray, secondFilter, genericMzArray.Min(), genericMzArray.Max(), 50, true);

            Array.Sort(genericMzArray, genericIntensityArray);

            //Scale the intensities

            const int rangeEnd = 75; //from J. Proteome Res., 2018, 17 (11), pp 3644–3656

            double[] scaledIntensities = new double[genericIntensityArray.Length];
            for (int i = 0; i < genericIntensityArray.Length; i++)
            {
                double scaleValue = 0;

                int low = Math.Max(0, i - rangeEnd);
                int high = Math.Min(genericIntensityArray.Length - 1, i + rangeEnd);
                int denominator = high - low + 1;

                for (int j = low; j <= high; j++)
                {
                    scaleValue += genericIntensityArray[j];
                }
                scaledIntensities[i] = Math.Max(0, genericIntensityArray[i] - 1d / (denominator) * scaleValue);
            }
            genericIntensityArray = scaledIntensities;

            intensities = genericIntensityArray;
            mArray = genericMzArray;
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