using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry;

namespace MassSpectrometry
{
    public static class WindowModeHelper
    {
        public static void Run(ref double[] intensities, ref double[] mArray,
            IFilteringParams filteringParams, double scanRangeMinMz, double scanRangeMaxMz, bool keepZeroPeaks = false)
        {
            // Backward-compatible thin wrapper: callers without a per-peak charge array
            // continue to use the original 5-parameter signature unchanged. The 6-parameter
            // overload below carries an optional `ref int[] chargeArray` so reader paths
            // that loaded a per-peak charge array (PSI-MS MS:1000516) can keep it in sync
            // with the windowing/filtering's reordering and pruning of mArray/intensities.
            int[] noChargeArray = null;
            Run(ref intensities, ref mArray, ref noChargeArray, filteringParams,
                scanRangeMinMz, scanRangeMaxMz, keepZeroPeaks);
        }

        /// <summary>
        /// Charge-array-aware overload. When <paramref name="chargeArray"/> is non-null and
        /// its length matches <paramref name="mArray"/>, the same selection / reordering /
        /// pruning that the windowing applies to mz+intensity is also applied to charges,
        /// preserving per-peak alignment. Pass <c>null</c> (or a length-mismatched array)
        /// to get the original behavior unchanged.
        /// </summary>
        public static void Run(ref double[] intensities, ref double[] mArray, ref int[] chargeArray,
            IFilteringParams filteringParams, double scanRangeMinMz, double scanRangeMaxMz, bool keepZeroPeaks = false)
        {
            // Charges only follow the selection if the caller provided a length-aligned
            // array. Mismatched lengths are ambiguous and would corrupt alignment further;
            // drop the array in that case so the caller sees an honest "no charge data".
            bool hasCharges = chargeArray != null && chargeArray.Length == mArray.Length;
            if (chargeArray != null && !hasCharges)
            {
                chargeArray = null;
            }

            if (hasCharges)
            {
                // Index-permutation sort: Array.Sort handles 2 parallel arrays at most;
                // for 3 parallel arrays we sort an index array by the intensity key and
                // then materialize all three in the new order. Locals are needed because
                // C# disallows capturing `ref` parameters inside the lambda comparator.
                double[] intLocal = intensities;
                double[] mzLocal = mArray;
                int[] chgLocal = chargeArray;
                int n = intLocal.Length;
                int[] perm = new int[n];
                for (int i = 0; i < n; i++) perm[i] = i;
                Array.Sort(perm, (a, b) => intLocal[a].CompareTo(intLocal[b]));
                var newInt = new double[n];
                var newMz = new double[n];
                var newChg = new int[n];
                for (int i = 0; i < n; i++)
                {
                    newInt[i] = intLocal[perm[i]];
                    newMz[i] = mzLocal[perm[i]];
                    newChg[i] = chgLocal[perm[i]];
                }
                intensities = newInt;
                mArray = newMz;
                chargeArray = newChg;
            }
            else
            {
                Array.Sort(intensities, mArray);
            }
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
            // Parallel reduced charge list — only populated when the caller passed a
            // length-aligned chargeArray. Indices into mArray/intensities are stored in
            // mzInRange[rangeIndex], so charges follow the same selection by indexing
            // into the (already-permuted) chargeArray with arrayIndex.
            List<int> reducedChargeList = hasCharges ? new List<int>() : null;

            foreach (int rangeIndex in mzInRange.Keys)
            {
                List<double> tempMzList = new List<double>();
                List<int> tempChargeList = hasCharges ? new List<int>() : null;
                foreach (int arrayIndex in mzInRange[rangeIndex])
                {
                    tempMzList.Add(mArray[arrayIndex]);
                    if (hasCharges) tempChargeList.Add(chargeArray[arrayIndex]);
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
                    int takeCount = Math.Min(tempMzList.Count, countOfPeaksToKeepPerWindow);
                    reducedMzList.AddRange(tempMzList.GetRange(0, takeCount));
                    reducedIntensityList.AddRange(mzRangeIntensities[rangeIndex].GetRange(0, Math.Min(mzRangeIntensities[rangeIndex].Count, countOfPeaksToKeepPerWindow)));
                    if (hasCharges)
                        reducedChargeList.AddRange(tempChargeList.GetRange(0, takeCount));
                }
            }

            intensities = reducedIntensityList.ToArray();
            mArray = reducedMzList.ToArray();
            if (hasCharges)
            {
                // Final m/z sort must permute charges in lockstep — repeat the
                // index-permutation idiom from the top of the method. Locals are needed
                // because C# disallows capturing `ref` parameters inside the lambda.
                chargeArray = reducedChargeList.ToArray();
                double[] mzLocal = mArray;
                double[] intLocal = intensities;
                int[] chgLocal = chargeArray;
                int n = mzLocal.Length;
                int[] perm = new int[n];
                for (int i = 0; i < n; i++) perm[i] = i;
                Array.Sort(perm, (a, b) => mzLocal[a].CompareTo(mzLocal[b]));
                var newMz = new double[n];
                var newInt = new double[n];
                var newChg = new int[n];
                for (int i = 0; i < n; i++)
                {
                    newMz[i] = mzLocal[perm[i]];
                    newInt[i] = intLocal[perm[i]];
                    newChg[i] = chgLocal[perm[i]];
                }
                mArray = newMz;
                intensities = newInt;
                chargeArray = newChg;
            }
            else
            {
                Array.Sort(mArray, intensities);
            }
        }
    }
}
