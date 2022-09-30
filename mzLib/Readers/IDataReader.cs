using MzLibUtil;
using MassSpectrometry; 
namespace Readers
{
    // only loads file.
    public interface IDataReader
    {
        public string FilePath { get; }
        public void LoadAllStaticData(FilteringParams filteringParams = null, int maxThreads = 1);

        public SourceFile GetSourceFile(); 
    }

    public static class MsDataFileHelpers
    {
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
    }
}