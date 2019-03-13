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

namespace MassSpectrometry
{
    public class FilteringParams : IFilteringParams
    {
        //Num: the number of windows used to filer; testSize: for comparing the amount of topN is used on
        public FilteringParams(int? numberOfPeaksToKeepPerWindow = null, double? minimumAllowedIntensityRatioToBasePeak = null, double? windowWidthDaltons = 0, int? numberOfWindows = 0, double? windowMaxNormalizationValue = null, bool applyTrimmingToMs1 = true, bool applyTrimmingToMsMs = true)
        {
            NumberOfPeaksToKeepPerWindow = numberOfPeaksToKeepPerWindow;
            MinimumAllowedIntensityRatioToBasePeakM = minimumAllowedIntensityRatioToBasePeak;

            //one can define only one, either window width or window number, not both. window width takes precendent
            if(windowWidthDaltons != null && windowWidthDaltons > 0)
            {
                WindowWidthDaltons = windowWidthDaltons.Value;
                NumberOfWindows = 0;
            }
            else if (numberOfWindows != null && numberOfWindows > 0)
            {
                WindowWidthDaltons = 0;
                NumberOfWindows = numberOfWindows.Value;
            }
            WindowMaxNormalizationValue = windowMaxNormalizationValue;
            ApplyTrimmingToMs1 = applyTrimmingToMs1;
            ApplyTrimmingToMsMs = applyTrimmingToMsMs;
        }

        public double? MinimumAllowedIntensityRatioToBasePeakM { get; }
        public int? NumberOfPeaksToKeepPerWindow { get; }
        public double? WindowWidthDaltons { get; }
        public int? NumberOfWindows { get; }
        public double? WindowMaxNormalizationValue { get; }
        public bool ApplyTrimmingToMs1 { get; }
        public bool ApplyTrimmingToMsMs { get; }


    }
}