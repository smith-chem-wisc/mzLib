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
        #region Public Constructors

        //Num: the number of windows used to filer; testSize: for comparing the amount of topN is used on
        public FilteringParams(int? top = null, double? ratio = null, int? windowNum = null, bool trimMs1Peaks = true, bool trimMsMsPeaks = true)
        {
            this.topNpeaks = top;
            this.minRatio = ratio;
            this.windowNum = windowNum;
            this.trimMs1Peaks = trimMs1Peaks;
            this.trimMsMsPeaks = trimMsMsPeaks;
        }

        #endregion Public Constructors

        #region Public Properties

        public double? minRatio { get; }
        public int? topNpeaks { get; }
        public int? windowNum { get; }
        public bool trimMs1Peaks { get; }
        public bool trimMsMsPeaks { get; }

        #endregion Public Properties
    }
}