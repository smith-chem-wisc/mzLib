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
    public interface IFilteringParams
    {
        double? MinimumAllowedIntensityRatioToBasePeakM { get; }
        int? NumberOfPeaksToKeepPerWindow { get; }
        double? WindowWidthThomsons { get; }
        int? NumberOfWindows { get; }
        bool NormalizePeaksAcrossAllWindows { get; }
        bool ApplyTrimmingToMs1 { get; }
        bool ApplyTrimmingToMsMs { get; }
    }
}