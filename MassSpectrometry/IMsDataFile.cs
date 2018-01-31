// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (IMsDataFile.cs) is part of MassSpectrometry.
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

using System;
using System.Collections.Generic;

namespace MassSpectrometry
{
    public interface IMsDataFile<out TScan> : IEnumerable<TScan>
        where TScan : IMsDataScan<IMzSpectrum<IMzPeak>>
    {
        #region Public Properties

        int NumSpectra { get; }
        SourceFile SourceFile { get; }

        #endregion Public Properties

        #region Public Methods

        TScan GetOneBasedScan(int oneBasedScanNumber);

        IEnumerable<TScan> GetMS1Scans();

        IEnumerable<DeconvolutionFeatureWithMassesAndScans> Deconvolute(int? minScan, int? maxScan, int minAssumedChargeState, int maxAssumedChargeState, double deconvolutionTolerancePpm, double intensityRatio, double aggregationTolerancePpm, Func<TScan, bool> scanFilterFunc);

        IEnumerable<TScan> GetMsScansInTimeRange(double firstRT, double lastRT);

        int GetClosestOneBasedSpectrumNumber(double retentionTime);

        #endregion Public Methods
    }
}