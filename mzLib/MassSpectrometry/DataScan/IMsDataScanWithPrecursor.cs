// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (IMsDataScan.cs) is part of MassSpectrometry.
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

namespace MassSpectrometry
{
    public interface IMsDataScanWithPrecursor<out TSpectrum> : IMsDataScan<TSpectrum>
        where TSpectrum : IMzSpectrum<IMzPeak>
    {

        #region Public Properties

        int OneBasedPrecursorScanNumber { get; }
        int? SelectedIonGuessChargeStateGuess { get; }
        double? SelectedIonGuessMonoisotopicIntensity { get; }
        double? SelectedIonGuessMonoisotopicMZ { get; }
        double? SelectedIonGuessIntensity { get; }
        double? SelectedIonGuessMZ { get; }
        DissociationType DissociationType { get; }
        double? IsolationWidth { get; }
        double IsolationMz { get; }
        MzRange IsolationRange { get; }

        #endregion Public Properties

        #region Public Methods

        void RecomputeChargeState(IMzSpectrum<IMzPeak> precursorSpectrum, double tolHere, int maxCharge);

        void RecomputeSelectedPeak(IMzSpectrum<IMzPeak> precursorSpectrum);

        void RecomputeMonoisotopicPeak(IMzSpectrum<IMzPeak> precursorSpectrum, double tolHere, double intensityFractionNeeded);

        void TranformByApplyingFunctionsToSpectraAndReplacingPrecursorMZs(Func<IMzPeak, double> convertorForSpectrum, double selectedIonGuessMZ, double selectedIonGuessMonoisotopicMZ);

        #endregion Public Methods

    }
}