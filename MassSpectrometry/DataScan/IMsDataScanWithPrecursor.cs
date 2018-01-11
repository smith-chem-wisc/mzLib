﻿// Copyright 2012, 2013, 2014 Derek J. Bailey
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
using Spectra;
using System;
using System.Collections.Generic;

namespace MassSpectrometry
{
    public interface IMsDataScanWithPrecursor<out TSpectrum> : IMsDataScan<TSpectrum>
        where TSpectrum : IMzSpectrum<IMzPeak>
    {
        #region Public Properties

        int? OneBasedPrecursorScanNumber { get; }

        double SelectedIonMZ { get; }
        double? SelectedIonIntensity { get; }
        int? SelectedIonChargeStateGuess { get; }
        double? SelectedIonMonoisotopicGuessMz { get; }
        double? SelectedIonMonoisotopicGuessIntensity { get; }

        DissociationType DissociationType { get; }

        double? IsolationMz { get; }
        MzRange IsolationRange { get; }

        #endregion Public Properties

        #region Public Methods

        /// <summary>
        /// Use to set value of SelectedIonIntensity based on SelectedIonMZ
        /// </summary>
        /// <param name="precursorSpectrum"></param>
        void ComputeSelectedPeakIntensity(IMzSpectrum<IMzPeak> precursorSpectrum);

        /// <summary>
        /// Used to set refine value of SelectedIonMZ and SelectedIonIntensity
        /// </summary>
        /// <param name="precursorSpectrum"></param>
        void RefineSelectedMzAndIntensity(IMzSpectrum<IMzPeak> precursorSpectrum);

        IEnumerable<IsotopicEnvelope> GetIsolatedMassesAndCharges(IMzSpectrum<IMzPeak> precursorSpectrum, int minAssumedChargeState, int maxAssumedChargeState, double deconvolutionTolerancePpm, double intensityRatio);

        void ComputeMonoisotopicPeakIntensity(IMzSpectrum<IMzPeak> precursorSpectrum);

        void TransformMzs(Func<IPeak, double> convertorForSpectrum, Func<IPeak, double> convertorForPrecursor);

        #endregion Public Methods
    }
}