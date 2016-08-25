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

using Spectra;
using System;

namespace MassSpectrometry
{

    public interface IMsDataScan<out TSpectrum>
        where TSpectrum : IMzSpectrum<MzPeak>
    {
        TSpectrum MassSpectrum { get; }
        int ScanNumber { get; }
        int MsnOrder { get; }
        double RetentionTime { get; }
        MzRange ScanWindowRange { get; }
        string ScanFilter { get; }
        string id { get; }
        bool isCentroid { get; }
        double InjectionTime { get; }
        double TotalIonCurrent { get; }
        Polarity Polarity { get; }
        MZAnalyzerType MzAnalyzer { get; }
        bool TryGetPrecursorScanNumber(out int precursorScanNumber);
        bool TryGetPrecursorID(out string PrecursorID);
        bool TryGetSelectedIonGuessChargeStateGuess(out int SelectedIonGuessChargeStateGuess);
        bool TryGetSelectedIonGuessMonoisotopicIntensity(out double SelectedIonGuessMonoisotopicIntensity);
        bool TryGetSelectedIonGuessMonoisotopicMZ(out double SelectedIonGuessMonoisotopicMZ);
        bool TryGetSelectedIonGuessIntensity(out double SelectedIonGuessIntensity);
        bool TryGetSelectedIonGuessMZ(out double SelectedIonGuessMZ);
        bool TryGetDissociationType(out DissociationType DissociationType);
        bool TryGetIsolationWidth(out double IsolationWidth);
        bool TryGetIsolationMZ(out double IsolationMZ);
        bool TryGetIsolationRange(out MzRange IsolationRange);
        void tranformByApplyingFunctionsToSpectraAndReplacingPrecursorMZs(Func<MzPeak, double> convertorForSpectrum, double newPrecursorMZ, double selectedIonGuessMonoisotopicMZ);
        void attemptToRefinePrecursorMonoisotopicPeak(ISpectrum<MzPeak> ms1Spectrum, double worstA = 0.0005, double worstA2 = 0.0005, double distBetweenPeaks = 1.003, double factorOfErrorAllowed = 3, double intensityDecreaseAllowed = 0.2);
    }
}