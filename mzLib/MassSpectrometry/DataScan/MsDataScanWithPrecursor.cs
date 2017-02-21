// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (MsDataScan.cs) is part of MassSpectrometry.
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
using System.Linq;

namespace MassSpectrometry
{
    public abstract class MsDataScanWithPrecursor<TSpectrum> : MsDataScan<TSpectrum>, IMsDataScanWithPrecursor<TSpectrum>
        where TSpectrum : IMzSpectrum<IMzPeak>
    {

        #region Private Fields

        //TODO: Validate these values
        private static readonly double[] mms = new double[] { 1.0025, 2.005, 3.0075, 4.010 };

        #endregion Private Fields

        #region Protected Constructors

        protected MsDataScanWithPrecursor(int ScanNumber, int MsnOrder, bool isCentroid, Polarity Polarity, double RetentionTime, MzRange MzRange, string ScanFilter, MZAnalyzerType MzAnalyzer, double TotalIonCurrent, double? selectedIonGuessMZ, int? selectedIonGuessChargeStateGuess, double? selectedIonGuessIntensity, double isolationMZ, double? isolationWidth, DissociationType dissociationType, int oneBasedPrecursorScanNumber, double? selectedIonGuessMonoisotopicMZ)
                                : base(ScanNumber, MsnOrder, isCentroid, Polarity, RetentionTime, MzRange, ScanFilter, MzAnalyzer, TotalIonCurrent)
        {
            this.IsolationMz = isolationMZ;
            this.SelectedIonGuessChargeStateGuess = selectedIonGuessChargeStateGuess;
            this.SelectedIonGuessIntensity = selectedIonGuessIntensity;
            this.SelectedIonGuessMZ = selectedIonGuessMZ;
            this.DissociationType = dissociationType;
            this.IsolationWidth = isolationWidth;
            this.OneBasedPrecursorScanNumber = oneBasedPrecursorScanNumber;
            this.SelectedIonGuessMonoisotopicMZ = selectedIonGuessMonoisotopicMZ;
        }

        #endregion Protected Constructors

        #region Public Properties

        public double IsolationMz { get; private set; }
        public int? SelectedIonGuessChargeStateGuess { get; private set; }
        public double? SelectedIonGuessIntensity { get; private set; }
        public double? SelectedIonGuessMZ { get; private set; }
        public DissociationType DissociationType { get; private set; }
        public double? IsolationWidth { get; private set; }
        public int OneBasedPrecursorScanNumber { get; private set; }
        public double? SelectedIonGuessMonoisotopicIntensity { get; private set; }
        public double? SelectedIonGuessMonoisotopicMZ { get; private set; }

        public MzRange IsolationRange
        {
            get
            {
                if (IsolationWidth.HasValue)
                    return new MzRange(IsolationMz - IsolationWidth.Value / 2, IsolationMz + IsolationWidth.Value / 2);
                return null;
            }
        }

        #endregion Public Properties

        #region Public Methods

        public void TranformByApplyingFunctionsToSpectraAndReplacingPrecursorMZs(Func<IMzPeak, double> convertorForSpectrum, double selectedIonGuessMZ, double selectedIonGuessMonoisotopicMZ)
        {
            MassSpectrum.ReplaceXbyApplyingFunction(convertorForSpectrum);
            this.SelectedIonGuessMZ = selectedIonGuessMZ;
            this.SelectedIonGuessMonoisotopicMZ = selectedIonGuessMonoisotopicMZ;
        }

        public void RecomputeChargeState(IMzSpectrum<IMzPeak> precursorSpectrum, double tolHere, int maxCharge)
        {
            var peaksCloseToIsolated = precursorSpectrum.Extract(IsolationMz - 2.1, IsolationMz + 2.1).ToList();
            int[] chargeCount = new int[maxCharge]; // charges 1,2,3,4
            for (int i = 0; i < peaksCloseToIsolated.Count; i++)
                for (int j = i + 1; j < peaksCloseToIsolated.Count; j++)
                    for (int charge = 1; charge <= maxCharge; charge++)
                        for (int isotope = 0; isotope < maxCharge; isotope++)
                            if (Math.Abs(peaksCloseToIsolated[j].X - peaksCloseToIsolated[i].X - mms[isotope] / charge) < tolHere)
                                chargeCount[charge - 1]++;
            SelectedIonGuessChargeStateGuess = Array.IndexOf(chargeCount, chargeCount.Max()) + 1;
        }

        public void RecomputeSelectedPeak(IMzSpectrum<IMzPeak> precursorSpectrum)
        {
            var thePeak = precursorSpectrum.GetClosestPeak(IsolationMz);
            SelectedIonGuessIntensity = thePeak.Intensity;
            SelectedIonGuessMZ = thePeak.Mz;
        }

        public void RecomputeMonoisotopicPeak(IMzSpectrum<IMzPeak> precursorSpectrum, double tolHere, double intensityFractionNeeded)
        {
            if (!SelectedIonGuessChargeStateGuess.HasValue)
                throw new Exception("Need charge state!");
            if (!SelectedIonGuessIntensity.HasValue || !SelectedIonGuessMZ.HasValue)
                RecomputeSelectedPeak(precursorSpectrum);
            IMzPeak mPeak = null;
            foreach (var ok in mms)
            {
                var closestPeak = precursorSpectrum.GetClosestPeak(IsolationMz - ok / SelectedIonGuessChargeStateGuess.Value);
                if ((Math.Abs(closestPeak.Mz
                    - (IsolationMz - ok / SelectedIonGuessChargeStateGuess.Value))
                    < tolHere)
                    && closestPeak.Intensity > SelectedIonGuessIntensity.Value * intensityFractionNeeded)
                {
                    mPeak = closestPeak;
                }
                else
                    break;
            }
            if (mPeak != null)
            {
                SelectedIonGuessMonoisotopicIntensity = mPeak.Intensity;
                SelectedIonGuessMonoisotopicMZ = mPeak.Mz;
            }
            else
            {
                SelectedIonGuessMonoisotopicIntensity = SelectedIonGuessIntensity;
                SelectedIonGuessMonoisotopicMZ = SelectedIonGuessMZ;
            }
        }

        #endregion Public Methods

    }
}