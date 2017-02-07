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

namespace MassSpectrometry
{
    public abstract class MsDataScanWithPrecursor<TSpectrum, TPeak> : MsDataScan<TSpectrum, TPeak>, IMsDataScanWithPrecursor<TSpectrum, TPeak>
        where TPeak : IMzPeak
        where TSpectrum : IMzSpectrum<TPeak>
    {

        #region Protected Constructors

        protected MsDataScanWithPrecursor(int ScanNumber, string id, int MsnOrder, bool isCentroid, Polarity Polarity, double RetentionTime, MzRange MzRange, string ScanFilter, MZAnalyzerType MzAnalyzer, double InjectionTime, double TotalIonCurrent, string precursorID, double selectedIonGuessMZ, int? selectedIonGuessChargeStateGuess, double selectedIonGuessIntensity, double isolationMZ, double isolationWidth, DissociationType dissociationType, int oneBasedPrecursorScanNumber, double selectedIonGuessMonoisotopicIntensity, double selectedIonGuessMonoisotopicMZ)
                        : base(ScanNumber, id, MsnOrder, isCentroid, Polarity, RetentionTime, MzRange, ScanFilter, MzAnalyzer, InjectionTime, TotalIonCurrent)
        {
            this.IsolationMz = isolationMZ;
            this.PrecursorID = precursorID;
            this.SelectedIonGuessChargeStateGuess = selectedIonGuessChargeStateGuess;
            this.SelectedIonGuessIntensity = selectedIonGuessIntensity;
            this.SelectedIonGuessMZ = selectedIonGuessMZ;
            this.DissociationType = dissociationType;
            this.IsolationWidth = isolationWidth;
            this.OneBasedPrecursorScanNumber = oneBasedPrecursorScanNumber;
            this.SelectedIonGuessMonoisotopicIntensity = selectedIonGuessMonoisotopicIntensity;
            this.SelectedIonGuessMonoisotopicMZ = selectedIonGuessMonoisotopicMZ;
        }

        #endregion Protected Constructors

        #region Public Properties

        public double IsolationMz { get; private set; }
        public string PrecursorID { get; private set; }
        public int? SelectedIonGuessChargeStateGuess { get; private set; }
        public double SelectedIonGuessIntensity { get; private set; }
        public double SelectedIonGuessMZ { get; private set; }
        public DissociationType DissociationType { get; private set; }
        public double IsolationWidth { get; private set; }
        public int OneBasedPrecursorScanNumber { get; private set; }
        public double SelectedIonGuessMonoisotopicIntensity { get; private set; }
        public double SelectedIonGuessMonoisotopicMZ { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            return string.Format("Scan #{0}", OneBasedScanNumber);
        }

        public MzRange IsolationRange
        {
            get
            {
                return new MzRange(IsolationMz - IsolationWidth / 2, IsolationMz + IsolationWidth / 2);
            }
        }

        public void TranformByApplyingFunctionsToSpectraAndReplacingPrecursorMZs(Func<IMzPeak, double> convertorForSpectrum, double selectedIonGuessMZ, double selectedIonGuessMonoisotopicMZ)
        {
            MassSpectrum.ReplaceXbyApplyingFunction(convertorForSpectrum);
            this.SelectedIonGuessMZ = selectedIonGuessMZ;
            this.SelectedIonGuessMonoisotopicMZ = selectedIonGuessMonoisotopicMZ;
        }

        #endregion Public Methods

    }
}