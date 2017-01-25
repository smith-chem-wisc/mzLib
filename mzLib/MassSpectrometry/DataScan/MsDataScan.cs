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

using Spectra;
using System;

namespace MassSpectrometry
{
    public class MsDataScan<TSpectrum> : IMsDataScan<TSpectrum>
        where TSpectrum : IMzSpectrum<MzPeak>
    {

        #region Private Fields

        private readonly double isolationMZ;
        private readonly string precursorID;
        private readonly int? selectedIonGuessChargeStateGuess;
        private readonly double selectedIonGuessIntensity;
        private double selectedIonGuessMZ;
        private readonly DissociationType dissociationType;
        private readonly double isolationWidth;
        private readonly int oneBasedPrecursorScanNumber;
        private readonly double selectedIonGuessMonoisotopicIntensity;
        private double selectedIonGuessMonoisotopicMZ;

        #endregion Private Fields

        #region Public Constructors

        public MsDataScan(int oneBasedScanNumber, TSpectrum massSpectrum, string id, int msnOrder, bool isCentroid, Polarity polarity, double retentionTime, MzRange scanWindowRange, string scanFilter, MZAnalyzerType mzAnalyzer, double injectionTime, double totalIonCurrent)
        {
            OneBasedScanNumber = oneBasedScanNumber;
            MassSpectrum = massSpectrum;
            Id = id;
            MsnOrder = msnOrder;
            IsCentroid = isCentroid;
            Polarity = polarity;
            RetentionTime = retentionTime;
            ScanWindowRange = scanWindowRange;
            ScanFilter = scanFilter;
            MzAnalyzer = mzAnalyzer;
			InjectionTime = injectionTime;
            TotalIonCurrent = totalIonCurrent;
        }

        public MsDataScan(int ScanNumber, TSpectrum MassSpectrum, string id, int MsnOrder, bool isCentroid, Polarity Polarity, double RetentionTime, MzRange MzRange, string ScanFilter, MZAnalyzerType MzAnalyzer, double InjectionTime, double TotalIonCurrent, string precursorID, double selectedIonGuessMZ, int? selectedIonGuessChargeStateGuess, double selectedIonGuessIntensity, double isolationMZ, double isolationWidth, DissociationType dissociationType, int oneBasedPrecursorScanNumber, double selectedIonGuessMonoisotopicIntensity, double selectedIonGuessMonoisotopicMZ)
                    : this(ScanNumber, MassSpectrum, id, MsnOrder, isCentroid, Polarity, RetentionTime, MzRange, ScanFilter, MzAnalyzer, InjectionTime, TotalIonCurrent)
        {
            this.isolationMZ = isolationMZ;
            this.precursorID = precursorID;
            this.selectedIonGuessChargeStateGuess = selectedIonGuessChargeStateGuess;
            this.selectedIonGuessIntensity = selectedIonGuessIntensity;
            this.selectedIonGuessMZ = selectedIonGuessMZ;
            this.dissociationType = dissociationType;
            this.isolationWidth = isolationWidth;
            this.oneBasedPrecursorScanNumber = oneBasedPrecursorScanNumber;
            this.selectedIonGuessMonoisotopicIntensity = selectedIonGuessMonoisotopicIntensity;
            this.selectedIonGuessMonoisotopicMZ = selectedIonGuessMonoisotopicMZ;
        }

        #endregion Public Constructors

        #region Public Properties

        /// <summary>
        /// The mass spectrum associated with the scan
        /// </summary>
        public TSpectrum MassSpectrum { get; private set; }

        public int OneBasedScanNumber { get; private set; }

        public int MsnOrder { get; private set; }

        public double RetentionTime { get; private set; }

        public Polarity Polarity { get; private set; }

        public MZAnalyzerType MzAnalyzer { get; private set; }

        public MzRange ScanWindowRange { get; private set; }

        public string ScanFilter { get; private set; }

        public bool IsCentroid { get; private set; }

        public string Id { get; private set; }

        public double InjectionTime { get; private set; }

        public double TotalIonCurrent { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            return string.Format("Scan #{0}", OneBasedScanNumber);
        }

        public bool TryGetPrecursorID(out string PrecursorID)
        {
            if (MsnOrder == 1)
            {
                PrecursorID = null;
                return false;
            }
            PrecursorID = precursorID;
            return true;
        }

        public bool TryGetSelectedIonGuessChargeStateGuess(out int? SelectedIonGuessChargeStateGuess)
        {
            if (MsnOrder == 1)
            {
                SelectedIonGuessChargeStateGuess = 0;
                return false;
            }
            SelectedIonGuessChargeStateGuess = selectedIonGuessChargeStateGuess;
            return true;
        }

        public bool TryGetSelectedIonGuessIntensity(out double SelectedIonGuessIntensity)
        {
            if (MsnOrder == 1)
            {
                SelectedIonGuessIntensity = double.NaN;
                return false;
            }
            SelectedIonGuessIntensity = selectedIonGuessIntensity;
            return true;
        }

        public bool TryGetSelectedIonGuessMZ(out double SelectedIonGuessMZ)
        {
            if (MsnOrder == 1)
            {
                SelectedIonGuessMZ = double.NaN;
                return false;
            }
            SelectedIonGuessMZ = selectedIonGuessMZ;
            return true;
        }

        public bool TryGetDissociationType(out DissociationType DissociationType)
        {
            if (MsnOrder == 1)
            {
                DissociationType = DissociationType.Unknown;
                return false;
            }
            DissociationType = dissociationType;
            return true;
        }

        public bool TryGetIsolationWidth(out double IsolationWidth)
        {
            if (MsnOrder == 1)
            {
                IsolationWidth = double.NaN;
                return false;
            }
            IsolationWidth = isolationWidth;
            return true;
        }

        public bool TryGetIsolationMZ(out double IsolationMZ)
        {
            if (MsnOrder == 1)
            {
                IsolationMZ = double.NaN;
                return false;
            }
            IsolationMZ = isolationMZ;
            return true;
        }

        public bool TryGetIsolationRange(out MzRange IsolationRange)
        {
            IsolationRange = null;
            if (MsnOrder == 1)
                return false;

            double isolationMzHere;
            TryGetIsolationMZ(out isolationMzHere);
            double isolationWidthHere;
            TryGetIsolationWidth(out isolationWidthHere);
            IsolationRange = new MzRange(isolationMzHere - isolationWidthHere / 2, isolationMzHere + isolationWidthHere / 2);

            return true;
        }

        public bool TryGetPrecursorOneBasedScanNumber(out int OneBasedPrecursorScanNumber)
        {
            if (MsnOrder == 1)
            {
                OneBasedPrecursorScanNumber = 0;
                return false;
            }
            OneBasedPrecursorScanNumber = oneBasedPrecursorScanNumber;
            return true;
        }

        public void TranformByApplyingFunctionsToSpectraAndReplacingPrecursorMZs(Func<MzPeak, double> convertorForSpectrum, double selectedIonGuessMZ, double selectedIonGuessMonoisotopicMZ)
        {
            MassSpectrum.ReplaceXbyApplyingFunction(convertorForSpectrum);
            this.selectedIonGuessMZ = selectedIonGuessMZ;
            this.selectedIonGuessMonoisotopicMZ = selectedIonGuessMonoisotopicMZ;
        }

        public bool TryGetSelectedIonGuessMonoisotopicIntensity(out double SelectedIonGuessMonoisotopicIntensity)
        {
            if (MsnOrder == 1)
            {
                SelectedIonGuessMonoisotopicIntensity = -1;
                return false;
            }
            SelectedIonGuessMonoisotopicIntensity = selectedIonGuessMonoisotopicIntensity;
            return true;
        }

        public bool TryGetSelectedIonGuessMonoisotopicMZ(out double SelectedIonGuessMonoisotopicMZ)
        {
            if (MsnOrder == 1)
            {
                SelectedIonGuessMonoisotopicMZ = -1;
                return false;
            }
            SelectedIonGuessMonoisotopicMZ = selectedIonGuessMonoisotopicMZ;
            return true;
        }

        #endregion Public Methods

    }
}
