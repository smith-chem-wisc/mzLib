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
using Spectra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    public abstract class MsDataScanWithPrecursor<TSpectrum> : MsDataScan<TSpectrum>, IMsDataScanWithPrecursor<TSpectrum>
        where TSpectrum : IMzSpectrum<IMzPeak>
    {
        #region Private Fields

        private MzRange isolationRange;

        #endregion Private Fields

        #region Protected Constructors

        protected MsDataScanWithPrecursor(TSpectrum massSpectrum, int ScanNumber, int MsnOrder, bool isCentroid, Polarity Polarity, double RetentionTime, MzRange MzRange, string ScanFilter, MZAnalyzerType MzAnalyzer, double TotalIonCurrent, double selectedIonMZ, int? selectedIonChargeStateGuess, double? selectedIonIntensity, double? isolationMZ, double? isolationWidth, DissociationType dissociationType, int? oneBasedPrecursorScanNumber, double? selectedIonMonoisotopicGuessMz, double? injectionTime, double[,] noiseData, string nativeId)
                                                        : base(massSpectrum, ScanNumber, MsnOrder, isCentroid, Polarity, RetentionTime, MzRange, ScanFilter, MzAnalyzer, TotalIonCurrent, injectionTime, noiseData, nativeId)
        {
            this.OneBasedPrecursorScanNumber = oneBasedPrecursorScanNumber;

            this.IsolationMz = isolationMZ;
            this.IsolationWidth = isolationWidth;

            this.DissociationType = dissociationType;

            this.SelectedIonMZ = selectedIonMZ;
            this.SelectedIonIntensity = selectedIonIntensity;
            this.SelectedIonChargeStateGuess = selectedIonChargeStateGuess;
            this.SelectedIonMonoisotopicGuessMz = selectedIonMonoisotopicGuessMz;
        }

        #endregion Protected Constructors

        #region Public Properties

        public double? IsolationMz { get; private set; } // May be adjusted by calibration
        public int? SelectedIonChargeStateGuess { get; }
        public double? SelectedIonIntensity { get; private set; } // May be refined
        public double SelectedIonMZ { get; private set; } // May be adjusted by calibration
        public DissociationType DissociationType { get; }
        public double? IsolationWidth { get; }
        public int? OneBasedPrecursorScanNumber { get; }
        public double? SelectedIonMonoisotopicGuessIntensity { get; private set; } // May be refined
        public double? SelectedIonMonoisotopicGuessMz { get; private set; } // May be refined

        public MzRange IsolationRange
        {
            get
            {
                if (isolationRange != null)
                    return isolationRange;
                if (IsolationWidth.HasValue && IsolationMz.HasValue)
                    isolationRange = new MzRange(IsolationMz.Value - IsolationWidth.Value / 2, IsolationMz.Value + IsolationWidth.Value / 2);
                return isolationRange;
            }
        }

        #endregion Public Properties

        #region Public Methods

        public IEnumerable<IsotopicEnvelope> GetIsolatedMassesAndCharges(IMzSpectrum<IMzPeak> precursorSpectrum, int minAssumedChargeState, int maxAssumedChargeState, double deconvolutionTolerancePpm, double intensityRatio)
        {
            if (IsolationRange == null)
                yield break;

            foreach (var haha in precursorSpectrum.Deconvolute(new MzRange(IsolationRange.Minimum - 8.5, IsolationRange.Maximum + 8.5), minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatio)
                                                  .Where(b => b.peaks.Any(cc => isolationRange.Contains(cc.mz))))
                yield return haha;
        }

        public void TransformMzs(Func<IPeak, double> convertorForSpectrum, Func<IPeak, double> convertorForPrecursor)
        {
            MassSpectrum.ReplaceXbyApplyingFunction(convertorForSpectrum);
            this.SelectedIonMZ = convertorForPrecursor(new MzPeak(SelectedIonMZ, SelectedIonIntensity.Value));
            if (SelectedIonMonoisotopicGuessMz.HasValue)
                this.SelectedIonMonoisotopicGuessMz = convertorForPrecursor(new MzPeak(SelectedIonMonoisotopicGuessMz.Value, SelectedIonMonoisotopicGuessIntensity.Value));
            if (IsolationMz.HasValue)
                this.IsolationMz = convertorForPrecursor(new MzPeak(IsolationMz.Value, SelectedIonIntensity.Value));

            // Will need to recompute this...
            isolationRange = null;
        }

        public void RefineSelectedMzAndIntensity(IMzSpectrum<IMzPeak> precursorSpectrum)
        {
            if (IsolationMz.HasValue)
            {
                var thePeak = precursorSpectrum.GetClosestPeakIndex(IsolationMz.Value);
                SelectedIonIntensity = precursorSpectrum.YArray[thePeak];
                SelectedIonMZ = precursorSpectrum.XArray[thePeak];
            }
        }

        public void ComputeSelectedPeakIntensity(IMzSpectrum<IMzPeak> precursorSpectrum)
        {
            var thePeak = precursorSpectrum.GetClosestPeakIndex(SelectedIonMZ);
            SelectedIonIntensity = precursorSpectrum.YArray[thePeak];
            SelectedIonMZ = precursorSpectrum.XArray[thePeak];
        }

        public void ComputeMonoisotopicPeakIntensity(IMzSpectrum<IMzPeak> precursorSpectrum)
        {
            var thePeak = precursorSpectrum.GetClosestPeakIndex(SelectedIonMonoisotopicGuessMz.Value);
            SelectedIonMonoisotopicGuessIntensity = precursorSpectrum.YArray[thePeak];
            SelectedIonMonoisotopicGuessMz = precursorSpectrum.XArray[thePeak];
        }

        public IEnumerable<Tuple<List<IMzPeak>, int>> GetIsolatedMassesAndChargesOld(IMzSpectrum<IMzPeak> precursorSpectrum, int maxAssumedChargeState, Tolerance massTolerance, double intensityRatio)
        {
            if (IsolationRange == null)
                yield break;

            foreach (var haha in precursorSpectrum.DeconvoluteOld(new MzRange(IsolationRange.Minimum - 8.5, IsolationRange.Maximum + 8.5), maxAssumedChargeState, massTolerance, intensityRatio)
                                                  .Where(b => b.Item1.Any(cc => isolationRange.Contains(cc.Mz))))
                yield return haha;
        }

        #endregion Public Methods
    }
}