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

using Chemistry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    public abstract class MsDataScanWithPrecursor<TSpectrum> : MsDataScan<TSpectrum>, IMsDataScanWithPrecursor<TSpectrum>
        where TSpectrum : IMzSpectrum<IMzPeak>
    {

        #region Private Fields

        private static readonly double[] mms = new double[] { 1.0029, 2.0052, 3.0077, 4.01, 5.012, 6.0139, 7.0154, 8.0164 };

        private MzRange isolationRange;

        #endregion Private Fields

        #region Protected Constructors

        protected MsDataScanWithPrecursor(int ScanNumber, int MsnOrder, bool isCentroid, Polarity Polarity, double RetentionTime, MzRange MzRange, string ScanFilter, MZAnalyzerType MzAnalyzer, double TotalIonCurrent, double selectedIonMZ, int? selectedIonChargeStateGuess, double? selectedIonIntensity, double isolationMZ, double? isolationWidth, DissociationType dissociationType, int oneBasedPrecursorScanNumber, double? selectedIonMonoisotopicMzGuess, double? injectionTime, double[,] noiseData)
                                                : base(ScanNumber, MsnOrder, isCentroid, Polarity, RetentionTime, MzRange, ScanFilter, MzAnalyzer, TotalIonCurrent, injectionTime, noiseData)
        {
            this.OneBasedPrecursorScanNumber = oneBasedPrecursorScanNumber;

            this.IsolationMz = isolationMZ;
            this.IsolationWidth = isolationWidth;

            this.DissociationType = dissociationType;

            this.SelectedIonMZ = selectedIonMZ;
            this.SelectedIonIntensity = selectedIonIntensity;
            this.SelectedIonChargeStateGuess = selectedIonChargeStateGuess;
            this.SelectedIonMonoisotopicMzGuess = selectedIonMonoisotopicMzGuess;
        }

        #endregion Protected Constructors

        #region Public Properties

        public double IsolationMz { get; private set; } // May be adjusted by calibration
        public int? SelectedIonChargeStateGuess { get; }
        public double? SelectedIonIntensity { get; private set; } // May be refined
        public double SelectedIonMZ { get; private set; } // May be adjusted by calibration
        public DissociationType DissociationType { get; }
        public double? IsolationWidth { get; }
        public int OneBasedPrecursorScanNumber { get; }
        public double? SelectedIonMonoisotopicIntensityGuess { get; private set; } // May be refined
        public double? SelectedIonMonoisotopicMzGuess { get; private set; } // May be refined

        public MzRange IsolationRange
        {
            get
            {
                if (isolationRange != null)
                    return isolationRange;
                if (IsolationWidth.HasValue)
                    isolationRange = new MzRange(IsolationMz - IsolationWidth.Value / 2, IsolationMz + IsolationWidth.Value / 2);
                return isolationRange;
            }
        }

        #endregion Public Properties

        #region Public Methods

        public List<double> GetIsolatedMasses(IMzSpectrum<IMzPeak> precursorSpectrum, int maxAssumedChargeState, Tolerance massTolerance, int maxMms)
        {
            if (IsolationRange == null)
                return null;

            var isolatedMasses = new List<double>();

            HashSet<IMzPeak> excluded = new HashSet<IMzPeak>();

            foreach (var peak in precursorSpectrum.Extract(new DoubleRange(IsolationRange.Minimum - maxMms - 0.5, IsolationRange.Maximum + maxMms + 0.5)))
            {
                if (excluded.Contains(peak))
                    continue;

                // This is a monoisotopic peak!

                var bestListOfPeaks = new List<IMzPeak>();
                int bestChargeState = 1;
                bool withinIsolationWindow = IsolationRange.Contains(peak.Mz);
                for (int chargeState = 1; chargeState <= maxAssumedChargeState; chargeState++)
                {
                    var mMass = peak.Mz.ToMass(chargeState);
                    var listOfAdditionalPeaksForThisChargeState = new List<IMzPeak>();
                    for (int mm = 1; mm <= maxMms; mm++)
                    {
                        double diffToNextMmPeak = mms[mm - 1];
                        double theorMass = mMass + diffToNextMmPeak;
                        var closestpeak = precursorSpectrum.GetClosestPeak(theorMass.ToMz(chargeState));
                        if (massTolerance.Within(closestpeak.Mz.ToMass(chargeState), theorMass))
                        {
                            // Found a match to an isotope peak for this charge state!
                            listOfAdditionalPeaksForThisChargeState.Add(closestpeak);
                        }
                        else
                            break;
                    }
                    if (listOfAdditionalPeaksForThisChargeState.Count > bestListOfPeaks.Count)
                    {
                        bestListOfPeaks = listOfAdditionalPeaksForThisChargeState;
                        if (!withinIsolationWindow)
                            withinIsolationWindow = listOfAdditionalPeaksForThisChargeState.Any(b => IsolationRange.Contains(b.Mz));
                        bestChargeState = chargeState;
                    }
                }

                foreach (var peakToExclude in bestListOfPeaks)
                    excluded.Add(peakToExclude);
                if (withinIsolationWindow)
                    isolatedMasses.Add(peak.Mz.ToMass(bestChargeState));
            }
            return isolatedMasses;
        }

        public void TransforMzs(Func<IMzPeak, double> convertorForSpectrum, Func<IMzPeak, double> convertorForPrecursor)
        {
            MassSpectrum.ReplaceXbyApplyingFunction(convertorForSpectrum);
            this.SelectedIonMZ = convertorForPrecursor(new MzPeak(SelectedIonMZ, SelectedIonIntensity.Value));
            if (SelectedIonMonoisotopicMzGuess.HasValue)
                this.SelectedIonMonoisotopicMzGuess = convertorForPrecursor(new MzPeak(SelectedIonMonoisotopicMzGuess.Value, SelectedIonMonoisotopicIntensityGuess.Value));
            this.IsolationMz = convertorForPrecursor(new MzPeak(IsolationMz, SelectedIonIntensity.Value));

            // Will need to recompute this...
            isolationRange = null;
        }

        public void RefineSelectedMzAndIntensity(IMzSpectrum<IMzPeak> precursorSpectrum)
        {
            var thePeak = precursorSpectrum.GetClosestPeak(IsolationMz);
            SelectedIonIntensity = thePeak.Intensity;
            SelectedIonMZ = thePeak.Mz;
        }

        public void ComputeSelectedPeakIntensity(IMzSpectrum<IMzPeak> precursorSpectrum)
        {
            var thePeak = precursorSpectrum.GetClosestPeak(SelectedIonMZ);
            SelectedIonIntensity = thePeak.Intensity;
            SelectedIonMZ = thePeak.Mz;
        }

        public void ComputeMonoisotopicPeakIntensity(IMzSpectrum<IMzPeak> precursorSpectrum)
        {
            var thePeak = precursorSpectrum.GetClosestPeak(SelectedIonMonoisotopicMzGuess.Value);
            SelectedIonMonoisotopicIntensityGuess = thePeak.Intensity;
            SelectedIonMonoisotopicMzGuess = thePeak.Mz;
        }

        #endregion Public Methods

    }
}