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
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;

namespace MassSpectrometry
{
    public class MsDataScan
    {
        public MsDataScan(MzSpectrum massSpectrum, 
            int oneBasedScanNumber, 
            int msnOrder, 
            bool isCentroid, 
            Polarity polarity, 
            double retentionTime, 
            MzRange scanWindowRange, 
            string scanFilter, 
            MZAnalyzerType mzAnalyzer,
            double totalIonCurrent, 
            double? injectionTime, 
            double[,] noiseData, 
            string nativeId, 
            double? selectedIonMz = null, 
            int? selectedIonChargeStateGuess = null, 
            double? selectedIonIntensity = null, 
            double? isolationMZ = null,
            double? isolationWidth = null, 
            DissociationType? dissociationType = null, 
            int? oneBasedPrecursorScanNumber = null, 
            double? selectedIonMonoisotopicGuessMz = null, 
            string hcdEnergy = null,
            string scanDescription = null)
        {
            OneBasedScanNumber = oneBasedScanNumber;
            MsnOrder = msnOrder;
            IsCentroid = isCentroid;
            Polarity = polarity;
            RetentionTime = retentionTime;
            ScanWindowRange = scanWindowRange;
            ScanFilter = scanFilter;
            MzAnalyzer = mzAnalyzer;
            TotalIonCurrent = totalIonCurrent;
            InjectionTime = injectionTime;
            NoiseData = noiseData;
            MassSpectrum = massSpectrum;
            NativeId = nativeId;
            OneBasedPrecursorScanNumber = oneBasedPrecursorScanNumber;
            IsolationMz = isolationMZ;
            IsolationWidth = isolationWidth;
            DissociationType = dissociationType;
            SelectedIonMZ = selectedIonMz;
            SelectedIonIntensity = selectedIonIntensity;
            SelectedIonMonoisotopicGuessMz = selectedIonMonoisotopicGuessMz;
            HcdEnergy = hcdEnergy;
            ScanDescription = scanDescription;

            // Ensure the charge of the selected ion matches the polarity of the scan 
            SelectedIonChargeStateGuess = Polarity switch
            {
                Polarity.Negative when selectedIonChargeStateGuess is > 0 => -selectedIonChargeStateGuess,
                Polarity.Positive when selectedIonChargeStateGuess is < 0 =>  Math.Abs(selectedIonChargeStateGuess.Value),
                _ => selectedIonChargeStateGuess
            };
        }

        /// <summary>
        /// The mass spectrum associated with the scan
        /// </summary>
        public MzSpectrum MassSpectrum { get; protected set; }

        public int OneBasedScanNumber { get; protected set; }
        public int MsnOrder { get; private set; }
        public double RetentionTime { get; }
        public Polarity Polarity { get; }
        public MZAnalyzerType MzAnalyzer { get; }
        public MzRange ScanWindowRange { get; }
        public string ScanFilter { get; }
        public string NativeId { get; private set; }
        public bool IsCentroid { get; }
        public double TotalIonCurrent { get; protected set; }
        public double? InjectionTime { get; }
        public double[,] NoiseData { get; }

        //MSn properties, all are nullable for MS1s, but MS1s are checked by evaluating if MsnOrder==1
        public double? IsolationMz { get; private set; } // May be adjusted by calibration

        public int? SelectedIonChargeStateGuess { get; }
        public double? SelectedIonIntensity { get; private set; } // May be refined
        public double? SelectedIonMZ { get; private set; } // May be adjusted by calibration
        public DissociationType? DissociationType { get; }
        public double? IsolationWidth { get; }
        public int? OneBasedPrecursorScanNumber { get; protected set; }
        public double? SelectedIonMonoisotopicGuessIntensity { get; private set; } // May be refined
        public double? SelectedIonMonoisotopicGuessMz { get; private set; } // May be refined
        public string HcdEnergy { get; private set; }
        public string ScanDescription { get; private set; }

        private MzRange isolationRange;

        public MzRange IsolationRange
        {
            get
            {
                if (isolationRange != null)
                {
                    return isolationRange;
                }
                if (IsolationWidth.HasValue && IsolationMz.HasValue)
                {
                    isolationRange = new MzRange(IsolationMz.Value - IsolationWidth.Value / 2, IsolationMz.Value + IsolationWidth.Value / 2);
                }
                return isolationRange;
            }
        }

        public override string ToString()
        {
            return string.Format("Scan #{0}", OneBasedScanNumber);
        }

        public byte[] Get64BitNoiseDataMass()
        {
            return MzSpectrum.Get64Bitarray(GetNoiseDataMass(NoiseData));
        }

        public byte[] Get64BitNoiseDataNoise()
        {
            return MzSpectrum.Get64Bitarray(GetNoiseDataNoise(NoiseData));
        }

        public byte[] Get64BitNoiseDataBaseline()
        {
            return MzSpectrum.Get64Bitarray(GetNoiseDataBaseline(NoiseData));
        }


        /// <summary>
        /// Get deconvoluted isotopic envelopes with peaks within the isolation range
        /// </summary>
        /// <param name="precursorSpectrum"> precursor spectrum</param>
        /// <param name="deconParameters">deconvolution parameters</param>
        /// <remarks>
        /// +- 8.5 are magic numbers from the original implementation of Classic Deconvolution
        /// it is believe that they represent the maximum mz space an isotopic envelopes can exist within
        /// This ensure that a peak is not cut in half by the mz isolation range
        /// </remarks>
        /// <returns></returns>
        public IEnumerable<IsotopicEnvelope> GetIsolatedMassesAndCharges(MzSpectrum precursorSpectrum,
            DeconvolutionParameters deconParameters)
        {
            return IsolationRange == null
                ? new List<IsotopicEnvelope>()
                : Deconvoluter.Deconvolute(precursorSpectrum, deconParameters,
                    new MzRange(IsolationRange.Minimum - 8.5, IsolationRange.Maximum + 8.5))
                    .Where(b => b.Peaks.Any(cc => isolationRange.Contains(cc.mz)));
        }

        /// <summary>
        /// Get deconvoluted isotopic envelopes with peaks within the isolation range
        /// </summary>
        /// <param name="precursorScan"> precursor scan</param>
        /// <param name="deconParameters">deconvolution parameters</param>
        /// <remarks>
        /// +- 8.5 are magic numbers from the original implementation of Classic Deconvolution
        /// it is believe that they represent the maximum mz space an isotopic envelopes can exist within
        /// This ensure that a peak is not cut in half by the mz isolation range
        /// </remarks>
        /// <returns></returns>
        public IEnumerable<IsotopicEnvelope> GetIsolatedMassesAndCharges(MsDataScan precursorScan,
            DeconvolutionParameters deconParameters)
        {
            return IsolationRange == null
                ? new List<IsotopicEnvelope>()
                : Deconvoluter.Deconvolute(precursorScan, deconParameters,
                    new MzRange(IsolationRange.Minimum - 8.5, IsolationRange.Maximum + 8.5))
                    .Where(b => b.Peaks.Any(cc => isolationRange.Contains(cc.mz)));
        }
        

        [Obsolete]
        public IEnumerable<IsotopicEnvelope> GetIsolatedMassesAndCharges(MzSpectrum precursorSpectrum, int minAssumedChargeState,
            int maxAssumedChargeState, double deconvolutionTolerancePpm, double intensityRatio)
        {
            if (IsolationRange == null)
            {
                yield break;
            }
            foreach (var haha in precursorSpectrum.Deconvolute(new MzRange(IsolationRange.Minimum - 8.5, IsolationRange.Maximum + 8.5),
                minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatio)
                                                  .Where(b => b.Peaks.Any(cc => isolationRange.Contains(cc.mz))))
            {
                yield return haha;
            }
        }

        public void TransformMzs(Func<MzPeak, double> convertorForSpectrum, Func<MzPeak, double> convertorForPrecursor)
        {
            MassSpectrum.ReplaceXbyApplyingFunction(convertorForSpectrum);
            SelectedIonMZ = convertorForPrecursor(new MzPeak(SelectedIonMZ.Value, SelectedIonIntensity.Value));
            if (SelectedIonMonoisotopicGuessMz.HasValue)
            {
                SelectedIonMonoisotopicGuessMz = convertorForPrecursor(new MzPeak(SelectedIonMonoisotopicGuessMz.Value, SelectedIonMonoisotopicGuessIntensity.Value));
            }
            if (IsolationMz.HasValue)
            {
                IsolationMz = convertorForPrecursor(new MzPeak(IsolationMz.Value, SelectedIonIntensity.Value));
            }

            // Will need to recompute this...
            isolationRange = null;
        }

        public void RefineSelectedMzAndIntensity(MzSpectrum precursorSpectrum)
        {
            if (!IsolationMz.HasValue)
            {
                throw new MzLibException("Could not define precursor ion because the isolation m/z window is undefined in the spectra file");
            }
            if (precursorSpectrum.Size == 0)
            {
                throw new MzLibException("Could not define precursor ion because the precursor scan contains no peaks");
            }
            int thePeak = precursorSpectrum.GetClosestPeakIndex(IsolationMz.Value);
            SelectedIonIntensity = precursorSpectrum.YArray[thePeak];
            SelectedIonMZ = precursorSpectrum.XArray[thePeak];
        }

        public void ComputeSelectedPeakIntensity(MzSpectrum precursorSpectrum)
        {
            if (precursorSpectrum.Size == 0)
            {
                throw new MzLibException("Could not compute selected peak intensity because the precursor scan contains no peaks");
            }
            int thePeak = precursorSpectrum.GetClosestPeakIndex(SelectedIonMZ.Value);
            SelectedIonIntensity = precursorSpectrum.YArray[thePeak];
            SelectedIonMZ = precursorSpectrum.XArray[thePeak];
        }

        public void ComputeMonoisotopicPeakIntensity(MzSpectrum precursorSpectrum)
        {
            if (precursorSpectrum.Size == 0)
            {
                throw new MzLibException("Could not compute monoisotopic peak intensity because the precursor scan contains no peaks");
            }
            int thePeak = precursorSpectrum.GetClosestPeakIndex(SelectedIonMonoisotopicGuessMz.Value);
            SelectedIonMonoisotopicGuessIntensity = precursorSpectrum.YArray[thePeak];
            SelectedIonMonoisotopicGuessMz = precursorSpectrum.XArray[thePeak];
        }

        public void SetOneBasedPrecursorScanNumber(int value)
        {
            this.OneBasedPrecursorScanNumber = value;
        }

        public void SetOneBasedScanNumber(int value)
        {
            this.OneBasedScanNumber = value;
        }

        public void SetMsnOrder(int value)
        {
            this.MsnOrder = value;
        }

        public void SetNativeID(string value)
        {
            this.NativeId = value;
        }

        public void SetIsolationMz(double value)
        {
            this.IsolationMz = value;
        }

        public void SetIsolationRange(double min, double max)
        {
            this.isolationRange = new MzRange(min, max);
        }

        private IEnumerable<double> GetNoiseDataMass(double[,] noiseData)
        {
            for (int i = 0; i < noiseData.Length / 3; i++)
            {
                yield return noiseData[0, i];
            }
        }

        private IEnumerable<double> GetNoiseDataNoise(double[,] noiseData)
        {
            for (int i = 0; i < noiseData.Length / 3; i++)
            {
                yield return noiseData[1, i];
            }
        }

        private IEnumerable<double> GetNoiseDataBaseline(double[,] noiseData)
        {
            for (int i = 0; i < noiseData.Length / 3; i++)
            {
                yield return noiseData[2, i];
            }
        }
    }
}