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
using System.Collections.Generic;
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
            string scanDescription = null,
            double? compensationVoltage = null,
            int[] chargeArray = null)
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
            CompensationVoltage = compensationVoltage;
            ChargeArray = chargeArray;

            // Ensure the charge of the selected ion matches the polarity of the scan
            SelectedIonChargeStateGuess = Polarity switch
            {
                Polarity.Negative when selectedIonChargeStateGuess is > 0 => -selectedIonChargeStateGuess,
                Polarity.Positive when selectedIonChargeStateGuess is < 0 =>  Math.Abs(selectedIonChargeStateGuess.Value),
                _ => selectedIonChargeStateGuess
            };
        }

        /// <summary>
        /// Per-peak charge state, parallel to <see cref="MassSpectrum"/>.XArray.
        /// Length must equal MassSpectrum.Size when set; null means "no charge data captured".
        /// Encoded in mzML as a third &lt;binaryDataArray&gt; with cvParam
        /// MS:1000516 ("charge array"), 32-bit float, no compression
        /// (charge arrays are small enough that the zlib overhead isn't worth it; could
        /// be made configurable in a follow-up if a downstream consumer needs it).
        ///
        /// Only set by deisotopers (e.g. YADA's annotate mode); vendor-converted mzML
        /// rarely includes it. Readers that don't understand MS:1000516 silently ignore it.
        /// </summary>
        public int[] ChargeArray { get; protected set; }

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
        public double? CompensationVoltage { get; private set; }

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
        /// Builds a <see cref="NoiseData"/> table from parallel arrays: row 0 the m/z of each sampled point, row 1 the
        /// noise level there and row 2 the baseline there, as the instrument reports them, ordered by m/z. Thermo reports
        /// such a table, a few dozen points per Orbitrap scan, and its per-peak noise is this table interpolated
        /// linearly at the peak (<see cref="TryGetNoiseAndBaseline"/>). Null when there are no points.
        /// </summary>
        /// <exception cref="ArgumentNullException">An array is null.</exception>
        /// <exception cref="ArgumentException">The arrays differ in length.</exception>
        public static double[,] CreateNoiseData(double[] mzs, double[] noise, double[] baseline)
        {
            ArgumentNullException.ThrowIfNull(mzs);
            ArgumentNullException.ThrowIfNull(noise);
            ArgumentNullException.ThrowIfNull(baseline);
            if (noise.Length != mzs.Length || baseline.Length != mzs.Length)
            {
                throw new ArgumentException($"The noise table needs one noise and one baseline value per m/z; got {mzs.Length} m/z, " +
                    $"{noise.Length} noise and {baseline.Length} baseline values.");
            }
            if (mzs.Length == 0)
            {
                return null;
            }

            int[] order = Enumerable.Range(0, mzs.Length).OrderBy(i => mzs[i]).ToArray();
            var table = new double[3, mzs.Length];
            for (int i = 0; i < order.Length; i++)
            {
                table[0, i] = mzs[order[i]];
                table[1, i] = noise[order[i]];
                table[2, i] = baseline[order[i]];
            }
            return table;
        }

        /// <summary>
        /// The noise level and baseline at <paramref name="mz"/>, interpolated linearly between the two points of
        /// <see cref="NoiseData"/> around it; beyond the first or last point, that point's values. False, with both
        /// outputs NaN, when this scan has no noise table: <see cref="NoiseData"/> is null or not three rows of at least
        /// one point (ion-trap scans, and files whose reader does not give noise).
        /// </summary>
        public bool TryGetNoiseAndBaseline(double mz, out double noise, out double baseline)
        {
            noise = double.NaN;
            baseline = double.NaN;
            if (NoiseData == null || NoiseData.GetLength(0) != 3 || NoiseData.GetLength(1) < 1)
            {
                return false;
            }

            int count = NoiseData.GetLength(1);
            if (mz <= NoiseData[0, 0])
            {
                noise = NoiseData[1, 0];
                baseline = NoiseData[2, 0];
                return true;
            }
            if (mz >= NoiseData[0, count - 1])
            {
                noise = NoiseData[1, count - 1];
                baseline = NoiseData[2, count - 1];
                return true;
            }

            // The last point at or below mz; the table is ordered by m/z, so the next point is above it.
            int low = 0;
            int high = count - 1;
            while (high - low > 1)
            {
                int middle = (low + high) / 2;
                if (NoiseData[0, middle] <= mz) low = middle; else high = middle;
            }
            double fraction = (mz - NoiseData[0, low]) / (NoiseData[0, high] - NoiseData[0, low]);
            noise = NoiseData[1, low] + fraction * (NoiseData[1, high] - NoiseData[1, low]);
            baseline = NoiseData[2, low] + fraction * (NoiseData[2, high] - NoiseData[2, low]);
            return true;
        }

        /// <summary>
        /// The signal-to-noise of a peak of <paramref name="intensity"/> at <paramref name="mz"/>, as Thermo reports it:
        /// (intensity - baseline) / (noise - baseline), with noise and baseline from <see cref="TryGetNoiseAndBaseline"/>.
        /// NaN when this scan has no noise table, or when noise does not exceed baseline there.
        /// </summary>
        public double GetSignalToNoise(double mz, double intensity)
        {
            if (!TryGetNoiseAndBaseline(mz, out double noise, out double baseline) || noise <= baseline)
            {
                return double.NaN;
            }
            return (intensity - baseline) / (noise - baseline);
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
            const double isolationPadding = 8.5;

            if (IsolationRange is null)
                return new List<IsotopicEnvelope>();

            var range = new MzRange(IsolationRange.Minimum - isolationPadding, IsolationRange.Maximum + isolationPadding);

            return Deconvoluter.Deconvolute(precursorSpectrum, deconParameters, range)
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
            const double isolationPadding = 8.5;

            if (IsolationRange is null)
                return new List<IsotopicEnvelope>();

            var range = new MzRange(IsolationRange.Minimum - isolationPadding, IsolationRange.Maximum + isolationPadding);

            return Deconvoluter.Deconvolute(precursorScan, deconParameters, range)
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