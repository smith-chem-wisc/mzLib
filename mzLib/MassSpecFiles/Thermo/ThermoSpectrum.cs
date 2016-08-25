// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
// 
// This file (ThermoSpectrum.cs) is part of MassSpecFiles.
// 
// MassSpecFiles is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// MassSpecFiles is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with MassSpecFiles. If not, see <http://www.gnu.org/licenses/>.

using Spectra;
using System;

namespace IO.Thermo
{
    /// <summary>
    /// A high resolution spectra from a Thermo raw file
    /// </summary>
    [Serializable]
    public sealed class ThermoSpectrum : MzSpectrum<ThermoMzPeak>
    {

        private readonly double[] _noises;
        private readonly double[] _resolutions;
        private readonly int[] _charges;

        internal ThermoSpectrum(double[,] peakData)
            : base(peakData)
        {
            int arrayLength = peakData.GetLength(1);
            int depthLength = peakData.GetLength(0);
            if (depthLength <= 2)
                return;
            _noises = new double[Count];
            _resolutions = new double[Count];
            var charges = new double[Count];

            Buffer.BlockCopy(peakData, sizeof(double) * arrayLength * (int)ThermoRawFile.RawLabelDataColumn.Resolution, _resolutions, 0, sizeof(double) * Count);
            Buffer.BlockCopy(peakData, sizeof(double) * arrayLength * (int)ThermoRawFile.RawLabelDataColumn.NoiseLevel, _noises, 0, sizeof(double) * Count);
            Buffer.BlockCopy(peakData, sizeof(double) * arrayLength * (int)ThermoRawFile.RawLabelDataColumn.Charge, charges, 0, sizeof(double) * Count);

            _charges = new int[Count];
            for (int i = 0; i < Count; i++)
            {
                _charges[i] = (int)charges[i];
            }
        }

        public ThermoSpectrum(double[] mz, double[] intensity, double[] noise, int[] charge, double[] resolutions, bool shouldCopy = true)
            : base(mz, intensity, shouldCopy)
        {
            if (!shouldCopy)
            {
                _noises = noise;
                _resolutions = resolutions;
                _charges = charge;
            }
            else
            {
                if (noise != null)
                {
                    _noises = new double[noise.Length];
                    Array.Copy(noise, _noises, noise.Length);
                }
                if (resolutions != null)
                {
                    _resolutions = new double[resolutions.Length];
                    Array.Copy(resolutions, _resolutions, resolutions.Length);
                }
                if (charge != null)
                {
                    _charges = new int[charge.Length];
                    Array.Copy(charge, _charges, charge.Length);
                }
            }
        }

        public ThermoSpectrum(ThermoSpectrum thermoSpectrum)
            : this(thermoSpectrum.xArray, thermoSpectrum.yArray, thermoSpectrum._noises, thermoSpectrum._charges, thermoSpectrum._resolutions)
        {

        }

        public double GetSignalToNoise(int index)
        {
            if (_noises == null)
                return double.NaN;
            double noise = _noises[index];
            return yArray[index] / noise;
        }

        public double[] GetNoises()
        {
            return _noises;
        }

        public double[] GetResolutions()
        {
            return _resolutions;
        }

        public int[] GetCharges()
        {
            return _charges;
        }

        public override ThermoMzPeak this[int index]
        {
            get
            {
                if (peakList[index] == null)
                    peakList[index] = _charges == null ?
                        new ThermoMzPeak(xArray[index], yArray[index]) :
                new ThermoMzPeak(xArray[index], yArray[index], _charges[index], _noises[index], _resolutions[index]);
                return peakList[index];
            }
        }

        public override double[,] CopyTo2DArray()
        {
            double[,] data = new double[5, Count];
            const int size = sizeof(double);
            int bytesToCopy = size * Count;
            Buffer.BlockCopy(xArray, 0, data, 0, bytesToCopy);
            Buffer.BlockCopy(yArray, 0, data, bytesToCopy, bytesToCopy);
            Buffer.BlockCopy(_resolutions, 0, data, 2 * bytesToCopy, bytesToCopy);
            Buffer.BlockCopy(_noises, 0, data, 3 * bytesToCopy, bytesToCopy);

            double[] charges = new double[Count];
            for (int i = 0; i < Count; i++)
                charges[i] = _charges[i];

            Buffer.BlockCopy(charges, 0, data, 4 * bytesToCopy, bytesToCopy);
            return data;
        }

        public new ThermoSpectrum newSpectrumExtract(double minMZ, double maxMZ)
        {

            int index = GetClosestPeakIndex(minMZ);
            if (this[index].X < minMZ)
                index++;
            int index2 = GetClosestPeakIndex(maxMZ);
            if (this[index2].X > maxMZ)
                index2--;

            int count = 1 + index2 - index;
            if (count <= 0)
                return new ThermoSpectrum(new double[5, 0]);

            double[] mz = new double[count];
            double[] intensity = new double[count];
            int[] charges = new int[count];
            double[] noises = new double[count];
            double[] resolutions = new double[count];
            int j = 0;

            while (index < Count && xArray[index] <= maxMZ)
            {
                mz[j] = xArray[index];
                intensity[j] = yArray[index];
                if (_charges != null)
                    charges[j] = _charges[index];
                if (_noises != null)
                    noises[j] = _noises[index];
                if (_resolutions != null)
                    resolutions[j] = _resolutions[index];
                index++;
                j++;
            }

            Array.Resize(ref mz, j);
            Array.Resize(ref intensity, j);
            Array.Resize(ref charges, j);
            Array.Resize(ref noises, j);
            Array.Resize(ref resolutions, j);
            return new ThermoSpectrum(mz, intensity, _noises == null ? null : noises, _charges == null ? null : charges, _resolutions == null ? null : resolutions, false);
        }


        public new ThermoSpectrum newSpectrumFilterByY(double minIntensity = 0, double maxIntensity = double.MaxValue)
        {

            int count = Count;
            double[] mz = new double[count];
            double[] intensities = new double[count];
            double[] resolutions = new double[count];
            double[] noises = new double[count];
            int[] charges = new int[count];

            int j = 0;
            for (int i = 0; i < count; i++)
            {
                double intensity = yArray[i];
                if (intensity >= minIntensity && intensity < maxIntensity)
                {
                    mz[j] = xArray[i];
                    intensities[j] = intensity;
                    if (_resolutions != null)
                        resolutions[j] = _resolutions[i];
                    if (_charges != null)
                        charges[j] = _charges[i];
                    if (_noises != null)
                        noises[j] = _noises[i];
                    j++;
                }
            }


            if (j != count)
            {
                Array.Resize(ref mz, j);
                Array.Resize(ref intensities, j);
                Array.Resize(ref resolutions, j);
                Array.Resize(ref noises, j);
                Array.Resize(ref charges, j);
            }

            return new ThermoSpectrum(mz, intensities, _noises == null ? null : noises, _charges == null ? null : charges, _resolutions == null ? null : resolutions, false);
        }

    }
}