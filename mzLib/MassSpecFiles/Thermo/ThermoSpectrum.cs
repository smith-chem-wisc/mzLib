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

using MassSpectrometry;
using System;

namespace IO.Thermo
{
    /// <summary>
    /// A high resolution spectra from a Thermo raw file
    /// </summary>
    [Serializable]
    public sealed class ThermoSpectrum : MzSpectrum<ThermoMzPeak>
    {
        #region Private Fields

        private readonly double[] _noises;
        private readonly double[] _resolutions;
        private readonly int[] _charges;

        #endregion Private Fields

        #region Public Constructors

        public ThermoSpectrum(double[] mz, double[] intensity, double[] noise, int[] charge, double[] resolutions, bool shouldCopy)
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
            : this(thermoSpectrum.XArray, thermoSpectrum.YArray, thermoSpectrum._noises, thermoSpectrum._charges, thermoSpectrum._resolutions, true)
        {
        }

        #endregion Public Constructors

        #region Internal Constructors

        internal ThermoSpectrum(double[,] peakData)
                            : base(peakData)
        {
            int arrayLength = peakData.GetLength(1);
            int depthLength = peakData.GetLength(0);
            if (depthLength <= 2)
                return;
            _noises = new double[Size];
            _resolutions = new double[Size];
            var charges = new double[Size];

            Buffer.BlockCopy(peakData, sizeof(double) * arrayLength * (int)ThermoRawFile.RawLabelDataColumn.Resolution, _resolutions, 0, sizeof(double) * Size);
            Buffer.BlockCopy(peakData, sizeof(double) * arrayLength * (int)ThermoRawFile.RawLabelDataColumn.NoiseLevel, _noises, 0, sizeof(double) * Size);
            Buffer.BlockCopy(peakData, sizeof(double) * arrayLength * (int)ThermoRawFile.RawLabelDataColumn.Charge, charges, 0, sizeof(double) * Size);

            _charges = new int[Size];
            for (int i = 0; i < Size; i++)
            {
                _charges[i] = (int)charges[i];
            }
        }

        #endregion Internal Constructors

        #region Public Indexers

        public override ThermoMzPeak this[int index]
        {
            get
            {
                if (peakList[index] == null)
                    peakList[index] = _charges == null ?
                        new ThermoMzPeak(XArray[index], YArray[index]) :
                new ThermoMzPeak(XArray[index], YArray[index], _charges[index], _noises[index], _resolutions[index]);
                return peakList[index];
            }
        }

        #endregion Public Indexers

        #region Public Methods

        public double GetSignalToNoise(int index)
        {
            if (_noises == null)
                return double.NaN;
            double noise = _noises[index];
            return YArray[index] / noise;
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

        public override double[,] CopyTo2DArray()
        {
            double[,] data = new double[5, Size];
            const int size = sizeof(double);
            int bytesToCopy = size * Size;
            Buffer.BlockCopy(XArray, 0, data, 0, bytesToCopy);
            Buffer.BlockCopy(YArray, 0, data, bytesToCopy, bytesToCopy);
            Buffer.BlockCopy(_resolutions, 0, data, 2 * bytesToCopy, bytesToCopy);
            Buffer.BlockCopy(_noises, 0, data, 3 * bytesToCopy, bytesToCopy);

            double[] charges = new double[Size];
            for (int i = 0; i < Size; i++)
                charges[i] = _charges[i];

            Buffer.BlockCopy(charges, 0, data, 4 * bytesToCopy, bytesToCopy);
            return data;
        }

        #endregion Public Methods

        #region Protected Methods

        protected override ThermoMzPeak GeneratePeak(int index)
        {
            return peakList[index] = _charges == null ?
                new ThermoMzPeak(XArray[index], YArray[index]) :
                new ThermoMzPeak(XArray[index], YArray[index], _charges[index], _noises[index], _resolutions[index]);
        }

        #endregion Protected Methods
    }
}