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
    public sealed class ThermoSpectrum : MzSpectrum<MzPeak>
    {
        #region Public Constructors

        public ThermoSpectrum(double[] mz, double[] intensity, bool shouldCopy)
            : base(mz, intensity, shouldCopy)
        {
        }

        public ThermoSpectrum(ThermoSpectrum thermoSpectrum)
            : this(thermoSpectrum.XArray, thermoSpectrum.YArray, true)
        {
        }

        #endregion Public Constructors

        #region Internal Constructors

        internal ThermoSpectrum(double[,] peakData)
                            : base(peakData)
        {
        }

        #endregion Internal Constructors

        #region Protected Methods

        protected override MzPeak GeneratePeak(int index)
        {
            return new MzPeak(XArray[index], YArray[index]);
        }

        #endregion Protected Methods
    }
}