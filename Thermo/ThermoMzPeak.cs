﻿// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (ThermoMzPeak.cs) is part of MassSpecFiles.
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

namespace IO.Thermo
{
    public class ThermoMzPeak : MzPeak
    {
        #region Public Constructors

        public ThermoMzPeak(double mz, double intensity)
            : base(mz, intensity)
        {
        }

        public ThermoMzPeak(double mz, double intensity, int charge, double noise, double resolution)
            : base(mz, intensity)
        {
            Charge = charge;
            Noise = noise;
            Resolution = resolution;
        }

        #endregion Public Constructors

        #region Public Properties

        public int Charge { get; private set; }

        public double Noise { get; private set; }

        public double Resolution { get; private set; }

        public double SignalToNoise
        {
            get
            {
                return Intensity / Noise;
            }
        }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            return string.Format("{0} z = {1:+#;-#;?} SN = {2:F2}", base.ToString(), Charge, SignalToNoise);
        }

        #endregion Public Methods
    }
}