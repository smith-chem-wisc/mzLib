// Copyright 2016 Stefan Solntsev
// 
// This file (DefaultSpectrum.cs) is part of MassSpectrometry.
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

namespace Spectra
{
    public class DefaultSpectrum : Spectrum<DefaultPeak>
    {
        public DefaultSpectrum(ISpectrum<Peak> spectrum) : base(spectrum)
        {
        }

        public DefaultSpectrum(double[] x, double[] y, bool shouldCopy) : base(x, y, shouldCopy)
        {
        }


        public DefaultSpectrum(double[,] xy)
            : base(xy)
        {
        }

    }
}
