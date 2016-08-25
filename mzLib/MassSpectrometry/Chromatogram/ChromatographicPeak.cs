// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
// 
// This file (ChromatographicPeak.cs) is part of MassSpectrometry.
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

namespace MassSpectrometry
{
    public sealed class ChromatographicPeak : Peak
    {
        public double Intensity
        {
            get
            {
                return Y;
            }
            private set
            {
                Y = value;
            }
        }

        public double Time
        {
            get
            {
                return X;
            }
            private set
            {
                X = value;
            }
        }


        public ChromatographicPeak(double time, double intensity)
        {
            Time = time;
            Intensity = intensity;
        }

        public override string ToString()
        {
            return string.Format("({0:G4}, {1:G4})", Time, Intensity);
        }
    }
}