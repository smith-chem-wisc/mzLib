// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
// 
// This file (ChromatogramElutionProfile.cs) is part of MassSpectrometry.
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
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    public class ChromatographicElutionProfile<T> where T : Peak
    {
        public T StartPeak
        {
            get { return _peaks[0]; }
        }

        public T EndPeak
        {
            get { return _peaks[Count - 1]; }
        }

        public DoubleRange TimeRange { get; private set; }

        public int Count { get; private set; }

        public double SummedArea { get; private set; }

        private readonly T[] _peaks;

        public ChromatographicElutionProfile(ICollection<T> peaks)
        {
            Count = peaks.Count;
            if (Count == 0)
            {
                return;
            }

            _peaks = peaks.ToArray();
            SummedArea = _peaks.Sum(p => p.Y);
            TimeRange = new DoubleRange(_peaks[0].X, _peaks[Count - 1].X);
        }

        public double TrapezoidalArea()
        {
            double area = 0;
            for (int i = 0; i < Count - 1; i++)
            {
                T peak1 = _peaks[i];
                T peak2 = _peaks[i + 1];
                area += (peak2.X - peak1.X) * (peak2.Y + peak1.Y);
            }
            return area / 2.0;
        }





    }
}
