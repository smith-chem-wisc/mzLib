// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (Tolerance.cs) is part of MassSpectrometry.
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

using System;

namespace MzLibUtil
{
    /// <summary>
    /// The tolerance, or error, of two points
    /// </summary>
    public class PpmTolerance : Tolerance
    {
        private readonly double _factor;

        /// <summary>
        /// Creates a new tolerance given value
        /// </summary>
        /// <param name="value">The numerical value of the tolerance</param>
        public PpmTolerance(double value)
            : base(value)
        {
            _factor = value / 1e6;
        }

        public override string ToString() => $"\u00b1{Value.ToString("f4", System.Globalization.CultureInfo.InvariantCulture)} PPM";

        public override DoubleRange GetRange(double mean)
        {
            double tol = _factor * mean;
            return new DoubleRange(mean - tol, mean + tol);
        }

        public override double GetMinimumValue(double mean) => mean * (1 - _factor);

        public override double GetMaximumValue(double mean) => mean * (1 + _factor);

        public override bool Within(double experimental, double theoretical)
        {
            double diff = experimental - theoretical;
            double scaledTolerance = theoretical * _factor;
            return -scaledTolerance <= diff && diff <= scaledTolerance;
        }
    }
}