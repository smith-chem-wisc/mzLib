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
using System.Text.RegularExpressions;

namespace MzLibUtil
{
    /// <summary>
    /// The tolerance, or error, of two points
    /// </summary>
    public class Tolerance
    {

        #region Private Fields

        /// <summary>
        /// A regex for parsing a string representation of a tolerance
        /// <para>
        /// i.e., "10 PPM", "-+10 PPM", "5 AbsoluteUnits", etc...
        /// </para>
        /// </summary>
        private static readonly Regex StringRegex = new Regex(@"(\+-|-\+|±)?\s*([\d.]+)\s*(PPM|AbsoluteUnits)", RegexOptions.Compiled | RegexOptions.IgnoreCase);

        #endregion Private Fields

        #region Public Constructors

        /// <summary>
        /// Creates a new tolerance given a unit, value, and whether the tolerance is ±
        /// </summary>
        /// <param name="unit">The units for this tolerance</param>
        /// <param name="value">The numerical value of the tolerance</param>
        public Tolerance(ToleranceUnit unit, double value)
        {
            Unit = unit;
            Value = Math.Abs(value);
        }

        /// <summary>
        /// Creates a new tolerance given a unit, two points (one experimental and one theoretical), and whether the tolerance is ±
        /// </summary>
        /// <param name="unit">The units for this tolerance</param>
        /// <param name="experimental">The experimental value</param>
        /// <param name="theoretical">The theoretical value</param>
        public Tolerance(ToleranceUnit unit, double experimental, double theoretical)
            : this(unit, GetTolerance(experimental, theoretical, unit))
        {
        }

        /// <summary>
        /// Calculates a tolerance from the string representation
        /// <para>
        /// i.e., "10 PPM", "-+10 PPM", "5 AbsoluteUnits", etc...
        /// </para>
        /// </summary>
        /// <param name="s"></param>
        public Tolerance(string s)
        {
            Match m = StringRegex.Match(s);
            Value = Math.Abs(double.Parse(m.Groups[2].Value));
            ToleranceUnit type;
            Enum.TryParse(m.Groups[3].Value, true, out type);
            Unit = type;
        }

        #endregion Public Constructors

        #region Public Properties

        /// <summary>
        /// The tolerance unit type
        /// </summary>
        public ToleranceUnit Unit { get; set; }

        /// <summary>
        /// The value of the tolerance
        /// </summary>
        public double Value { get; set; }

        #endregion Public Properties

        #region Public Methods

        public static double GetTolerance(double experimental, double theoretical, ToleranceUnit type)
        {
            switch (type)
            {
                case ToleranceUnit.PPM:
                    return Math.Abs((experimental - theoretical) / theoretical * 1e6);

                default:
                    return Math.Abs(experimental - theoretical);
            }
        }

        public static Tolerance FromPpm(double value)
        {
            return new Tolerance(ToleranceUnit.PPM, value);
        }

        public static Tolerance FromAbsolute(double value)
        {
            return new Tolerance(ToleranceUnit.Absolute, value);
        }

        /// <summary>
        /// Gets the range of values encompassed by this tolerance
        /// </summary>
        /// <param name="mean">The mean value</param>
        /// <returns></returns>
        public DoubleRange GetRange(double mean)
        {
            double value = Value * 2;

            double tol;
            switch (Unit)
            {
                case ToleranceUnit.PPM:
                    tol = value * mean / 2e6;
                    break;

                default:
                    tol = value / 2.0;
                    break;
            }
            return new DoubleRange(mean - tol, mean + tol);
        }

        /// <summary>
        /// Gets the minimum value that is still within this tolerance
        /// </summary>
        /// <param name="mean"></param>
        /// <returns></returns>
        public double GetMinimumValue(double mean)
        {
            double value = Value;

            switch (Unit)
            {
                case ToleranceUnit.PPM:
                    return mean * (1 - (value / 2e6));

                default:
					return mean - value;
            }
        }

        /// <summary>
        /// Gets the maximum value that is still within this tolerance
        /// </summary>
        /// <param name="mean"></param>
        /// <returns></returns>
        public double GetMaximumValue(double mean)
        {
            double value = Value;

            switch (Unit)
            {
                case ToleranceUnit.PPM:
                    return mean * (1 + (value / 2e6));

                default:
                    return mean + value ;
            }
        }

        /// <summary>
        /// Indicates if the two values provided are within this tolerance
        /// </summary>
        /// <param name="experimental">The experimental value</param>
        /// <param name="theoretical">The theoretical value</param>
        /// <returns>Returns true if the value is within this tolerance  </returns>
        public bool Within(double experimental, double theoretical)
        {
            double tolerance = Math.Abs(GetTolerance(experimental, theoretical, Unit));
            return tolerance <= Value;
        }

        public override string ToString()
        {
            return string.Format("{0}{1:f4} {2}", "±", Value, Enum.GetName(typeof(ToleranceUnit), Unit));
        }

        #endregion Public Methods

    }
}