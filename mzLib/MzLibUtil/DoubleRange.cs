// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (DoubleRange.cs) is part of MassSpectrometry.
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

namespace MzLibUtil
{
    public class DoubleRange
    {
        /// <summary>
        /// Creates a range from the minimum to maximum values
        /// </summary>
        /// <param name="minimum">The minimum value of the range</param>
        /// <param name="maximum">The maximum value of the range</param>
        public DoubleRange(double minimum, double maximum)
        {
            Minimum = minimum;
            Maximum = maximum;
        }

        /// <summary>
        /// Creates a range from another double range. This is the
        /// clone constructor.
        /// </summary>
        /// <param name="range">The other range to copy</param>
        public DoubleRange(DoubleRange range)
            : this(range.Minimum, range.Maximum)
        {
        }

        /// <summary>
        /// The maximum value of the range
        /// </summary>
        public double Maximum { get; protected set; }

        /// <summary>
        /// The minimum value of the range
        /// </summary>
        public double Minimum { get; protected set; }

        /// <summary>
        /// The mean value of this range:
        /// (Max + Min) / 2
        /// </summary>
        public double Mean
        {
            get { return (Maximum + Minimum) / 2.0; }
        }

        /// <summary>
        /// The width of this range:
        /// (Max - Min)
        /// </summary>
        public double Width
        {
            get { return Maximum - Minimum; }
        }

        public override string ToString()
        {
            return ToString("G9");
        }

        public virtual string ToString(string format)
        {
            return $"[{Minimum.ToString(format, System.Globalization.CultureInfo.InvariantCulture)};{Maximum.ToString(format, System.Globalization.CultureInfo.InvariantCulture)}]";
        }

        public int CompareTo(double item)
        {
            if (Minimum.CompareTo(item) > 0)
                return -1;
            if (Maximum.CompareTo(item) < 0)
                return 1;
            return 0;
        }

        /// <summary>
        /// Checks to see if this range is a proper super range of another range (inclusive)
        /// </summary>
        /// <param name="other">The other range to compare to</param>
        /// <returns>True if this range is fully encloses the other range, false otherwise</returns>
        public bool IsSuperRange(DoubleRange other)
        {
            return Maximum.CompareTo(other.Maximum) >= 0 && Minimum.CompareTo(other.Minimum) <= 0;
        }

        /// <summary>
        /// Checks to see if this range is a proper sub range of another range (inclusive)
        /// </summary>
        /// <param name="other">The other range to compare to</param>
        /// <returns>True if this range is fully enclosed by the other range, false otherwise</returns>
        public bool IsSubRange(DoubleRange other)
        {
            return Maximum.CompareTo(other.Maximum) <= 0 && Minimum.CompareTo(other.Minimum) >= 0;
        }

        /// <summary>
        /// Checks to see if this range overlaps another range (inclusive)
        /// </summary>
        /// <param name="other">The other range to compare to</param>
        /// <returns>True if the other range in any way overlaps this range, false otherwise</returns>
        public bool IsOverlapping(DoubleRange other)
        {
            return Maximum.CompareTo(other.Minimum) >= 0 && Minimum.CompareTo(other.Maximum) <= 0;
        }

        /// <summary>
        /// Determines if the item is contained within a range of values
        /// </summary>
        /// <param name="item">The item to compare against the range</param>
        /// <returns>True if the item is within the range (inclusive), false otherwise</returns>
        public bool Contains(double item)
        {
            return CompareTo(item).Equals(0);
        }
    }
}