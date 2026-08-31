// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (MzRange.cs) is part of MassSpectrometry.
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
    public class MzRange : DoubleRange
    {
        public MzRange(double minMZ, double maxMZ)
            : base(minMZ, maxMZ)
        {
        }

        public override string ToString(string format)
        {
            return $"[{Minimum.ToString(format, System.Globalization.CultureInfo.InvariantCulture)} to {Maximum.ToString(format, System.Globalization.CultureInfo.InvariantCulture)}] m/z";
        }
    }

    /// <summary>
    /// MzRange augmented with a retention-time window. Used by deconvolution algorithms
    /// (e.g. <c>FromFileDeconvolutionAlgorithm</c>) that need to filter results by both
    /// m/z and RT simultaneously. Inherited <see cref="MzRange"/> methods still filter
    /// by m/z only.
    /// </summary>
    public class MzRtRange : MzRange
    {
        public double MinimumRt { get; }
        public double MaximumRt { get; }

        public MzRtRange(double minMZ, double maxMZ, double minRt, double maxRt)
            : base(minMZ, maxMZ)
        {
            MinimumRt = minRt;
            MaximumRt = maxRt;
        }

        public MzRtRange(DoubleRange mzRange, double retentionTime, double rtTolerance = 0.001)
            : base(mzRange.Minimum, mzRange.Maximum)
        {
            MinimumRt = retentionTime - rtTolerance;
            MaximumRt = retentionTime + rtTolerance;
        }

        public double MinimumMZ => Minimum;
        public double MaximumMZ => Maximum;

        public double MeanMZ => (Minimum + Maximum) / 2;
        public double MeanRt => (MinimumRt + MaximumRt) / 2;

        public double MzRangeWidth => Maximum - Minimum;
        public double RtRangeWidth => MaximumRt - MinimumRt;

        public int CompareTo(double mz, double rt)
        {
            if (mz < Minimum) return 1;
            if (mz > Maximum) return -1;
            if (rt < MinimumRt) return 1;
            if (rt > MaximumRt) return -1;
            return 0;
        }

        public bool Contains(double mz, double rt)
        {
            return mz >= Minimum && mz <= Maximum && rt >= MinimumRt && rt <= MaximumRt;
        }

        public override string ToString(string format)
        {
            return $"[{Minimum.ToString(format, System.Globalization.CultureInfo.InvariantCulture)} to {Maximum.ToString(format, System.Globalization.CultureInfo.InvariantCulture)}] m/z, [{MinimumRt.ToString(format, System.Globalization.CultureInfo.InvariantCulture)} to {MaximumRt.ToString(format, System.Globalization.CultureInfo.InvariantCulture)}] RT";
        }
    }
}