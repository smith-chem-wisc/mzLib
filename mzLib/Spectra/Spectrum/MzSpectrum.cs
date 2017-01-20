// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (MzSpectrum.cs) is part of MassSpectrometry.
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
using System.Collections.Generic;

namespace Spectra
{
    public abstract class MzSpectrum<TPeak> : Spectrum<TPeak>, IMzSpectrum<TPeak>
        where TPeak : MzPeak
    {

        #region Protected Constructors

        protected MzSpectrum(double[,] mzintensities) : base(mzintensities)
        {
        }

        protected MzSpectrum(ISpectrum<Peak> mZSpectrum) : base(mZSpectrum)
        {
        }

        protected MzSpectrum(double[] mz, double[] intensities, bool shouldCopy) : base(mz, intensities, shouldCopy)
        {
        }

        #endregion Protected Constructors

        #region Public Properties

        new public MzRange Range
        {
            get
            {
                return new MzRange(FirstX, LastX);
            }
        }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            return string.Format("{0} (Peaks {1})", Range, Count);
        }

        public new IMzSpectrum<MzPeak> NewSpectrumFilterByNumberOfMostIntense(int topNPeaks)
        {
            var ok = FilterByNumberOfMostIntense(topNPeaks);
            return new DefaultMzSpectrum(ok.Item1, ok.Item2, false);
        }

        public new IMzSpectrum<MzPeak> NewSpectrumExtract(DoubleRange xRange)
        {
            return NewSpectrumExtract(xRange.Minimum, xRange.Maximum);
        }

        public new IMzSpectrum<MzPeak> NewSpectrumExtract(double minX, double maxX)
        {
            var ok = Extract(minX, maxX);
            return new DefaultMzSpectrum(ok.Item1, ok.Item2, false);
        }

        public new IMzSpectrum<MzPeak> NewSpectrumWithRangesRemoved(IEnumerable<DoubleRange> xRanges)
        {
            var ok = WithRangesRemoved(xRanges);
            return new DefaultMzSpectrum(ok.Item1, ok.Item2, false);
        }

        public new IMzSpectrum<MzPeak> NewSpectrumWithRangeRemoved(DoubleRange xRange)
        {
            return NewSpectrumWithRangeRemoved(xRange.Minimum, xRange.Maximum);
        }

        public new IMzSpectrum<MzPeak> NewSpectrumWithRangeRemoved(double minX, double maxX)
        {
            var ok = WithRangeRemoved(minX, maxX);
            return new DefaultMzSpectrum(ok.Item1, ok.Item2, false);
        }

        public new IMzSpectrum<MzPeak> NewSpectrumFilterByY(double minY, double maxY)
        {
            var ok = FilterByY(minY, maxY);
            return new DefaultMzSpectrum(ok.Item1, ok.Item2, false);
        }

        public new IMzSpectrum<MzPeak> NewSpectrumFilterByY(DoubleRange yRange)
        {
            return NewSpectrumFilterByY(yRange.Minimum, yRange.Maximum);
        }

        public new IMzSpectrum<MzPeak> NewSpectrumApplyFunctionToX(Func<double, double> convertor)
        {
            var ok = ApplyFunctionToX(convertor);
            return new DefaultMzSpectrum(ok.Item1, ok.Item2, false);
        }

        #endregion Public Methods

    }
}