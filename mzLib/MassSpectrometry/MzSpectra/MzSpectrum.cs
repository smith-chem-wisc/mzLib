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

using MzLibUtil;
using Spectra;
using System;
using System.Collections.Generic;

namespace MassSpectrometry
{
    public abstract class MzSpectrum<TPeak> : Spectrum<TPeak>, IMzSpectrum<TPeak>
        where TPeak : IMzPeak
    {

        #region Protected Constructors

        protected MzSpectrum(double[,] mzintensities) : base(mzintensities)
        {
        }

        protected MzSpectrum(ISpectrum<IMzPeak> mZSpectrum) : base(mZSpectrum)
        {
        }

        protected MzSpectrum(double[] mz, double[] intensities, bool shouldCopy) : base(mz, intensities, shouldCopy)
        {
        }

        public void ReplaceXbyApplyingFunction(Func<IMzPeak, double> convertor)
        {
            for (int i = 0; i < Count; i++)
                XArray[i] = convertor(this[i]);
            ResetSpectrum();
        }
        private void ResetSpectrum()
        {
            peakList = new TPeak[Count];
            yofPeakWithHighestY = double.NaN;
            sumOfAllY = double.NaN;
            peakWithHighestY = default(TPeak);
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

        public new IMzSpectrum<TPeak> NewSpectrumFilterByNumberOfMostIntense(int topNPeaks)
        {
            var ok = FilterByNumberOfMostIntense(topNPeaks);
            return GetMzSpectrumFromTwoArrays(ok.Item1, ok.Item2, false);
        }

        public new IMzSpectrum<TPeak> NewSpectrumExtract(DoubleRange xRange)
        {
            return NewSpectrumExtract(xRange.Minimum, xRange.Maximum);
        }

        public new IMzSpectrum<TPeak> NewSpectrumExtract(double minX, double maxX)
        {
            var ok = Extract(minX, maxX);
            return GetMzSpectrumFromTwoArrays(ok.Item1, ok.Item2, false);
        }

        public new IMzSpectrum<TPeak> NewSpectrumWithRangesRemoved(IEnumerable<DoubleRange> xRanges)
        {
            var ok = WithRangesRemoved(xRanges);
            return GetMzSpectrumFromTwoArrays(ok.Item1, ok.Item2, false);
        }

        public new IMzSpectrum<TPeak> NewSpectrumWithRangeRemoved(DoubleRange xRange)
        {
            return NewSpectrumWithRangeRemoved(xRange.Minimum, xRange.Maximum);
        }

        public new IMzSpectrum<TPeak> NewSpectrumWithRangeRemoved(double minX, double maxX)
        {
            var ok = WithRangeRemoved(minX, maxX);
            return GetMzSpectrumFromTwoArrays(ok.Item1, ok.Item2, false);
        }

        public new IMzSpectrum<TPeak> NewSpectrumFilterByY(double minY, double maxY)
        {
            var ok = FilterByY(minY, maxY);
            return GetMzSpectrumFromTwoArrays(ok.Item1, ok.Item2, false);
        }

        public abstract MzSpectrum<TPeak> GetMzSpectrumFromTwoArrays(double[] item1, double[] item2, bool v);

        public new IMzSpectrum<TPeak> NewSpectrumFilterByY(DoubleRange yRange)
        {
            return NewSpectrumFilterByY(yRange.Minimum, yRange.Maximum);
        }

        public new IMzSpectrum<TPeak> NewSpectrumApplyFunctionToX(Func<double, double> convertor)
        {
            var ok = ApplyFunctionToX(convertor);
            return GetMzSpectrumFromTwoArrays(ok.Item1, ok.Item2, false);
        }

        #endregion Public Methods

    }
}