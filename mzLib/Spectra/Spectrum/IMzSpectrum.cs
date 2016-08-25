// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (IMzSpectrum.cs) is part of MassSpectrometry.
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
    public interface IMzSpectrum<out TPeak> : ISpectrum<TPeak>
        where TPeak : MzPeak
    {
        new MzRange Range { get; }
        new IMzSpectrum<MzPeak> newSpectrumFilterByNumberOfMostIntense(int topNPeaks);
        new IMzSpectrum<MzPeak> newSpectrumExtract(DoubleRange xRange);
        new IMzSpectrum<MzPeak> newSpectrumExtract(double minX, double maxX);
        new IMzSpectrum<MzPeak> newSpectrumWithRangesRemoved(IEnumerable<DoubleRange> xRanges);
        new IMzSpectrum<MzPeak> newSpectrumWithRangeRemoved(DoubleRange xRange);
        new IMzSpectrum<MzPeak> newSpectrumWithRangeRemoved(double minX, double maxX);
        new IMzSpectrum<MzPeak> newSpectrumFilterByY(double minY, double maxY);
        new IMzSpectrum<MzPeak> newSpectrumFilterByY(DoubleRange yRange);
        new IMzSpectrum<MzPeak> newSpectrumApplyFunctionToX(Func<double, double> convertor);
    }
}