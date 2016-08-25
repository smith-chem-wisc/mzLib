// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (ISpectrum.cs) is part of MassSpectrometry.
//
// MassSpectrometry is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MassSpectrometry is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY ors
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
// License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with MassSpectrometry. If not, see <http://www.gnu.org/licenses/>.

using System;
using System.Collections.Generic;

namespace Spectra
{
    public interface ISpectrum<out TPeak> : IEnumerable<TPeak>
        where TPeak : Peak
    {
        double[] xArray { get; }
        double[] yArray { get; }
        double FirstX { get; }
        double LastX { get; }
        int Count { get; }
        double YofPeakWithHighestY { get; }
        double SumOfAllY { get; }
        DoubleRange Range { get; }
        TPeak PeakWithHighestY { get; }
        TPeak this[int index] { get; }
        double[,] CopyTo2DArray();
        void replaceXbyApplyingFunction(Func<TPeak, double> convertor);
        int NumPeaksWithinRange(double minX, double maxX);
        TPeak GetClosestPeak(double x);
        double GetClosestPeakXvalue(double x);
        ISpectrum<Peak> newSpectrumFilterByNumberOfMostIntense(int topNPeaks);
        ISpectrum<Peak> newSpectrumExtract(DoubleRange xRange);
        ISpectrum<Peak> newSpectrumExtract(double minX, double maxX);
        ISpectrum<Peak> newSpectrumWithRangesRemoved(IEnumerable<DoubleRange> xRanges);
        ISpectrum<Peak> newSpectrumWithRangeRemoved(DoubleRange xRange);
        ISpectrum<Peak> newSpectrumWithRangeRemoved(double minX, double maxX);
        ISpectrum<Peak> newSpectrumFilterByY(double minY, double maxY);
        ISpectrum<Peak> newSpectrumFilterByY(DoubleRange yRange);
        ISpectrum<Peak> newSpectrumApplyFunctionToX(Func<double, double> convertor);
    }
}