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

using MzLibUtil;
using System;
using System.Collections.Generic;

namespace Spectra
{
    /// <summary>
    /// Spectrum that has at least an X array and a Y array
    /// </summary>
    /// <typeparam name="TPeak"></typeparam>
    public interface ISpectrum<out TPeak> : IEnumerable<TPeak>
        where TPeak : IPeak
    {

        #region Public Properties
        
        double FirstX { get; }
        double LastX { get; }
        int Size { get; }
        double YofPeakWithHighestY { get; }
        double SumOfAllY { get; }
        DoubleRange Range { get; }
        TPeak PeakWithHighestY { get; }

        #endregion Public Properties

        #region Public Indexers

        TPeak this[int index] { get; }

        #endregion Public Indexers

        #region Public Methods

        double[,] CopyTo2DArray();

        int NumPeaksWithinRange(double minX, double maxX);

        TPeak GetClosestPeak(double x);

        double GetClosestPeakXvalue(double x);

        IEnumerable<TPeak> FilterByNumberOfMostIntense(int topNPeaks);

        IEnumerable<TPeak> Extract(DoubleRange xRange);

        IEnumerable<TPeak> Extract(double minX, double maxX);

        IEnumerable<TPeak> FilterByY(double minY, double maxY);

        IEnumerable<TPeak> FilterByY(DoubleRange yRange);

        #endregion Public Methods

    }
}