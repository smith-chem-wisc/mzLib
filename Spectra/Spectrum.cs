﻿// Copyright 2016 Stefan Solnts// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work copyright 2016 Stefan Solntsev
//
// This file (Spectrum.cs) is part of MassSpectrometry.
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

using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Spectra
{
    /// <summary>
    /// Spectrum that is defined by its arrays
    /// </summary>
    /// <typeparam name="TPeak"></typeparam>
    public abstract class Spectrum<TPeak> : ISpectrum<TPeak>
        where TPeak : IPeak
    {
        #region Private Fields

        private TPeak[] peakList;
        private int? indexOfpeakWithHighestY;
        private double? sumOfAllY;

        #endregion Private Fields

        #region Protected Constructors

        protected Spectrum(double[] x, double[] y, bool shouldCopy)
        {
            if (shouldCopy)
            {
                XArray = new double[x.Length];
                YArray = new double[y.Length];
                Array.Copy(x, XArray, x.Length);
                Array.Copy(y, YArray, y.Length);
            }
            else
            {
                XArray = x;
                YArray = y;
            }
            peakList = new TPeak[Size];
        }

        protected Spectrum(double[,] xy)
        {
            var count = xy.GetLength(1);

            XArray = new double[count];
            YArray = new double[count];
            Buffer.BlockCopy(xy, 0, XArray, 0, sizeof(double) * count);
            Buffer.BlockCopy(xy, sizeof(double) * count, YArray, 0, sizeof(double) * count);
            peakList = new TPeak[Size];
        }

        #endregion Protected Constructors

        #region Public Properties

        public double[] XArray { get; private set; }
        public double[] YArray { get; private set; }

        public double? FirstX
        {
            get
            {
                if (Size == 0)
                    return null;
                return XArray[0];
            }
        }

        public double? LastX
        {
            get
            {
                if (Size == 0)
                    return null;
                return XArray[Size - 1];
            }
        }

        public int Size { get { return XArray.Length; } }

        public int? IndexOfPeakWithHighesetY
        {
            get
            {
                if (Size == 0)
                    return null;
                if (!indexOfpeakWithHighestY.HasValue)
                    indexOfpeakWithHighestY = Array.IndexOf(YArray, YArray.Max());
                return indexOfpeakWithHighestY.Value;
            }
        }

        public double? YofPeakWithHighestY
        {
            get
            {
                if (Size == 0)
                    return null;
                return YArray[IndexOfPeakWithHighesetY.Value];
            }
        }

        public double? XofPeakWithHighestY
        {
            get
            {
                if (Size == 0)
                    return null;
                return XArray[IndexOfPeakWithHighesetY.Value];
            }
        }

        public double SumOfAllY
        {
            get
            {
                if (!sumOfAllY.HasValue)
                    sumOfAllY = YArray.Sum();
                return sumOfAllY.Value;
            }
        }

        public DoubleRange Range
        {
            get
            {
                if (Size == 0)
                    return null;
                return new DoubleRange(FirstX.Value, LastX.Value);
            }
        }

        #endregion Public Properties

        #region Public Methods

        public void ReplaceXbyApplyingFunction(Func<IPeak, double> convertor)
        {
            for (int i = 0; i < Size; i++)
                XArray[i] = convertor(GetPeak(i));
            peakList = new TPeak[Size];
        }

        public virtual double[,] CopyTo2DArray()
        {
            double[,] data = new double[2, Size];
            const int size = sizeof(double);
            Buffer.BlockCopy(XArray, 0, data, 0, size * Size);
            Buffer.BlockCopy(YArray, 0, data, size * Size, size * Size);
            return data;
        }

        public int? GetClosestPeakIndex(double x)
        {
            if (Size == 0)
                return null;
            int index = Array.BinarySearch(XArray, x);
            if (index >= 0)
                return index;
            index = ~index;

            if (index >= Size)
                return index - 1;
            if (index == 0)
                return index;

            if (x - XArray[index - 1] > XArray[index] - x)
                return index;
            return index - 1;
        }

        public double? GetClosestPeakXvalue(double x)
        {
            if (Size == 0)
                return null;
            return XArray[GetClosestPeakIndex(x).Value];
        }

        public int NumPeaksWithinRange(double minX, double maxX)
        {
            int startingIndex = Array.BinarySearch(XArray, minX);
            if (startingIndex < 0)
                startingIndex = ~startingIndex;
            if (startingIndex >= Size)
                return 0;
            int endIndex = Array.BinarySearch(XArray, maxX);
            if (endIndex < 0)
                endIndex = ~endIndex;
            if (endIndex == 0)
                return 0;

            return endIndex - startingIndex;
        }

        public IEnumerable<TPeak> FilterByNumberOfMostIntense(int topNPeaks)
        {
            var quantile = 1.0 - (double)topNPeaks / Size;
            quantile = Math.Max(0, quantile);
            quantile = Math.Min(1, quantile);
            double cutoffYvalue = YArray.Quantile(quantile);

            for (int i = 0; i < Size; i++)
                if (YArray[i] >= cutoffYvalue)
                    yield return GetPeak(i);
        }

        public IEnumerable<TPeak> Extract(DoubleRange xRange)
        {
            return Extract(xRange.Minimum, xRange.Maximum);
        }

        public IEnumerable<TPeak> Extract(double minX, double maxX)
        {
            int ind = Array.BinarySearch(XArray, minX);
            if (ind < 0)
                ind = ~ind;
            while (ind < Size && XArray[ind] <= maxX)
            {
                yield return GetPeak(ind);
                ind++;
            }
        }

        public IEnumerable<int> ExtractIndices(double minX, double maxX)
        {
            int ind = Array.BinarySearch(XArray, minX);
            if (ind < 0)
                ind = ~ind;
            while (ind < Size && XArray[ind] <= maxX)
            {
                yield return ind;
                ind++;
            }
        }

        public IEnumerable<TPeak> FilterByY(double minY, double maxY)
        {
            for (int i = 0; i < Size; i++)
                if (YArray[i] >= minY && YArray[i] <= maxY)
                    yield return GetPeak(i);
        }

        public IEnumerable<TPeak> FilterByY(DoubleRange yRange)
        {
            return FilterByY(yRange.Minimum, yRange.Maximum);
        }

        #endregion Public Methods

        #region Protected Methods

        protected abstract TPeak GeneratePeak(int index);

        protected TPeak GetPeak(int index)
        {
            if (peakList[index] == null)
                peakList[index] = GeneratePeak(index);
            return peakList[index];
        }

        #endregion Protected Methods
    }
}