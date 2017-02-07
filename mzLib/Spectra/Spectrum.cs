// Copyright 2016 Stefan Solnts// Copyright 2012, 2013, 2014 Derek J. Bailey
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
using System.Collections;
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

        protected readonly double[] XArray;
        protected readonly double[] YArray;

        private double? yofPeakWithHighestY;
        private double? sumOfAllY;
        protected TPeak[] peakList;
        protected TPeak peakWithHighestY;

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
            int length = xy.GetLength(1);

            XArray = new double[count];
            YArray = new double[count];
            Buffer.BlockCopy(xy, 0, XArray, 0, sizeof(double) * count);
            Buffer.BlockCopy(xy, sizeof(double) * length, YArray, 0, sizeof(double) * count);
            peakList = new TPeak[Size];
        }

        #endregion Protected Constructors

        #region Public Properties

        public double FirstX { get { return XArray[0]; } }

        public double LastX { get { return XArray[Size - 1]; } }

        public int Size { get { return XArray.Length; } }

        public double YofPeakWithHighestY
        {
            get
            {
                if (!yofPeakWithHighestY.HasValue)
                    yofPeakWithHighestY = YArray.Max();
                return yofPeakWithHighestY.Value;
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
                return new DoubleRange(FirstX, LastX);
            }
        }

        public TPeak PeakWithHighestY
        {
            get
            {
                if (EqualityComparer<TPeak>.Default.Equals(peakWithHighestY, default(TPeak)))
                    peakWithHighestY = this[Array.IndexOf(YArray, YArray.Max())];
                return peakWithHighestY;
            }
        }

        #endregion Public Properties

        #region Public Indexers

        public virtual TPeak this[int index]
        {
            get
            {
                if (peakList[index] == null)
                    peakList[index] = GeneratePeak(index);
                return peakList[index];
            }
        }

        #endregion Public Indexers

        #region Public Methods

        public override string ToString()
        {
            return string.Format("{0} (Peaks {1})", Range, Size);
        }

        public virtual double[,] CopyTo2DArray()
        {
            double[,] data = new double[2, Size];
            const int size = sizeof(double);
            Buffer.BlockCopy(XArray, 0, data, 0, size * Size);
            Buffer.BlockCopy(YArray, 0, data, size * Size, size * Size);
            return data;
        }

        public TPeak GetClosestPeak(double x)
        {
            return this[GetClosestPeakIndex(x)];
        }

        public double GetClosestPeakXvalue(double x)
        {
            return XArray[GetClosestPeakIndex(x)];
        }

        public IEnumerator<TPeak> GetEnumerator()
        {
            for (int i = 0; i < Size; i++)
                yield return this[i];
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return GetEnumerator();
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
            double cutoffYvalue = YArray.Quantile(1.0 - (double)topNPeaks / Size);

            for (int i = 0; i < Size; i++)
                if (YArray[i] >= cutoffYvalue)
                    yield return this[i];
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
            while (this[ind].X <= maxX && ind < Size)
            {
                yield return this[ind];
                ind++;
            }
        }

        public IEnumerable<TPeak> FilterByY(double minY, double maxY)
        {
            for (int i = 0; i < Size; i++)
                if (YArray[i] >= minY && YArray[i] <= maxY)
                    yield return this[i];
        }

        public IEnumerable<TPeak> FilterByY(DoubleRange yRange)
        {
            return FilterByY(yRange.Minimum, yRange.Maximum);
        }

        #endregion Public Methods

        #region Protected Methods

        protected abstract TPeak GeneratePeak(int index);


        #endregion Protected Methods

        #region Private Methods

        private int GetClosestPeakIndex(double targetX)
        {
            int index = Array.BinarySearch(XArray, targetX);
            if (index >= 0)
                return index;
            index = ~index;

            if (index >= Size)
                return index - 1;
            if (index == 0)
                return index;

            if (targetX - XArray[index - 1] > XArray[index] - targetX)
                return index;
            return index - 1;
        }

        #endregion Private Methods

    }
}