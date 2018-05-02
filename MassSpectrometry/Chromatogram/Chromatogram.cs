﻿// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (Chromatogram.cs) is part of MassSpectrometry.
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
using System.Linq;

namespace MassSpectrometry
{
    public class Chromatogram : Chromatogram<ChromatographicPeak>
    {
        #region Public Constructors

        public Chromatogram(double[] times, double[] intensities, bool shouldCopy)
            : base(times, intensities, shouldCopy)
        {
        }

        public Chromatogram(double[,] timeintensities)
            : base(timeintensities)
        {
        }

        public Chromatogram(Chromatogram other)
            : base(other)
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public Chromatogram CreateSmoothChromatogram(SmoothingType smoothing, int points)
        {
            switch (smoothing)
            {
                case SmoothingType.BoxCar:
                    double[] newTimes = XArray.BoxCarSmooth(points);
                    double[] newIntensities = YArray.BoxCarSmooth(points);
                    return new Chromatogram(newTimes, newIntensities, false);

                default:
                    return new Chromatogram(this);
            }
        }

        #endregion Public Methods

        #region Protected Methods

        protected override ChromatographicPeak GeneratePeak(int index)
        {
            return new ChromatographicPeak(XArray[index], YArray[index]);
        }

        #endregion Protected Methods
    }

    public abstract class Chromatogram<TPeak> : Spectrum<TPeak>
        where TPeak : IPeak
    {
        #region Protected Constructors

        protected Chromatogram(double[] times, double[] intensities, bool shouldCopy) : base(times, intensities, shouldCopy)
        {
        }

        protected Chromatogram(double[,] timeintensities) : base(timeintensities)
        {
        }

        protected Chromatogram(Chromatogram<TPeak> other)
            : this(other.XArray, other.YArray, true)
        {
        }

        #endregion Protected Constructors

        #region Public Properties

        public double FirstTime
        {
            get { return XArray[0]; }
        }

        public double LastTime
        {
            get { return XArray[Size - 1]; }
        }

        #endregion Public Properties

        #region Public Methods

        public double[] GetTimes()
        {
            double[] times = new double[Size];
            Buffer.BlockCopy(XArray, 0, times, 0, sizeof(double) * Size);
            return times;
        }

        public double[] GetIntensities()
        {
            double[] intensities = new double[Size];
            Buffer.BlockCopy(YArray, 0, intensities, 0, sizeof(double) * Size);
            return intensities;
        }

        public double GetTime(int index)
        {
            return XArray[index];
        }

        public double GetIntensity(int index)
        {
            return YArray[index];
        }

        public virtual TPeak GetApex(DoubleRange timeRange)
        {
            return GetApex(timeRange.Minimum, timeRange.Maximum);
        }

        public virtual TPeak GetApex(double mintime, double maxTime)
        {
            int index = Array.BinarySearch(XArray, mintime);
            if (index < 0)
                index = ~index;

            if (index >= Size)
            {
                return GetPeak(Size - 1);
            }

            double maxvalue = -1; // double.negative infinity?
            int apexIndex = index;
            while (index < Size && XArray[index] <= maxTime)
            {
                double intensity = YArray[index];
                if (intensity > maxvalue)
                {
                    apexIndex = index;
                    maxvalue = intensity;
                }
                index++;
            }
            return GetPeak(apexIndex);
        }

        public virtual ChromatographicElutionProfile<TPeak> GetElutionProfile(DoubleRange timeRange)
        {
            return GetElutionProfile(timeRange.Minimum, timeRange.Maximum);
        }

        public virtual ChromatographicElutionProfile<TPeak> GetElutionProfile(double mintime, double maxTime)
        {
            int index = Array.BinarySearch(XArray, mintime);
            if (index < 0)
                index = ~index;

            List<TPeak> peaks = new List<TPeak>();
            while (index < Size && XArray[index] <= maxTime)
            {
                peaks.Add(GetPeak(index));
                index++;
            }
            return new ChromatographicElutionProfile<TPeak>(peaks);
        }

        public virtual TPeak GetApex()
        {
            if (Size == 0)
                return default(TPeak);
            return GetPeak(IndexOfPeakWithHighesetY.Value);
        }

        public TPeak FindNearestApex(double rt, int skipablePts)
        {
            int index = Array.BinarySearch(XArray, rt);
            if (index < 0)
                index = ~index;

            if (index >= Size)
                index--;

            int bestApex = index;
            double apexValue = YArray[bestApex];

            int i = index - 1;
            int count = 0;
            while (i >= 0)
            {
                if (YArray[i] > apexValue)
                {
                    bestApex = i;
                    apexValue = YArray[bestApex];
                    count = 0;
                }
                else
                {
                    count++;
                    if (count >= skipablePts)
                        break;
                }
                i--;
            }

            i = index + 1;
            count = 0;
            while (i < Size)
            {
                if (YArray[i] > apexValue)
                {
                    bestApex = i;
                    apexValue = YArray[bestApex];
                    count = 0;
                }
                else
                {
                    count++;
                    if (count >= skipablePts)
                        break;
                }
                i++;
            }

            return GetPeak(bestApex);
        }

        public override string ToString()
        {
            return string.Format("Count = {0:N0} TIC = {1:G4}", Size, YArray.Sum());
        }

        #endregion Public Methods
    }
}