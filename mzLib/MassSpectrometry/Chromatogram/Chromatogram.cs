// Copyright 2012, 2013, 2014 Derek J. Bailey
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

using Spectra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    public class Chromatogram : Chromatogram<ChromatographicPeak>
    {
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

        public Chromatogram CreateSmoothChromatogram(SmoothingType smoothing, int points)
        {
            switch (smoothing)
            {
                case SmoothingType.BoxCar:
                    double[] newTimes = xArray.BoxCarSmooth(points);
                    double[] newIntensities = yArray.BoxCarSmooth(points);
                    return new Chromatogram(newTimes, newIntensities, false);
                default:
                    return new Chromatogram(this);
            }
        }

        public override ChromatographicPeak GetPeak(int index)
        {
            return new ChromatographicPeak(xArray[index], yArray[index]);
        }

    }

    public abstract class Chromatogram<TPeak> : Spectrum<TPeak>
        where TPeak : Peak
    {

        public double FirstTime
        {
            get { return xArray[0]; }
        }

        public double LastTime
        {
            get { return xArray[Count - 1]; }
        }

        protected Chromatogram(double[] times, double[] intensities, bool shouldCopy = true) : base(times, intensities, shouldCopy)
        {
        }

        protected Chromatogram(double[,] timeintensities) : base(timeintensities)
        {
        }

        protected Chromatogram(Chromatogram<TPeak> other)
            : this(other.xArray, other.yArray)
        {
        }

        public abstract TPeak GetPeak(int index);

        public double[] GetTimes()
        {
            double[] times = new double[Count];
            Buffer.BlockCopy(xArray, 0, times, 0, sizeof(double) * Count);
            return times;
        }

        public double[] GetIntensities()
        {
            double[] intensities = new double[Count];
            Buffer.BlockCopy(yArray, 0, intensities, 0, sizeof(double) * Count);
            return intensities;
        }

        public double GetTime(int index)
        {
            return xArray[index];
        }

        public double GetIntensity(int index)
        {
            return yArray[index];
        }

        public virtual TPeak GetApex(DoubleRange timeRange)
        {
            return GetApex(timeRange.Minimum, timeRange.Maximum);
        }

        public virtual TPeak GetApex(double mintime, double maxTime)
        {
            int index = Array.BinarySearch(xArray, mintime);
            if (index < 0)
                index = ~index;

            if (index >= Count)
            {
                return GetPeak(Count - 1);
            }

            double maxvalue = -1; // double.negative infinity?
            int apexIndex = index;
            while (index < Count && xArray[index] <= maxTime)
            {
                double intensity = yArray[index];
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
            int index = Array.BinarySearch(xArray, mintime);
            if (index < 0)
                index = ~index;

            List<TPeak> peaks = new List<TPeak>();
            while (index < Count && xArray[index] <= maxTime)
            {
                peaks.Add(GetPeak(index));
                index++;
            }
            return new ChromatographicElutionProfile<TPeak>(peaks);
        }

        public virtual TPeak GetApex()
        {
            return PeakWithHighestY;
        }

        public TPeak FindNearestApex(double rt, int skipablePts = 1)
        {
            int index = Array.BinarySearch(xArray, rt);
            if (index < 0)
                index = ~index;

            if (index >= Count)
                index--;

            int bestApex = index;
            double apexValue = yArray[bestApex];

            int i = index - 1;
            int count = 0;
            while (i >= 0)
            {
                if (yArray[i] > apexValue)
                {
                    bestApex = i;
                    apexValue = yArray[bestApex];
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
            while (i < Count)
            {
                if (yArray[i] > apexValue)
                {
                    bestApex = i;
                    apexValue = yArray[bestApex];
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

        public DoubleRange GetPeakWidth(double time, double fraction = 0.1, int upPts = 3, double upPrecent = 1.4, double minValue = 0)
        {
            int index = Array.BinarySearch(xArray, time);
            if (index < 0)
                index = ~index;

            if (index == xArray.Length)
                index--;

            double maxTime = xArray[index];
            double minTime = maxTime;
            double threshold = Math.Max(yArray[index] * fraction, minValue);

            int count = 0;
            double localMin = yArray[index];
            for (int i = index + 1; i < Count; i++)
            {
                double peakIntensity = yArray[i];

                if (peakIntensity > localMin * upPrecent)
                {
                    // Going up
                    count++;
                    if (count > upPts)
                        break;
                    continue;
                }

                maxTime = xArray[i];

                if (peakIntensity < localMin)
                    localMin = peakIntensity;

                count = 0;

                if (peakIntensity <= threshold)
                    break;
            }

            localMin = yArray[index];
            count = 0;
            for (int i = index - 1; i >= 0; i--)
            {
                double peakIntensity = yArray[i];

                if (peakIntensity > localMin * upPrecent)
                {
                    // Going up
                    count++;
                    if (count > upPts)
                        break;

                    continue;
                }

                minTime = xArray[i];

                if (peakIntensity < localMin)
                    localMin = peakIntensity;

                count = 0;

                if (peakIntensity < threshold)
                    break;
            }

            return new DoubleRange(minTime, maxTime);
        }

        public override string ToString()
        {
            return string.Format("Count = {0:N0} TIC = {1:G4}", Count, yArray.Sum());
        }

    }
}