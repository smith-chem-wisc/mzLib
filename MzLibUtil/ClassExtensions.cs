// Copyright 2012, 2013, 2014 Derek J. Bailey
// Modified work Copyright 2016 Stefan Solntsev
//
// This file (ClassExtensions.cs) is part of MassSpectrometry.
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
using System.Linq;

namespace MzLibUtil
{
    public static class ClassExtensions
    {
        public static double[] BoxCarSmooth(this double[] data, int points)
        {
            // Force to be odd
            points = points - (1 - points % 2);

            int count = data.Length;

            int newCount = count - points + 1;

            double[] smoothedData = new double[newCount];

            for (int i = 0; i < newCount; i++)
            {
                double value = 0;

                for (int j = i; j < i + points; j++)
                {
                    value += data[j];
                }

                smoothedData[i] = value / points;
            }
            return smoothedData;
        }

        public static T[] SubArray<T>(this T[] data, int index, int length)
        {
            T[] result = new T[length];
            Array.Copy(data, index, result, 0, length);
            return result;
        }

        /// <summary>
        /// Checks if two collections are equivalent, regardless of the order of their contents
        /// </summary>
        public static bool ScrambledEquals<T>(this IEnumerable<T> list1, IEnumerable<T> list2)
        {
            var cnt = new Dictionary<T, int>();
            foreach (T s in list1)
            {
                if (cnt.ContainsKey(s))
                    cnt[s]++;
                else
                    cnt.Add(s, 1);
            }
            foreach (T s in list2)
            {
                if (cnt.ContainsKey(s))
                    cnt[s]--;
                else
                    return false;
            }
            return cnt.Values.All(c => c == 0);
        }
    }
}