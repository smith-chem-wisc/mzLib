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
using System.Text;

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

        /// <summary>
        /// Determines if all items in collection are equal
        /// </summary>
        /// <typeparam name="T">type to check</typeparam>
        /// <param name="list">collection to check</param>
        /// <returns></returns>
        public static bool AllSame<T>(this IEnumerable<T> list)
        {
            var enumerable = list.ToList();
            T comparand = enumerable.First();

            bool first = true;
            foreach (T item in enumerable)
            {
                if (first) comparand = item;
                else if (!item.Equals(comparand)) return false;
                first = false;
            }

            return true;
        }
        
        /// <summary>
        /// Transcribes a DNA sequence into an RNA sequence
        /// </summary>
        /// <param name="dna">The input dna sequence</param>
        /// <param name="isCodingStrand">True if the input sequence is the coding strand, False if the input sequence is the template strand</param>
        /// <returns></returns>
        public static string Transcribe(this string dna, bool isCodingStrand = true)
        {
            var sb = new StringBuilder();
            foreach (var t in dna)
            {
                if (isCodingStrand)
                {
                    sb.Append(t == 'T' ? 'U' : t);
                }
                else
                {
                    switch (t)
                    {
                        case 'A':
                            sb.Append('U');
                            break;
                        case 'T':
                            sb.Append('A');
                            break;
                        case 'C':
                            sb.Append('G');
                            break;
                        case 'G':
                            sb.Append('C');
                            break;
                        default:
                            sb.Append(t);
                            break;
                    }
                }
            }
            return sb.ToString();
        }
    }
}