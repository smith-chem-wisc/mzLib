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
using System.Text.RegularExpressions;

namespace MzLibUtil
{
    public static class ClassExtensions
    {
        /// <summary>
        /// Applies a boxcar smoothing algorithm to the input data.
        /// </summary>
        /// <param name="data">The input data.</param>
        /// <param name="points">The number of points to use for smoothing. Must be an odd number.</param>
        /// <returns>The smoothed data.</returns>
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

        /// <summary>
        /// Returns a subarray of the input array.
        /// </summary>
        /// <typeparam name="T">The type of the array elements.</typeparam>
        /// <param name="data">The input array.</param>
        /// <param name="index">The starting index of the subarray.</param>
        /// <param name="length">The length of the subarray.</param>
        /// <returns>The subarray.</returns>
        public static T[] SubArray<T>(this T[] data, int index, int length)
        {
            T[] result = new T[length];
            Array.Copy(data, index, result, 0, length);
            return result;
        }

        public static bool ToEnum<T>(this int modeInt, out T result) where T : Enum
        {
            Type enumType = typeof(T);
            if (!Enum.IsDefined(enumType, modeInt))
            {
                result = default(T);
                return false;
            }
            result = (T)Enum.ToObject(enumType, modeInt);
            return true;
        }

        /// <summary>
        /// Checks if two collections are equivalent, regardless of the order of their contents.
        /// </summary>
        /// <typeparam name="T">The type of the collection elements.</typeparam>
        /// <param name="list1">The first collection.</param>
        /// <param name="list2">The second collection.</param>
        /// <returns>True if the collections are equivalent, false otherwise.</returns>
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
        /// Determines if all items in the collection are equal.
        /// </summary>
        /// <typeparam name="T">The type of the collection elements.</typeparam>
        /// <param name="list">The collection to check.</param>
        /// <returns>True if all items in the collection are equal, false otherwise.</returns>
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
        /// Finds the index of all instances of a specified substring within the source string.
        /// </summary>
        /// <param name="sourceString">The source string to be searched.</param>
        /// <param name="subString">The substring to be located.</param>
        /// <returns>An enumerable of the indices of the substring.</returns>
        public static IEnumerable<int> IndexOfAll(this string sourceString, string subString)
        {
            return Regex.Matches(sourceString, subString).Cast<Match>().Select(m => m.Index);
        }

        /// <summary>
        /// Gets the filename without extension in a period-tolerant manner.
        /// </summary>
        /// <param name="filePath">The file path.</param>
        /// <returns>The filename without extension.</returns>
        public static string GetPeriodTolerantFilenameWithoutExtension(this string filePath)
        {
            return PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(filePath);
        }

        /// <summary>
        /// Determines whether the collection is null or contains no elements.
        /// </summary>
        /// <typeparam name="T">The IEnumerable type.</typeparam>
        /// <param name="enumerable">The enumerable, which may be null or empty.</param>
        /// <returns>
        ///     <c>true</c> if the IEnumerable is null or empty; otherwise, <c>false</c>.
        /// </returns>
        public static bool IsNullOrEmpty<T>(this IEnumerable<T> enumerable)
        {
            if (enumerable == null)
            {
                return true;
            }
            /* If this is a list, use the Count property.
             * The Count property is O(1) while IEnumerable.Count() is O(N). */
            if (enumerable is ICollection<T> collection)
            {
                return collection.Count < 1;
            }
            return !enumerable.Any();
        }

        /// <summary>
        /// Determines whether the dictionary is null or contains no elements.
        /// </summary>
        /// <param name="dictionary">The dictionary, which may be null or empty.</param>
        /// <returns>
        ///     <c>true</c> if the dictionary is null or empty; otherwise, <c>false</c>.
        /// </returns>
        public static bool IsNullOrEmpty<TKey, TValue>(this IDictionary<TKey, TValue> dictionary)
        {
            return dictionary == null || dictionary.Count < 1;
        }
        /// <summary>
        /// Converts a string to a nullable double.
        /// </summary>
        /// <param name="value">The string value.</param>
        /// <returns>The nullable double value, or null if the conversion fails.</returns>
        public static double? ToNullableDouble(this string value)
        {
            if (double.TryParse(value, out var result))
            {
                return result;
            }
            return null;
        }

        /// <summary>
        /// Converts a string to a nullable integer.
        /// </summary>
        /// <param name="value">The string value.</param>
        /// <returns>The nullable integer value, or null if the conversion fails.</returns>
        public static int? ToNullableInt(this string value)
        {
            if (int.TryParse(value, out var result))
            {
                return result;
            }
            return null;
        }
    }
}