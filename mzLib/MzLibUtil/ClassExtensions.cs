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
using System.Runtime.CompilerServices;
using System.Runtime.ConstrainedExecution;
using System.Text.RegularExpressions;

namespace MzLibUtil
{
    public static class ClassExtensions
    {
        /// <summary>
        /// Parses the full sequence to identify mods. Note: This method has been updated to NOT handle ambiguous mods on a given position (e.g. M[modA]|[modB]).
        /// If ambiguity exists, generate a separate full sequence for each mod and parse each separately.
        /// </summary>
        /// <param name="fullSequence"> Full sequence of the peptide in question.</param>
        /// <param name="ignoreTerminusMod"> If true, terminal modifications will be ignored.</param>
        /// <returns> Dictionary with the key being the amino acid position of the mod and the value being the string representing the mod</returns>
        public static Dictionary<int, string> ParseModifications(this string fullSequence, bool ignoreTerminusMod = false)
        {
            // use a regex to get modifications
            string modPattern = @"-?\[(.+?)\](?<!\[I+\])"; //The "look-behind" condition prevents matching ] for metal ion modifications
            Regex modRegex = new(modPattern);

            var fullSeq = fullSequence;
            Dictionary<int, string> modDict = new();

            MatchCollection matches = modRegex.Matches(fullSeq);
            int captureLengthSum = 0;
            int positionToAddToDict = 0;
            foreach (Match match in matches)
            {
                GroupCollection group = match.Groups;
                string rawModString = group[0].Value;
                string mod = group[1].Value;
                int startIndex = group[0].Index;
                int captureLength = group[0].Length;

                // The position of the amino acids is tracked by the positionToAddToDict variable. It takes the 
                // startIndex of the modification Match and removes the cumulative length of the modifications
                // found (including the brackets). The difference will be the number of nonmodification characters, 
                // or the number of amino acids prior to the startIndex in the sequence. 
                positionToAddToDict = startIndex - captureLengthSum;

                if (((positionToAddToDict == 0) || rawModString.StartsWith("-")) && ignoreTerminusMod) // ignore terminal mods
                {
                    captureLengthSum += captureLength;
                    continue;
                }

                if (rawModString.StartsWith("-"))
                {
                    positionToAddToDict++;
                }

                modDict.Add(positionToAddToDict, mod);
                captureLengthSum += captureLength;
            }
            return modDict;
        }

        // This method is a WIP. It is not currently used, and may be removed in the future depending on how/if we want to handle ambiguity here.
        public static Dictionary<int, string> ParseModificationsWithAmbiguity(this string ambiguousFullSequences, bool ignoreTerminusMod = false)
        {
            var modDicts = ambiguousFullSequences.Split('|').Select(fullSeq => fullSeq.ParseModifications(ignoreTerminusMod)).ToList();

            if (modDicts.Count == 1) { return modDicts[0]; }
            else
            {
                var modDict = modDicts.First();

                foreach (var md in modDicts.Skip(1))
                {
                    foreach (var mod in md)
                    {
                        if (modDict.ContainsKey(mod.Key))
                        {
                            if (!modDict[mod.Key].Split('|').Contains(mod.Value))
                            {
                                modDict[mod.Key] += "|" + mod.Value;
                            }
                        }
                        else
                        {
                            modDict.Add(mod.Key, mod.Value);
                        }
                    }
                }
                return modDict;
            }
        }

        /// <summary>
        /// Fixes an issue where the | appears and throws off the numbering if there are multiple mods on a single amino acid.
        /// </summary>
        /// <param name="fullSequence"></param>
        /// <param name="replacement"></param>
        /// <param name="specialCharacter"></param>
        /// <returns></returns>
        public static void RemoveSpecialCharacters(ref string fullSequence, string replacement = @"", string specialCharacter = @"\|")
        {
            // next regex is used in the event that multiple modifications are on a missed cleavage Lysine (K)
            Regex regexSpecialChar = new(specialCharacter);
            fullSequence = regexSpecialChar.Replace(fullSequence, replacement);
        }

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
        /// Finds the index of all instances of a specified substring within the source string.
        /// The index returned is the position of the first character of the substring within the source tring
        /// </summary>
        /// <param name="sourceString">Haystack: string to be searched</param>
        /// <param name="subString">Needle: substring to be located</param>
        public static IEnumerable<int> IndexOfAll(this string sourceString, string subString)
        {
            return Regex.Matches(sourceString, subString).Cast<Match>().Select(m => m.Index);
        }

        /// <summary>
        /// Extension method to invoke the GetPeriodTolerantFileNameWithoutExtension method
        /// </summary>
        /// <param name="filePath"></param>
        /// <returns></returns>
        public static string GetPeriodTolerantFilenameWithoutExtension(this string filePath)
        {
            return PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(filePath);
        }

    }
}