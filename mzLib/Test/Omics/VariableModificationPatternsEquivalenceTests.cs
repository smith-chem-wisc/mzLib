using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Omics;
using Omics.Digestion;
using Omics.Modifications;

namespace Test.Omics
{
    /// <summary>
    /// GetVariableModificationPatterns sizes each pattern's dictionary to the modifications it will hold and reads each
    /// position's modifications from an array instead of SortedSet.ElementAt. These tests hold it to the previous
    /// implementation, kept verbatim below, on random inputs: the same patterns in the same order, each with the same
    /// keys in the same order holding the same modification objects - also after fixed modifications are appended.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class VariableModificationPatternsEquivalenceTests
    {
        [Test]
        public static void PatternsAreThoseThePreviousCodeMade()
        {
            var random = new Random(1015);
            List<Modification> pool = Enumerable.Range(0, 12)
                .Select(i => NewModification("mod" + i, "ACDEFGHIKLMNPQRSTVWY"[i].ToString(), 1 + i)).ToList();
            long patterns = 0, modifications = 0;

            for (int trial = 0; trial < 1000; trial++)
            {
                int length = 1 + random.Next(30);
                var possible = new Dictionary<int, SortedSet<Modification>>();
                // terminal sets first and sometimes empty, as PopulateVariableModifications makes them
                possible.Add(1, RandomSet(random, pool, random.Next(3) == 0 ? 0 : random.Next(3)));
                possible.Add(length + 2, RandomSet(random, pool, random.Next(3) == 0 ? 0 : random.Next(3)));
                foreach (int position in Enumerable.Range(2, length).Where(_ => random.Next(3) == 0).OrderBy(_ => random.Next()))
                {
                    possible.Add(position, RandomSet(random, pool, 1 + random.Next(3)));
                }
                if (random.Next(20) == 0)
                {
                    possible.Clear();
                }
                int maxMods = random.Next(5);
                var fixedMods = new Dictionary<int, Modification>();
                foreach (int position in Enumerable.Range(1, length + 2).Where(_ => random.Next(6) == 0))
                {
                    fixedMods.Add(position, pool[random.Next(pool.Count)]);
                }

                List<Dictionary<int, Modification>> expected = GetVariableModificationPatternsBefore(possible, maxMods, length).ToList();
                List<Dictionary<int, Modification>> actual = Access.Patterns(possible, maxMods, length, fixedMods.Count).ToList();

                Assert.That(actual.Count, Is.EqualTo(expected.Count), $"trial {trial}");
                for (int p = 0; p < expected.Count; p++)
                {
                    AssertSame(expected[p], actual[p]);
                    Access.Append(fixedMods, expected[p], out int expectedFixed);
                    Access.Append(fixedMods, actual[p], out int actualFixed);
                    Assert.That(actualFixed, Is.EqualTo(expectedFixed));
                    AssertSame(expected[p], actual[p]);
                    modifications += actual[p].Count;
                }
                patterns += expected.Count;
            }

            Assert.That(patterns, Is.GreaterThan(30000));
            Assert.That(modifications, Is.GreaterThan(30000));
        }

        [Test]
        public static void TheOverloadWithoutCapacityIsUnchanged()
        {
            List<Modification> pool = Enumerable.Range(0, 3).Select(i => NewModification("m" + i, "K", 1 + i)).ToList();
            var possible = new Dictionary<int, SortedSet<Modification>>
            {
                { 1, new SortedSet<Modification>() },
                { 6, new SortedSet<Modification>() },
                { 3, new SortedSet<Modification>(pool) },
                { 4, new SortedSet<Modification>(pool.Take(2)) },
            };
            var expected = GetVariableModificationPatternsBefore(possible, 2, 4).ToList();
            var actual = Access.Patterns(possible, 2, 4).ToList();
            Assert.That(actual.Count, Is.EqualTo(expected.Count));
            for (int p = 0; p < expected.Count; p++)
            {
                AssertSame(expected[p], actual[p]);
            }
        }

        private static void AssertSame(Dictionary<int, Modification> expected, Dictionary<int, Modification> actual)
        {
            Assert.That(actual.Keys.ToList(), Is.EqualTo(expected.Keys.ToList()));
            foreach (int key in expected.Keys)
            {
                Assert.That(actual[key], Is.SameAs(expected[key]));
            }
        }

        private static SortedSet<Modification> RandomSet(Random random, List<Modification> pool, int count)
        {
            var set = new SortedSet<Modification>();
            while (set.Count < count)
            {
                set.Add(pool[random.Next(pool.Count)]);
            }
            return set;
        }

        private static Modification NewModification(string id, string motif, double mass)
        {
            ModificationMotif.TryGetMotif(motif, out ModificationMotif target);
            return new Modification(_originalId: id, _modificationType: "Test", _target: target, _locationRestriction: "Anywhere.", _monoisotopicMass: mass);
        }

        /// <summary>Reaches DigestionProduct's protected helpers; never instantiated.</summary>
        private abstract class Access : DigestionProduct
        {
            private Access() : base(null, 0, 0, 0, CleavageSpecificity.Full) { }

            internal static IEnumerable<Dictionary<int, Modification>> Patterns(Dictionary<int, SortedSet<Modification>> possible, int maxMods, int length, int additionalCapacity)
                => GetVariableModificationPatterns(possible, maxMods, length, additionalCapacity);

            internal static IEnumerable<Dictionary<int, Modification>> Patterns(Dictionary<int, SortedSet<Modification>> possible, int maxMods, int length)
                => GetVariableModificationPatterns(possible, maxMods, length);

            // AppendFixedModificationsToVariable, which is an instance method that reads no instance state
            internal static void Append(Dictionary<int, Modification> fixedModDict, Dictionary<int, Modification> variableModPattern, out int numFixedMods)
            {
                numFixedMods = 0;
                foreach (var fixedModPattern in fixedModDict)
                {
                    if (variableModPattern.ContainsKey(fixedModPattern.Key))
                        continue;
                    numFixedMods++;
                    variableModPattern.Add(fixedModPattern.Key, fixedModPattern.Value);
                }
            }
        }

        #region the previous implementation, verbatim

        private static IEnumerable<Dictionary<int, Modification>> GetVariableModificationPatternsBefore(Dictionary<int, SortedSet<Modification>> possibleVariableModifications, int maxModsForPeptide, int peptideLength)
        {
            if (possibleVariableModifications.Count <= 0)
                yield break;

            int[] baseVariableModificationPattern = new int[peptideLength + 4];
            int totalAvailableMods = possibleVariableModifications.Values.Sum(modList => modList?.Count ?? 0);
            int maxVariableMods = Math.Min(totalAvailableMods, maxModsForPeptide);
            var variableModKvpList = possibleVariableModifications.ToList();

            for (int variable_modifications = 0; variable_modifications <= maxVariableMods; variable_modifications++)
            {
                foreach (int[] variable_modification_pattern in GetVariableModificationPatternsRecursiveBefore(variableModKvpList,
                             possibleVariableModifications.Count - variable_modifications, baseVariableModificationPattern, 0))
                {
                    // use modification pattern to construct a dictionary of modifications for the peptide
                    var modificationPattern = new Dictionary<int, Modification>(possibleVariableModifications.Count);

                    foreach (var variableModSet in possibleVariableModifications)
                    {
                        int modIndex = variable_modification_pattern[variableModSet.Key] - 1;
                        if (modIndex >= 0)
                        {
                            modificationPattern.Add(variableModSet.Key, variableModSet.Value.ElementAt(modIndex));
                        }
                    }

                    yield return modificationPattern;
                }
            }
        }

        private static IEnumerable<int[]> GetVariableModificationPatternsRecursiveBefore(List<KeyValuePair<int, SortedSet<Modification>>> possibleVariableModifications,
            int unmodifiedResiduesDesired, int[] variableModificationPattern, int index)
        {
            if (index < possibleVariableModifications.Count - 1)
            {
                if (unmodifiedResiduesDesired > 0)
                {
                    variableModificationPattern[possibleVariableModifications[index].Key] = 0;
                    foreach (int[] new_variable_modification_pattern in GetVariableModificationPatternsRecursiveBefore(possibleVariableModifications,
                        unmodifiedResiduesDesired - 1, variableModificationPattern, index + 1))
                    {
                        yield return new_variable_modification_pattern;
                    }
                }
                if (unmodifiedResiduesDesired < possibleVariableModifications.Count - index)
                {
                    for (int i = 1; i <= possibleVariableModifications[index].Value.Count; i++)
                    {
                        variableModificationPattern[possibleVariableModifications[index].Key] = i;
                        foreach (int[] new_variable_modification_pattern in GetVariableModificationPatternsRecursiveBefore(possibleVariableModifications,
                            unmodifiedResiduesDesired, variableModificationPattern, index + 1))
                        {
                            yield return new_variable_modification_pattern;
                        }
                    }
                }
            }
            else
            {
                if (unmodifiedResiduesDesired > 0)
                {
                    variableModificationPattern[possibleVariableModifications[index].Key] = 0;
                    yield return variableModificationPattern;
                }
                else
                {
                    for (int i = 1; i <= possibleVariableModifications[index].Value.Count; i++)
                    {
                        variableModificationPattern[possibleVariableModifications[index].Key] = i;
                        yield return variableModificationPattern;
                    }
                }
            }
        }

        #endregion
    }
}
