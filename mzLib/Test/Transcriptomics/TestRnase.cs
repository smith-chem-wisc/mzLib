using NUnit.Framework;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Proteomics.ProteolyticDigestion;
using Transcriptomics;
using Transcriptomics.Digestion;

namespace Test.Transcriptomics
{
    [ExcludeFromCodeCoverage]
    public class TestRnase
    {
        [Test]
        public void TestRnaseDictionaryLoading()
        {
            // Verify the dictionary loads correctly from embedded resource
            Assert.That(RnaseDictionary.Dictionary.Count, Is.GreaterThan(0));

            // Verify expected RNases are present
            Assert.That(RnaseDictionary.Dictionary.ContainsKey("RNase T1"));
            Assert.That(RnaseDictionary.Dictionary.ContainsKey("RNase A"));
            Assert.That(RnaseDictionary.Dictionary.ContainsKey("RNase 4"));
            Assert.That(RnaseDictionary.Dictionary.ContainsKey("top-down"));
        }

        [Test]
        public void TestRnase4_CleavesAfterUridineOnlyBeforePurines()
        {
            var rnase4 = RnaseDictionary.Dictionary["RNase 4"];

            var products = rnase4.GetUnmodifiedOligos(new RNA("AUGCUGA"), 0, 1, int.MaxValue)
                .ToArray();

            var distinctBaseSequences = products.Select(p => p.BaseSequence).Distinct().ToArray();
            Assert.That(distinctBaseSequences, Is.EqualTo(new[] { "AU", "GCU", "GA" }));

            var grouped = products.GroupBy(p => p.BaseSequence).ToDictionary(g => g.Key, g => g.Select(o => (o.ThreePrimeTerminus, o.FivePrimeTerminus)));
            foreach (var group in grouped)
            {
                var distinctBaseSeqs = group.Value.Select(v => v).Distinct().ToArray();
                Assert.That(distinctBaseSeqs.Length, Is.EqualTo(group.Value.Count()), $"Base sequence {group.Key} has multiple distinct termini: {string.Join(", ", distinctBaseSeqs.Select(v => $"({v.ThreePrimeTerminus}, {v.FivePrimeTerminus})"))}");
            }
        }

        [Test]
        public void TestRnase4_DoesNotCleaveAfterUridineBeforePyrimidinesOrUridine()
        {
            var rnase4 = RnaseDictionary.Dictionary["RNase 4"];

            var products = rnase4.GetUnmodifiedOligos(new RNA("AUUCUGA"), 0, 1, int.MaxValue)
                .ToArray();

            var distinctBaseSequences = products.Select(p => p.BaseSequence).Distinct().ToArray();
            Assert.That(distinctBaseSequences, Is.EqualTo(new[] { "AUUCU", "GA" }));

            var grouped = products.GroupBy(p => p.BaseSequence).ToDictionary(g => g.Key, g => g.Select(o => (o.ThreePrimeTerminus, o.FivePrimeTerminus)));
            foreach (var group in grouped)
            {
                var distinctBaseSeqs = group.Value.Select(v => v).Distinct().ToArray();
                Assert.That(distinctBaseSeqs.Length, Is.EqualTo(group.Value.Count()), $"Base sequence {group.Key} has multiple distinct termini: {string.Join(", ", distinctBaseSeqs.Select(v => $"({v.ThreePrimeTerminus}, {v.FivePrimeTerminus})"))}");
            }
        }

        /// <summary>
        /// Cusativin motif C[C]|, U|A, X|U produces three cleavage sites in GUACUG:
        ///   X|U fires after G (pos 1)
        ///   U|A fires after U (pos 2)
        ///   C[C]| and X|U both fire after C (pos 4)
        /// Resulting in four distinct base fragments: G, U, AC, UG.
        /// </summary>
        [Test]
        public void TestCusativin_CleavagePattern()
        {
            var cusativin = RnaseDictionary.Dictionary["Cusativin"];

            var products = cusativin.GetUnmodifiedOligos(new RNA("GUACUG"), 0, 1, int.MaxValue)
                .ToArray();

            var distinctBaseSequences = products.Select(p => p.BaseSequence).Distinct().ToArray();
            Assert.That(distinctBaseSequences, Is.EqualTo(new[] { "G", "U", "AC", "UG" }));
        }

        /// <summary>
        /// Cusativin does not cleave C when followed by C (C[C]| motif prevention).
        /// </summary>
        [Test]
        public void TestCusativin_DoesNotCleaveBeforeCytidine()
        {
            var cusativin = RnaseDictionary.Dictionary["Cusativin"];

            var products = cusativin.GetUnmodifiedOligos(new RNA("GCCUGA"), 0, 1, int.MaxValue)
                .ToArray();

            // C at pos 2 is followed by C at pos 3 → no cleavage there
            // X|U fires after C (pos 3) → pos 3→4, C→U
            // U|A is absent (no U followed by A)
            // C[C]| fires after C at pos 3 is NOT followed by C - wait, pos3=C, pos4=U → C[C]| fires after pos 3
            // Actually pos2=C followed by pos3=C → C[C]| does NOT fire (C followed by C = prevented)
            var distinctBaseSequences = products.Select(p => p.BaseSequence).Distinct().ToArray();
            Assert.That(distinctBaseSequences, Does.Not.Contain("GC"),
                "C followed by C should not be cleaved");
        }

        [Test]
        public void TestRnaseDictionaryCustomLoadAndMerge()
        {
            int originalCount = RnaseDictionary.Dictionary.Count;
            string tempPath = Path.Combine(Path.GetTempPath(), "custom_rnases.tsv");
            try
            {
                File.WriteAllText(tempPath, "Name\tMotif\tSpecificity\nCustomRNase\tA|\tfull\n");

                var result = RnaseDictionary.LoadAndMergeCustomRnases(tempPath);

                // CustomRNase is a new name — must be added, not skipped
                Assert.That(result.Added.Count, Is.EqualTo(1));
                Assert.That(result.Added[0], Is.EqualTo("CustomRNase"));
                Assert.That(result.Skipped, Is.Empty);
                Assert.That(RnaseDictionary.Dictionary.ContainsKey("CustomRNase"));
                Assert.That(RnaseDictionary.Dictionary.Count, Is.EqualTo(originalCount + 1));
            }
            finally
            {
                // Remove all custom entries that were added so other tests are not affected
                if (RnaseDictionary.Dictionary.ContainsKey("CustomRNase"))
                    RnaseDictionary.Dictionary.Remove("CustomRNase");

                if (File.Exists(tempPath))
                    File.Delete(tempPath);
            }
        }

        [Test]
        public void TestRnaseDictionaryCustomLoadAndMerge_CollisionWithEmbedded_IsSkipped()
        {
            int originalCount = RnaseDictionary.Dictionary.Count;
            string tempPath = Path.Combine(Path.GetTempPath(), "collision_rnases.tsv");
            try
            {
                // "RNase T1" exists in the embedded resource — must be skipped, not overwritten
                File.WriteAllText(tempPath, "Name\tMotif\tSpecificity\nRNase T1\tA|\tfull\nCustomRNase2\tA|\tfull\n");

                var originalT1 = RnaseDictionary.Dictionary["RNase T1"];
                var result = RnaseDictionary.LoadAndMergeCustomRnases(tempPath);

                Assert.That(result.Skipped, Contains.Item("RNase T1"),
                    "Embedded RNase name must appear in Skipped");
                Assert.That(result.Added, Does.Not.Contain("RNase T1"),
                    "Embedded RNase name must not appear in Added");
                Assert.That(result.Added, Contains.Item("CustomRNase2"),
                    "New RNase name must appear in Added");

                // Embedded definition must be unchanged
                Assert.That(ReferenceEquals(RnaseDictionary.Dictionary["RNase T1"], originalT1), Is.True,
                    "The embedded RNase T1 object must not have been replaced");

                // Count increases by exactly 1 (only CustomRNase2 was added)
                Assert.That(RnaseDictionary.Dictionary.Count, Is.EqualTo(originalCount + 1));
            }
            finally
            {
                RnaseDictionary.Dictionary.Remove("CustomRNase2");

                if (File.Exists(tempPath))
                    File.Delete(tempPath);
            }
        }

        [Test]
        public void TestRnaseEqualityProperties()
        {
            Rnase t1 = RnaseDictionary.Dictionary["RNase T1"];
            Rnase t1Duplicate = RnaseDictionary.Dictionary["RNase T1"];
            Rnase t2 = RnaseDictionary.Dictionary["RNase T2"];

            Assert.That(t1.ToString(), Is.EqualTo("RNase T1"));
            Assert.That(t1.Equals(t1Duplicate));
            Assert.That(t1.Equals(t1));
            Assert.That(!t1.Equals(t2));
            Assert.That(!t1.Equals(null));
            Assert.That(t1.GetHashCode(), Is.EqualTo(t1Duplicate.GetHashCode()));
            Assert.That(t1.GetHashCode(), Is.Not.EqualTo(t2.GetHashCode()));
            Assert.That(t1.Equals((object)t1Duplicate));
            Assert.That(t1.Equals((object)t1));
            Assert.That(!t1.Equals((object)t2));
            Assert.That(!t1.Equals((object)null));
            // ReSharper disable once SuspiciousTypeConversion.Global
            Assert.That(!t1.Equals((object)ProteaseDictionary.Dictionary["top-down"]));
        }
    }
}
