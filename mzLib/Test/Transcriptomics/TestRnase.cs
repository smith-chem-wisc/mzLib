using NUnit.Framework;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using Proteomics.ProteolyticDigestion;
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
            Assert.That(RnaseDictionary.Dictionary.ContainsKey("top-down"));
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