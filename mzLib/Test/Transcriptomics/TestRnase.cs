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
            // Reset to defaults first
            RnaseDictionary.ResetToDefaults();
            int originalCount = RnaseDictionary.Dictionary.Count;

            // Create a temporary custom RNase file
            string tempPath = Path.Combine(Path.GetTempPath(), "custom_rnases.tsv");
            try
            {
                File.WriteAllText(tempPath, "Name\tMotif\tSpecificity\nCustomRNase\tA|\tfull\n");
                
                var addedOrUpdated = RnaseDictionary.LoadAndMergeCustomRnases(tempPath);
                
                Assert.That(addedOrUpdated.Count, Is.EqualTo(1));
                Assert.That(addedOrUpdated[0], Is.EqualTo("CustomRNase"));
                Assert.That(RnaseDictionary.Dictionary.ContainsKey("CustomRNase"));
                Assert.That(RnaseDictionary.Dictionary.Count, Is.EqualTo(originalCount + 1));

                // Reset and verify custom RNase is gone
                RnaseDictionary.ResetToDefaults();
                Assert.That(RnaseDictionary.Dictionary.ContainsKey("CustomRNase"), Is.False);
                Assert.That(RnaseDictionary.Dictionary.Count, Is.EqualTo(originalCount));
            }
            finally
            {
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
