using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using System;
using System.IO;
using System.Linq;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// The "Cleavage Requirement" column of <c>proteases.tsv</c>, and the compact subsite syntax it
    /// carries. The column is the only way a glycoprotease can be configured from data rather than code,
    /// so it has to survive both a malformed value and a file written before it existed.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class CleavageRequirementParsingTests
    {
        [Test]
        [TestCase("P2:O-glycan", false, 2, GlycosylationClass.OLinked, TestName = "Parse StcE's P2 O-glycan")]
        [TestCase("P1':O-glycan", true, 1, GlycosylationClass.OLinked, TestName = "Parse the OgpA family's P1-prime O-glycan")]
        [TestCase("P4:O-glycan", false, 4, GlycosylationClass.OLinked, TestName = "Parse ZmpC's distal P4 requirement")]
        [TestCase("P1:N-glycan", false, 1, GlycosylationClass.NLinked, TestName = "Parse an N-glycan requirement")]
        [TestCase("  p1' : o-glycan ", true, 1, GlycosylationClass.OLinked, TestName = "Parse tolerates spacing and case")]
        public static void Parse_ReadsTheSubsiteAddressAsTheLiteratureWritesIt(
            string text, bool expectPrime, int expectSubsite, GlycosylationClass expectClass)
        {
            CleavageRequirement requirement = CleavageRequirementParser.Parse(text);

            Assert.IsNotNull(requirement);
            Assert.AreEqual(expectPrime, requirement.IsPrimeSide);
            Assert.AreEqual(expectSubsite, requirement.Subsite);
            Assert.AreEqual(expectClass, requirement.RequiredClass);
        }

        [Test]
        [TestCase(null)]
        [TestCase("")]
        [TestCase("   ")]
        public static void Parse_ReturnsNullForNoRequirement(string text)
        {
            // Nearly every protease cleaves on sequence alone, so an empty cell is the common case and
            // must not be an error.
            Assert.IsNull(CleavageRequirementParser.Parse(text));
        }

        [Test]
        [TestCase("O-glycan", TestName = "Reject a requirement with no subsite")]
        [TestCase("P2", TestName = "Reject a subsite with no modification class")]
        [TestCase("X2:O-glycan", TestName = "Reject a subsite that does not start with P")]
        [TestCase("P0:O-glycan", TestName = "Reject subsite zero, since subsites count from one")]
        [TestCase("Pn:O-glycan", TestName = "Reject a non-numeric subsite")]
        [TestCase("P2:phospho", TestName = "Reject a modification class that cannot be resolved")]
        public static void Parse_ThrowsOnAMalformedRequirement(string text)
        {
            // Loud, not silent. A requirement that failed to parse and was ignored would let a
            // glycoprotease digest as though it had none -- which is precisely the over-digestion this
            // whole feature exists to stop, and it would look like a data problem rather than a typo.
            Assert.Throws<MzLibException>(() => CleavageRequirementParser.Parse(text));
        }

        [Test]
        public static void TheShippedGlycoproteasesDeclareTheirRequirement()
        {
            foreach ((string name, bool prime, int subsite) in new[]
                     {
                         ("StcE", false, 2),
                         ("StcE-trypsin", false, 2),
                         ("OpeRATOR", true, 1),
                         ("IMPa", true, 1),
                         ("SmE", true, 1),
                     })
            {
                Assert.IsTrue(ProteaseDictionary.Dictionary.ContainsKey(name), name + " is missing from proteases.tsv");
                Protease protease = ProteaseDictionary.Dictionary[name];

                Assert.IsTrue(protease.HasCleavageRequirement, name + " must declare a cleavage requirement");

                CleavageRequirement requirement = protease.DigestionMotifs
                    .SelectMany(m => m.CleavageRequirements).First(r => !r.IsForbidden);

                Assert.AreEqual(prime, requirement.IsPrimeSide, name + " subsite side");
                Assert.AreEqual(subsite, requirement.Subsite, name + " subsite number");
                Assert.AreEqual(GlycosylationClass.OLinked, requirement.RequiredClass, name + " required class");
            }
        }

        [Test]
        public static void InACompositeTheRequirementAttachesOnlyToTheMotifsThatHaveThatSubsite()
        {
            // StcE-trypsin is the case worth pinning: its four StcE motifs need the glycan at P2, and its
            // two tryptic motifs must NOT inherit it. If they did, a tryptic cut would need a glycan and
            // trypsin would effectively be switched off inside the composite.
            Protease composite = ProteaseDictionary.Dictionary["StcE-trypsin"];

            var withRequirement = composite.DigestionMotifs.Where(m => m.HasCleavageRequirement).ToList();
            var withoutRequirement = composite.DigestionMotifs.Where(m => !m.HasCleavageRequirement).ToList();

            Assert.AreEqual(4, withRequirement.Count, "the four StcE motifs must carry the requirement");
            Assert.AreEqual(2, withoutRequirement.Count, "the two tryptic motifs must not");
            CollectionAssert.AreEquivalent(new[] { "K", "R" },
                withoutRequirement.Select(m => m.InducingCleavage).ToList(),
                "K| and R| are the motifs that must stay requirement-free");
        }

        [Test]
        public static void AProteaseFileWrittenBeforeTheColumnExistedStillLoads()
        {
            // The trap this test exists for: making the column REQUIRED would break every custom
            // proteases.tsv already in the wild, because none of them has it. Read through GetFieldValue,
            // a missing column is simply absent and the protease declares no requirement.
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "legacy-proteases-" + Guid.NewGuid().ToString("N") + ".tsv");

            File.WriteAllLines(path, new[]
            {
                "Name\tSequences Inducing Cleavage\tCleavage Specificity\tPSI-MS Accession Number\tPSI-MS Name",
                "legacy-no-requirement-column\tK|\tfull\t\tTrypsin",
            });

            try
            {
                var loaded = ProteaseDictionary.LoadProteaseDictionary(path);

                Assert.IsTrue(loaded.ContainsKey("legacy-no-requirement-column"));
                Protease protease = loaded["legacy-no-requirement-column"];
                Assert.IsFalse(protease.HasCleavageRequirement,
                    "a file with no requirement column must produce proteases that require nothing");
                Assert.IsTrue(protease.DigestionMotifs.All(m => !m.HasCleavageRequirement));
            }
            finally
            {
                File.Delete(path);
            }
        }
    }
}
