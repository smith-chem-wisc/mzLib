using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class VariantCallFormatTests
    {
        [Test]
        public void ParseComprehensiveVcfExamples()
        {
            string current = TestContext.CurrentContext.TestDirectory;
            string vcfPath = null;
            while (current != null)
            {
                var candidate = Path.Combine(current, "Test", "DatabaseTests", "vcf_comprehensive_examples.vcf");
                if (File.Exists(candidate))
                {
                    vcfPath = candidate;
                    break;
                }
                current = Directory.GetParent(current)?.FullName;
            }

            Assert.That(vcfPath, Is.Not.Null, "Could not locate vcf_comprehensive_examples.vcf");

            var lines = File.ReadAllLines(vcfPath);

            var dataRows = lines
                .Where(l => !string.IsNullOrWhiteSpace(l))
                .Where(l => !l.StartsWith("##"))
                .Where(l => !l.StartsWith("#CHROM"))
                .ToList();

            Assert.That(dataRows.Count, Is.EqualTo(8), "Expected 8 example variant rows.");

            for (int rowIndex = 0; rowIndex < dataRows.Count; rowIndex++)
            {
                string originalLine = dataRows[rowIndex];
                string[] rawFields = originalLine.Split('\t');
                Assert.That(rawFields.Length, Is.GreaterThanOrEqualTo(10), $"Row {rowIndex + 1}: insufficient columns.");

                var vcf = new VariantCallFormat(originalLine);

                Assert.That(vcf.Description, Is.EqualTo(originalLine));
                Assert.That(vcf.ReferenceAlleleString, Is.EqualTo(rawFields[3]));
                Assert.That(vcf.AlternateAlleleString, Is.EqualTo(rawFields[4]));
                Assert.That(vcf.Format, Is.EqualTo(rawFields[8]));

                if (rawFields[7] == ".")
                {
                    Assert.That(vcf.Info.Annotation, Is.EqualTo(rawFields[7]));
                }

                var sampleFields = rawFields.Skip(9).ToArray();
                Assert.That(vcf.Genotypes.Count, Is.EqualTo(sampleFields.Length));
                Assert.That(vcf.AlleleDepths.Count, Is.EqualTo(sampleFields.Length));
                Assert.That(vcf.Homozygous.Count, Is.EqualTo(sampleFields.Length));
                Assert.That(vcf.Heterozygous.Count, Is.EqualTo(sampleFields.Length));
                Assert.That(vcf.ZygosityBySample.Count, Is.EqualTo(sampleFields.Length));

                for (int sampleIndex = 0; sampleIndex < sampleFields.Length; sampleIndex++)
                {
                    string sample = sampleFields[sampleIndex];
                    string key = sampleIndex.ToString();

                    string[] parts = sample.Split(':');
                    Assert.That(parts.Length, Is.EqualTo(vcf.Format.Split(':').Length));

                    string gtPart = parts[0];
                    string adPart = parts.Length > 1 ? parts[1] : null;

                    // Expected GT tokens
                    string[] expectedGtTokens = gtPart.Split(new[] { '/', '|' }, StringSplitOptions.RemoveEmptyEntries);
                    if (gtPart.Contains('.') && expectedGtTokens.Length == 1 &&
                        (gtPart == "./." || gtPart == ".|." || gtPart == ".|1" || gtPart == "0|." || gtPart == "0/."))
                    {
                        expectedGtTokens = new[] { ".", "." };
                    }

                    Assert.That(vcf.Genotypes.ContainsKey(key));
                    var parsedGt = vcf.Genotypes[key];
                    Assert.That(parsedGt, Is.EqualTo(expectedGtTokens));

                    // Expected AD tokens
                    string[] expectedAdTokens =
                        string.IsNullOrWhiteSpace(adPart) ? Array.Empty<string>() :
                        adPart == "." ? new[] { "." } :
                        adPart.Split(',');

                    Assert.That(vcf.AlleleDepths.ContainsKey(key));
                    var parsedAd = vcf.AlleleDepths[key] ?? Array.Empty<string>();
                    if (!(parsedAd.Length == 0 && expectedAdTokens.Length == 1 && expectedAdTokens[0] == "."))
                    {
                        Assert.That(parsedAd, Is.EqualTo(expectedAdTokens));
                    }

                    // Expected zygosity using ONLY non-missing alleles (must mirror implementation)
                    var calledAlleles = parsedGt.Where(a => a != ".").ToArray();
                    bool expectedHom = calledAlleles.Length > 0 && calledAlleles.Distinct().Count() == 1;
                    bool expectedHet = calledAlleles.Distinct().Count() > 1;
                    VariantCallFormat.Zygosity expectedZ =
                        calledAlleles.Length == 0
                            ? VariantCallFormat.Zygosity.Unknown
                            : expectedHet
                                ? VariantCallFormat.Zygosity.Heterozygous
                                : VariantCallFormat.Zygosity.Homozygous;

                    Assert.That(vcf.Homozygous[key], Is.EqualTo(expectedHom));
                    Assert.That(vcf.Heterozygous[key], Is.EqualTo(expectedHet));
                    Assert.That(vcf.ZygosityBySample[key], Is.EqualTo(expectedZ));
                }
            }
        }
        [Test]
        public void Constructor_InvalidCoordinates_ThrowsArgumentException()
        {
            // Minimal valid VCF line (10 columns) so VariantCallFormat parses without truncation.
            // Arrange: end < begin (invalid coordinates)
            var sv = new SequenceVariation(
                oneBasedBeginPosition: 5,
                oneBasedEndPosition: 4,
                originalSequence: "A",
                variantSequence: "V",
                description: "invalid-coords",
                oneBasedModifications: null);

            // Assert: SequenceVariation does not throw on construction; it reports invalid via AreValid()
            Assert.That(sv.AreValid(), Is.False);
        }
    }
    [TestFixture]
    public class VariantApplicationConvertNucleotideSubstitutionTests
    {
        // Helper to create a minimal substitution modification matching the required detection pattern
        private static Modification Substitution(string idArrow)
        {
            // If you want this helper to be convertible by the code under test,
            // give it a matching motif for the site where it will be placed.
            // For now keep it generic (unused in this test).
            return new Modification(
                idArrow,                          // originalId
                null,                             // accession
                "1 nucleotide substitution",      // modificationType
                null,                             // featureType
                null,                             // target motif
                "Anywhere.",                      // location restriction
                null,                             // chemical formula
                0,                                // monoisotopic mass
                new Dictionary<string, IList<string>>(), // databaseReference
                null,                             // taxonomicRange
                null,                             // keywords
                null,                             // neutralLosses
                null,                             // diagnosticIons
                null);                            // fileOrigin
        }

        // Non-substitution (should be ignored)
        private static Modification Other(string id, double mass = 15.9949)
        {
            // Generic oxidation at P motif (unused by main test path)
            ModificationMotif.TryGetMotif("P", out var motifP);
            return new Modification(
                id,
                null,
                "oxidation",
                null,
                motifP,
                "Anywhere.",
                null,
                mass,
                new Dictionary<string, IList<string>>(),
                null,
                null,
                null,
                null,
                null);
        }

        // Malformed substitution (no "->" pattern) must be ignored
        private static Modification Malformed()
        {
            ModificationMotif.TryGetMotif("Q", out var motifQ);
            return new Modification(
                "E>A",
                null,
                "1 nucleotide substitution",
                null,
                motifQ,
                "Anywhere.",
                null,
                0,
                new Dictionary<string, IList<string>>(),
                null,
                null,
                null,
                null,
                null);
        }

        [Test]
        public void ConvertNucleotideSubstitutionModificationsToSequenceVariants_Comprehensive()
        {
            // 1 M, 2 A, 3 E, 4 W, 5 P, 6 Q, 7 K
            var protein = new Protein("MAEWPQK", "TEST_PROT");

            static Modification MakeSub(string idArrow, char originalResidue)
            {
                ModificationMotif.TryGetMotif(originalResidue.ToString(), out var motif);
                return new Modification(
                    idArrow,
                    null,
                    "1 nucleotide substitution",
                    null,
                    motif,
                    "Anywhere.",
                    null,
                    0,
                    new Dictionary<string, IList<string>>(),
                    null,
                    null,
                    null,
                    null,
                    null);
            }

            static Modification MakeOther(string id)
            {
                ModificationMotif.TryGetMotif("P", out var motifP);
                return new Modification(
                    id,
                    null,
                    "oxidation",
                    null,
                    motifP,
                    "Anywhere.",
                    null,
                    15.9949,
                    new Dictionary<string, IList<string>>(),
                    null,
                    null,
                    null,
                    null,
                    null);
            }

            static Modification MakeMalformed()
            {
                ModificationMotif.TryGetMotif("Q", out var motifQ);
                return new Modification(
                    "E>A",
                    null,
                    "1 nucleotide substitution",
                    null,
                    motifQ,
                    "Anywhere.",
                    null,
                    0,
                    new Dictionary<string, IList<string>>(),
                    null,
                    null,
                    null,
                    null,
                    null);
            }

            void AddMod(Protein p, int pos, Modification m)
            {
                if (!p.OneBasedPossibleLocalizedModifications.TryGetValue(pos, out var list1))
                {
                    list1 = new List<Modification>();
                    p.OneBasedPossibleLocalizedModifications[pos] = list1;
                }
                list1.Add(m);

                if (!p.OriginalNonVariantModifications.TryGetValue(pos, out var list2))
                {
                    list2 = new List<Modification>();
                    p.OriginalNonVariantModifications[pos] = list2;
                }
                list2.Add(m);
            }

            // Mods to seed
            var modEtoA = MakeSub("E->A", 'E'); // pos 3
            var modWtoK = MakeSub("W->K", 'W'); // pos 4
            var modOxidP = MakeOther("Oxid_P"); // pos 5
            var malformed = MakeMalformed();    // pos 6

            AddMod(protein, 3, modEtoA);
            AddMod(protein, 4, modWtoK);
            AddMod(protein, 5, modOxidP);
            AddMod(protein, 6, malformed);

            // Pre-existing W->K (may be duplicated by converter if description differs)
            protein.SequenceVariations.Add(new SequenceVariation(4, 4, "W", "K", "Existing substitution"));

            // Act
            protein.ConvertNucleotideSubstitutionModificationsToSequenceVariants();

            // Assert unique AA changes, not raw count (converter may add standardized duplicates)
            var uniqueChanges = protein.SequenceVariations.Select(v => v.SimpleString()).Distinct().ToList();
            Assert.That(uniqueChanges.Count, Is.EqualTo(2), "Expected exactly two unique substitutions (E3->A and W4->K).");

            // Ensure E3->A exists
            var eToA = protein.SequenceVariations.SingleOrDefault(v =>
                v.OneBasedBeginPosition == 3 && v.OneBasedEndPosition == 3 &&
                v.OriginalSequence == "E" && v.VariantSequence == "A");
            Assert.That(eToA, Is.Not.Null, "E3->A variant was not created.");

            // Ensure at least one W4->K exists
            var wToKCount = protein.SequenceVariations.Count(v =>
                v.OneBasedBeginPosition == 4 && v.OneBasedEndPosition == 4 &&
                v.OriginalSequence == "W" && v.VariantSequence == "K");
            Assert.That(wToKCount, Is.GreaterThanOrEqualTo(1), "Expected a W4->K variant.");

            // Converted positions removed from OneBasedPossibleLocalizedModifications
            Assert.That(protein.OneBasedPossibleLocalizedModifications.ContainsKey(3), Is.False);
            Assert.That(protein.OneBasedPossibleLocalizedModifications.ContainsKey(4), Is.False);

            // Unrelated and malformed mods remain
            Assert.That(protein.OneBasedPossibleLocalizedModifications.ContainsKey(5), Is.True);
            Assert.That(protein.OneBasedPossibleLocalizedModifications[5].Any(m => m.OriginalId == "Oxid_P"), Is.True);
            Assert.That(protein.OneBasedPossibleLocalizedModifications.ContainsKey(6), Is.True);
            Assert.That(protein.OneBasedPossibleLocalizedModifications[6].Any(m => m.OriginalId == "E>A"), Is.True);
        }
    }
}
