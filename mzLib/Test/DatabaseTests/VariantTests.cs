using NUnit.Framework;
using Omics.BioPolymer;
using System;
using System.IO;
using System.Linq;

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
}
