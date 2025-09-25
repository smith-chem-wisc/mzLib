using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymer;

namespace Test
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

                Assert.That(vcf.Description, Is.EqualTo(originalLine), $"Row {rowIndex + 1}: Description mismatch.");
                Assert.That(vcf.ReferenceAlleleString, Is.EqualTo(rawFields[3]), $"Row {rowIndex + 1}: REF mismatch.");
                Assert.That(vcf.AlternateAlleleString, Is.EqualTo(rawFields[4]), $"Row {rowIndex + 1}: ALT mismatch.");
                Assert.That(vcf.Format, Is.EqualTo(rawFields[8]), $"Row {rowIndex + 1}: FORMAT mismatch.");

                if (rawFields[7] == ".")
                {
                    Assert.That(vcf.Info.Annotation, Is.EqualTo(rawFields[7]), $"Row {rowIndex + 1}: INFO mismatch.");
                }

                var sampleFields = rawFields.Skip(9).ToArray();
                Assert.That(vcf.Genotypes.Count, Is.EqualTo(sampleFields.Length), $"Row {rowIndex + 1}: genotype count mismatch.");
                Assert.That(vcf.AlleleDepths.Count, Is.EqualTo(sampleFields.Length), $"Row {rowIndex + 1}: AD count mismatch.");
                Assert.That(vcf.Homozygous.Count, Is.EqualTo(sampleFields.Length), $"Row {rowIndex + 1}: Homozygous count mismatch.");
                Assert.That(vcf.Heterozygous.Count, Is.EqualTo(sampleFields.Length), $"Row {rowIndex + 1}: Heterozygous count mismatch.");

                for (int sampleIndex = 0; sampleIndex < sampleFields.Length; sampleIndex++)
                {
                    string sample = sampleFields[sampleIndex];
                    string key = sampleIndex.ToString();

                    string[] parts = sample.Split(':');
                    Assert.That(parts.Length, Is.EqualTo(vcf.Format.Split(':').Length),
                        $"Row {rowIndex + 1}, Sample {sampleIndex + 1}: FORMAT parts mismatch.");

                    string gtPart = parts[0];
                    string adPart = parts.Length > 1 ? parts[1] : null;

                    string[] expectedGtTokens = gtPart.Split(new[] { '/', '|' }, StringSplitOptions.RemoveEmptyEntries);

                    if (gtPart.Contains('.') && expectedGtTokens.Length == 1 &&
                        (gtPart == "./." || gtPart == ".|." || gtPart == ".|1" || gtPart == "0|." || gtPart == "0/."))
                    {
                        expectedGtTokens = new[] { ".", "." };
                    }

                    Assert.That(vcf.Genotypes.ContainsKey(key), $"Row {rowIndex + 1}, Sample {sampleIndex + 1}: genotype key missing.");
                    var parsedGt = vcf.Genotypes[key];
                    Assert.That(parsedGt, Is.EqualTo(expectedGtTokens),
                        $"Row {rowIndex + 1}, Sample {sampleIndex + 1}: GT mismatch.");

                    string[] expectedAdTokens;
                    if (string.IsNullOrWhiteSpace(adPart))
                    {
                        expectedAdTokens = Array.Empty<string>();
                    }
                    else if (adPart == ".")
                    {
                        expectedAdTokens = new[] { "." };
                    }
                    else
                    {
                        expectedAdTokens = adPart.Split(',');
                    }

                    Assert.That(vcf.AlleleDepths.ContainsKey(key), $"Row {rowIndex + 1}, Sample {sampleIndex + 1}: AD key missing.");
                    var parsedAd = vcf.AlleleDepths[key] ?? Array.Empty<string>();

                    if (!(parsedAd.Length == 0 && expectedAdTokens.Length == 1 && expectedAdTokens[0] == "."))
                    {
                        Assert.That(parsedAd, Is.EqualTo(expectedAdTokens),
                            $"Row {rowIndex + 1}, Sample {sampleIndex + 1}: AD mismatch.");
                    }

                    int distinctAlleles = parsedGt.Distinct().Count();
                    bool expectedHom = distinctAlleles == 1;
                    bool expectedHet = distinctAlleles > 1;

                    Assert.That(vcf.Homozygous[key], Is.EqualTo(expectedHom),
                        $"Row {rowIndex + 1}, Sample {sampleIndex + 1}: Homozygous flag mismatch.");
                    Assert.That(vcf.Heterozygous[key], Is.EqualTo(expectedHet),
                        $"Row {rowIndex + 1}, Sample {sampleIndex + 1}: Heterozygous flag mismatch.");
                }
            }
        }
    }
}