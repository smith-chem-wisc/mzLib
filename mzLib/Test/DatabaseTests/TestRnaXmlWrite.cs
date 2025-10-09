using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Xml;
using NUnit.Framework.Legacy;
using Transcriptomics;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal class TestRnaXmlWrite
    {
        [Test]
        public void RnaSequenceVariantDescription_Fallbacks()
        {
            // RNA: A U G C; apply U2C (position 2)
            var rna = new RNA(
                sequence: "AUGC",
                accession: "RNA0001",
                oneBasedPossibleModifications: null,
                fivePrimeTerminus: null,
                threePrimeTerminus: null,
                name: "rna1",
                organism: "Test organism",
                databaseFilePath: null,
                isContaminant: false,
                isDecoy: false,
                geneNames: new List<System.Tuple<string, string>> { System.Tuple.Create("primary", "GENE1") },
                databaseAdditionalFields: null,
                truncationProducts: null,
                sequenceVariations: new List<SequenceVariation>
                {
                    // Empty description + no VCF ? writer must synthesize a fallback (SimpleString "U2C")
                    new SequenceVariation(
                        oneBasedPosition: 2,
                        originalSequence: "U",
                        variantSequence: "C",
                        description: string.Empty,
                        variantCallFormatDataString: null,
                        oneBasedModifications: null)
                },
                appliedSequenceVariations: null,
                sampleNameForVariants: null,
                fullName: "full rna name");

            string outPath = Path.Combine(TestContext.CurrentContext.WorkDirectory, "rna_variant_write.xml");
            try
            {
                var newModRes = ProteinDbWriter.WriteXmlDatabase(
                    additionalModsToAddToNucleicAcids: new Dictionary<string, HashSet<System.Tuple<int, Modification>>>(),
                    nucleicAcidList: new List<RNA> { rna },
                    outputFileName: outPath,
                    updateTimeStamp: false);

                FileAssert.Exists(outPath, "RNA XML was not written.");

                // Parse XML and find sequence variant feature
                var doc = new XmlDocument();
                doc.Load(outPath);

                var featureNodes = doc.GetElementsByTagName("feature")
                    .Cast<XmlElement>()
                    .Where(e => e.HasAttribute("type") && e.GetAttribute("type") == "sequence variant")
                    .ToList();

                Assert.That(featureNodes, Is.Not.Empty, "No RNA sequence variant feature found in XML.");
                Assert.That(featureNodes, Has.Count.EqualTo(1), "Expected exactly one RNA sequence variant feature.");

                // There is exactly one, and its description should be "U2C" (fallback from SimpleString)
                var desc = featureNodes[0].GetAttribute("description");
                Assert.That(desc, Does.Match(@".*\S.*"), "RNA variant description should not be empty.");
                Assert.That(desc, Is.EqualTo("U2C"), "RNA variant description fallback mismatch.");
            }
            finally
            {
                if (File.Exists(outPath))
                    File.Delete(outPath);
            }
        }
    }
}