using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Xml.Linq;
using NUnit.Framework;
using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class ProteinDbWriterSequenceVariantFeatureTests
    {
        // Creates a modification guaranteeing a non-null IdWithMotif (needed by ProteinDbWriter)
        private static Modification CreateModWithId(string id)
        {
            var mod = new Modification(_originalId: id, _modificationType: "TestType");
            // If the implementation exposes IdWithMotif privately, try to set it via reflection
            var prop = mod.GetType().GetProperty("IdWithMotif", BindingFlags.Instance | BindingFlags.Public | BindingFlags.NonPublic);
            if (prop != null && prop.CanWrite)
            {
                prop.SetValue(mod, id, null);
            }
            // Fallback: some code paths may derive IdWithMotif from another property (OriginalId already set)
            return mod;
        }

        private static Protein MakeBaseProtein(string accession, string sequence = "MPEPTIDESEQ")
        {
            var attrs = new UniProtSequenceAttributes(
                length: sequence.Length,
                mass: 1234,
                checkSum: "CHK",
                entryModified: new DateTime(2024, 1, 1),
                sequenceVersion: 1,
                isPrecursor: true,
                fragment: UniProtSequenceAttributes.FragmentType.single);

            return new Protein(
                sequence: sequence,
                accession: accession,
                organism: "TestOrg",
                geneNames: new List<Tuple<string,string>> { Tuple.Create("primary","GENE") },
                oneBasedModifications: null,
                proteolysisProducts: new List<TruncationProduct>(),
                name: "ProtName",
                fullName: "Protein Full Name",
                isDecoy: false,
                isContaminant: false,
                databaseReferences: new List<DatabaseReference>(),
                sequenceVariations: new List<SequenceVariation>(),
                disulfideBonds: new List<DisulfideBond>(),
                spliceSites: new List<SpliceSite>(),
                databaseFilePath: null,
                uniProtSequenceAttributes: attrs,
                appliedSequenceVariations: new List<SequenceVariation>(),
                sampleNameForVariants: null);
        }

        private static Protein GetConsensusCarrier(Protein baseProtein) =>
            (baseProtein.ConsensusVariant as Protein) ?? baseProtein;

        private static XDocument WriteAndLoad(Protein baseProtein,
            string testName,
            Dictionary<string, HashSet<Tuple<int, Modification>>> extraMods = null)
        {
            var path = Path.Combine(Path.GetTempPath(),
                $"ProteinVariantWriter_{testName}_{Guid.NewGuid():N}.xml");

            ProteinDbWriter.WriteXmlDatabase(extraMods, new List<Protein> { baseProtein }, path);
            return XDocument.Parse(File.ReadAllText(path));
        }

        private static IEnumerable<XElement> VariantFeatures(XDocument doc) =>
            doc
                .Descendants()
                .Where(f => f.Name.LocalName == "feature"
                            && string.Equals((string)f.Attribute("type"), "sequence variant", StringComparison.Ordinal));

        private static XElement AssertSingleVariantFeature(XDocument doc)
        {
            var feats = VariantFeatures(doc).ToList();
            Assert.That(feats.Count, Is.EqualTo(1),
                $"Expected exactly 1 sequence variant feature, found {feats.Count}. Raw XML:\n{doc}");
            return feats[0];
        }

        private static XElement FirstChild(XElement parent, string localName) =>
            parent.Elements().FirstOrDefault(e => e.Name.LocalName == localName);

        [Test]
        public void NoSequenceVariations_ProducesNoSequenceVariantFeatures()
        {
            var prot = MakeBaseProtein("ACC_NO_VAR");
            GetConsensusCarrier(prot); // ensure access
            var doc = WriteAndLoad(prot, nameof(NoSequenceVariations_ProducesNoSequenceVariantFeatures));
            Assert.That(VariantFeatures(doc), Is.Empty);
        }

        [Test]
        public void Variation_WithExplicitDescription_UsesDescriptionUnchanged()
        {
            var prot = MakeBaseProtein("ACC_EXPLICIT");
            var carrier = GetConsensusCarrier(prot);
            carrier.SequenceVariations.Add(new SequenceVariation(3, 3, "E", "K", "ExpDesc_E3K", variantCallFormatDataString: null));

            var doc = WriteAndLoad(prot, nameof(Variation_WithExplicitDescription_UsesDescriptionUnchanged));
            Assert.That((string)AssertSingleVariantFeature(doc).Attribute("description"),
                Is.EqualTo("ExpDesc_E3K"));
        }

        [Test]
        public void Variation_NullDescription_UsesVcfDescription()
        {
            var prot = MakeBaseProtein("ACC_VCF_DESC");
            var carrier = GetConsensusCarrier(prot);
            carrier.SequenceVariations.Add(new SequenceVariation(5, 5, "T", "A", null,
                variantCallFormatDataString:
                    "1\t100\t.\tT\tA\t.\tPASS\tANN=A|missense_variant\tGT:AD:DP\t0/1:5,6:11"));

            var doc = WriteAndLoad(prot, nameof(Variation_NullDescription_UsesVcfDescription));
            var desc = (string)AssertSingleVariantFeature(doc).Attribute("description");
            Assert.That(desc, Does.Contain("1\t100\t.\tT\tA\t"));
        }

        [Test]
        public void Variation_WhitespaceDescription_PointSubstitution_SynthesizesPointCode()
        {
            var prot = MakeBaseProtein("ACC_POINT_SYN");
            var carrier = GetConsensusCarrier(prot);
            carrier.SequenceVariations.Add(new SequenceVariation(2, 2, "P", "A", "  ", variantCallFormatDataString: null));

            var doc = WriteAndLoad(prot, nameof(Variation_WhitespaceDescription_PointSubstitution_SynthesizesPointCode));
            Assert.That((string)AssertSingleVariantFeature(doc).Attribute("description"), Is.EqualTo("P2A"));
        }

        [Test]
        public void Variation_WhitespaceDescription_MultiResidueRange_SynthesizesRangeCode()
        {
            var prot = MakeBaseProtein("ACC_RANGE_SYN");
            var carrier = GetConsensusCarrier(prot);
            carrier.SequenceVariations.Add(new SequenceVariation(4, 6, "PTI", "KAA", " \t ", variantCallFormatDataString: null));

            var doc = WriteAndLoad(prot, nameof(Variation_WhitespaceDescription_MultiResidueRange_SynthesizesRangeCode));
            Assert.That((string)AssertSingleVariantFeature(doc).Attribute("description"), Is.EqualTo("PTI4-6KAA"));
        }

        [Test]
        public void Variation_Deletion_SynthesizesFallbackSequenceVariant()
        {
            var prot = MakeBaseProtein("ACC_DEL");
            var carrier = GetConsensusCarrier(prot);
            carrier.SequenceVariations.Add(new SequenceVariation(3, 3, "E", "", "  ", variantCallFormatDataString: null));

            var doc = WriteAndLoad(prot, nameof(Variation_Deletion_SynthesizesFallbackSequenceVariant));
            Assert.That((string)AssertSingleVariantFeature(doc).Attribute("description"),
                Is.EqualTo("sequence variant"));
        }

        [Test]
        public void Variation_Insertion_SynthesizesFallbackSequenceVariant()
        {
            var prot = MakeBaseProtein("ACC_INS");
            var carrier = GetConsensusCarrier(prot);
            carrier.SequenceVariations.Add(new SequenceVariation(5, (string)null, "AA", "  ", variantCallFormatDataString: null));

            var doc = WriteAndLoad(prot, nameof(Variation_Insertion_SynthesizesFallbackSequenceVariant));
            Assert.That((string)AssertSingleVariantFeature(doc).Attribute("description"),
                Is.EqualTo("sequence variant"));
        }
            
        [Test]
        public void MultipleVariants_AreOrdered_ByBeginThenVariantSequence()
        {
            var prot = MakeBaseProtein("ACC_ORDER");
            var carrier = GetConsensusCarrier(prot);
            carrier.SequenceVariations.Add(new SequenceVariation(7, 7, "S", "R", "Z", variantCallFormatDataString: null));
            carrier.SequenceVariations.Add(new SequenceVariation(3, 3, "E", "K", "DescK", variantCallFormatDataString: null));
            carrier.SequenceVariations.Add(new SequenceVariation(3, 3, "E", "A", "DescA", variantCallFormatDataString: null));

            var doc = WriteAndLoad(prot, nameof(MultipleVariants_AreOrdered_ByBeginThenVariantSequence));
            var ordered = VariantFeatures(doc)
                .Select(f =>
                {
                    var loc = FirstChild(f, "location");
                    var posNode = loc.Elements().First(e => e.Name.LocalName == "position" || e.Name.LocalName == "begin");
                    int pos = int.Parse(posNode.Attribute("position").Value, CultureInfo.InvariantCulture);
                    string variation = FirstChild(f, "variation")?.Value ?? "";
                    return (pos, variation);
                })
                .ToList();

            Assert.That(ordered.Count, Is.EqualTo(3));
            Assert.That(ordered[0].pos, Is.EqualTo(3));
            Assert.That(ordered[1].pos, Is.EqualTo(3));
            Assert.That(ordered[2].pos, Is.EqualTo(7));
            Assert.That(ordered[0].variation, Is.EqualTo("A"));
            Assert.That(ordered[1].variation, Is.EqualTo("K"));
            Assert.That(ordered[2].variation, Is.EqualTo("R"));
        }

        [Test]
        public void VariantSpecificModifications_WrittenAsSubfeatures()
        {
            var prot = MakeBaseProtein("ACC_VAR_MOD");
            var carrier = GetConsensusCarrier(prot);

            var varMods = new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification>{ CreateModWithId("VarModX") } }
            };

            carrier.SequenceVariations.Add(new SequenceVariation(1, 1, "M", "K", "  ",
                variantCallFormatDataString: null,
                oneBasedModifications: varMods));

            var doc = WriteAndLoad(prot, nameof(VariantSpecificModifications_WrittenAsSubfeatures));
            var feature = AssertSingleVariantFeature(doc);
            var subfeatures = feature
                .Descendants()
                .Where(sf => sf.Name.LocalName == "subfeature"
                             && string.Equals((string)sf.Attribute("type"), "modified residue", StringComparison.Ordinal))
                .ToList();

            Assert.That(subfeatures.Count, Is.EqualTo(1), "Expected exactly one modified residue subfeature.");
            var desc = (string)subfeatures[0].Attribute("description");
            Assert.That(desc, Is.EqualTo("VarModX"), "Subfeature description should use IdWithMotif (VarModX).");
            Assert.That(subfeatures[0]
                .Descendants()
                .Any(sp => sp.Name.LocalName == "subposition"
                           && (string)sp.Attribute("subposition") == "1"), Is.True);
        }

        [Test]
        public void AdditionalExternallySuppliedMods_DoNotAffectDescriptionLogic()
        {
            var prot = MakeBaseProtein("ACC_EXTRA_MOD");
            var carrier = GetConsensusCarrier(prot);
            carrier.SequenceVariations.Add(new SequenceVariation(2, 2, "P", "A", "   ", variantCallFormatDataString: null));

            var externalMod = CreateModWithId("ExtraMod1");

            var extraMods = new Dictionary<string, HashSet<Tuple<int, Modification>>>
            {
                { carrier.Accession, new HashSet<Tuple<int, Modification>>
                    {
                        Tuple.Create(2, externalMod)
                    }
                }
            };

            var doc = WriteAndLoad(prot, nameof(AdditionalExternallySuppliedMods_DoNotAffectDescriptionLogic), extraMods);
            Assert.That((string)AssertSingleVariantFeature(doc).Attribute("description"),
                Is.EqualTo("P2A"), "External mods must not alter synthesized variant description.");
        }
    }
}