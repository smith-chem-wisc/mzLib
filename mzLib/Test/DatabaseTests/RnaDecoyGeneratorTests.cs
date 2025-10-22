using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using UsefulProteomicsDatabases;
using Transcriptomics;
using Omics.Modifications;
using Omics.BioPolymer;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class RnaDecoyGeneratorTests
    {
        // Build per-position modifications whose motif matches the nucleotide at that position
        private static Dictionary<int, List<Modification>> BuildModsForSequence(string sequence, params int[] positions)
        {
            var dict = new Dictionary<int, List<Modification>>();
            foreach (var pos in positions.Distinct())
            {
                if (pos < 1 || pos > sequence.Length)
                    throw new ArgumentOutOfRangeException(nameof(positions), $"Position {pos} out of range for length {sequence.Length}");
                char baseChar = sequence[pos - 1];
                if (!ModificationMotif.TryGetMotif(baseChar.ToString(), out var motif))
                {
                    ModificationMotif.TryGetMotif(char.ToUpperInvariant(baseChar).ToString(), out motif);
                }

                var mod = new Modification(
                    _originalId: $"Mod_{pos}_{baseChar}",
                    _modificationType: "TestType",
                    _target: motif,
                    _locationRestriction: "Anywhere.");

                dict[pos] = new List<Modification> { mod };
            }
            return dict;
        }

        private static SequenceVariation MakeVariant(string seq,
                                                     int begin,
                                                     int end,
                                                     string original,
                                                     string variant,
                                                     string description,
                                                     Dictionary<int, List<Modification>> variantSiteMods = null,
                                                     string vcf = null)
        {
            if (variantSiteMods != null)
            {
                var rebuilt = new Dictionary<int, List<Modification>>();
                foreach (var kvp in variantSiteMods)
                {
                    int pos = kvp.Key;
                    char baseChar = (pos >= begin && pos <= end && variant.Length > 0)
                        ? variant[Math.Min(variant.Length - 1, pos - begin)]
                        : (pos - 1 < seq.Length ? seq[pos - 1] : 'A');

                    if (!ModificationMotif.TryGetMotif(baseChar.ToString(), out var motif))
                    {
                        ModificationMotif.TryGetMotif("A", out motif);
                    }

                    rebuilt[pos] = kvp.Value.Select(v =>
                        new Modification(
                            _originalId: v.OriginalId,
                            _modificationType: v.ModificationType,
                            _target: motif,
                            _locationRestriction: "Anywhere.")).ToList();
                }
                variantSiteMods = rebuilt;
            }

            return new SequenceVariation(begin, end, original, variant, description, vcf, variantSiteMods);
        }

        private static RNA MakeSimpleRna(string accession, string sequence = "AUGCUA")
        {
            var mods = BuildModsForSequence(sequence, 2, 5);
            return new RNA(sequence, accession,
                oneBasedPossibleModifications: mods,
                fivePrimeTerminus: null, threePrimeTerminus: null,
                name: accession + "_NAME",
                organism: "TestOrg",
                databaseFilePath: "inMemory",
                isContaminant: false,
                isDecoy: false,
                geneNames: null,
                databaseAdditionalFields: null,
                truncationProducts: new List<TruncationProduct>(),
                sequenceVariations: new List<SequenceVariation>(),
                appliedSequenceVariations: new List<SequenceVariation>(),
                sampleNameForVariants: null,
                fullName: accession + "_FULL");
        }

        private static RNA MakeComplexRnaWithVariants(string accession)
        {
            string seq = "AUGCGAUCGU";
            var baseMods = BuildModsForSequence(seq, 1, 4, 10);
            string vcf = "1\t100\t.\tA\tG\t.\tPASS\tANN=G|.\tGT:AD:DP\t0/1:5,7:12";

            var varSiteMods = BuildModsForSequence(seq, 3, 4, 5);
            var baseVar = MakeVariant(seq, 3, 5, "GCG", "AAA", "BaseVar", varSiteMods, vcf);

            var appliedSiteMods = BuildModsForSequence(seq, 6, 7);
            var appliedVar = MakeVariant(seq, 6, 7, "AU", "GG", "AppVar", appliedSiteMods, vcf);

            var trunc = new TruncationProduct(2, 8, "Internal");

            return new RNA(
                sequence: seq,
                accession: accession,
                oneBasedPossibleModifications: baseMods,
                fivePrimeTerminus: null,
                threePrimeTerminus: null,
                name: accession + "_Name",
                organism: "TestOrg",
                databaseFilePath: "inMemory",
                isContaminant: false,
                isDecoy: false,
                geneNames: null,
                databaseAdditionalFields: null,
                truncationProducts: new List<TruncationProduct> { trunc },
                sequenceVariations: new List<SequenceVariation> { baseVar },
                appliedSequenceVariations: new List<SequenceVariation> { appliedVar },
                sampleNameForVariants: null,
                fullName: accession + "_Full");
        }

        private static Dictionary<int, int> IndexMapping(string seq) =>
            Enumerable.Range(1, seq.Length).ToDictionary(i => i, i => seq.Length - i + 1);

        private static void AssertBaseModsReversed(RNA original, RNA decoy)
        {
            var map = IndexMapping(original.BaseSequence);

            var expected = original.OneBasedPossibleLocalizedModifications
                .SelectMany(kvp => kvp.Value.Select(m => (newPos: map[kvp.Key], m.OriginalId)))
                .GroupBy(x => x.newPos)
                .ToDictionary(g => g.Key, g => g.Select(x => x.OriginalId).OrderBy(s => s).ToList());

            var actual = decoy.OneBasedPossibleLocalizedModifications
                .SelectMany(kvp => kvp.Value.Select(m => (pos: kvp.Key, m.OriginalId)))
                .GroupBy(x => x.pos)
                .ToDictionary(g => g.Key, g => g.Select(x => x.OriginalId).OrderBy(s => s).ToList());

            Assert.That(actual.Keys.OrderBy(i => i), Is.EquivalentTo(expected.Keys.OrderBy(i => i)),
                "Reversed modification site positions mismatch");
            foreach (var kv in expected)
            {
                Assert.That(actual[kv.Key], Is.EqualTo(kv.Value), $"Mismatch at reversed site {kv.Key}");
            }
        }

        [Test]
        public void GenerateDecoys_None_ReturnsEmpty_OriginalUnchanged()
        {
            var rna = MakeSimpleRna("ACC_NONE");
            var originalHash = rna.BaseSequence.GetHashCode();
            var decoys = RnaDecoyGenerator.GenerateDecoys(new List<RNA> { rna }, DecoyType.None, 1, "D");
            Assert.That(decoys, Is.Empty);
            Assert.That(rna.BaseSequence.GetHashCode(), Is.EqualTo(originalHash));
        }

        [Test]
        public void GenerateDecoys_Reverse_Simple_ModificationsMoveWithBases()
        {
            var rna = MakeSimpleRna("ACC_SIMPLE", "AUGCUA");
            var decoys = RnaDecoyGenerator.GenerateDecoys(new List<RNA> { rna }, DecoyType.Reverse, 1, "REV");
            Assert.That(decoys.Count, Is.EqualTo(1));
            var rev = decoys[0];

            Assert.That(rev.BaseSequence, Is.EqualTo(new string(rna.BaseSequence.Reverse().ToArray())));
            AssertBaseModsReversed(rna, rev);
        }

        [Test]
        public void GenerateDecoys_Reverse_Complex_WithVariantsAndTruncations()
        {
            var rna = MakeComplexRnaWithVariants("ACC_COMPLEX");
            int L = rna.BaseSequence.Length;

            var decoys = RnaDecoyGenerator.GenerateDecoys(new List<RNA> { rna }, DecoyType.Reverse, 1, "REV");
            Assert.That(decoys.Count, Is.EqualTo(1));
            var rev = decoys[0] as RNA;
            Assert.That(rev, Is.Not.Null);

            Assert.That(rev.BaseSequence, Is.EqualTo(new string(rna.BaseSequence.Reverse().ToArray())));
            AssertBaseModsReversed(rna, rev);

            var baseVarOrig = rna.SequenceVariations.Single();
            var baseVarRev = rev.SequenceVariations.Single(v => v.Description == baseVarOrig.Description);
            Assert.That(baseVarRev.OneBasedBeginPosition, Is.EqualTo(L - baseVarOrig.OneBasedEndPosition + 1));
            Assert.That(baseVarRev.OneBasedEndPosition, Is.EqualTo(L - baseVarOrig.OneBasedBeginPosition + 1));

            var expectedVarModSites = baseVarOrig.OneBasedModifications.Keys
                .Select(k => L - k + 1)
                .OrderBy(i => i)
                .ToArray();
            var actualVarModSites = baseVarRev.OneBasedModifications.Keys.OrderBy(i => i).ToArray();
            Assert.That(actualVarModSites, Is.EquivalentTo(expectedVarModSites));

            var appliedOrig = rna.AppliedSequenceVariations.Single();
            var appliedRev = rev.AppliedSequenceVariations.Single(v => v.Description == appliedOrig.Description);
            Assert.That(appliedRev.OneBasedBeginPosition, Is.EqualTo(L - appliedOrig.OneBasedEndPosition + 1));
            Assert.That(appliedRev.OneBasedEndPosition, Is.EqualTo(L - appliedOrig.OneBasedBeginPosition + 1));

            var expectedAppliedModSites = appliedOrig.OneBasedModifications.Keys
                .Select(k => L - k + 1)
                .OrderBy(i => i)
                .ToArray();
            var actualAppliedModSites = appliedRev.OneBasedModifications.Keys.OrderBy(i => i).ToArray();
            Assert.That(actualAppliedModSites, Is.EquivalentTo(expectedAppliedModSites));

            var truncOrig = rna.TruncationProducts.Single();
            var truncRev = rev.TruncationProducts.Single();
            Assert.That(truncRev.OneBasedBeginPosition, Is.EqualTo(L - truncOrig.OneBasedEndPosition!.Value + 1));
            Assert.That(truncRev.OneBasedEndPosition, Is.EqualTo(L - truncOrig.OneBasedBeginPosition!.Value + 1));
        }

        [Test]
        public void GenerateDecoys_Reverse_InsertionVariant_PointMappingPreserved()
        {
            string seq = "AUGCUGCA";
            var mods = BuildModsForSequence(seq, 1, 8);

            var insVarMods = BuildModsForSequence(seq, 5);
            var insertionVar = new SequenceVariation(
                oneBasedPosition: 5,
                originalSequence: null,
                variantSequence: "GG",
                description: "InsGG",
                variantCallFormatDataString: "1\t50\t.\t.\tGG\t.\tPASS\tANN=GG|.\tGT:AD:DP\t0/1:4,6:10",
                oneBasedModifications: insVarMods);

            var rna = new RNA(seq, "ACC_INS",
                oneBasedPossibleModifications: mods,
                fivePrimeTerminus: null,
                threePrimeTerminus: null,
                name: "ACC_INS_Name",
                organism: "TestOrg",
                databaseFilePath: "mem",
                isContaminant: false,
                isDecoy: false,
                geneNames: null,
                databaseAdditionalFields: null,
                truncationProducts: new List<TruncationProduct>(),
                sequenceVariations: new List<SequenceVariation> { insertionVar },
                appliedSequenceVariations: new List<SequenceVariation>(),
                sampleNameForVariants: null,
                fullName: "ACC_INS_Full");

            var decoys = RnaDecoyGenerator.GenerateDecoys(new List<RNA> { rna }, DecoyType.Reverse, 1, "REV");
            Assert.That(decoys.Count, Is.EqualTo(1));
            var rev = decoys[0];
            int L = seq.Length;

            AssertBaseModsReversed(rna, rev);

            var insRev = rev.SequenceVariations.Single();
            int expectedPoint = L - 5 + 1;
            Assert.That(insRev.OneBasedBeginPosition, Is.EqualTo(expectedPoint));
            Assert.That(insRev.OneBasedEndPosition, Is.EqualTo(expectedPoint));
            Assert.That(insRev.OneBasedModifications.Keys.Single(), Is.EqualTo(expectedPoint));
        }

        [Test]
        public void GenerateDecoys_Reverse_MultipleTruncations_CorrectlyMapped()
        {
            string seq = "AUGCGAUCGU";
            var rna = new RNA(seq, "ACC_TRUNC",
                oneBasedPossibleModifications: new Dictionary<int, List<Modification>>(),
                fivePrimeTerminus: null,
                threePrimeTerminus: null,
                name: "ACC_TRUNC_Name",
                organism: "TestOrg",
                databaseFilePath: "mem",
                isContaminant: false,
                isDecoy: false,
                geneNames: null,
                databaseAdditionalFields: null,
                truncationProducts: new List<TruncationProduct>
                {
                    new TruncationProduct(1,5,"FragA"),
                    new TruncationProduct(3,8,"FragB"),
                    new TruncationProduct(9,10,"FragC")
                },
                sequenceVariations: new List<SequenceVariation>(),
                appliedSequenceVariations: new List<SequenceVariation>(),
                sampleNameForVariants: null,
                fullName: "ACC_TRUNC_Full");

            var decoys = RnaDecoyGenerator.GenerateDecoys(new List<RNA> { rna }, DecoyType.Reverse, 1, "REV");
            Assert.That(decoys.Count, Is.EqualTo(1));
            var rev = decoys[0];
            int L = seq.Length;

            (int b, int e, string type) Map(int begin, int end, string t)
                => (L - end + 1, L - begin + 1, t);

            var expected = rna.TruncationProducts
                .Select(t => Map(t.OneBasedBeginPosition!.Value, t.OneBasedEndPosition!.Value, t.Type))
                .Select(t => (begin: Math.Min(t.b, t.e), end: Math.Max(t.b, t.e), t.type))
                .ToList();

            var actual = rev.TruncationProducts
                .Select(t => (t.OneBasedBeginPosition!.Value, t.OneBasedEndPosition!.Value, t.Type))
                .ToList();

            Assert.That(actual.Count, Is.EqualTo(expected.Count));
            foreach (var exp in expected)
            {
                Assert.That(actual.Any(a => a.ValueTupleEquals(exp) && a.Item3.Contains("REV")),
                    Is.True, $"Missing reversed truncation {exp}");
            }
        }

        [Test]
        public void GenerateDecoys_Reverse_PalindromicSequence_ModsSymmetricallyRemapped()
        {
            string seq = "AUGUA";
            var mods = BuildModsForSequence(seq, 1, 3, 5);

            var rna = new RNA(seq, "ACC_PAL",
                oneBasedPossibleModifications: mods,
                fivePrimeTerminus: null,
                threePrimeTerminus: null,
                name: "ACC_PAL_Name",
                organism: "TestOrg",
                databaseFilePath: "mem",
                isContaminant: false,
                isDecoy: false,
                geneNames: null,
                databaseAdditionalFields: null,
                truncationProducts: new List<TruncationProduct>(),
                sequenceVariations: new List<SequenceVariation>(),
                appliedSequenceVariations: new List<SequenceVariation>(),
                sampleNameForVariants: null,
                fullName: "ACC_PAL_Full");

            var decoys = RnaDecoyGenerator.GenerateDecoys(new List<RNA> { rna }, DecoyType.Reverse, 1, "REV");
            Assert.That(decoys.Count, Is.EqualTo(1));
            var rev = decoys[0];

            Assert.That(rev.BaseSequence, Is.EqualTo(seq));
            AssertBaseModsReversed(rna, rev);
        }
    }

    internal static class TupleExtensions
    {
        public static bool ValueTupleEquals(this (int, int, string) a, (int begin, int end, string type) b) =>
            a.Item1 == b.begin && a.Item2 == b.end && a.Item3.Contains(b.type);
    }
}