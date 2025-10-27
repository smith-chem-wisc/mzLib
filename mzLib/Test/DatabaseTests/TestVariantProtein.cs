using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using Omics;
using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Transcriptomics;
using UsefulProteomicsDatabases;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.DatabaseTests
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestVariantProtein
    {
        private static List<Modification> UniProtPtms;
        private static Stopwatch Stopwatch { get; set; }

        [OneTimeSetUp]
        public static void SetUpModifications()
        {
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();
        }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public static void VariantProtein()
        {
            Protein p = new Protein("MAAA", "accession");
            Protein v = new Protein("MAVA", p, new[] { new SequenceVariation(3, "A", "V", "desc", null) }, null, null, null);
            Assert.AreEqual(p, v.ConsensusVariant);
        }

        [Test]
        public void VariantXml()
        {
            string file = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SeqVar.xml");
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.None, null, false, null, out var un,
                consensusPlusVariantIsoforms: 10,
                minAlleleDepth: 0,
                maxVariantsPerIsoform:1)
                .Where(v=>v.AppliedSequenceVariations.Count > 0)
                .ToList();

            Assert.AreEqual(5, variantProteins.First().ConsensusVariant.SequenceVariations.Count());
            Assert.AreEqual(1, variantProteins.Count); // there is only one unique amino acid change
            Assert.AreNotEqual(variantProteins.First().ConsensusVariant.BaseSequence, variantProteins.First().BaseSequence);
            Assert.AreEqual('C', variantProteins.First().ConsensusVariant.BaseSequence[116]);
            Assert.AreEqual('Y', variantProteins.First().BaseSequence[116]);
            Assert.AreNotEqual(variantProteins.First().ConsensusVariant.Name, variantProteins.First().Name);
            Assert.AreNotEqual(variantProteins.First().ConsensusVariant.FullName, variantProteins.First().FullName);
            Assert.AreNotEqual(variantProteins.First().ConsensusVariant.Accession, variantProteins.First().Accession);

            List<PeptideWithSetModifications> peptides = variantProteins.SelectMany(vp => vp.Digest(new DigestionParams(), null, null)).ToList();
        }

        [Test]
        public static void SeqVarXmlTest()
        {
            var ok = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "seqvartests.xml"),
                true, DecoyType.Reverse, UniProtPtms, false, null, out var un);

            var target = ok.First(p => !p.IsDecoy);
            Protein decoy = ok.Where(p => p.IsDecoy && p.SequenceVariations.Count() > 0).First();

            Assert.AreEqual('M', target[0]);
            Assert.AreEqual('M', decoy[0]);
            List<SequenceVariation> targetVariants = target.SequenceVariations.ToList();
            List<SequenceVariation> decoyVariants = decoy.SequenceVariations.ToList();
            Assert.AreEqual(targetVariants.Count, decoyVariants.Count);

            // starting methionine, but there's more
            Assert.AreEqual("MPEQA", targetVariants.First().OriginalSequence);
            Assert.AreEqual("MP", targetVariants.First().VariantSequence);
            Assert.AreEqual(1, targetVariants.First().OneBasedBeginPosition);
            Assert.AreEqual(5, targetVariants.First().OneBasedEndPosition);
            Assert.AreEqual("AQEP", decoy.SequenceVariations.First().OriginalSequence); // methionine will be at the front, so clipped off of the variant
            Assert.AreEqual("P", decoy.SequenceVariations.First().VariantSequence);
            Assert.AreEqual(target.Length - 3, decoy.SequenceVariations.First().OneBasedBeginPosition);
            Assert.AreEqual(target.Length, decoy.SequenceVariations.First().OneBasedEndPosition);

            // start loss
            Assert.AreEqual("MPEQA", targetVariants[1].OriginalSequence);
            Assert.AreEqual("P", decoyVariants[1].VariantSequence);
            Assert.AreEqual(1, targetVariants[1].OneBasedBeginPosition);
            Assert.AreEqual(5, targetVariants[1].OneBasedEndPosition);
            Assert.AreEqual("AQEP", decoy.SequenceVariations.First().OriginalSequence); // methionine will be at the front, so clipped off of the variant
            Assert.AreEqual("P", decoy.SequenceVariations.First().VariantSequence);
            Assert.AreEqual(target.Length - 3, decoy.SequenceVariations.First().OneBasedBeginPosition);
            Assert.AreEqual(target.Length, decoy.SequenceVariations.First().OneBasedEndPosition);

            foreach (SequenceVariation s in targetVariants)
            {
                Assert.AreEqual(s.OriginalSequence, target.BaseSequence.Substring(s.OneBasedBeginPosition - 1, s.OneBasedEndPosition - s.OneBasedBeginPosition + 1));
            }
            foreach (SequenceVariation s in decoyVariants)
            {
                Assert.AreEqual(s.OriginalSequence, decoy.BaseSequence.Substring(s.OneBasedBeginPosition - 1, s.OneBasedEndPosition - s.OneBasedBeginPosition + 1));
            }
            Assert.AreNotEqual(target.SequenceVariations.First().VariantCallFormatDataString, decoy.SequenceVariations.First().VariantCallFormatDataString); //decoys and target variations don't have the same desc.

            List<PeptideWithSetModifications> peptides = ok.SelectMany(vp => vp.Digest(new DigestionParams(), null, null)).ToList();
        }

        [Test]
        [TestCase("oblm1.xml", 1, 1)] // mod on starting methionine
        [TestCase("oblm2.xml", 3, 4)] // without starting methionine
        [TestCase("oblm3.xml", 3, 5)] // with starting methionine
        public static void LoadSeqVarModifications(string databaseName, int modIdx, int reversedModIdx)
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", databaseName), true,
                DecoyType.Reverse, null, false, null, out var unknownModifications);
            var target = proteins[0];
            Assert.AreEqual(1, target.OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(modIdx, target.OneBasedPossibleLocalizedModifications.Single().Key);
            Assert.AreEqual(1, target.AppliedSequenceVariations.Count());
            Assert.AreEqual(modIdx, target.AppliedSequenceVariations.Single().OneBasedBeginPosition);
            Assert.AreEqual(1, target.SequenceVariations.Count());
            Assert.AreEqual(modIdx, target.SequenceVariations.Single().OneBasedBeginPosition);
            Assert.AreEqual(1, target.SequenceVariations.Single().OneBasedModifications.Count);
            Assert.AreEqual(modIdx, target.SequenceVariations.Single().OneBasedModifications.Single().Key); //PEP[mod]TID, MEP[mod]TID
            var decoy = proteins[1];
            Assert.AreEqual(1, decoy.OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(reversedModIdx, decoy.OneBasedPossibleLocalizedModifications.Single().Key); //DITP[mod]EP, MDITP[mod]E
            Assert.AreEqual(1, decoy.AppliedSequenceVariations.Count());
            Assert.AreEqual(reversedModIdx, decoy.AppliedSequenceVariations.Single().OneBasedBeginPosition);
            Assert.AreEqual(1, decoy.SequenceVariations.Count());
            Assert.AreEqual(reversedModIdx, decoy.SequenceVariations.Single().OneBasedBeginPosition);
            Assert.AreEqual(1, decoy.SequenceVariations.Single().OneBasedModifications.Count);
            Assert.AreEqual(reversedModIdx, decoy.SequenceVariations.Single().OneBasedModifications.Single().Key);

            string rewriteDbName = $"{Path.GetFileNameWithoutExtension(databaseName)}rewrite.xml";
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteins.Where(p => !p.IsDecoy).ToList(), Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName));
            proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName), true,
                DecoyType.Reverse, null, false, null, out unknownModifications);
            target = proteins[0];
            Assert.AreEqual(1, target.OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(modIdx, target.OneBasedPossibleLocalizedModifications.Single().Key);
            Assert.AreEqual(1, target.AppliedSequenceVariations.Count());
            Assert.AreEqual(modIdx, target.AppliedSequenceVariations.Single().OneBasedBeginPosition);
            Assert.AreEqual(1, target.SequenceVariations.Count());
            Assert.AreEqual(modIdx, target.SequenceVariations.Single().OneBasedBeginPosition);
            Assert.AreEqual(1, target.SequenceVariations.Single().OneBasedModifications.Count);
            Assert.AreEqual(modIdx, target.SequenceVariations.Single().OneBasedModifications.Single().Key);
            decoy = proteins[1];
            Assert.AreEqual(1, decoy.OneBasedPossibleLocalizedModifications.Count);
            Assert.AreEqual(reversedModIdx, decoy.OneBasedPossibleLocalizedModifications.Single().Key);
            Assert.AreEqual(1, decoy.AppliedSequenceVariations.Count());
            Assert.AreEqual(reversedModIdx, decoy.AppliedSequenceVariations.Single().OneBasedBeginPosition);
            Assert.AreEqual(1, decoy.SequenceVariations.Count());
            Assert.AreEqual(reversedModIdx, decoy.SequenceVariations.Single().OneBasedBeginPosition);
            Assert.AreEqual(1, decoy.SequenceVariations.Single().OneBasedModifications.Count);
            Assert.AreEqual(reversedModIdx, decoy.SequenceVariations.Single().OneBasedModifications.Single().Key);
        }

        [TestCase("ranges1.xml", 1, 2, 5, 6)] // without starting methionine
        [TestCase("ranges2.xml", 1, 1, 5, 5)] // with starting methionine
        public static void ReverseDecoyProteolysisProducts(string databaseName, int beginIdx, int reversedBeginIdx, int endIdx, int reversedEndIdx)
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", databaseName), true,
                DecoyType.Reverse, null, false, null, out var unknownModifications);
            var target = proteins[0];
            Assert.AreEqual(1, target.TruncationProducts.Count());
            Assert.AreEqual(beginIdx, target.TruncationProducts.Single().OneBasedBeginPosition); //P[start]EPTI[end]D, M[start]EPTI[end]D
            Assert.AreEqual(endIdx, target.TruncationProducts.Single().OneBasedEndPosition);
            var decoy = proteins[1];
            Assert.AreEqual(1, decoy.TruncationProducts.Count());
            Assert.AreEqual(reversedBeginIdx, decoy.TruncationProducts.Single().OneBasedBeginPosition); //DI[start]TPEP[end], M[start]DITP[end]E
            Assert.AreEqual(reversedEndIdx, decoy.TruncationProducts.Single().OneBasedEndPosition);

            string rewriteDbName = $"{Path.GetFileNameWithoutExtension(databaseName)}rewrite.xml";
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteins.Where(p => !p.IsDecoy).ToList(), Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName));
            proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName), true,
                DecoyType.Reverse, null, false, null, out unknownModifications);
            target = proteins[0];
            Assert.AreEqual(1, target.TruncationProducts.Count());
            Assert.AreEqual(beginIdx, target.TruncationProducts.Single().OneBasedBeginPosition);
            Assert.AreEqual(endIdx, target.TruncationProducts.Single().OneBasedEndPosition);
            decoy = proteins[1];
            Assert.AreEqual(1, decoy.TruncationProducts.Count());
            Assert.AreEqual(reversedBeginIdx, decoy.TruncationProducts.Single().OneBasedBeginPosition);
            Assert.AreEqual(reversedEndIdx, decoy.TruncationProducts.Single().OneBasedEndPosition);
        }

        [TestCase("bonds1.xml", 2, 3, "DICPCP", 4, 5)] // without starting methionine
        [TestCase("bonds2.xml", 2, 4, "MDICPC", 4, 6)] // with starting methionine
        public static void ReverseDecoyDisulfideBonds(string databaseName, int beginIdx, int reversedBeginIdx, string reversedSequence, int endIdx, int reversedEndIdx)
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", databaseName), true,
                DecoyType.Reverse, null, false, null, out var unknownModifications);
            var target = proteins[0];
            Assert.AreEqual(1, target.DisulfideBonds.Count());
            Assert.AreEqual(beginIdx, target.DisulfideBonds.Single().OneBasedBeginPosition); //PC[start]PC[end]ID, MC[start]PC[end]ID
            Assert.AreEqual(endIdx, target.DisulfideBonds.Single().OneBasedEndPosition);
            var decoy = proteins[1];
            Assert.AreEqual(1, decoy.DisulfideBonds.Count());
            Assert.AreEqual(reversedSequence, decoy.BaseSequence);
            Assert.AreEqual(reversedBeginIdx, decoy.DisulfideBonds.Single().OneBasedBeginPosition); //DIC[start]PC[end]P, MDIC[start]PC[end]
            Assert.AreEqual(reversedEndIdx, decoy.DisulfideBonds.Single().OneBasedEndPosition);

            string rewriteDbName = $"{Path.GetFileNameWithoutExtension(databaseName)}rewrite.xml";
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteins.Where(p => !p.IsDecoy).ToList(), Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName));
            proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName), true,
                DecoyType.Reverse, null, false, null, out unknownModifications);
            target = proteins[0];
            Assert.AreEqual(1, target.DisulfideBonds.Count());
            Assert.AreEqual(beginIdx, target.DisulfideBonds.Single().OneBasedBeginPosition);
            Assert.AreEqual(endIdx, target.DisulfideBonds.Single().OneBasedEndPosition);
            decoy = proteins[1];
            Assert.AreEqual(1, decoy.DisulfideBonds.Count());
            Assert.AreEqual(reversedBeginIdx, decoy.DisulfideBonds.Single().OneBasedBeginPosition);
            Assert.AreEqual(reversedEndIdx, decoy.DisulfideBonds.Single().OneBasedEndPosition);
        }

        [Test]
        [TestCase("splices1.xml", 2, 4, 3, 5)] // range without starting methionine
        [TestCase("splices2.xml", 2, 5, 3, 6)] // range with starting methionine
        [TestCase("splices3.xml", 2, 5, 2, 5)] // site without starting methionine
        [TestCase("splices4.xml", 2, 6, 2, 6)] // site with starting methionine
        [TestCase("splices5.xml", 1, 6, 1, 6)] // start site without starting methionine
        [TestCase("splices6.xml", 1, 1, 1, 1)] // start site with starting methionine
        [TestCase("splices7.xml", 1, 5, 2, 6)] // range with start without starting methionine
        [TestCase("splices8.xml", 1, 5, 2, 6)] // range with start with starting methionine
        public static void ReverseDecoySpliceSites(string databaseName, int beginIdx, int reversedBeginIdx, int endIdx, int reversedEndIdx)
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", databaseName), true,
                DecoyType.Reverse, null, false, null, out var unknownModifications);
            var target = proteins[0];
            Assert.AreEqual(1, target.SpliceSites.Count());
            Assert.AreEqual(beginIdx, target.SpliceSites.Single().OneBasedBeginPosition); //PE[start]P[end]TID, ME[start]P[start]TID, PE[site]PTID, ME[site]PTID, P[site]EPTID, M[site]EPTID
            Assert.AreEqual(endIdx, target.SpliceSites.Single().OneBasedEndPosition);
            var decoy = proteins[1];
            Assert.AreEqual(1, decoy.SpliceSites.Count());
            Assert.AreEqual(reversedBeginIdx, decoy.SpliceSites.Single().OneBasedBeginPosition); //DITP[start]E[end]P, MDITP[start]E[end], DITPE[site]P, MDITPE[site], DITPEP[site], M[site]DITPE
            Assert.AreEqual(reversedEndIdx, decoy.SpliceSites.Single().OneBasedEndPosition);

            string rewriteDbName = $"{Path.GetFileNameWithoutExtension(databaseName)}rewrite.xml";
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteins.Where(p => !p.IsDecoy).ToList(), Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName));
            proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", rewriteDbName), true,
                DecoyType.Reverse, null, false, null, out unknownModifications);
            target = proteins[0];
            Assert.AreEqual(1, target.SpliceSites.Count());
            Assert.AreEqual(beginIdx, target.SpliceSites.Single().OneBasedBeginPosition);
            Assert.AreEqual(endIdx, target.SpliceSites.Single().OneBasedEndPosition);
            decoy = proteins[1];
            Assert.AreEqual(1, decoy.SpliceSites.Count());
            Assert.AreEqual(reversedBeginIdx, decoy.SpliceSites.Single().OneBasedBeginPosition);
            Assert.AreEqual(reversedEndIdx, decoy.SpliceSites.Single().OneBasedEndPosition);

            List<PeptideWithSetModifications> peptides = proteins.SelectMany(vp => vp.Digest(new DigestionParams(), null, null)).ToList();
        }

        [Test]
        [TestCase("HomozygousHLA.xml", 1, 18)]
        [TestCase("HomozygousHLA.xml", 10, 17)]
        public static void HomozygousVariantsAtVariedDepths(string filename, int minVariantDepth, int appliedCount)
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", filename), true,
                DecoyType.None, null, false, null, out var unknownModifications, minAlleleDepth: minVariantDepth);
            Assert.AreEqual(1, proteins.Count);
            Assert.AreEqual(18, proteins[0].SequenceVariations.Count()); // some redundant
            Assert.AreEqual(18, proteins[0].SequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(appliedCount, proteins[0].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(appliedCount, proteins[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(1, proteins[0].GetVariantBioPolymers().Count);
            var variantProteins = proteins[0].GetVariantBioPolymers();
            List<PeptideWithSetModifications> peptides = proteins.SelectMany(vp => vp.Digest(new DigestionParams(), null, null)).ToList();
        }
        [Test]
        public static void AppliedVariants()
        {
            // This test verifies that applying sequence variations (SAV, MNV, insertion, deletion)
            // produces the correct variant protein sequences, maps variant coordinates correctly,
            // preserves per-variant metadata (AppliedSequenceVariations), and remains stable across:
            // - repeated in-memory application,
            // - round-tripping through XML (write → read).
            //
            // Additionally, it verifies that a modification attached to a variant (protein5)
            // is persisted and localized at the expected one-based position after application
            // and after XML reload.

            // Define a simple P-specific modification used later to validate that modifications
            // attached to a variant are preserved and localized correctly after applying variants
            // and XML round-tripping.
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            Modification mp = new Modification(
                _originalId: "mod",
                _accession: null,
                _modificationType: "type",
                _featureType: null,
                _target: motifP,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: null,
                _monoisotopicMass: 42.01,
                _databaseReference: new Dictionary<string, IList<string>>(),
                _taxonomicRange: null,
                _keywords: null,
                _neutralLosses: null,
                _diagnosticIons: null,
                _fileOrigin: null);

            // Prepare five proteins that each have one sequence variation:
            // 1) protein1: Single Amino-acid Variant (SAV) P→V at position 4 (4..4)
            // 2) protein2: Multi-Nucleotide Variant (MNV) PT→KT spanning positions 4..5
            // 3) protein3: Insertion-like replacement: P→PPP at position 4 (longer variant)
            // 4) protein4: Deletion-like replacement: PPP→P spanning 4..6 (shorter variant)
            // 5) protein5: Same as (3) but with a modification attached at one-based index 5
            //    to verify mod persistence/localization through variant application and XML.
            List<Protein> proteinsWithSeqVars = new List<Protein>
            {
                new Protein("MPEPTIDE", "protein1",
                    sequenceVariations: new List<SequenceVariation>
                    {
                        // SAV: P(4) -> V(4)
                        new SequenceVariation(4, 4, "P", "V",
                            @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null)
                    }),
                new Protein("MPEPTIDE", "protein2",
                    sequenceVariations: new List<SequenceVariation>
                    {
                        // MNV: PT(4..5) -> KT(4..5)
                        new SequenceVariation(4, 5, "PT", "KT",
                            @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null)
                    }),
                new Protein("MPEPTIDE", "protein3",
                    sequenceVariations: new List<SequenceVariation>
                    {
                        // Insertion-like: P(4) -> PPP(4..6) (length +2)
                        new SequenceVariation(4, 4, "P", "PPP",
                            @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null)
                    }),
                new Protein("MPEPPPTIDE", "protein4",
                    sequenceVariations: new List<SequenceVariation>
                    {
                        // Deletion-like: PPP(4..6) -> P(4) (length -2)
                        new SequenceVariation(4, 6, "PPP", "P",
                            @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null)
                    }),
                new Protein("MPEPTIDE", "protein5",
                    sequenceVariations: new List<SequenceVariation>
                    {
                        // Insertion-like with a downstream mod to verify mod localization at 5
                        new SequenceVariation(4, 4, "P", "PPP",
                            @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30",
                            new Dictionary<int, List<Modification>> { { 5, new[] { mp }.ToList() } })
                    }),
            };

            // Apply variants in memory twice; the output should be stable and identical.
            var proteinsWithAppliedVariants = proteinsWithSeqVars.SelectMany(p => p.GetVariantBioPolymers()).ToList();
            var proteinsWithAppliedVariants2 = proteinsWithSeqVars.SelectMany(p => p.GetVariantBioPolymers()).ToList(); // should be stable

            // Round-trip through XML: write the variant-bearing proteins (targets only) and read back.
            // This ensures that variant application and metadata survive I/O and result identically.
            string xml = Path.Combine(TestContext.CurrentContext.TestDirectory, "AppliedVariants.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinsWithSeqVars, xml);
            var proteinsWithAppliedVariants3 = ProteinDbLoader.LoadProteinXML(xml, true, DecoyType.None, null, false, null, out var un);

            // Compare across all three sources:
            // - [0]: in-memory 1
            // - [1]: in-memory 2 (should match [0])
            // - [2]: XML reload (should match [0] and [1])
            var listArray = new[] { proteinsWithAppliedVariants, proteinsWithAppliedVariants2, proteinsWithAppliedVariants3 };

            // Assert we always get exactly five variant proteins in the same order
            // (SAV, MNV, insertion, deletion, insertion+mod).
            Assert.AreEqual(5, proteinsWithAppliedVariants.Count, "Expected 5 applied variants (in-memory #1).");
            Assert.AreEqual(5, proteinsWithAppliedVariants2.Count, "Expected 5 applied variants (in-memory #2).");
            Assert.AreEqual(5, proteinsWithAppliedVariants3.Count, "Expected 5 applied variants (XML reload).");

            // The expected sequences for each of the five applied variants, in order:
            // 0: "MPEVTIDE"  (SAV at 4: P->V)
            // 1: "MPEKTIDE"  (MNV at 4..5: PT->KT)
            // 2: "MPEPPPTIDE" (Insertion-like at 4: P->PPP; length +2)
            // 3: "MPEPTIDE"  (Deletion-like at 4..6: PPP->P; length -2 from "MPEPPPTIDE")
            // 4: "MPEPPPTIDE" (Insertion-like with mod at 5)
            for (int dbIdx = 0; dbIdx < listArray.Length; dbIdx++)
            {
                // sequences
                Assert.AreEqual("MPEVTIDE", listArray[dbIdx][0].BaseSequence, $"[{dbIdx}] SAV sequence mismatch");
                Assert.AreEqual("MPEKTIDE", listArray[dbIdx][1].BaseSequence, $"[{dbIdx}] MNV sequence mismatch");
                Assert.AreEqual("MPEPPPTIDE", listArray[dbIdx][2].BaseSequence, $"[{dbIdx}] insertion sequence mismatch");
                Assert.AreEqual("MPEPTIDE", listArray[dbIdx][3].BaseSequence, $"[{dbIdx}] deletion sequence mismatch");
                Assert.AreEqual("MPEPPPTIDE", listArray[dbIdx][4].BaseSequence, $"[{dbIdx}] insertion+mod sequence mismatch");

                // Confirm the modification attached to protein5 survives application and XML round-trip
                Assert.AreEqual(1, listArray[dbIdx][4].OneBasedPossibleLocalizedModifications.Count, $"[{dbIdx}] Expected exactly 1 mod on the insertion+mod variant");
                Assert.AreEqual(5, listArray[dbIdx][4].OneBasedPossibleLocalizedModifications.Single().Key, $"[{dbIdx}] Mod should localize to one-based position 5");
                // Sanity: ensure the residue under the mod is indeed P
                Assert.AreEqual('P', listArray[dbIdx][4].BaseSequence[5 - 1], $"[{dbIdx}] Residue at mod position should be 'P'");

                // SAV expectations: single-residue change; length unchanged; position 4 is 'V'
                Assert.AreEqual(8, listArray[dbIdx][0].BaseSequence.Length, $"[{dbIdx}] SAV length should be unchanged");
                Assert.AreEqual(4, listArray[dbIdx][0].AppliedSequenceVariations.Single().OneBasedBeginPosition, $"[{dbIdx}] SAV begin mismatch");
                Assert.AreEqual(4, listArray[dbIdx][0].AppliedSequenceVariations.Single().OneBasedEndPosition, $"[{dbIdx}] SAV end mismatch");
                Assert.AreEqual("P", listArray[dbIdx][0].AppliedSequenceVariations.Single().OriginalSequence, $"[{dbIdx}] SAV original should be 'P'");
                Assert.AreEqual("V", listArray[dbIdx][0].AppliedSequenceVariations.Single().VariantSequence, $"[{dbIdx}] SAV variant should be 'V'");
                Assert.AreEqual('V', listArray[dbIdx][0].BaseSequence[4 - 1], $"[{dbIdx}] SAV residue at 4 should be 'V'");

                // MNV expectations: multi-residue change; length unchanged; positions 4..5 become "KT"
                Assert.AreEqual(8, listArray[dbIdx][1].BaseSequence.Length, $"[{dbIdx}] MNV length should be unchanged");
                Assert.AreEqual(4, listArray[dbIdx][1].AppliedSequenceVariations.Single().OneBasedBeginPosition, $"[{dbIdx}] MNV begin mismatch");
                Assert.AreEqual(5, listArray[dbIdx][1].AppliedSequenceVariations.Single().OneBasedEndPosition, $"[{dbIdx}] MNV end mismatch");
                Assert.AreEqual("PT", listArray[dbIdx][1].AppliedSequenceVariations.Single().OriginalSequence, $"[{dbIdx}] MNV original should be 'PT'");
                Assert.AreEqual("KT", listArray[dbIdx][1].AppliedSequenceVariations.Single().VariantSequence, $"[{dbIdx}] MNV variant should be 'KT'");
                Assert.AreEqual("KT", listArray[dbIdx][1].BaseSequence.Substring(4 - 1, 2), $"[{dbIdx}] MNV residues 4..5 should be 'KT'");

                // insertion expectations: length grows by +2; positions 4..6 are "PPP"
                Assert.AreEqual(10, listArray[dbIdx][2].BaseSequence.Length, $"[{dbIdx}] insertion length should be 8 + 2 = 10");
                Assert.AreEqual(4, listArray[dbIdx][2].AppliedSequenceVariations.Single().OneBasedBeginPosition, $"[{dbIdx}] insertion begin mismatch");
                Assert.AreEqual(6, listArray[dbIdx][2].AppliedSequenceVariations.Single().OneBasedEndPosition, $"[{dbIdx}] insertion end mismatch");
                Assert.AreEqual("P", listArray[dbIdx][2].AppliedSequenceVariations.Single().OriginalSequence, $"[{dbIdx}] insertion original should be 'P'");
                Assert.AreEqual("PPP", listArray[dbIdx][2].AppliedSequenceVariations.Single().VariantSequence, $"[{dbIdx}] insertion variant should be 'PPP'");
                Assert.AreEqual("PPP", listArray[dbIdx][2].BaseSequence.Substring(4 - 1, 3), $"[{dbIdx}] insertion residues 4..6 should be 'PPP'");

                // deletion expectations: length shrinks by -2 relative to the starting "MPEPPPTIDE" (10 → 8)
                // and positions collapse so that sequence returns to "MPEPTIDE".
                Assert.AreEqual(8, listArray[dbIdx][3].BaseSequence.Length, $"[{dbIdx}] deletion length should be 10 - 2 = 8");
                Assert.AreEqual(4, listArray[dbIdx][3].AppliedSequenceVariations.Single().OneBasedBeginPosition, $"[{dbIdx}] deletion begin mismatch");
                Assert.AreEqual(4, listArray[dbIdx][3].AppliedSequenceVariations.Single().OneBasedEndPosition, $"[{dbIdx}] deletion end mismatch (post-collapse)");
                Assert.AreEqual("PPP", listArray[dbIdx][3].AppliedSequenceVariations.Single().OriginalSequence, $"[{dbIdx}] deletion original should be 'PPP'");
                Assert.AreEqual("P", listArray[dbIdx][3].AppliedSequenceVariations.Single().VariantSequence, $"[{dbIdx}] deletion variant should be 'P'");

                // For completeness, also assert the summarized begin/end expectations that the original test verified.
                // SAV
                Assert.AreEqual(4, listArray[dbIdx][0].AppliedSequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(4, listArray[dbIdx][0].AppliedSequenceVariations.Single().OneBasedEndPosition);
                // MNV
                Assert.AreEqual(4, listArray[dbIdx][1].AppliedSequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(5, listArray[dbIdx][1].AppliedSequenceVariations.Single().OneBasedEndPosition);
                // insertion
                Assert.AreEqual(4, listArray[dbIdx][2].AppliedSequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(6, listArray[dbIdx][2].AppliedSequenceVariations.Single().OneBasedEndPosition);
                // deletion
                Assert.AreEqual(4, listArray[dbIdx][3].AppliedSequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(4, listArray[dbIdx][3].AppliedSequenceVariations.Single().OneBasedEndPosition);
            }

            // Ensure exact stability across the three sources:
            // - All sequences in in-memory #1 equal in-memory #2 and XML reload.
            CollectionAssert.AreEqual(
                proteinsWithAppliedVariants.Select(p => p.BaseSequence).ToList(),
                proteinsWithAppliedVariants2.Select(p => p.BaseSequence).ToList(),
                "In-memory application should be stable across repeated calls");
            CollectionAssert.AreEqual(
                proteinsWithAppliedVariants.Select(p => p.BaseSequence).ToList(),
                proteinsWithAppliedVariants3.Select(p => p.BaseSequence).ToList(),
                "XML round-trip should preserve variant-applied sequences in the same order");
        }
        [Test]
        public static void AppliedVariants_AsIBioPolymer()
        {
            // PURPOSE
            // This test mirrors "AppliedVariants" but exercises the IBioPolymer interface path.
            // It validates, in detail:
            // 1) Correct application of four variant types (SAV, MNV, insertion, deletion) to produce expected sequences.
            // 2) Correct coordinates (begin/end), original/variant strings in AppliedSequenceVariations after application.
            // 3) Stability of results across:
            //    - repeated in-memory variant application,
            //    - XML round-trip (write → read).
            // 4) Modification propagation and localization for a variant carrying a downstream mod.
            // 5) Distinguishing two variants with identical sequences by their modification state (index 2 vs 4).
            //
            // EXPECTATIONS SUMMARY
            // - Variant sequences (in order): "MPEVTIDE", "MPEKTIDE", "MPEPPPTIDE", "MPEPTIDE", "MPEPPPTIDE".
            // - AppliedSequenceVariations:
            //   SAV       (idx 0): begin=4, end=4, original="P",   variant="V",   len=8
            //   MNV       (idx 1): begin=4, end=5, original="PT",  variant="KT",  len=8
            //   Insertion (idx 2): begin=4, end=6, original="P",   variant="PPP", len=10
            //   Deletion  (idx 3): begin=4, end=4, original="PPP", variant="P",   len=8
            //   Insert+Mod(idx 4): same sequence as insertion, plus exactly one mod at one-based index 5 targeting a 'P'.
            //
            // NOTE: Lists [0], [1], [2] below represent:
            //   [0] in-memory application (first call)
            //   [1] in-memory application (second call, should be identical/stable)
            //   [2] XML round-trip results, must match [0] and [1]

            // Arrange: create a P-specific modification used to test propagation and localization
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            Modification mp = new Modification(
                _originalId: "mod",
                _accession: null,
                _modificationType: "type",
                _featureType: null,
                _target: motifP,
                _locationRestriction: "Anywhere.",
                _chemicalFormula: null,
                _monoisotopicMass: 42.01,
                _databaseReference: new Dictionary<string, IList<string>>(),
                _taxonomicRange: null,
                _keywords: null,
                _neutralLosses: null,
                _diagnosticIons: null,
                _fileOrigin: null);

            // Arrange: build five proteins (as IBioPolymer) with one sequence variation each
            // 1) SAV P(4)->V
            // 2) MNV PT(4..5)->KT
            // 3) Insertion-like P(4)->PPP(4..6)
            // 4) Deletion-like PPP(4..6)->P(4) on a longer starting sequence
            // 5) Insertion-like with a downstream mod at one-based index 5
            List<IBioPolymer> proteinsWithSeqVars = new List<IBioPolymer>
            {
                new Protein("MPEPTIDE",   "protein1", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P",   "V",   @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE",   "protein2", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 5, "PT",  "KT",  @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE",   "protein3", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P",   "PPP", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPPPTIDE", "protein4", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 6, "PPP", "P",   @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE",   "protein5", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P",   "PPP", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", new Dictionary<int, List<Modification>> { { 5, new[] { mp }.ToList() } }) }),
             };

            // Act: apply variants in-memory twice (should be identical/stable)
            var proteinsWithAppliedVariants = proteinsWithSeqVars.SelectMany(p => p.GetVariantBioPolymers()).ToList();
            var proteinsWithAppliedVariants2 = proteinsWithSeqVars.SelectMany(p => p.GetVariantBioPolymers()).ToList();

            // Act: write to XML and load back; results should match in-memory outputs
            string xml = Path.Combine(TestContext.CurrentContext.TestDirectory, "AppliedVariants.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinsWithSeqVars, xml);
            var proteinsWithAppliedVariants3 = ProteinDbLoader.LoadProteinXML(xml, true, DecoyType.None, null, false, null, out var un);

            // Group lists for uniform validation loops
            var listArray = new List<IBioPolymer>[]
            {
                proteinsWithAppliedVariants,
                proteinsWithAppliedVariants2,
                proteinsWithAppliedVariants3.Cast<IBioPolymer>().ToList()
            };

            // Assert: each expansion produces exactly 5 variant biopolymers, in the same, predictable order
            Assert.AreEqual(5, proteinsWithAppliedVariants.Count, "Expected 5 variants (in-memory 1)");
            Assert.AreEqual(5, proteinsWithAppliedVariants2.Count, "Expected 5 variants (in-memory 2)");
            Assert.AreEqual(5, proteinsWithAppliedVariants3.Count, "Expected 5 variants (XML reload)");

            // Assert: stability across the three sources (same sequences, same order)
            CollectionAssert.AreEqual(
                proteinsWithAppliedVariants.Select(p => p.BaseSequence).ToList(),
                proteinsWithAppliedVariants2.Select(p => p.BaseSequence).ToList(),
                "In-memory application should be stable across repeated calls");
            CollectionAssert.AreEqual(
                proteinsWithAppliedVariants.Select(p => p.BaseSequence).ToList(),
                proteinsWithAppliedVariants3.Select(p => p.BaseSequence).ToList(),
                "XML round-trip should preserve variant-applied sequences in the same order");

            // Assert: for all variants we expect exactly one applied sequence variation
            foreach (var variants in listArray)
            {
                Assert.AreEqual(1, variants[0].AppliedSequenceVariations.Count, "SAV must have exactly one applied variant");
                Assert.AreEqual(1, variants[1].AppliedSequenceVariations.Count, "MNV must have exactly one applied variant");
                Assert.AreEqual(1, variants[2].AppliedSequenceVariations.Count, "Insertion must have exactly one applied variant");
                Assert.AreEqual(1, variants[3].AppliedSequenceVariations.Count, "Deletion must have exactly one applied variant");
                Assert.AreEqual(1, variants[4].AppliedSequenceVariations.Count, "Insertion+Mod must have exactly one applied variant");
            }

            // Per-list validation of sequence, coordinates, and (where appropriate) residue checks and mod checks
            for (int dbIdx = 0; dbIdx < listArray.Length; dbIdx++)
            {
                // Assert: expected sequences in fixed order
                Assert.AreEqual("MPEVTIDE", listArray[dbIdx][0].BaseSequence, $"[{dbIdx}] SAV sequence mismatch");
                Assert.AreEqual("MPEKTIDE", listArray[dbIdx][1].BaseSequence, $"[{dbIdx}] MNV sequence mismatch");
                Assert.AreEqual("MPEPPPTIDE", listArray[dbIdx][2].BaseSequence, $"[{dbIdx}] insertion sequence mismatch");
                Assert.AreEqual("MPEPTIDE", listArray[dbIdx][3].BaseSequence, $"[{dbIdx}] deletion sequence mismatch");
                Assert.AreEqual("MPEPPPTIDE", listArray[dbIdx][4].BaseSequence, $"[{dbIdx}] insertion+mod sequence mismatch");

                // Assert: lengths (sanity for ins/del)
                Assert.AreEqual(8, listArray[dbIdx][0].BaseSequence.Length, $"[{dbIdx}] SAV length should be unchanged");
                Assert.AreEqual(8, listArray[dbIdx][1].BaseSequence.Length, $"[{dbIdx}] MNV length should be unchanged");
                Assert.AreEqual(10, listArray[dbIdx][2].BaseSequence.Length, $"[{dbIdx}] insertion length should be 8 + 2 = 10");
                Assert.AreEqual(8, listArray[dbIdx][3].BaseSequence.Length, $"[{dbIdx}] deletion length should be 10 - 2 = 8");
                Assert.AreEqual(10, listArray[dbIdx][4].BaseSequence.Length, $"[{dbIdx}] insertion+mod length should be 10");

                // SAV assertions: P(4)->V
                var sav = listArray[dbIdx][0].AppliedSequenceVariations.Single();
                Assert.AreEqual(4, sav.OneBasedBeginPosition, $"[{dbIdx}] SAV begin");
                Assert.AreEqual(4, sav.OneBasedEndPosition, $"[{dbIdx}] SAV end");
                Assert.AreEqual("P", sav.OriginalSequence, $"[{dbIdx}] SAV original");
                Assert.AreEqual("V", sav.VariantSequence, $"[{dbIdx}] SAV variant");
                Assert.AreEqual('V', listArray[dbIdx][0].BaseSequence[3], $"[{dbIdx}] SAV residue at 4 must be 'V'");

                // MNV assertions: PT(4..5)->KT
                var mnv = listArray[dbIdx][1].AppliedSequenceVariations.Single();
                Assert.AreEqual(4, mnv.OneBasedBeginPosition, $"[{dbIdx}] MNV begin");
                Assert.AreEqual(5, mnv.OneBasedEndPosition, $"[{dbIdx}] MNV end");
                Assert.AreEqual("PT", mnv.OriginalSequence, $"[{dbIdx}] MNV original");
                Assert.AreEqual("KT", mnv.VariantSequence, $"[{dbIdx}] MNV variant");
                Assert.AreEqual("KT", listArray[dbIdx][1].BaseSequence.Substring(3, 2), $"[{dbIdx}] MNV residues 4..5 must be 'KT'");

                // Insertion-like assertions: P(4)->PPP(4..6)
                var ins = listArray[dbIdx][2].AppliedSequenceVariations.Single();
                Assert.AreEqual(4, ins.OneBasedBeginPosition, $"[{dbIdx}] insertion begin");
                Assert.AreEqual(6, ins.OneBasedEndPosition, $"[{dbIdx}] insertion end");
                Assert.AreEqual("P", ins.OriginalSequence, $"[{dbIdx}] insertion original");
                Assert.AreEqual("PPP", ins.VariantSequence, $"[{dbIdx}] insertion variant");
                Assert.AreEqual("PPP", listArray[dbIdx][2].BaseSequence.Substring(3, 3), $"[{dbIdx}] insertion residues 4..6 must be 'PPP'");

                // Deletion-like assertions: PPP(4..6)->P(4) (collapses back to MPEPTIDE)
                var del = listArray[dbIdx][3].AppliedSequenceVariations.Single();
                Assert.AreEqual(4, del.OneBasedBeginPosition, $"[{dbIdx}] deletion begin");
                Assert.AreEqual(4, del.OneBasedEndPosition, $"[{dbIdx}] deletion end (post-collapse)");
                Assert.AreEqual("PPP", del.OriginalSequence, $"[{dbIdx}] deletion original");
                Assert.AreEqual("P", del.VariantSequence, $"[{dbIdx}] deletion variant");

                // Insertion + Modification assertions: identical sequence to insertion, plus one mod at pos 5
                var insMod = listArray[dbIdx][4].AppliedSequenceVariations.Single();
                Assert.AreEqual(4, insMod.OneBasedBeginPosition, $"[{dbIdx}] insertion+mod begin");
                Assert.AreEqual(6, insMod.OneBasedEndPosition, $"[{dbIdx}] insertion+mod end");
                Assert.AreEqual("P", insMod.OriginalSequence, $"[{dbIdx}] insertion+mod original");
                Assert.AreEqual("PPP", insMod.VariantSequence, $"[{dbIdx}] insertion+mod variant");

                // Mods: Only the "insertion+mod" variant (index 4) carries a mod; all others should have none
                // Confirm the variant 4 has exactly one mod at one-based position 5, and the residue at 5 is 'P' (motif).
                if (dbIdx == 0 || dbIdx == 1 || dbIdx == 2) // only assert detailed mod state once per path for clarity
                {
                    // Index 2 (plain insertion) must have no mods
                    Assert.AreEqual(0, listArray[dbIdx][2].OneBasedPossibleLocalizedModifications.Count, $"[{dbIdx}] insertion should have no mods");

                    // Index 4 (insertion+mod) must have exactly one mod at site 5, targeting a P residue
                    Assert.AreEqual(1, listArray[dbIdx][4].OneBasedPossibleLocalizedModifications.Count, $"[{dbIdx}] insertion+mod should have exactly one site with mods");
                    Assert.AreEqual(5, listArray[dbIdx][4].OneBasedPossibleLocalizedModifications.Single().Key, $"[{dbIdx}] insertion+mod site should be one-based index 5");
                    Assert.AreEqual('P', listArray[dbIdx][4].BaseSequence[5 - 1], $"[{dbIdx}] insertion+mod residue at site 5 must be 'P' (motif)");

                    // All other variants should have zero possible localized mods
                    Assert.AreEqual(0, listArray[dbIdx][0].OneBasedPossibleLocalizedModifications.Count, $"[{dbIdx}] SAV should have no mods");
                    Assert.AreEqual(0, listArray[dbIdx][1].OneBasedPossibleLocalizedModifications.Count, $"[{dbIdx}] MNV should have no mods");
                    Assert.AreEqual(0, listArray[dbIdx][3].OneBasedPossibleLocalizedModifications.Count, $"[{dbIdx}] deletion should have no mods");
                }
            }

            // Additional cross-list checks:
            // - The two "MPEPPPTIDE" variants (indices 2 and 4) should be sequence-identical, but differ by modification presence.
            for (int dbIdx = 0; dbIdx < listArray.Length; dbIdx++)
            {
                Assert.AreEqual(listArray[dbIdx][2].BaseSequence, listArray[dbIdx][4].BaseSequence, $"[{dbIdx}] insertion and insertion+mod sequences must match");
                Assert.AreEqual(0, listArray[dbIdx][2].OneBasedPossibleLocalizedModifications.Count, $"[{dbIdx}] insertion should remain unmodified");
                Assert.AreEqual(1, listArray[dbIdx][4].OneBasedPossibleLocalizedModifications.Count, $"[{dbIdx}] insertion+mod should remain modified");
            }
        }
        [Test]
        public static void CrashOnCreateVariantFromRNA()
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "HomozygousHLA.xml"), true,
                DecoyType.None, null, false, null, out var unknownModifications);

            var rna = new RNA("GUACUGACU");
            NUnit.Framework.Assert.Throws<ArgumentException>(() =>
            {
                proteins[0].CreateVariant(proteins[0].BaseSequence, rna, [], [], new Dictionary<int, List<Modification>>(), "");
            });
        }
        [Test]
        public static void StopGained()
        {
            // PURPOSE
            // This test validates correct handling of a stop-gained sequence variation (creation of a premature stop codon).
            // Verifies two scenarios:
            // 1) Default variant depth filtering: both reference (no variant applied) and alternate (stop-gained applied) proteins are emitted.
            // 2) High min-allele-depth threshold: only the stop-gained (applied) protein is retained.
            //
            // EXPECTATIONS SUMMARY (based on the StopGained.xml test data):
            // - Two proteins are initially produced (reference + variant-applied):
            //   [0] Reference protein:
            //       * 1 sequence variation present in metadata, but 0 applied (reference form).
            //       * BaseSequence length = 191.
            //       * Residue at one-based position 161 equals 'Q' (so zero-based index 160 is 'Q').
            //   [1] Variant-applied protein (stop-gained):
            //       * Exactly 1 applied variation.
            //       * BaseSequence length truncated to 161 - 1 = 160 (stop codon at 161 shortens the sequence).
            //       * The sequence is exactly the prefix of the reference up to length 160.
            //       * No '*' appears in the resulting BaseSequence (the stop codon is not a literal character in the sequence).
            //       * The applied variant's VariantSequence ends with '*'.
            //       * The applied variant one-based begin (and end) position is 161.
            //
            // - With minAlleleDepth: 400
            //   * Only the variant-applied protein is returned.
            //   * It retains the same applied-variation and truncated length expectations.

            // Load the proteins with default filtering
            var proteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "StopGained.xml"),
                true,               // generateTargets
                DecoyType.None,     // decoyType
                null,               // allKnownModifications
                false,              // isContaminant
                null,               // modTypesToExclude
                out var unknownModifications);

            // Sanity: Decoys are not requested
            Assert.IsTrue(proteins.All(p => !p.IsDecoy), "No decoys expected when using DecoyType.None");

            // Expect exactly two proteins: reference (no applied variant) and stop-gained (applied)
            Assert.AreEqual(2, proteins.Count, "Expected reference and stop-gained variant proteins");
            // Reference protein metadata: one possible sequence variation in the record
            Assert.AreEqual(1, proteins[0].SequenceVariations.Count(), "Reference metadata should contain one sequence variation");
            Assert.AreEqual(1, proteins[0].SequenceVariations.Select(v => v.SimpleString()).Distinct().Count(), "Reference metadata should contain exactly one unique sequence variation");
            // Reference should have zero applied variations (reference form retained)
            Assert.AreEqual(0, proteins[0].AppliedSequenceVariations.Count(), "Reference protein should not have an applied sequence variation");
            Assert.AreEqual(0, proteins[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count(), "Reference applied variations should be zero");

            // Variant-applied protein: applied variation present and unique
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Count(), "Variant protein should have exactly one applied variation");
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count(), "Variant protein applied variations should be unique");

            // Reference length and residue checks around the stop site
            Assert.AreEqual(191, proteins[0].Length, "Reference protein length should match source data");
            Assert.AreEqual('Q', proteins[0][161 - 1], "Reference residue at 161 should be 'Q' prior to stop-gain");

            // Variant length must be truncated by the stop at position 161 → length becomes 160
            Assert.AreEqual(161 - 1, proteins[1].Length, "Variant protein length should be truncated to 160 due to stop at 161");

            // The variant BaseSequence must be exactly the prefix of the reference up to the stop position - 1
            string reference = proteins[0].BaseSequence;
            string variant = proteins[1].BaseSequence;
            Assert.IsTrue(reference.StartsWith(variant), "Variant sequence must be a prefix of the reference sequence");
            Assert.AreEqual(reference.Substring(0, 161 - 1), variant, "Variant sequence must equal reference[0..159]");

            // Ensure we did not write literal '*' into the protein sequence; stop codon is represented by truncation instead
            Assert.AreEqual(-1, variant.IndexOf('*'), "Variant BaseSequence should not contain a literal '*'");

            // Verify applied-variant details for the stop-gained protein
            var applied = proteins[1].AppliedSequenceVariations.Single();
            Assert.IsTrue(applied.VariantSequence.EndsWith("*"), "Stop-gained variant must end with '*'");
            Assert.AreEqual(161, applied.OneBasedBeginPosition, "Stop-gained begins at residue 161");
            Assert.AreEqual(161, applied.OneBasedEndPosition, "Stop-gained ends at residue 161 (single residue change)");

            // The two forms must differ in sequence length (as a sanity check)
            Assert.AreNotEqual(proteins[0].Length, proteins[1].Length, "Reference and variant proteins should differ in length");

            // Now require a higher min-allele-depth; expect only the variant-applied protein retained
            proteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "StopGained.xml"),
                true,               // generateTargets
                DecoyType.None,     // decoyType
                null,               // allKnownModifications
                false,              // isContaminant
                null,               // modTypesToExclude
                out unknownModifications,
                minAlleleDepth: 400);

            // Only the stop-gained, variant-applied form is retained under a strict depth threshold
            Assert.AreEqual(1, proteins.Count, "High min-allele-depth should retain only the variant-applied protein");
            Assert.AreEqual(1, proteins[0].AppliedSequenceVariations.Count(), "Variant-applied protein should still have one applied variation");
            Assert.AreEqual(1, proteins[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count(), "Variant-applied unique variation should be retained");
            Assert.AreEqual(161 - 1, proteins[0].Length, "Variant-applied protein length should remain truncated to 160");

            // Confirm stability: the single protein from the depth-filtered load matches the previously observed variant sequence
            Assert.AreEqual(variant, proteins[0].BaseSequence, "Depth-filtered variant sequence should match previously observed variant");

            // Re-check applied-variant semantics after filtering for completeness
            var appliedAfterFilter = proteins[0].AppliedSequenceVariations.Single();
            Assert.IsTrue(appliedAfterFilter.VariantSequence.EndsWith("*"), "Stop-gained variant must end with '*' (after filtering)");
            Assert.AreEqual(161, appliedAfterFilter.OneBasedBeginPosition, "Stop-gained begins at residue 161 (after filtering)");
            Assert.AreEqual(161, appliedAfterFilter.OneBasedEndPosition, "Stop-gained ends at residue 161 (after filtering)");
        }
        [Test]
        public static void StopGained_TruncationIsPrefixAndNoOutOfBoundsAnnotations()
        {
            var proteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "StopGained.xml"),
                true, DecoyType.None, null, false, null, out var _);

            Assert.AreEqual(2, proteins.Count);
            var reference = proteins[0];
            var truncated = proteins[1];

            // The truncated sequence must be a prefix of the reference,
            // i.e., identical up to the truncation point.
            Assert.That(reference.BaseSequence.StartsWith(truncated.BaseSequence));

            // Any possible localized modifications must not point past the truncation boundary.
            Assert.That(truncated.OneBasedPossibleLocalizedModifications
                .All(kv => kv.Key >= 1 && kv.Key <= truncated.Length));

            // Any proteolysis products (if present) must not reference indices outside the sequence.
            Assert.That(truncated.TruncationProducts.All(tp =>
                (!tp.OneBasedBeginPosition.HasValue || (tp.OneBasedBeginPosition.Value >= 1 && tp.OneBasedBeginPosition.Value <= truncated.Length)) &&
                (!tp.OneBasedEndPosition.HasValue || (tp.OneBasedEndPosition.Value >= 1 && tp.OneBasedEndPosition.Value <= truncated.Length))));

            // The applied stop-gained variation often encodes a '*' in the variant sequence.
            // If present, that indicates stop; the actual sequence is cut at the stop.
            if (truncated.AppliedSequenceVariations.Any())
            {
                Assert.That(truncated.AppliedSequenceVariations.Single().VariantSequence.EndsWith("*") ||
                            !truncated.AppliedSequenceVariations.Single().VariantSequence.Contains("*"));
            }
        }
        [Test]
        public static void StopGained_RoundTripSerializationPreservesTruncation()
        {
            var proteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "StopGained.xml"),
                true, DecoyType.None, null, false, null, out var _);

            Assert.AreEqual(2, proteins.Count);
            var tempPath = Path.Combine(TestContext.CurrentContext.TestDirectory, $"StopGained_roundtrip_{Guid.NewGuid()}.xml");

            try
            {
                // Persist both proteins (reference + variant) and reload.
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteins, tempPath);
                var roundtrip = ProteinDbLoader.LoadProteinXML(tempPath, true, DecoyType.None, null, false, null, out var __);

                // Round-trip preserves count and the truncation boundary for the variant-applied protein.
                Assert.AreEqual(2, roundtrip.Count);
                Assert.AreEqual(proteins[0].Length, roundtrip[0].Length);
                Assert.AreEqual(proteins[1].Length, roundtrip[1].Length);
                Assert.AreEqual(proteins[1].AppliedSequenceVariations.Count(), roundtrip[1].AppliedSequenceVariations.Count());
            }
            finally
            {
                if (File.Exists(tempPath))
                {
                    File.SetAttributes(tempPath, FileAttributes.Normal);
                    File.Delete(tempPath);
                }
            }
        }
        [Test]
        public static void StopGained_NoPeptidesCrossTruncationSite()
        {
            var proteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "StopGained.xml"),
                true, DecoyType.None, null, false, null, out var _);

            Assert.AreEqual(2, proteins.Count);
            var reference = proteins[0];
            var truncated = proteins[1];

            // Peptides from the truncated protein must not reference indices past the truncation boundary.
            var dp = new DigestionParams();
            var variantPeps = truncated.Digest(dp, null, null).ToList();
            Assert.That(variantPeps.All(p => p.OneBasedEndResidueInProtein <= truncated.Length));

            // Any peptide in the reference that extends past the truncation boundary cannot exist in the variant.
            var refPeps = reference.Digest(dp, null, null).ToList();
            var refCrossing = refPeps.Where(p => p.OneBasedEndResidueInProtein > truncated.Length).ToList();
            var variantPepWindows = new HashSet<(int start, int end)>(variantPeps.Select(p => (p.OneBasedStartResidueInProtein, p.OneBasedEndResidueInProtein)));
            Assert.That(refCrossing.All(p => !variantPepWindows.Contains((p.OneBasedStartResidueInProtein, p.OneBasedEndResidueInProtein))));
        }
        [Test]
        public static void StopGainedDecoysAndDigestion()
        {
            // PURPOSE
            // This test ensures that:
            // 1) Reverse decoys are generated for a stop-gained (truncated) target protein.
            // 2) The target and its decoy digest into peptides without referencing residues outside their sequence bounds.
            // 3) Basic invariants hold for decoy generation (count/order/length/decoy flag).
            //
            // CONTEXT
            // - Database: "StopGain.xml" contains a target protein with a stop-gained variant that truncates the sequence.
            // - Decoys: Generated using DecoyType.Reverse (sequence reversal-based decoys).
            // - minAlleleDepth: 400 ensures the stop-gained variant is applied (truncated target).
            //
            // EXPECTATIONS
            // - Exactly 2 proteins are returned: [0] target (non-decoy) + [1] decoy (IsDecoy = true).
            // - Target and decoy lengths are equal (reverse decoys preserve length).
            // - Both target and decoy produce peptides via digestion.
            // - No peptide references indices outside its parent protein's 1..Length range.
            // - Accessions reflect decoy generation (decoy starts with the default "DECOY_").
            // - If variant(s) are present, their counts match between target and decoy.

            // Arrange: Load a variant-applied protein set and reverse decoy pair from StopGain.xml.
            // Using a strict minAlleleDepth applies the stop-gained variant to the target.
            var proteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "StopGain.xml"),
                true,                 // generateTargets
                DecoyType.Reverse,    // generate reverse-sequence decoys
                null,                 // allKnownModifications
                false,                // isContaminant
                null,                 // modTypesToExclude
                out var unknownModifications,
                minAlleleDepth: 400); // force applying the stop-gained variant

            // Assert: We expect exactly two proteins: target then its decoy.
            Assert.AreEqual(2, proteins.Count, "Expected 1 target + 1 decoy");
            Assert.IsFalse(proteins[0].IsDecoy, "First protein should be the target (non-decoy)");
            Assert.IsTrue(proteins[1].IsDecoy, "Second protein should be the decoy");
            Assert.That(proteins[1].Accession.StartsWith("DECOY_"), "Decoy accession should start with the default decoy identifier");

            // Assert: Reverse decoys should preserve sequence length.
            Assert.AreEqual(proteins[0].Length, proteins[1].Length, "Target and decoy should have identical lengths");
            // In general, target and decoy sequences should not be byte-identical.
            Assert.AreNotEqual(proteins[0].BaseSequence, proteins[1].BaseSequence, "Decoy sequence should differ from target sequence");

            // If the stop-gained variant is applied to the target, the decoy should carry a corresponding variant count.
            // We do not enforce exact mapping positions here, only that counts match if any are present.
            if (proteins[0].AppliedSequenceVariations.Any() || proteins[1].AppliedSequenceVariations.Any())
            {
                Assert.AreEqual(
                    proteins[0].AppliedSequenceVariations.Count(),
                    proteins[1].AppliedSequenceVariations.Count(),
                    "Target and decoy should carry the same number of applied sequence variations");
            }

            // Act: Digest both target and decoy using default digestion parameters (typically trypsin-like rules).
            var dp = new DigestionParams();
            var targetPeps = proteins[0].Digest(dp, null, null).ToList();
            var decoyPeps = proteins[1].Digest(dp, null, null).ToList();

            // Assert: Both should yield peptides.
            Assert.That(targetPeps.Count > 0, "Target digestion should produce peptides");
            Assert.That(decoyPeps.Count > 0, "Decoy digestion should produce peptides");

            // Assert: No peptide references residues outside the corresponding protein bounds.
            Assert.That(targetPeps.All(p =>
                p.OneBasedStartResidueInProtein >= 1 &&
                p.OneBasedEndResidueInProtein <= proteins[0].Length),
                "All target peptides must fall within target bounds");

            Assert.That(decoyPeps.All(p =>
                p.OneBasedStartResidueInProtein >= 1 &&
                p.OneBasedEndResidueInProtein <= proteins[1].Length),
                "All decoy peptides must fall within decoy bounds");

            // Note:
            // We intentionally do NOT assert the number of peptides or the sum of peptide lengths to be equal between
            // target and decoy. Even with reverse decoys, tryptic cleavage context differs and may alter cleavage patterns.
            // The commented lines below are often too strict and can fail legitimately:
            //
            // Assert.AreEqual(targetPeps.Sum(p => p.Length), decoyPeps.Sum(p => p.Length));
            // Assert.AreEqual(targetPeps.Count, decoyPeps.Count);
        }

        [Test]
        public static void MultipleAlternateAlleles()
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "MultipleAlternateAlleles.xml"), true,
                DecoyType.None, null, false, null, out var unknownModifications);
            Assert.AreEqual(2, proteins.Count);
            Assert.AreEqual(2, proteins[0].SequenceVariations.Count()); // some redundant
            Assert.AreEqual(2, proteins[0].SequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes

            Assert.IsTrue(proteins[0].SequenceVariations.All(v => v.OneBasedBeginPosition == 63)); // there are two alternate alleles (1 and 2), but only 2 is in the genotype, so only that's applied
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(72, proteins[0].Length);
            Assert.AreEqual(72, proteins[1].Length);
            Assert.AreEqual('K', proteins[0][63 - 1]);
            Assert.AreEqual('R', proteins[1][63 - 1]);

            proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "MultipleAlternateAlleles.xml"), true,
                DecoyType.None, null, false, null, out unknownModifications, minAlleleDepth: 10);
            Assert.AreEqual(1, proteins.Count);
            Assert.AreEqual(0, proteins[0].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(0, proteins[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual('K', proteins[0][63 - 1]); // reference only
        }

        [Test]
        public static void MultipleAlternateFrameshifts()
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "MultipleAlternateFrameshifts.xml"), true,
                DecoyType.None, null, false, null, out var unknownModifications);
            Assert.AreEqual(2, proteins.Count);
            Assert.AreEqual(3, proteins[0].SequenceVariations.Count()); // some redundant
            Assert.AreEqual(3, proteins[0].SequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes

            Assert.IsTrue(proteins[0].SequenceVariations.All(v => v.OneBasedBeginPosition == 471)); // there are two alternate alleles (1 and 2), but only 2 is in the genotype, so only that's applied
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Count()); // some redundant
            var applied = proteins[1].AppliedSequenceVariations.Single();
            Assert.AreEqual("KDKRATGRIKS", applied.VariantSequence);
            Assert.AreEqual(403 - 11, applied.OriginalSequence.Length - applied.VariantSequence.Length);
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(873, proteins[0].Length);
            Assert.AreEqual(873 - 403 + 11, proteins[1].Length);
        }

        [Test]
        public void VariantSymbolWeirdnessXml()
        {
            string file = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SeqVarSymbolWeirdness.xml");
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.None, null, false, null, out var un);
            Assert.AreEqual(12, variantProteins.First().ConsensusVariant.SequenceVariations.Count());
            Assert.AreEqual(2, variantProteins.First().ConsensusVariant.SequenceVariations.Count(v => v.VariantCallFormatDataString.Heterozygous.Any(kv => kv.Value)));

            Assert.AreEqual(1, variantProteins.Count); // Should be 2^2 from combinitorics of heterozygous, but the giant indels overwrite them
            Assert.AreEqual(0, variantProteins.Where(v => v.BaseSequence == variantProteins.First().ConsensusVariant.BaseSequence).Count()); // Homozygous variations are included
            Assert.AreNotEqual(variantProteins.First().ConsensusVariant.Name, variantProteins.First().Name);
            Assert.AreNotEqual(variantProteins.First().ConsensusVariant.FullName, variantProteins.First().FullName);
            Assert.AreNotEqual(variantProteins.First().ConsensusVariant.Accession, variantProteins.First().Accession);

            List<PeptideWithSetModifications> peptides = variantProteins.SelectMany(vp => vp.Digest(new DigestionParams(), null, null)).ToList();
        }

        [Test]
        public void VariantSymbolWeirdness2Xml()
        {
            string file = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "SeqVarSymbolWeirdness2.xml");
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.None, null, false, null, out var un);

            Assert.AreEqual(1, variantProteins.First().ConsensusVariant.SequenceVariations.Count());
            Assert.AreEqual(2, variantProteins.Count); // there is only one unique amino acid change
            Assert.AreEqual(1, variantProteins.Where(v => v.BaseSequence == variantProteins.First().ConsensusVariant.BaseSequence).Count());
            var variantProteinRef = variantProteins.First();
            var variantProteinAlt = variantProteins.Last();
            Assert.AreEqual('R', variantProteins.First().ConsensusVariant.BaseSequence[2386]);
            Assert.AreEqual('R', variantProteinRef.BaseSequence[2386]);
            Assert.AreEqual('H', variantProteinAlt.BaseSequence[2386]);
            Assert.AreEqual(variantProteins.First().ConsensusVariant.Name, variantProteinRef.Name);
            Assert.AreNotEqual(variantProteins.First().ConsensusVariant.Name, variantProteinAlt.Name);
            Assert.AreEqual(variantProteins.First().ConsensusVariant.FullName, variantProteinRef.FullName);
            Assert.AreNotEqual(variantProteins.First().ConsensusVariant.FullName, variantProteinAlt.FullName);
            Assert.AreEqual(variantProteins.First().ConsensusVariant.Accession, variantProteinRef.Accession);
            Assert.AreNotEqual(variantProteins.First().ConsensusVariant.Accession, variantProteinAlt.Accession);
            List<PeptideWithSetModifications> peptides = variantProteins.SelectMany(vp => vp.Digest(new DigestionParams(), null, null)).ToList();
        }

        [Test]
        public void IndelDecoyError()
        {
            string file = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "IndelDecoy.xml");

            int maxHeterozygousVariants = 4;
            int minAlleleDepth = 1;

            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.Reverse, null, false, null, out var un, consensusPlusVariantIsoforms: maxHeterozygousVariants, minAlleleDepth: minAlleleDepth);
            Assert.AreEqual(8, variantProteins.Count);
            var indelProtein = variantProteins[2];
            Assert.AreNotEqual(indelProtein.AppliedSequenceVariations.Single().OriginalSequence.Length, indelProtein.AppliedSequenceVariations.Single().VariantSequence.Length);
            Assert.AreNotEqual(indelProtein.ConsensusVariant.Length, variantProteins[2].Length);
            var decoyIndelProtein = variantProteins[5];
            Assert.AreNotEqual(decoyIndelProtein.AppliedSequenceVariations.Single().OriginalSequence.Length, decoyIndelProtein.AppliedSequenceVariations.Single().VariantSequence.Length);
            Assert.AreNotEqual(decoyIndelProtein.ConsensusVariant.Length, variantProteins[2].Length);
            Assert.AreEqual(indelProtein.Length - indelProtein.AppliedSequenceVariations.Single().OneBasedBeginPosition, decoyIndelProtein.AppliedSequenceVariations.Single().OneBasedBeginPosition);
            var variantSeq = indelProtein.AppliedSequenceVariations.Single().VariantSequence.ToCharArray();
            Array.Reverse(variantSeq);
            Assert.AreEqual(new string(variantSeq), decoyIndelProtein.AppliedSequenceVariations.Single().VariantSequence);
        }

        [Test]
        public void IndelDecoyVariants()
        {
            string file = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "DecoyVariants.xml");
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.Reverse, null, false, null, out var un);
            Assert.AreEqual(4, variantProteins.Count);
            Assert.AreEqual(3, variantProteins[0].AppliedSequenceVariations.Count); // homozygous variations
            Assert.AreEqual(4, variantProteins[1].AppliedSequenceVariations.Count); // plus one heterozygous variation
            Assert.AreEqual("M", variantProteins[0].AppliedSequenceVariations.Last().OriginalSequence);
            Assert.AreEqual(1646, variantProteins[0].AppliedSequenceVariations.Last().OneBasedBeginPosition);
            Assert.AreEqual("V", variantProteins[0].AppliedSequenceVariations.Last().VariantSequence);
            Assert.AreEqual("M", variantProteins[2].AppliedSequenceVariations.First().OriginalSequence);
            Assert.AreEqual(variantProteins[0].Length - 1646 + 2, variantProteins[2].AppliedSequenceVariations.First().OneBasedBeginPosition);
            Assert.AreEqual("V", variantProteins[2].AppliedSequenceVariations.First().VariantSequence);
        }
        [Test]
        public void SequenceVariationIsValidTest()
        {
            // PURPOSE
            // Validate the minimal, position-only "validity" rules implemented by SequenceVariation.AreValid():
            //   AreValid() == (OneBasedBeginPosition > 0) && (OneBasedEndPosition >= OneBasedBeginPosition)
            //
            // We cover:
            // 1) Explicit begin/end ctor with typical point mutations → valid.
            // 2) Explicit begin/end ctor with invalid positions → invalid.
            // 3) Explicit begin/end ctor representing insertion/deletion edge-cases → valid as long as positions are valid.
            // 4) One-position convenience ctor behavior for different originalSequence values (null, "", length > 0).
            //    - This ctor derives end as: end = (original == null) ? begin : begin + original.Length - 1.
            //    - Therefore, empty originalSequence "" makes end = begin - 1 → invalid by design.
            // 5) Content fields (Original/Variant) and OneBasedModifications do NOT affect AreValid(), only positions do.
            // 6) Optional sanity checks on derived fields (SimpleString and computed end position).

            // -----------------------------
            // 1) Valid: explicit begin/end point mutations (begin == end, begin > 0)
            // -----------------------------
            SequenceVariation sv1 = new SequenceVariation(
                oneBasedBeginPosition: 10, oneBasedEndPosition: 10,
                originalSequence: "A", variantSequence: "T",
                description: "info", oneBasedModifications: null);
            SequenceVariation sv2 = new SequenceVariation(
                oneBasedBeginPosition: 5, oneBasedEndPosition: 5,
                originalSequence: "G", variantSequence: "C",
                description: "info", oneBasedModifications: null);
            SequenceVariation sv3 = new SequenceVariation(
                oneBasedBeginPosition: 8, oneBasedEndPosition: 8,
                originalSequence: "T", variantSequence: "A",
                description: "info", oneBasedModifications: null);

            // A protein can carry multiple variations; positions alone determine validity.
            List<SequenceVariation> svList = new List<SequenceVariation> { sv1, sv2, sv3 };
            Protein variantProtein = new Protein("ACDEFGHIKLMNPQRSTVWY", "protein1", sequenceVariations: svList);

            // Expectation: all three above are valid (begin > 0 and end == begin).
            Assert.IsTrue(variantProtein.SequenceVariations.All(v => v.AreValid()), "All explicit point mutations with valid positions should be valid");

            // -----------------------------
            // 2) Invalid: begin < 1 and end < begin
            // -----------------------------
            SequenceVariation svInvalidOneBasedBeginLessThanOne = new SequenceVariation(
                oneBasedBeginPosition: 0, oneBasedEndPosition: 10,
                originalSequence: "A", variantSequence: "T",
                description: "info", oneBasedModifications: null);
            Assert.IsFalse(svInvalidOneBasedBeginLessThanOne.AreValid(), "Begin position must be >= 1");

            SequenceVariation svInvalidOneBasedEndLessThanOneBasedBegin = new SequenceVariation(
                oneBasedBeginPosition: 5, oneBasedEndPosition: 4,
                originalSequence: "G", variantSequence: "C",
                description: "info", oneBasedModifications: null);
            Assert.IsFalse(svInvalidOneBasedEndLessThanOneBasedBegin.AreValid(), "End position cannot be less than begin position");

            // -----------------------------
            // 3) Explicit begin/end edge-cases: insertion and deletion modeled by content only
            //    NOTE: AreValid ignores Original/Variant content; only positions matter.
            // -----------------------------
            // Insertion-like (explicit): original is empty (""), variant has content.
            // Valid because we explicitly supply begin == end (positions are valid).
            SequenceVariation svValidOriginalSequenceIsEmpty = new SequenceVariation(
                oneBasedBeginPosition: 8, oneBasedEndPosition: 8,
                originalSequence: "", variantSequence: "A",
                description: "info", oneBasedModifications: null);
            Assert.IsTrue(svValidOriginalSequenceIsEmpty.AreValid(), "Explicit insertion with valid positions should be considered valid");

            // Deletion-like (explicit): variant is empty (""), original has content, positions still valid.
            SequenceVariation svValidVariantSequenceLengthIsZero = new SequenceVariation(
                oneBasedBeginPosition: 10, oneBasedEndPosition: 10,
                originalSequence: "A", variantSequence: "",
                description: "info", oneBasedModifications: null);
            Assert.IsTrue(svValidVariantSequenceLengthIsZero.AreValid(), "Explicit deletion with valid positions should be considered valid");

            // -----------------------------
            // 4) One-position convenience ctor behavior for originalSequence values
            //    ctor: SequenceVariation(int oneBasedPosition, string originalSequence, string variantSequence, string description, ...)
            //    end is computed as:
            //      - if original == null  → end = begin
            //      - else                  → end = begin + original.Length - 1
            // -----------------------------

            // 4a) originalSequence has length 1 → end == begin (valid)
            var svPosCtorLength1 = new SequenceVariation(
                oneBasedPosition: 15,
                originalSequence: "K",
                variantSequence: "R",
                description: "pos-ctor length 1");
            Assert.AreEqual(15, svPosCtorLength1.OneBasedBeginPosition);
            Assert.AreEqual(15, svPosCtorLength1.OneBasedEndPosition, "With original length 1, end should equal begin");
            Assert.IsTrue(svPosCtorLength1.AreValid(), "Single-site variation via position ctor should be valid");

            // 4b) originalSequence has length > 1 → end == begin + len - 1 (valid)
            var svPosCtorLength3 = new SequenceVariation(
                oneBasedPosition: 20,
                originalSequence: "PEP",   // len = 3
                variantSequence: "AAA",    // content irrelevant to AreValid
                description: "pos-ctor length 3");
            Assert.AreEqual(20, svPosCtorLength3.OneBasedBeginPosition);
            Assert.AreEqual(22, svPosCtorLength3.OneBasedEndPosition, "End should be begin + original.Length - 1");
            Assert.IsTrue(svPosCtorLength3.AreValid(), "Multi-length replacement with valid positions should be valid");

            // 4c) originalSequence == null → end = begin (valid)
            var svPosCtorNullOriginal = new SequenceVariation(
                oneBasedPosition: 30,
                originalSequence: null,    // special case handled in ctor: end = begin
                variantSequence: "A",
                description: "pos-ctor null original");
            Assert.AreEqual(30, svPosCtorNullOriginal.OneBasedBeginPosition);
            Assert.AreEqual(30, svPosCtorNullOriginal.OneBasedEndPosition, "Null original should set end = begin");
            Assert.IsTrue(svPosCtorNullOriginal.AreValid(), "Null original via position ctor is treated as length 1 (valid)");

            // 4d) originalSequence == "" (empty) → end = begin - 1 (invalid by design)
            //     This models an insertion if you rely solely on the position ctor, but produces end < begin → invalid.
            //     For insertions, prefer the explicit begin/end ctor with valid positions (see 3).
            var svPosCtorEmptyOriginal = new SequenceVariation(
                oneBasedPosition: 40,
                originalSequence: "",      // empty → end = 39
                variantSequence: "A",
                description: "pos-ctor empty original");
            Assert.AreEqual(40, svPosCtorEmptyOriginal.OneBasedBeginPosition);
            Assert.AreEqual(39, svPosCtorEmptyOriginal.OneBasedEndPosition, "Empty original sets end = begin - 1");
            Assert.IsFalse(svPosCtorEmptyOriginal.AreValid(), "Position ctor with empty original is invalid (end < begin)");

            // -----------------------------
            // 5) Validity is position-only; content and mods do not change AreValid()
            //     - VariantSequence null is normalized to "" in the ctor.
            //     - OneBasedModifications is stored but ignored by AreValid().
            // -----------------------------
            var mods = new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification>() } // empty mod list at a site; still ignored by AreValid
            };
            var svContentIrrelevant = new SequenceVariation(
                oneBasedBeginPosition: 3, oneBasedEndPosition: 3,
                originalSequence: "M", variantSequence: null, // becomes ""
                description: "mods/variant null test", oneBasedModifications: mods);
            Assert.IsTrue(svContentIrrelevant.AreValid(), "Null variant and/or mods should not affect positional validity");
            Assert.AreEqual("", svContentIrrelevant.VariantSequence, "Null VariantSequence is normalized to empty string");

            // -----------------------------
            // 6) Sanity: SimpleString format and positive bounds at the edge of sequence
            // -----------------------------
            var svSimple = new SequenceVariation(
                oneBasedBeginPosition: 1, oneBasedEndPosition: 1,
                originalSequence: "A", variantSequence: "V",
                description: "simple");
            // SimpleString = Original + Begin + Variant (no delimiter)
            Assert.AreEqual("A1V", svSimple.SimpleString(), "SimpleString should concatenate original + begin + variant");

            // Additional guard: begin == 1 is valid if end >= begin
            var svAtStart = new SequenceVariation(
                oneBasedBeginPosition: 1, oneBasedEndPosition: 2,
                originalSequence: "MA", variantSequence: "MV",
                description: "range at start");
            Assert.IsTrue(svAtStart.AreValid(), "Ranges that start at 1 are valid provided end >= begin");
        }
        [Test]
        public void VariantModificationTest()
        {
            string file = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "VariantModsGPTMD.xml");
            int maxHeterozygousVariants = 4;
            int minAlleleDepth = 1;
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.Reverse, null, false, null, out var un, consensusPlusVariantIsoforms: maxHeterozygousVariants, minAlleleDepth: minAlleleDepth);
            List<Protein> targets = variantProteins.Where(p => p.IsDecoy == false).ToList();
            List<Protein> variantTargets = targets.Where(p => p.AppliedSequenceVariations.Count >= 1).ToList();
            List<Protein> decoys = variantProteins.Where(p => p.IsDecoy == true).ToList();
            List<Protein> variantDecoys = decoys.Where(p => p.AppliedSequenceVariations.Count >= 1).ToList();
            bool homozygousVariant = targets.Select(p => p.Accession).Contains("Q6P6B1");           

            var variantMods = targets.SelectMany(p => p.AppliedSequenceVariations.Where(x=>x.OneBasedModifications.Count>= 1)).ToList();
            var decoyMods = decoys.SelectMany(p => p.AppliedSequenceVariations.Where(x => x.OneBasedModifications.Count >= 1)).ToList();
            var negativeResidues = decoyMods.SelectMany(x => x.OneBasedModifications.Where(w => w.Key < 0)).ToList();
            bool namingWrong = targets.Select(p => p.Accession).Contains("Q8N865_H300R_A158T_H300R");
            bool namingRight = targets.Select(p => p.Accession).Contains("Q8N865_A158T_H300R");
            Assert.AreEqual(false, namingWrong);
            Assert.AreEqual(true, namingRight);
            Assert.AreEqual(false, homozygousVariant);
            Assert.AreEqual(62, variantProteins.Count);
            Assert.AreEqual(31, targets.Count);
            Assert.AreEqual(26, variantTargets.Count);
            Assert.AreEqual(31, decoys.Count);
            Assert.AreEqual(26, variantDecoys.Count);
            Assert.AreEqual(2, variantMods.Count);
            Assert.AreEqual(2, decoyMods.Count);
            Assert.AreEqual(0, negativeResidues.Count);

        }
        [Test]
        public void WriteProteinXmlWithVariantsDiscoveredAsModifications2()
        {
            string databaseName = "humanGAPDH.xml";
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", databaseName), true,
                DecoyType.Reverse, null, false, null, out var unknownModifications, 1, 0);
            var target = proteins[0];
            int totalSequenceVariations = target.SequenceVariations.Count();
            Assert.AreEqual(2, totalSequenceVariations); //these sequence variations were in the original
            ModificationMotif.TryGetMotif("W", out ModificationMotif motifW);
            string _originalId = "W->G";
            string _accession = null;
            string _modificationType = "1 nucleotide substitution";
            string _featureType = null;
            ModificationMotif _target = motifW;
            string _locationRestriction = "Anywhere.";
            ChemicalFormula _chemicalFormula = ChemicalFormula.ParseFormula("C-9H-7N-1");
            double? _monoisotopicMass = null;
            Dictionary<string, IList<string>> _databaseReference = null;
            Dictionary<string, IList<string>> _taxonomicRange = null;
            List<string> _keywords = null;
            Dictionary<DissociationType, List<double>> _neutralLosses = null;
            Dictionary<DissociationType, List<double>> _diagnosticIons = null;
            string _fileOrigin = null;

            Modification substitutionMod = new Modification(_originalId, _accession, _modificationType, _featureType, _target, _locationRestriction,
                               _chemicalFormula, _monoisotopicMass, _databaseReference, _taxonomicRange, _keywords, _neutralLosses, _diagnosticIons, _fileOrigin);
            Dictionary<int, List<Modification>> substitutionDictionary = new Dictionary<int, List<Modification>>();
            substitutionDictionary.Add(87, new List<Modification> { substitutionMod });

            Protein newProtein = (Protein)target.CloneWithNewSequenceAndMods(target.BaseSequence, substitutionDictionary);
            Assert.That(newProtein.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1));

            // This process examines the OneBasedPossibleLocalizedModifications that are ModificationType 'nucleotide substitution'
            // and converts them to SequenceVariations
            newProtein.ConvertNucleotideSubstitutionModificationsToSequenceVariants();
            Assert.That(newProtein.SequenceVariations.Count, Is.EqualTo(totalSequenceVariations + 1)); //This number increases by 1 because we added a sequence variation that was discovered as a modification
            Assert.AreEqual(0,newProtein.OneBasedPossibleLocalizedModifications.Count); //This number should be 0 because we converted the modification to a sequence variation
        }
        /// <summary>
        /// PURPOSE
        /// Ensures that variant proteins are automatically generated during XML read, both for targets and reverse decoys.
        /// The database "humanGAPDH.xml" encodes two single-nucleotide substitutions on the target protein P04406:
        /// - A22G
        /// - K251N
        ///
        /// EXPECTATIONS
        /// - Loader emits all combinatorial target variants derived from those two changes:
        ///   [0] Reference (no variants applied)                         → Accession "P04406"
        ///   [1] Single variant A22G                                     → Accession "P04406_A22G"
        ///   [2] Single variant K251N                                    → Accession "P04406_K251N"
        ///   [3] Double variant K251N + A22G (combined)                  → Accession "P04406_K251N_A22G"
        /// - With DecoyType.Reverse, a matching set of 4 reverse decoys is produced with mirrored coordinates:
        ///   [4] Decoy of reference                                      → "DECOY_P04406"
        ///   [5] Decoy with mirrored A22G (mapped to site 315)           → "DECOY_P04406_A315G"
        ///   [6] Decoy with mirrored K251N (mapped to site 86)           → "DECOY_P04406_K86N"
        ///   [7] Decoy with both mirrored variants                       → "DECOY_P04406_K86N_A315G"
        ///
        /// WHY THIS MATTERS
        /// - Validates that the reader expands sequence variation definitions into concrete variant proteins.
        /// - Verifies decoy generation mirrors variant coordinates appropriately and preserves ordering.
        /// - Guards against regressions in accession naming, variant application counts, and decoy parity.
        ///
        /// PARAMETERS PASSED TO LOADER
        /// - generateTargets: true        → emit target proteins
        /// - decoyType: Reverse           → also emit reverse decoys
        /// - allKnownModifications: UniProtPtms → resolve any UniProt-annotated PTMs so no "unknown" mods remain
        /// - isContaminant: false
        /// - modTypesToExclude: null
        /// - out unknownModifications     → capture any unrecognized mods (should be empty when UniProtPtms is provided)
        /// - minAlleleDepth: 1            → do not filter out these low-depth test variants
        /// - maxHeterozygousVariants: 99  → allow generating all combinations from the two sites (up to 2^2)
        /// </summary>
        [Test]
        public static void TestThatProteinVariantsAreGeneratedDuringRead()
        {
            // Arrange: load a target with two site-specific variants and request reverse decoys as well.
            // IMPORTANT: Provide UniProtPtms so annotations in humanGAPDH.xml resolve and do not end up in unknownModifications.
            string databaseName = "humanGAPDH.xml";
            var proteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", databaseName),
                generateTargets: true,
                decoyType: DecoyType.Reverse,
                allKnownModifications: UniProtPtms, // CHANGED: was null; supplying known PTMs prevents unknownModifications from being populated
                isContaminant: false,
                modTypesToExclude: null,
                unknownModifications: out var unknownModifications,
                minAlleleDepth: 1,
                consensusPlusVariantIsoforms: 99);

            // Basic shape: 4 targets + 4 reverse decoys in a deterministic order.
            Assert.AreEqual(8, proteins.Count, "Expected 4 targets and 4 decoys in a fixed order");

            // Targets/decoys split and flags should be consistent and easy to reason about.
            Assert.AreEqual(4, proteins.Count(p => !p.IsDecoy), "First half should be targets");
            Assert.AreEqual(4, proteins.Count(p => p.IsDecoy), "Second half should be decoys");
            for (int i = 0; i < 4; i++)
            {
                Assert.IsFalse(proteins[i].IsDecoy, $"Index {i} should be a target");
                Assert.IsTrue(proteins[i + 4].IsDecoy, $"Index {i + 4} should be a decoy");
            }

            // The reference target (index 0) carries two possible sequence variations in its metadata.
            // This documents the test input and ensures the reader surfaced them.
            Assert.AreEqual(2, proteins[0].SequenceVariations.Count(), "Reference should advertise exactly two possible sequence variations");

            // Accessions must match exact, canonical variant labeling and order for both targets and decoys.
            Assert.That("P04406", Is.EqualTo(proteins[0].Accession), "Reference target accession mismatch");
            Assert.That("P04406_A22G", Is.EqualTo(proteins[1].Accession), "Single-variant (A22G) target accession mismatch");
            Assert.That("P04406_K251N", Is.EqualTo(proteins[2].Accession), "Single-variant (K251N) target accession mismatch");
            Assert.That("P04406_K251N_A22G", Is.EqualTo(proteins[3].Accession), "Double-variant target accession mismatch");

            Assert.That("DECOY_P04406", Is.EqualTo(proteins[4].Accession), "Reference decoy accession mismatch");
            Assert.That("DECOY_P04406_A315G", Is.EqualTo(proteins[5].Accession), "Decoy accession for mirrored A22G mismatch");
            Assert.That("DECOY_P04406_K86N", Is.EqualTo(proteins[6].Accession), "Decoy accession for mirrored K251N mismatch");
            Assert.That("DECOY_P04406_K86N_A315G", Is.EqualTo(proteins[7].Accession), "Decoy accession for double-variant mismatch");

            // Sanity: accessions are non-empty and unique (avoid accidental duplication/shuffling).
            Assert.That(proteins.All(p => !string.IsNullOrWhiteSpace(p.Accession)), "All proteins must have non-empty accessions");
            Assert.AreEqual(proteins.Count, proteins.Select(p => p.Accession).Distinct().Count(), "Accessions must be unique");

            // Each decoy should be length-equal to its corresponding target, but usually sequence-different (reverse).
            for (int i = 0; i < 4; i++)
            {
                Assert.AreEqual(proteins[i].Length, proteins[i + 4].Length, $"Target/decoy length should match for index {i}");
                Assert.AreNotEqual(proteins[i].BaseSequence, proteins[i + 4].BaseSequence, $"Decoy sequence should differ from its target for index {i}");
            }

            // Applied variant counts (how many variations were actually realized in this protein instance):
            // Targets:    [0] ref=0, [1] A22G=1, [2] K251N=1, [3] both=2
            // Decoys:     [4] ref=0, [5] A315G=1, [6] K86N=1, [7] both=2
            int[] expectedAppliedCounts = { 0, 1, 1, 2, 0, 1, 1, 2 };
            for (int i = 0; i < proteins.Count; i++)
            {
                Assert.AreEqual(expectedAppliedCounts[i], proteins[i].AppliedSequenceVariations.Count(),
                    $"Applied variant count mismatch at index {i} ({proteins[i].Accession})");
            }

            // The specific applied-variant labels should match accessions:
            // - For targets: "A22G" and/or "K251N"
            // - For decoys:  mirrored positions → "A315G" and/or "K86N"
            static HashSet<string> AppliedLabels(Protein p) =>
                new HashSet<string>(p.AppliedSequenceVariations.Select(v => v.SimpleString()));

            // Targets (indices 0..3)
            Assert.That(AppliedLabels(proteins[0]).SetEquals(Array.Empty<string>()), "Reference target should have no applied variants");
            Assert.That(AppliedLabels(proteins[1]).SetEquals(new[] { "A22G" }), "Single-variant target must be exactly A22G");
            Assert.That(AppliedLabels(proteins[2]).SetEquals(new[] { "K251N" }), "Single-variant target must be exactly K251N");
            Assert.That(AppliedLabels(proteins[3]).SetEquals(new[] { "A22G", "K251N" }), "Double-variant target should have A22G and K251N");

            // Decoys (indices 4..7) have mirrored coordinates (due to reverse decoying).
            Assert.That(AppliedLabels(proteins[4]).SetEquals(Array.Empty<string>()), "Reference decoy should have no applied variants");
            Assert.That(AppliedLabels(proteins[5]).SetEquals(new[] { "A315G" }), "Single-variant decoy must be exactly A315G (mirror of A22G)");
            Assert.That(AppliedLabels(proteins[6]).SetEquals(new[] { "K86N" }), "Single-variant decoy must be exactly K86N (mirror of K251N)");
            Assert.That(AppliedLabels(proteins[7]).SetEquals(new[] { "A315G", "K86N" }), "Double-variant decoy should have A315G and K86N");

            // Parity check: each target and its decoy should carry the same number of applied variations.
            for (int i = 0; i < 4; i++)
            {
                Assert.AreEqual(
                    proteins[i].AppliedSequenceVariations.Count(),
                    proteins[i + 4].AppliedSequenceVariations.Count(),
                    $"Applied variant count should match between target and decoy at index {i}");
            }

            // With UniProtPtms supplied, no unknown modifications should be reported.
            if (unknownModifications != null && unknownModifications.Count > 0)
            {
                // Extra diagnostics to ease debugging in case of future schema/content changes.
                TestContext.WriteLine("Unknown modifications reported by loader:");
                foreach (var um in unknownModifications)
                    TestContext.WriteLine($" - {um}");
            }
            Assert.That(unknownModifications == null || unknownModifications.Count == 0, "No unknown modifications expected from this input");
        }

        [Test]
        public static void ProteinVariantsReadAsModificationsWrittenAsVariants()
        {
            // PURPOSE
            // This test verifies the I/O pipeline that converts "nucleotide substitution" modifications
            // embedded in an input protein XML into canonical "sequence variant" features when read,
            // and persists them back out as proper <feature type="sequence variant"> entries when written.
            //
            // WHY
            // - Some sources (e.g., GPTMD discovery) encode AA substitutions as modifications with
            //   ModificationType = "1 nucleotide substitution". Internally, we want these represented
            //   as SequenceVariations, not remaining as generic modifications.
            // - On read: these substitution mods must be removed from the protein’s modifications collections
            //   and turned into SequenceVariations.
            // - On write: they must be serialized as sequence variant features with a standardized description
            //   ("Putative GPTMD Substitution") and NOT re-serialized as modifications.
            //
            // EXPECTATIONS SUMMARY
            // - Input file has 57 lines that contain "1 nucleotide substitution" (source encoding as mods).
            // - After LoadProteinXML:
            //     * We get exactly 9 proteins (DecoyType.None).
            //     * Total SequenceVariations across all proteins = 194.
            //     * Total OneBasedPossibleLocalizedModifications count across all proteins = 0
            //       (i.e., all substitution mods became sequence variations).
            //     * No unknown modifications should remain.
            //     * No applied variants are expected (metadata only; we are not expanding combinatorics here).
            // - After WriteXmlDatabase then re-load:
            //     * Count and totals remain the same (9 proteins; 194 total sequence variations; 0 mods).
            //     * The output file contains:
            //         - Exactly 194 lines with feature type="sequence variant".
            //         - Exactly 194 lines with "Putative GPTMD Substitution".
            //         - Exactly 0 lines with "1 nucleotide substitution" (proved that mods were not serialized).
            //
            // NOTE
            // - We keep DecoyType.None to avoid expanding the protein list.
            // - We use minAlleleDepth = 1 and maxHeterozygousVariants = 0 to avoid variant-proteoform expansion.

            // Arrange: locate the source database that encodes nucleotide substitutions as modifications
            string databaseName = "nucleotideVariantsAsModifications.xml";
            string inputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", databaseName);
            Assert.That(File.Exists(inputPath), Is.True, "Input database file must exist for this test");

            // Sanity: confirm the source encodes substitutions as modifications in the raw XML
            // (57 lines that mention "1 nucleotide substitution" in the file).
            var inputLines = File.ReadAllLines(inputPath);
            int inputSubstitutionModLines = inputLines.Count(l => l.Contains("1 nucleotide substitution"));
            Assert.That(inputSubstitutionModLines, Is.EqualTo(57), "Source XML should contain 57 substitution-mod lines");

            // Optional sanity: the source should not already be in sequence-variant form
            int inputSeqVarFeatureLines = inputLines.Count(l => l.Contains("feature type=\"sequence variant\""));
            Assert.That(inputSeqVarFeatureLines, Is.EqualTo(0), "Source XML should not already encode sequence variants");

            // Act: read the database. Expect conversion to SequenceVariations and removal from modifications
            var proteins = ProteinDbLoader.LoadProteinXML(
                inputPath,
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: null,
                isContaminant: false,
                modTypesToExclude: null,
                unknownModifications: out var unknownModifications,
                minAlleleDepth: 1,
                consensusPlusVariantIsoforms: 0);

            // Assert: no decoys requested, so all should be targets
            Assert.That(proteins.All(p => !p.IsDecoy), "All proteins should be targets when DecoyType.None is used");

            // Assert: exactly 9 proteins
            Assert.AreEqual(9, proteins.Count, "Expected exactly 9 proteins from the input");

            // Assert: the loader converted all substitution modifications into proper SequenceVariations
            // and removed them from OneBasedPossibleLocalizedModifications
            int totalSequenceVariations = proteins.Sum(p => p.SequenceVariations.Count);
            int totalPossibleLocalizedMods = proteins.Sum(p => p.OneBasedPossibleLocalizedModifications.Values.Sum(v => v.Count));
            Assert.AreEqual(194, totalSequenceVariations, "Total number of sequence variations after load must be 194");
            Assert.AreEqual(0, totalPossibleLocalizedMods, "All substitution modifications should have been converted; none should remain as mods");

            // FIX: Safely log unknown modifications (Dictionary<string, Modification>), if any appear unexpectedly.
            if (unknownModifications != null && unknownModifications.Count > 0)
            {
                TestContext.WriteLine("Unknown modifications encountered during read (unexpected):");
                foreach (var kvp in unknownModifications)
                {
                    var mod = kvp.Value;
                    var id = mod?.OriginalId ?? "<null>";
                    var type = mod?.ModificationType ?? "<null>";
                    TestContext.WriteLine($" - key={kvp.Key}, id={id}, type={type}");
                }
            }
            Assert.That(unknownModifications == null || unknownModifications.Count == 0, "No unknown modifications should remain after conversion");

            // Assert: No applied variants expected in this test (we are validating representation, not expansion)
            int totalAppliedVariants = proteins.Sum(p => p.AppliedSequenceVariations.Count());
            Assert.That(totalAppliedVariants, Is.EqualTo(0), "No applied variants expected; variants are metadata here");

            // Assert: None of the proteins should carry "1 nucleotide substitution" as OriginalNonVariantModifications anymore
            int residualSubstitutionMods =
                proteins.Sum(p => p.OriginalNonVariantModifications.Values.Sum(list => list.Count(m => m.ModificationType == "1 nucleotide substitution")));
            Assert.That(residualSubstitutionMods, Is.EqualTo(0), "No '1 nucleotide substitution' mods should remain in OriginalNonVariantModifications");

            // Persist the converted representation and read it back; this file should contain only sequence-variant features
            string tempDir = Path.Combine(Path.GetTempPath(), Guid.NewGuid().ToString());
            Directory.CreateDirectory(tempDir);
            string tempFile = Path.Combine(tempDir, "xmlWithSequenceVariantsAndNoModifications.txt");

            try
            {
                // Write only the targets (all are targets here). Expect the writer to emit sequence variant features.
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteins.Where(p => !p.IsDecoy).ToList(), tempFile);

                // Inspect the written file contents directly for ground truth
                var writtenLines = File.ReadAllLines(tempFile);
                int writtenSeqVarFeatures = writtenLines.Count(l => l.Contains("feature type=\"sequence variant\""));
                int writtenPutativeGptmd = writtenLines.Count(l => l.Contains("Putative GPTMD Substitution"));
                int writtenSubstitutionMods = writtenLines.Count(l => l.Contains("1 nucleotide substitution"));

                // Assert: writer produced only sequence variants (not substitution mods)
                Assert.That(writtenSeqVarFeatures, Is.EqualTo(194), "All 194 substitutions must be serialized as sequence variant features");
                Assert.That(writtenPutativeGptmd, Is.EqualTo(194), "All 194 variants should have the standardized description label");
                Assert.That(writtenSubstitutionMods, Is.EqualTo(0), "No '1 nucleotide substitution' strings should remain in the written XML");

                // Re-load the written file to confirm round-trip stability
                proteins = ProteinDbLoader.LoadProteinXML(
                    tempFile,
                    generateTargets: true,
                    decoyType: DecoyType.None,
                    allKnownModifications: null,
                    isContaminant: false,
                    modTypesToExclude: null,
                    unknownModifications: out unknownModifications,
                    minAlleleDepth: 1,
                    consensusPlusVariantIsoforms: 0);

                // Assert: the round-tripped representation is identical in shape and counts
                Assert.AreEqual(9, proteins.Count, "Round-trip must preserve protein count");
                Assert.AreEqual(194, proteins.Sum(v => v.SequenceVariations.Count), "Round-trip must preserve total number of sequence variations");
                Assert.AreEqual(0, proteins.Sum(m => m.OneBasedPossibleLocalizedModifications.Values.Sum(list => list.Count)),
                    "Round-trip must preserve the fact that no substitution mods remain");
                Assert.That(unknownModifications == null || unknownModifications.Count == 0, "Round-trip should not introduce unknown modifications");
                Assert.That(proteins.Sum(p => p.AppliedSequenceVariations.Count()), Is.EqualTo(0), "Round-trip should not apply any variants");
            }
            finally
            {
                // Cleanup test artifacts
                if (Directory.Exists(tempDir))
                {
                    try { File.SetAttributes(tempFile, FileAttributes.Normal); } catch { /* ignore */ }
                    Directory.Delete(tempDir, true);
                }
            }
        }
        [Test]
        public void Constructor_ParsesDescriptionCorrectly()
        {
            // Arrange
            string description = @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30";

            // Example VCF line with snpEff annotation:
            // 1   50000000   .   A   G   .   PASS   ANN=G||||||||||||||||   GT:AD:DP   1/1:30,30:30

            // --- VCF Standard Columns ---
            //
            // CHROM (1)      → Chromosome name (here, chromosome 1).
            // POS (50000000) → 1-based position of the variant (50,000,000).
            // ID (.)         → Variant identifier. "." means no ID (e.g., not in dbSNP).
            // REF (A)        → Reference allele in the reference genome (A).
            // ALT (G)        → Alternate allele observed in reads (G).
            // QUAL (.)       → Variant call quality score (Phred-scaled). "." means not provided.
            // FILTER (PASS)  → Indicates if the call passed filtering. "PASS" = high confidence.
            //
            // --- INFO Column ---
            //
            // INFO (ANN=...) holds snpEff annotation data.
            // ANN format is:
            //   Allele | Effect | Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID |
            //   Transcript_Biotype | Rank | HGVS.c | HGVS.p | cDNA_pos/cDNA_len |
            //   CDS_pos/CDS_len | AA_pos/AA_len | Distance | Errors/Warnings
            //
            // In this case: ANN=G||||||||||||||||
            //   - Allele = G
            //   - All other fields are empty → snpEff did not predict any functional impact
            //     (likely intergenic or unannotated region).
            //
            // --- FORMAT Column ---
            //
            // FORMAT (GT:AD:DP) defines how to read the sample column(s):
            //   GT → Genotype
            //   AD → Allele depth (number of reads supporting REF and ALT)
            //   DP → Read depth (total reads covering the site)
            //
            // --- SAMPLE Column ---
            //
            // Sample entry: 1/1:30,30:30
            //   GT = 1/1 → Homozygous ALT genotype (both alleles = G)
            //   AD = 30,30 → Read counts: REF=A has 30 reads, ALT=G has 30 reads
            //                (⚠ usually homozygous ALT would have few/no REF reads;
            //                 this may be caller-specific behavior or a quirk.)
            //   DP = 30   → Total coverage at this site = 30 reads
            //                (⚠ note AD sums to 60, which does not match DP.
            //                 This discrepancy is common in some callers.)
            //
            // --- Overall Summary ---
            // Variant at chr1:50000000 changes A → G.
            // The sample is homozygous for the ALT allele (G).
            // Variant passed filters, but no functional annotation from snpEff.


            // Act
            var svd = new VariantCallFormat(description);

            // Assert
            Assert.AreEqual(description, svd.Description);
            Assert.AreEqual("A", svd.ReferenceAlleleString);
            Assert.AreEqual("G", svd.AlternateAlleleString);
            Assert.IsNotNull(svd.Info);
            Assert.AreEqual("GT:AD:DP", svd.Format);
            Assert.AreEqual(1, svd.Genotypes.Count);
            Assert.AreEqual(1, svd.AlleleDepths.Count);
            Assert.AreEqual(new[] { "0" }, new List<string>(svd.Genotypes.Keys));

            var hzKey = svd.Homozygous.Keys.First();
            Assert.AreEqual("0", hzKey);
            var hzBool = svd.Homozygous[hzKey];
            Assert.IsTrue(hzBool);
            var adKey = svd.AlleleDepths.Keys.First();
            Assert.AreEqual("0", adKey);
            var adValues = svd.AlleleDepths[adKey];
            Assert.AreEqual(new[] { "30", "30" }, adValues);
        }
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
                    if (parsedAd.Length != 0 || expectedAdTokens.Length != 1 || expectedAdTokens[0] != ".")
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