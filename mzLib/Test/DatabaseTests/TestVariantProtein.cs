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
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.None, null, false, null, out var un);

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
            Assert.AreNotEqual(target.SequenceVariations.First().Description, decoy.SequenceVariations.First().Description); //decoys and target variations don't have the same desc.

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
            // test decoys and digestion
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "StopGain.xml"), true,
                DecoyType.Reverse, null, false, null, out var unknownModifications, minAlleleDepth: 400);
            Assert.AreEqual(2, proteins.Count);
            var targetPeps = proteins[0].Digest(new DigestionParams(), null, null).ToList();
            var decoyPeps = proteins[1].Digest(new DigestionParams(), null, null).ToList();
            //Assert.AreEqual(targetPeps.Sum(p => p.Length), decoyPeps.Sum(p => p.Length));
            //Assert.AreEqual(targetPeps.Count, decoyPeps.Count);
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
            Assert.AreEqual(2, variantProteins.First().ConsensusVariant.SequenceVariations.Count(v => v.Description.Heterozygous.Any(kv => kv.Value)));

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
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.Reverse, null, false, null, out var un);
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
            SequenceVariation sv1 = new SequenceVariation(10, 10, "A", "T", "info", null);
            SequenceVariation sv2 = new SequenceVariation(5, 5, "G", "C", "info", null);
            SequenceVariation sv3 = new SequenceVariation(8, 8, "T", "A", "info", null);
            List<SequenceVariation> svList = new List<SequenceVariation> { sv1, sv2, sv3 };

            Protein variantProtein = new Protein("ACDEFGHIKLMNPQRSTVWY", "protein1", sequenceVariations: svList);
            Assert.IsTrue(variantProtein.SequenceVariations.All(v => v.AreValid()));
            SequenceVariation svInvalidOneBasedBeginLessThanOne = new SequenceVariation(0, 10, "A", "T", "info", null);
            SequenceVariation svInvalidOneBasedEndLessThanOneBasedBegin = new SequenceVariation(5, 4, "G", "C", "info", null);
            SequenceVariation svValidOriginalSequenceIsEmpty = new SequenceVariation(8, 8, "", "A", "info", null);
            SequenceVariation svValidVariantSequenceLenthIsZero = new SequenceVariation(10, 10, "A", "", "info", null);
            Assert.IsFalse(svInvalidOneBasedBeginLessThanOne.AreValid());
            Assert.IsFalse(svInvalidOneBasedEndLessThanOneBasedBegin.AreValid());
            Assert.IsTrue(svValidOriginalSequenceIsEmpty.AreValid()); //This is valid because it is an insertion
            Assert.IsTrue(svValidVariantSequenceLenthIsZero.AreValid()); // This is valid because it is a deletion
        }
        [Test]
        public void VariantModificationTest()
        {
            string file = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "VariantModsGPTMD.xml");
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.Reverse, null, false, null, out var un);
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

        [Test]
        public static void TestThatProteinVariantsAreGeneratedDuringRead()
        {
            string databaseName = "humanGAPDH.xml";
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", databaseName), true,
                DecoyType.Reverse, null, false, null, out var unknownModifications, 1, 99);
            Assert.AreEqual(8, proteins.Count); // 4 target + 4 decoy
            Assert.AreEqual(2, proteins[0].SequenceVariations.Count()); // these sequence variations were in the original
            Assert.That("P04406", Is.EqualTo(proteins[0].Accession));
            Assert.That("P04406_A22G", Is.EqualTo(proteins[1].Accession));
            Assert.That("P04406_K251N", Is.EqualTo(proteins[2].Accession));
            Assert.That("P04406_K251N_A22G", Is.EqualTo(proteins[3].Accession));
            Assert.That("DECOY_P04406", Is.EqualTo(proteins[4].Accession));
            Assert.That("DECOY_P04406_A315G", Is.EqualTo(proteins[5].Accession));
            Assert.That("DECOY_P04406_K86N", Is.EqualTo(proteins[6].Accession));
            Assert.That("DECOY_P04406_K86N_A315G", Is.EqualTo(proteins[7].Accession));
        }
        [Test]
        public static void ProteinVariantsReadAsModificationsWrittenAsVariants()
        {
            string databaseName = "nucleotideVariantsAsModifications.xml";

            Assert.That(File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", databaseName)).Count(l => l.Contains("1 nucleotide substitution")), Is.EqualTo(57));

            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", databaseName), true,
                DecoyType.None, null, false, null, out var unknownModifications, 1, 0);
            Assert.AreEqual(9, proteins.Count); // 1 target
            Assert.AreEqual(194, proteins.Select(v=>v.SequenceVariations.Count).Sum()); // there are no sequence variations in the original proteins
            Assert.AreEqual(0, proteins.Select(m => m.OneBasedPossibleLocalizedModifications.Sum(list=>list.Value.Count)).Sum()); // there are 194 sequence variants as modifications in the original proteins

            string tempDir = Path.Combine(Path.GetTempPath(), Guid.NewGuid().ToString());
            Directory.CreateDirectory(tempDir);
            string tempFile = Path.Combine(tempDir, "xmlWithSequenceVariantsAndNoModifications.txt");

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteins.Where(p => !p.IsDecoy).ToList(), tempFile);
            proteins = ProteinDbLoader.LoadProteinXML(tempFile, true,
                DecoyType.None, null, false, null, out unknownModifications, 1, 0);
            Assert.AreEqual(9, proteins.Count); // 1 target
            Assert.AreEqual(194, proteins.Select(v => v.SequenceVariations.Count).Sum()); // there are 194 sequence variations in the revised proteins
            Assert.AreEqual(0, proteins.Select(m => m.OneBasedPossibleLocalizedModifications.Sum(list => list.Value.Count)).Sum()); // there are 0 sequence variants as modifications in the original proteins
            
            Assert.That(File.ReadAllLines(tempFile).Count(l => l.Contains("feature type=\"sequence variant\"")), Is.EqualTo(194));
            Assert.That(File.ReadAllLines(tempFile).Count(l => l.Contains("Putative GPTMD Substitution")), Is.EqualTo(194));
            Assert.That(File.ReadAllLines(tempFile).Count(l => l.Contains("1 nucleotide substitution")), Is.EqualTo(0));
            if (Directory.Exists(tempDir)) Directory.Delete(tempDir, true);
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
    }
}