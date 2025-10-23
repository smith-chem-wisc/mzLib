using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Omics.BioPolymer;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases;
using Stopwatch = System.Diagnostics.Stopwatch;
using Omics;
using Transcriptomics;
using MassSpectrometry;
using Chemistry;

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
            // Build a tiny variant suite on the same base sequence to exercise:
            // - SAV (single amino-acid variant) P→V at position 4
            // - MNV (multi-residue substitution) PT→KT at positions 4–5
            // - Insertion: P→PPP at position 4 (sequence grows)
            // - Deletion: PPP→P at positions 4–6 (sequence shrinks)
            // - Variant-scoped PTM (on the inserted segment) to ensure PTMs map to post-variation coordinates
            //
            // Why this test exists:
            // - Verifies end-to-start application order of variants yields stable, deterministic results.
            // - Confirms index remapping for applied variants, proteoform sequence output, and PTM carry-forward.
            // - Ensures serialization round-trip (write + re-load) preserves exactly the same proteoforms.

            // Create a dummy modification to attach to a variant (variant-scoped PTM at residue 5)
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            Modification mp = new Modification(
                _originalId: "mod", _accession: null, _modificationType: "type", _featureType: null,
                _target: motifP, _locationRestriction: "Anywhere.", _chemicalFormula: null,
                _monoisotopicMass: 42.01, _databaseReference: new Dictionary<string, IList<string>>(),
                _taxonomicRange: null, _keywords: null, _neutralLosses: null, _diagnosticIons: null, _fileOrigin: null);

            // Build five independent proteins, each carrying a single sequence variation.
            // Note: here the 5th ctor arg is a free-text description; we intentionally pass the VCF-like text there (VariantCallFormatData is null),
            // which exercises the "non-VCF/mixed" path and forces combinatorial expansion.
            List<Protein> proteinsWithSeqVars = new List<Protein>
            {
                // protein1: SAV P(4)→V                  => "MPEVTIDE"
                new Protein("MPEPTIDE", "protein1",
                    sequenceVariations: new List<SequenceVariation>
                    {
                        new SequenceVariation(4, 4, "P", "V", "",
                            @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30",
                            null) // VariantCallFormatData omitted on purpose for mixed path
                    }),

                // protein2: MNV PT(4-5)→KT              => "MPEKTIDE"
                new Protein("MPEPTIDE", "protein2",
                    sequenceVariations: new List<SequenceVariation>
                    {
                        new SequenceVariation(4, 5, "PT", "KT","",
                            @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30",
                            null)
                    }),

                // protein3: insertion P(4)→PPP          => "MPEPPPTIDE"
                new Protein("MPEPTIDE", "protein3",
                    sequenceVariations: new List<SequenceVariation>
                    {
                        new SequenceVariation(4, 4, "P", "PPP","",
                            @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30",
                            null)
                    }),

                // protein4: deletion PPP(4-6)→P on base "MPEPPPTIDE" => "MPEPTIDE"
                new Protein("MPEPPPTIDE", "protein4",
                    sequenceVariations: new List<SequenceVariation>
                    {
                        new SequenceVariation(4, 6, "PPP", "P","",
                            @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30",
                            null)
                    }),

                // protein5: insertion plus variant-PTM at position 5 => "MPEPPPTIDE" and mods at site 5
                new Protein("MPEPTIDE", "protein5",
                    sequenceVariations: new List<SequenceVariation>
                    {
                        new SequenceVariation(4, 4, "P", "PPP","",
                            @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30",
                            new Dictionary<int, List<Modification>> { { 5, new List<Modification> { mp } } })
                    }),
            };

            // ——— DIAGNOSTICS: Show whether each input protein’s variants look VCF-backed and how many proteoforms expansion returns
            for (int i = 0; i < proteinsWithSeqVars.Count; i++)
            {
                var p = proteinsWithSeqVars[i];
                var vars = p.SequenceVariations.ToList();
                TestContext.WriteLine($"[Input {i}] Accession={p.Accession} Base='{p.BaseSequence}' Variations={vars.Count}");
                foreach (var v in vars)
                {
                    bool vcfNonNull = v.VariantCallFormatData != null;
                    int gtCount = v.VariantCallFormatData?.Genotypes?.Count ?? 0;
                    string fmt = v.VariantCallFormatData?.Format ?? "(null)";
                    string alleles = v.VariantCallFormatData == null
                        ? "(no VCF)"
                        : $"REF={v.VariantCallFormatData.ReferenceAlleleString ?? "null"} ALT={v.VariantCallFormatData.AlternateAlleleString ?? "null"}";
                    TestContext.WriteLine($"    SV [{v.OneBasedBeginPosition}-{v.OneBasedEndPosition}] '{v.OriginalSequence}'→'{v.VariantSequence}' " +
                                          $"VCF? {vcfNonNull} GTs={gtCount} FORMAT='{fmt}' {alleles}");
                }
                var expanded = p.GetVariantBioPolymers();
                int baseOnly = expanded.Count(ep => !ep.AppliedSequenceVariations.Any());
                int variantApplied = expanded.Count(ep => ep.AppliedSequenceVariations.Any());
                TestContext.WriteLine($"    Expansion count={expanded.Count} BaseOnly={baseOnly} VariantApplied={variantApplied}");
                for (int j = 0; j < expanded.Count; j++)
                {
                    var ep = expanded[j];
                    var applied = ep.AppliedSequenceVariations.ToList();
                    string appliedSummary = applied.Count == 0
                        ? "(none)"
                        : string.Join(";", applied.Select(a => $"{a.OriginalSequence}{a.OneBasedBeginPosition}{a.VariantSequence}"));
                    TestContext.WriteLine($"       -> [{j}] Len={ep.Length} AppliedCount={applied.Count} Applied={appliedSummary} Seq='{ep.BaseSequence}'");
                }
            }

            // Expand each protein into its proteoform(s). Because each input protein here has a single variant,
            // each expansion should yield exactly one applied-variant proteoform (no reference-only copy is expected here).
            var proteinsWithAppliedVariants = proteinsWithSeqVars.SelectMany(p => p.GetVariantBioPolymers()).ToList();

            // ——— DIAGNOSTICS: Summarize the flattened list
            int flatBaseOnly = proteinsWithAppliedVariants.Count(x => !x.AppliedSequenceVariations.Any());
            int flatVariantApplied = proteinsWithAppliedVariants.Count(x => x.AppliedSequenceVariations.Any());
            TestContext.WriteLine($"[Flattened] Total={proteinsWithAppliedVariants.Count} BaseOnly={flatBaseOnly} VariantApplied={flatVariantApplied}");
            for (int k = 0; k < proteinsWithAppliedVariants.Count; k++)
            {
                var ep = proteinsWithAppliedVariants[k];
                var applied = ep.AppliedSequenceVariations.ToList();
                string appliedSummary = applied.Count == 0
                    ? "(none)"
                    : string.Join(";", applied.Select(a => $"{a.OriginalSequence}{a.OneBasedBeginPosition}{a.VariantSequence}"));
                TestContext.WriteLine($"   Flat[{k}] From={ep.ConsensusVariant?.Accession} Len={ep.Length} AppliedCount={applied.Count} Applied={appliedSummary}");
            }

            // Invoke expansion again to assert determinism (should match 1:1 by sequence and variation semantics).
            var proteinsWithAppliedVariants2 = proteinsWithSeqVars.SelectMany(p => p.GetVariantBioPolymers()).ToList();

            // Persist the original variant-bearing proteins (not the expansions) and reload to validate round-trip stability.
            string xml = Path.Combine(TestContext.CurrentContext.TestDirectory, "AppliedVariants.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinsWithSeqVars, xml);
            var proteinsWithAppliedVariants3 = ProteinDbLoader.LoadProteinXML(xml, true, DecoyType.None, null, false, null, out var un);

            // ——— DIAGNOSTICS: Round-trip summary and why list[2] (reload) differs
            TestContext.WriteLine($"[Roundtrip] L1={proteinsWithAppliedVariants.Count} L2={proteinsWithAppliedVariants2.Count} L3={proteinsWithAppliedVariants3.Count}");
            int l3BaseOnly = proteinsWithAppliedVariants3.Count(x => !x.AppliedSequenceVariations.Any());
            int l3VariantApplied = proteinsWithAppliedVariants3.Count(x => x.AppliedSequenceVariations.Any());
            TestContext.WriteLine($"[Roundtrip] L3 BaseOnly={l3BaseOnly} VariantApplied={l3VariantApplied}");

            var l3Groups = proteinsWithAppliedVariants3.GroupBy(p => p.ConsensusVariant?.Accession ?? "(null)").ToList();
            foreach (var g in l3Groups)
            {
                TestContext.WriteLine($"[Roundtrip] Group Accession={g.Key} Count={g.Count()} " +
                                      $"RefOnly={g.Count(p => !p.AppliedSequenceVariations.Any())} " +
                                      $"Alt={g.Count(p => p.AppliedSequenceVariations.Any())}");
            }

            for (int idx = 0; idx < proteinsWithAppliedVariants3.Count; idx++)
            {
                var p3 = proteinsWithAppliedVariants3[idx];
                var applied = p3.AppliedSequenceVariations.ToList();
                string kind = applied.Count == 0 ? "REF" : "ALT";
                string appliedSummary = applied.Count == 0
                    ? "(none)"
                    : string.Join(";", applied.Select(a => $"[{a.OneBasedBeginPosition}-{a.OneBasedEndPosition}] {a.OriginalSequence}->{a.VariantSequence} VCF={(a.VariantCallFormatData != null ? "Y" : "N")} GTs={(a.VariantCallFormatData?.Genotypes?.Count ?? 0)}"));
                TestContext.WriteLine($"[Roundtrip] L3[{idx}] {kind} Acc={p3.Accession} From={p3.ConsensusVariant?.Accession} Len={p3.Length} AppliedCount={applied.Count} Applied={appliedSummary}");
            }

            // ——— COORDINATE AUDIT: compare original vs applied positions with computed diff regions
            static (int refStart, int refEnd, int altStart, int altEnd) DiffSpan(string refSeq, string altSeq)
            {
                int p = 0;
                int maxPrefix = Math.Min(refSeq.Length, altSeq.Length);
                while (p < maxPrefix && refSeq[p] == altSeq[p]) p++;

                int s = 0;
                int maxSuffix = Math.Min(refSeq.Length - p, altSeq.Length - p);
                while (s < maxSuffix && refSeq[refSeq.Length - 1 - s] == altSeq[altSeq.Length - 1 - s]) s++;

                int refStart = p + 1; // 1-based
                int refEnd = refSeq.Length - s;
                int altStart = p + 1; // 1-based
                int altEnd = altSeq.Length - s;
                if (refStart > refEnd) { refStart = refEnd = p; } // no diff (shouldn’t happen here)
                if (altStart > altEnd) { altStart = altEnd = p; }
                return (refStart, refEnd, altStart, altEnd);
            }

            for (int k = 0; k < proteinsWithAppliedVariants.Count; k++)
            {
                var ep = proteinsWithAppliedVariants[k];
                var svIn = proteinsWithSeqVars[k].SequenceVariations.Single();
                var svOut = ep.AppliedSequenceVariations.Single();

                string refSeq = ep.ConsensusVariant.BaseSequence;
                string altSeq = ep.BaseSequence;

                var (refStart, refEnd, altStart, altEnd) = DiffSpan(refSeq, altSeq);

                TestContext.WriteLine(
                    $"[Audit k={k}] RefLen={refSeq.Length} AltLen={altSeq.Length} " +
                    $"InputSV=[{svIn.OneBasedBeginPosition}-{svIn.OneBasedEndPosition}] '{svIn.OriginalSequence}'→'{svIn.VariantSequence}' " +
                    $"AppliedSV=[{svOut.OneBasedBeginPosition}-{svOut.OneBasedEndPosition}] '{svOut.OriginalSequence}'→'{svOut.VariantSequence}' " +
                    $"DiffSpan Ref=[{refStart}-{refEnd}] Alt=[{altStart}-{altEnd}]");

                // Show local 1-based windows (±2) around input vs applied begin sites
                int win = 2;
                int inBeg = svIn.OneBasedBeginPosition;
                int apBeg = svOut.OneBasedBeginPosition;

                string refWinIn = refSeq.Substring(Math.Max(0, inBeg - 1 - win), Math.Min(refSeq.Length - Math.Max(0, inBeg - 1 - win), 1 + 2 * win));
                string refWinAp = refSeq.Substring(Math.Max(0, apBeg - 1 - win), Math.Min(refSeq.Length - Math.Max(0, apBeg - 1 - win), 1 + 2 * win));
                string altWinAp = altSeq.Substring(Math.Max(0, apBeg - 1 - win), Math.Min(altSeq.Length - Math.Max(0, apBeg - 1 - win), 1 + 2 * win));

                TestContext.WriteLine($"[Audit k={k}] Around InputBegin={inBeg} RefWin='{refWinIn}'");
                TestContext.WriteLine($"[Audit k={k}] Around AppliedBegin={apBeg} RefWin='{refWinAp}' AltWin='{altWinAp}'");
            }

            // Convenience: aggregate all three lists to run the same assertions over them.
            var listArray = new[] { proteinsWithAppliedVariants, proteinsWithAppliedVariants2, proteinsWithAppliedVariants3 };

            // Sanity: each collection should have exactly the same five proteoforms in the same order
            Assert.AreEqual(5, proteinsWithAppliedVariants.Count);
            Assert.AreEqual(5, proteinsWithAppliedVariants2.Count);
            Assert.AreEqual(5, proteinsWithAppliedVariants3.Count);

            // Validate each expanded collection, index-by-index.
            for (int dbIdx = 0; dbIdx < listArray.Length; dbIdx++)
            {
                // ————————————————————————————————————————————————————————————————————————————————
                // [0] SAV P(4)→V
                // Expect a single residue change at position 4: M P E P T I D E → M P E V T I D E
                Assert.AreEqual("MPEVTIDE", listArray[dbIdx][0].BaseSequence);
                Assert.AreEqual(1, listArray[dbIdx][0].AppliedSequenceVariations.Count());
                Assert.AreEqual(4, listArray[dbIdx][0].AppliedSequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(4, listArray[dbIdx][0].AppliedSequenceVariations.Single().OneBasedEndPosition);
                Assert.AreEqual("P", listArray[dbIdx][0].AppliedSequenceVariations.Single().OriginalSequence);
                Assert.AreEqual("V", listArray[dbIdx][0].AppliedSequenceVariations.Single().VariantSequence);
                Assert.AreEqual("P4V", listArray[dbIdx][0].AppliedSequenceVariations.Single().SimpleString()); // begin-index encoding

                // ————————————————————————————————————————————————————————————————————————————————
                // [1] MNV PT(4-5)→KT
                // Expect a two-residue substitution: PT → KT across indices 4–5
                Assert.AreEqual("MPEKTIDE", listArray[dbIdx][1].BaseSequence);
                Assert.AreEqual(1, listArray[dbIdx][1].AppliedSequenceVariations.Count());
                Assert.AreEqual(4, listArray[dbIdx][1].AppliedSequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(5, listArray[dbIdx][1].AppliedSequenceVariations.Single().OneBasedEndPosition);
                Assert.AreEqual("PT", listArray[dbIdx][1].AppliedSequenceVariations.Single().OriginalSequence);
                Assert.AreEqual("KT", listArray[dbIdx][1].AppliedSequenceVariations.Single().VariantSequence);
                Assert.AreEqual("PT4KT", listArray[dbIdx][1].AppliedSequenceVariations.Single().SimpleString());

                // ————————————————————————————————————————————————————————————————————————————————
                // [2] Insertion P(4)→PPP
                // Expect sequence growth: single P becomes PPP at position 4
                Assert.AreEqual("MPEPPPTIDE", listArray[dbIdx][2].BaseSequence);
                Assert.AreEqual(1, listArray[dbIdx][2].AppliedSequenceVariations.Count());
                Assert.AreEqual(4, listArray[dbIdx][2].AppliedSequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(6, listArray[dbIdx][2].AppliedSequenceVariations.Single().OneBasedEndPosition); // insertion expands the span
                Assert.AreEqual("P", listArray[dbIdx][2].AppliedSequenceVariations.Single().OriginalSequence);
                Assert.AreEqual("PPP", listArray[dbIdx][2].AppliedSequenceVariations.Single().VariantSequence);
                Assert.AreEqual("P4PPP", listArray[dbIdx][2].AppliedSequenceVariations.Single().SimpleString());
                Assert.Greater(listArray[dbIdx][2].Length, proteinsWithSeqVars[0].Length); // length increased

                // ————————————————————————————————————————————————————————————————————————————————
                // [3] Deletion PPP(4-6)→P on base "MPEPPPTIDE"
                // Expect sequence shrinkage: PPP collapses to P (net -2 aa vs index [2])
                Assert.AreEqual("MPEPTIDE", listArray[dbIdx][3].BaseSequence);
                Assert.AreEqual(1, listArray[dbIdx][3].AppliedSequenceVariations.Count());
                Assert.AreEqual(4, listArray[dbIdx][3].AppliedSequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(4, listArray[dbIdx][3].AppliedSequenceVariations.Single().OneBasedEndPosition); // compressed span
                Assert.AreEqual("PPP", listArray[dbIdx][3].AppliedSequenceVariations.Single().OriginalSequence);
                Assert.AreEqual("P", listArray[dbIdx][3].AppliedSequenceVariations.Single().VariantSequence);
                Assert.AreEqual("PPP4P", listArray[dbIdx][3].AppliedSequenceVariations.Single().SimpleString());
                Assert.Less(listArray[dbIdx][3].Length, proteinsWithSeqVars[3].ConsensusVariant.Length); // length decreased vs its own consensus

                // ————————————————————————————————————————————————————————————————————————————————
                // [4] Insertion with variant-scoped PTM at residue 5
                // Same sequence as [2] but ensures the PTM dictionary is present and mapped to the correct 1-based site.
                Assert.AreEqual("MPEPPPTIDE", listArray[dbIdx][4].BaseSequence);
                Assert.AreEqual(1, listArray[dbIdx][4].AppliedSequenceVariations.Count());
                Assert.AreEqual(4, listArray[dbIdx][4].AppliedSequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(6, listArray[dbIdx][4].AppliedSequenceVariations.Single().OneBasedEndPosition);
                Assert.AreEqual("P", listArray[dbIdx][4].AppliedSequenceVariations.Single().OriginalSequence);
                Assert.AreEqual("PPP", listArray[dbIdx][4].AppliedSequenceVariations.Single().VariantSequence);

                // Variant-scoped PTM should project to a valid protein index (and persist across round-trip).
                Assert.AreEqual(1, listArray[dbIdx][4].OneBasedPossibleLocalizedModifications.Count);
                Assert.AreEqual(5, listArray[dbIdx][4].OneBasedPossibleLocalizedModifications.Single().Key);
                Assert.That(listArray[dbIdx][4].OneBasedPossibleLocalizedModifications[5], Is.Not.Null);
                Assert.That(listArray[dbIdx][4].OneBasedPossibleLocalizedModifications[5].Count, Is.GreaterThanOrEqualTo(1));
            }

            // Cross-collection consistency checks (index-by-index, sequences only).
            // Ensures determinism across two in-memory expansions and the serialization round-trip.
            for (int i = 0; i < proteinsWithAppliedVariants.Count; i++)
            {
                Assert.AreEqual(proteinsWithAppliedVariants[i].BaseSequence, proteinsWithAppliedVariants2[i].BaseSequence);
                Assert.AreEqual(proteinsWithAppliedVariants[i].BaseSequence, proteinsWithAppliedVariants3[i].BaseSequence);
            }

            // Additional invariants:
            // - Every proteoform here has exactly one applied variant (by construction).
            // - No applied-variant set is empty (i.e., no reference-only proteoforms slipped in).
            Assert.That(proteinsWithAppliedVariants.All(p => p.AppliedSequenceVariations.Count() == 1));
            Assert.That(proteinsWithAppliedVariants2.All(p => p.AppliedSequenceVariations.Count() == 1));
            Assert.That(proteinsWithAppliedVariants3.All(p => p.AppliedSequenceVariations.Count() == 1));
        }
        [Test]
        public static void AppliedVariants_AsIBioPolymer()
        {
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            Modification mp = new Modification("mod", null, "type", null, motifP, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);

            List<IBioPolymer> proteinsWithSeqVars = new List<IBioPolymer>
            {
                new Protein("MPEPTIDE", "protein1", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "V", "", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein2", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 5, "PT", "KT", "", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein3", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "PPP", "", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPPPTIDE", "protein4", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 6, "PPP", "P", "", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein5",
                    sequenceVariations: new List<SequenceVariation> {
                        new SequenceVariation(4, 4, "P", "PPP","",
                            @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30",
                            new Dictionary<int, List<Modification>> {{ 5, new[] { mp }.ToList() } })
                    }),
            };
            var proteinsWithAppliedVariants = proteinsWithSeqVars.SelectMany(p => p.GetVariantBioPolymers()).ToList();
            var proteinsWithAppliedVariants2 = proteinsWithSeqVars.SelectMany(p => p.GetVariantBioPolymers()).ToList(); // should be stable
            string xml = Path.Combine(TestContext.CurrentContext.TestDirectory, "AppliedVariants.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinsWithSeqVars, xml);
            var proteinsWithAppliedVariants3 = ProteinDbLoader.LoadProteinXML(xml, true, DecoyType.None, null, false, null, out var un);

            var listArray = new List<IBioPolymer>[]
            {
                proteinsWithAppliedVariants,
                proteinsWithAppliedVariants2,
                proteinsWithAppliedVariants3.Cast<IBioPolymer>().ToList()
            };

            for (int dbIdx = 0; dbIdx < listArray.Length; dbIdx++)
            {
                // sequences
                Assert.AreEqual("MPEVTIDE", listArray[dbIdx][0].BaseSequence);
                Assert.AreEqual("MPEKTIDE", listArray[dbIdx][1].BaseSequence);
                Assert.AreEqual("MPEPPPTIDE", listArray[dbIdx][2].BaseSequence);
                Assert.AreEqual("MPEPTIDE", listArray[dbIdx][3].BaseSequence);
                Assert.AreEqual("MPEPPPTIDE", listArray[dbIdx][4].BaseSequence);
                Assert.AreEqual(5, listArray[dbIdx][4].OneBasedPossibleLocalizedModifications.Single().Key);

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
            // Load a protein that contains a stop-gained sequence variation.
            // With default settings (minAlleleDepth default), the loader emits both:
            // - a reference (consensus) protein (no variant applied),
            // - and a variant-applied protein (truncated at the stop site).
            var proteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "StopGained.xml"),
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: null,
                isContaminant: false,
                modTypesToExclude: null,
                unknownModifications: out var unknownModifications);

            // Expect both reference and truncated variant proteins to be present.
            Assert.AreEqual(2, proteins.Count);

            // Reference protein should report exactly one sequence variation in the entry
            // (the stop-gained annotation), even though it is not applied to the consensus.
            Assert.AreEqual(1, proteins[0].SequenceVariations.Count()); // some redundant
            Assert.AreEqual(1, proteins[0].SequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes

            // On the consensus (reference) protein, no variations are applied by default,
            // so the reference sequence remains intact.
            Assert.AreEqual(0, proteins[0].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(0, proteins[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes

            // The second protein should be the variant-applied one (stop-gained),
            // so it should have exactly one applied variation (the truncating stop-gained).
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes

            // The reference protein retains its original full length.
            Assert.AreEqual(191, proteins[0].Length);

            // The reference residue at position 161 is 'Q' (which becomes a stop in the variant).
            Assert.AreEqual('Q', proteins[0][161 - 1]);

            // The variant-applied (stop-gained) protein is truncated to length 160,
            // i.e., everything at and after the stop is removed.
            Assert.AreEqual(161 - 1, proteins[1].Length);

            // Sanity: the reference and the variant-applied protein should differ in length.
            Assert.AreNotEqual(proteins[0].Length, proteins[1].Length);

            // Now re-load the same input but require very high minimum allele depth.
            // This threshold will cause only the deeply-supported stop-gained variant to be returned.
            proteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "StopGained.xml"),
                generateTargets: true,
                decoyType: DecoyType.None,
                allKnownModifications: null,
                isContaminant: false,
                modTypesToExclude: null,
                unknownModifications: out unknownModifications,
                minAlleleDepth: 400);

            // Only the variant-applied protein passes the depth filter, so only one protein remains.
            Assert.AreEqual(1, proteins.Count);

            // That remaining protein must have the stop-gained variation applied.
            Assert.AreEqual(1, proteins[0].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(1, proteins[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes

            // It must still be truncated at the stop site (length 160).
            Assert.AreEqual(161 - 1, proteins[0].Length);
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
            Assert.AreEqual(2, variantProteins.First().ConsensusVariant.SequenceVariations.Count(v => v.VariantCallFormatData.Heterozygous.Any(kv => kv.Value)));

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
            // This test documents the current validation semantics implemented in SequenceVariation.AreValid():
            // - A SequenceVariation is considered valid if and only if:
            //      OneBasedBeginPosition > 0  AND  OneBasedEndPosition >= OneBasedBeginPosition.
            // - The validation does NOT consider the contents of OriginalSequence or VariantSequence.
            //   Therefore:
            //      • Point substitutions (len = 1 → 1) are valid if indices are valid.
            //      • Insertions     (len = 0 → >0) are valid if indices are valid.
            //      • Deletions      (>0 → 0) are valid if indices are valid.
            // - Sequences may be empty strings. The constructor normalizes null strings to "".
            // - Indices are 1-based and inclusive.

            // Local helpers: replicate the library's internal spatial relations for use in tests.
            // These mirror Omics.BioPolymer.SequenceVariation's internal methods to avoid visibility issues across assemblies.
            static bool IntersectsPos(SequenceVariation sv, int pos) =>
                sv.OneBasedBeginPosition <= pos && pos <= sv.OneBasedEndPosition; // same as internal Intersects(int pos)

            static bool IntersectsVar(SequenceVariation a, SequenceVariation b) =>
                b.OneBasedEndPosition >= a.OneBasedBeginPosition && b.OneBasedBeginPosition <= a.OneBasedEndPosition; // same as internal Intersects(SequenceVariation segment)

            static bool IncludesVar(SequenceVariation a, SequenceVariation b) =>
                a.OneBasedBeginPosition <= b.OneBasedBeginPosition && a.OneBasedEndPosition >= b.OneBasedEndPosition; // same as internal Includes(SequenceVariation segment)

            // --- Construct three simple, valid point substitutions (begin == end). ---
            // sv1: at position 10, A → T
            SequenceVariation sv1 = new SequenceVariation(
                oneBasedBeginPosition: 10,    // must be > 0
                oneBasedEndPosition: 10,      // must be >= begin (equal is allowed; this encodes a single position)
                originalSequence: "A",        // original residue (not used by AreValid)
                variantSequence: "T",         // variant residue (not used by AreValid)
                description: "info",          // free text; wrapped in SequenceVariantDescription internally
                oneBasedModifications: null); // null is allowed; ctor normalizes to an empty dictionary

            // sv2: at position 5, G → C
            SequenceVariation sv2 = new SequenceVariation(5, 5, "G", "C", "", "info", null);

            // sv3: at position 8, T → A
            SequenceVariation sv3 = new SequenceVariation(8, 8, "T", "A", "", "info", null);

            // Group the above three variations. All should be valid because indices satisfy: begin > 0 and end >= begin.
            var svList = new List<SequenceVariation> { sv1, sv2, sv3 };

            // Construct a dummy protein that carries these three variations.
            // Note: AreValid() checks only the variant's indices; it does not check that the residues actually match
            // the protein sequence at those indices. This test documents that behavior.
            Protein variantProtein = new Protein(
                sequence: "ACDEFGHIKLMNPQRSTVWY", // 20 aa, indexes 1..20
                accession: "protein1",
                sequenceVariations: svList);

            // Assert: all three SequenceVariations are valid under the current AreValid() rule.
            Assert.IsTrue(variantProtein.SequenceVariations.All(v => v.AreValid()));

            // --- Construct two INVALID variations to document the exact failure modes. ---

            // Invalid because begin <= 0 (begin == 0). End is irrelevant once begin is invalid.
            SequenceVariation svInvalidOneBasedBeginLessThanOne = new SequenceVariation(
                oneBasedBeginPosition: 0,  // invalid (must be > 0)
                oneBasedEndPosition: 10,
                originalSequence: "A",
                variantSequence: "T",
                description: "info",
                oneBasedModifications: null);

            // Invalid because end < begin (5..4 is a negative-length interval).
            SequenceVariation svInvalidOneBasedEndLessThanOneBasedBegin = new SequenceVariation(
                oneBasedBeginPosition: 5,
                oneBasedEndPosition: 4, // invalid (must be >= begin)
                originalSequence: "G",
                variantSequence: "C",
                description: "info",
                oneBasedModifications: null);

            Assert.IsFalse(svInvalidOneBasedBeginLessThanOne.AreValid(), "Begin must be >= 1 for 1-based coordinates.");
            Assert.IsFalse(svInvalidOneBasedEndLessThanOneBasedBegin.AreValid(), "End must be >= begin (inclusive interval).");

            // --- Construct two VALID structural edits to document insertion and deletion semantics. ---

            // Insertion at position 8:
            // - OriginalSequence is empty ("") and VariantSequence has content ("A").
            // - AreValid() only looks at indices, so this is valid.
            SequenceVariation svValidOriginalSequenceIsEmpty = new SequenceVariation(
                oneBasedBeginPosition: 8,
                oneBasedEndPosition: 8, // for an insertion, begin == end represents the insertion site
                originalSequence: "",   // empty original implies insertion
                variantSequence: "A",   // inserted content
                description: "info",
                oneBasedModifications: null);

            // Deletion at position 10:
            // - OriginalSequence has content ("A") and VariantSequence is empty ("").
            // - AreValid() only looks at indices, so this is valid.
            SequenceVariation svValidVariantSequenceLenthIsZero = new SequenceVariation(
                oneBasedBeginPosition: 10,
                oneBasedEndPosition: 10, // single-residue deletion example
                originalSequence: "A",   // removed content
                variantSequence: "",     // empty variant implies deletion
                description: "info",
                oneBasedModifications: null);

            Assert.IsTrue(svValidOriginalSequenceIsEmpty.AreValid(), "Insertions (original == \"\") are valid if indices are valid.");
            Assert.IsTrue(svValidVariantSequenceLenthIsZero.AreValid(), "Deletions (variant == \"\") are valid if indices are valid.");

            // --- Additional explicit documentation asserts (non-functional but clarifying current behavior). ---

            // 1) One-based, inclusive interval rule at boundary:
            // Begin == 1 and End == 1 is valid (first residue).
            var svBoundary = new SequenceVariation(1, 1, "A", "T", "", "info", null);
            Assert.IsTrue(svBoundary.AreValid(), "Begin == 1 and End == 1 is valid (first position).");

            // 2) Constructor null handling:
            // Passing null for sequence strings is normalized to "" (empty) internally.
            var svNulls = new SequenceVariation(3, 3, null, null, "", "info", null);
            Assert.AreEqual("", svNulls.OriginalSequence, "Null originalSequence is normalized to empty string.");
            Assert.AreEqual("", svNulls.VariantSequence, "Null variantSequence is normalized to empty string.");

            // 3) Description wrapping:
            // SequenceVariation stores Description as a SequenceVariantDescription; the original text is accessible via .Description.
            Assert.IsNotNull(sv1.Description, "Description is always constructed.");
            Assert.AreEqual("info", sv1.Description, "Free-text description preserved in SequenceVariantDescription.");

            // 4) SimpleString() formatting:
            // SimpleString() returns: OriginalSequence + OneBasedBeginPosition + VariantSequence (no separators).
            // For a substitution (A@10 -> T) the simple string is "A10T".
            Assert.AreEqual("A10T", sv1.SimpleString(), "Substitution simple-string format: Original + Begin + Variant.");

            // For an insertion at position 8 ("" -> "A"), the simple string becomes "8A"
            // because OriginalSequence is empty.
            Assert.AreEqual("8A", svValidOriginalSequenceIsEmpty.SimpleString(), "Insertion simple-string format: '' + Begin + Variant.");

            // For a deletion at position 10 ("A" -> ""), the simple string becomes "A10"
            // because VariantSequence is empty.
            Assert.AreEqual("A10", svValidVariantSequenceLenthIsZero.SimpleString(), "Deletion simple-string format: Original + Begin + ''.");

            // 5) AreValid() ignores sequence content/lengths:
            //    - The method does not require Original and Variant to be the same length.
            //    - It does not check that Original matches the underlying protein sequence at that location.
            // The following examples demonstrate that length changes do not affect validity when indices are valid.
            var svInsertionLong = new SequenceVariation(7, 7, "", "PPP", "", "info", null);  // insertion of "PPP" at 7
            var svDeletionLong = new SequenceVariation(12, 14, "ABC", "", "", "info", null); // deletion of 3 residues starting at 12
            Assert.IsTrue(svInsertionLong.AreValid(), "Length-increasing edits are valid if indices are valid.");
            Assert.IsTrue(svDeletionLong.AreValid(), "Length-decreasing edits are valid if indices are valid.");

            // 6) Intersects/Includes are not part of AreValid(); they are separate spatial relations.
            //    Because those methods are internal in the library, we mirror their logic here via helpers.
            var p5 = new SequenceVariation(5, 5, "G", "C", "", "info", null);
            var p6 = new SequenceVariation(6, 6, "H", "R", "", "info", null);
            Assert.IsTrue(IntersectsPos(p5, 5), "Intersects(pos): returns true when pos lies within [begin..end].");
            Assert.IsFalse(IntersectsPos(p5, 6), "Intersects(pos): returns false when pos is outside [begin..end].");
            Assert.IsFalse(IntersectsVar(p5, p6), "Intersects(SequenceVariation): disjoint singletons at 5 and 6 do not intersect.");
            Assert.IsTrue(IncludesVar(p5, p5), "Includes(SequenceVariation): a segment includes itself.");

            // Summary: This test intentionally focuses on the exact current behavior of AreValid()
            // and related constructors/utilities, to provide unambiguous documentation of accepted inputs.
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
            // Purpose of this test:
            // - Load a small, known XML (humanGAPDH.xml) that contains:
            //     • a single target (reference) protein entry with two annotated sequence variations
            //     • and instruct the loader to generate reverse decoys alongside the targets.
            // - Validate that variant-expansion during read produces:
            //     • 4 target proteins: reference-only + each single-variant-applied + both-variants-applied
            //     • 4 decoys mirroring the same structure (reference-only + variants) but using decoy accessions
            // - Confirm accessions, variant counts, decoy flags, and basic invariants about applied vs. consensus variants.

            // Database file name (relative to the DatabaseTests directory)
            string databaseName = "humanGAPDH.xml";

            // Load the XML, asking for decoys (__Reverse__) and allowing up to 99 heterozygous combinations (not limiting here).
            // Parameters:
            //   - generateTargets: true        → read target entries
            //   - decoyType: Reverse           → generate reverse-sequence decoys for each target
            //   - allKnownModifications: null  → not needed in this test
            //   - isContaminant: false         → not contaminant
            //   - modTypesToExclude: null      → include all
            //   - unknownModifications: out _  → capture but not used here
            //   - minAlleleDepth: 1            → include all annotated variants (no depth filter)
            //   - maxHeterozygousVariants: 99 → ensure combinatorics not curtailed
            var proteins = ProteinDbLoader.LoadProteinXML(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", databaseName),
                generateTargets: true,
                decoyType: DecoyType.Reverse,
                allKnownModifications: null,
                isContaminant: false,
                modTypesToExclude: null,
                unknownModifications: out var unknownModifications,
                minAlleleDepth: 1,
                maxHeterozygousVariants: 99);

            // Expect 8 proteins total:
            //   - 4 targets (reference + 2 single-variant proteoforms + 1 double-variant proteoform)
            //   - 4 decoys mirroring the above set
            Assert.AreEqual(8, proteins.Count, "Expected: 4 target + 4 decoy entries.");

            // The first target entry is the consensus (reference-only) protein. It should carry two annotated sequence variations
            // in its entry (not applied to the consensus), i.e., the XML encodes 2 unique variants for GAPDH.
            Assert.AreEqual(2, proteins[0].SequenceVariations.Count(), "Reference entry should list exactly two annotated sequence variations.");

            // Accessions for targets are expected to be:
            //   [0] reference-only (no suffix)
            //   [1] A22G applied
            //   [2] K251N applied
            //   [3] K251N and A22G applied (order in name is stable and deterministic)
            Assert.That("P04406", Is.EqualTo(proteins[0].Accession));
            Assert.That("P04406_A22G", Is.EqualTo(proteins[1].Accession));
            Assert.That("P04406_K251N", Is.EqualTo(proteins[2].Accession));
            Assert.That("P04406_K251N_A22G", Is.EqualTo(proteins[3].Accession));

            // Accessions for decoys mirror the same shape, but prefixed with "DECOY_".
            // Note that the positions encoded in the decoy names reflect reverse-decoy coordinate transforms:
            //   A22G (target) becomes A315G (decoy), and K251N (target) becomes K86N (decoy).
            Assert.That("DECOY_P04406", Is.EqualTo(proteins[4].Accession));
            Assert.That("DECOY_P04406_A315G", Is.EqualTo(proteins[5].Accession));
            Assert.That("DECOY_P04406_K86N", Is.EqualTo(proteins[6].Accession));
            Assert.That("DECOY_P04406_K86N_A315G", Is.EqualTo(proteins[7].Accession));

            // ——— Additional explicit validations to capture the entire behavior ———

            // Targets first, then decoys (ordering guarantee used throughout this test).
            Assert.IsTrue(proteins.Take(4).All(p => !p.IsDecoy), "First 4 entries must be targets.");
            Assert.IsTrue(proteins.Skip(4).All(p => p.IsDecoy), "Last 4 entries must be decoys.");

            // Reference entries (target[0] and decoy[4]) should not have variants applied to the BaseSequence.
            Assert.AreEqual(0, proteins[0].AppliedSequenceVariations.Count(), "Target reference entry must have no applied variants.");
            Assert.AreEqual(0, proteins[4].AppliedSequenceVariations.Count(), "Decoy reference entry must have no applied variants.");

            // The next two target entries each carry exactly one applied variant.
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Count(), "Target A22G should have exactly one applied variant.");
            Assert.AreEqual(1, proteins[2].AppliedSequenceVariations.Count(), "Target K251N should have exactly one applied variant.");

            // The fourth target entry carries both variants applied.
            Assert.AreEqual(2, proteins[3].AppliedSequenceVariations.Count(), "Target double-variant should have exactly two applied variants.");

            // Decoys mirror the above shape: [5] one variant, [6] one variant, [7] two variants.
            Assert.AreEqual(1, proteins[5].AppliedSequenceVariations.Count(), "Decoy A315G should have exactly one applied variant.");
            Assert.AreEqual(1, proteins[6].AppliedSequenceVariations.Count(), "Decoy K86N should have exactly one applied variant.");
            Assert.AreEqual(2, proteins[7].AppliedSequenceVariations.Count(), "Decoy double-variant should have exactly two applied variants.");

            // Reference accessions act as consensus for their variant-applied proteoforms:
            // targets[1..3] should point back to targets[0], decoys[5..7] back to decoys[4].
            Assert.AreEqual(proteins[0].Accession, proteins[1].ConsensusVariant.Accession, "Variant targets must share the same consensus accession.");
            Assert.AreEqual(proteins[0].Accession, proteins[2].ConsensusVariant.Accession);
            Assert.AreEqual(proteins[0].Accession, proteins[3].ConsensusVariant.Accession);
            Assert.AreEqual(proteins[4].Accession, proteins[5].ConsensusVariant.Accession, "Variant decoys must share the same consensus accession.");
            Assert.AreEqual(proteins[4].Accession, proteins[6].ConsensusVariant.Accession);
            Assert.AreEqual(proteins[4].Accession, proteins[7].ConsensusVariant.Accession);

            // By construction:
            // - Reference entries (no applied variants) keep BaseSequence identical to their own ConsensusVariant.BaseSequence.
            // - Variant entries (>=1 applied) must differ in BaseSequence from their ConsensusVariant.BaseSequence.
            static bool SequenceEqualsConsensus(Protein p) => p.BaseSequence == p.ConsensusVariant.BaseSequence;
            Assert.IsTrue(SequenceEqualsConsensus(proteins[0]), "Target reference: BaseSequence must equal ConsensusVariant.BaseSequence.");
            Assert.IsTrue(SequenceEqualsConsensus(proteins[4]), "Decoy reference: BaseSequence must equal ConsensusVariant.BaseSequence.");
            Assert.IsTrue(proteins[1].BaseSequence != proteins[1].ConsensusVariant.BaseSequence, "A22G should differ from consensus.");
            Assert.IsTrue(proteins[2].BaseSequence != proteins[2].ConsensusVariant.BaseSequence, "K251N should differ from consensus.");
            Assert.IsTrue(proteins[3].BaseSequence != proteins[3].ConsensusVariant.BaseSequence, "Double-variant should differ from consensus.");
            Assert.IsTrue(proteins[5].BaseSequence != proteins[5].ConsensusVariant.BaseSequence, "Decoy A315G should differ from consensus.");
            Assert.IsTrue(proteins[6].BaseSequence != proteins[6].ConsensusVariant.BaseSequence, "Decoy K86N should differ from consensus.");
            Assert.IsTrue(proteins[7].BaseSequence != proteins[7].ConsensusVariant.BaseSequence, "Decoy double-variant should differ from consensus.");

            // Validate the simple-string encodings of the applied variants for targets:
            //   A22G   → "A22G"
            //   K251N  → "K251N"
            // The double-variant entry should contain both of those simple strings (order-insensitive check).
            Assert.That(proteins[1].AppliedSequenceVariations.Single().SimpleString(), Is.EqualTo("A22G"));
            Assert.That(proteins[2].AppliedSequenceVariations.Single().SimpleString(), Is.EqualTo("K251N"));
            var targetDouble = proteins[3].AppliedSequenceVariations.Select(v => v.SimpleString()).OrderBy(s => s).ToArray();
            Assert.That(targetDouble, Is.EqualTo(new[] { "A22G", "K251N" }.OrderBy(s => s).ToArray()));

            // Validate the simple-string encodings of applied variants for decoys:
            //   A22G(target) → A315G(decoy)
            //   K251N(target) → K86N(decoy)
            // Again, confirm both are present on the double-variant decoy entry (order-insensitive).
            Assert.That(proteins[5].AppliedSequenceVariations.Single().SimpleString(), Is.EqualTo("A315G"));
            Assert.That(proteins[6].AppliedSequenceVariations.Single().SimpleString(), Is.EqualTo("K86N"));
            var decoyDouble = proteins[7].AppliedSequenceVariations.Select(v => v.SimpleString()).OrderBy(s => s).ToArray();
            Assert.That(decoyDouble, Is.EqualTo(new[] { "A315G", "K86N" }.OrderBy(s => s).ToArray()));

            // Optional diagnostics to aid future debugging and to document the read-time expansion.
            TestContext.WriteLine($"Loaded {proteins.Count} proteins from {databaseName}.");
            for (int i = 0; i < proteins.Count; i++)
            {
                var p = proteins[i];
                string kind = p.IsDecoy ? "DECOY" : "TARGET";
                string appliedStr = p.AppliedSequenceVariations.Any()
                    ? string.Join(",", p.AppliedSequenceVariations.Select(v => v.SimpleString()))
                    : "(none)";
                TestContext.WriteLine($"[{i}] {kind} Acc={p.Accession} Consensus={p.ConsensusVariant?.Accession} AppliedCount={p.AppliedSequenceVariations.Count()} Applied={appliedStr}");
            }
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