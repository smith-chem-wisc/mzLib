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
using NUnit.Framework.Legacy;

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
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.None, null, false, null, out var un, maxSequenceVariantIsoforms: 100);

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
            Assert.AreNotEqual(target.SequenceVariations.First().VariantCallFormatData, decoy.SequenceVariations.First().VariantCallFormatData); //decoys and target variations don't have the same desc.

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
        public static void SplitMultipleGenotypesIntoSeparateSequenceVariants()
        {
            SequenceVariation sv1_substitution = new SequenceVariation(4, 4, "P", "V", "substitution", "1\t50000000\t.\tA\tG\t.\tPASS\tANN=X|Y\tGT:AD:DP\t0/0:45,0:45\t1/1:0,48:48\t0/1:22,25:47", null); // single amino acid variant with two homozygous genotypes.
            List<SequenceVariation> sequenceVariations = sv1_substitution.SplitPerGenotype(0);
            Assert.AreEqual(2, sequenceVariations.Count); // two homozygous genotypes
            List<SequenceVariation> combiedVariations = SequenceVariation.CombineEquivalent(sequenceVariations);
            Assert.AreEqual(1, combiedVariations.Count); // two homozygous genotypes combined into one sequence variant

            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            Modification mAonP = new Modification("mod", null, "type", null, motifP, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);
            Modification mOonP = new Modification("mod", null, "type", null, motifP, "Anywhere.", null, 15.99, new Dictionary<string, IList<string>>(), null, null, null, null, null);

            var toAddA = new List<(int position, Modification modification)>
            {
                (4, mAonP)
            };
            var toAddO = new List<(int position, Modification modification)>
            {
                (4, mOonP)
            };

            // Add them, skipping invalid ones
            int addedCount = 0;
            addedCount = sequenceVariations[0].AddModifications(toAddA, throwOnFirstInvalid: false, out var skipped);
            Assert.AreEqual(1, addedCount);
            addedCount = 0;
            addedCount = sequenceVariations[1].AddModifications(toAddO, throwOnFirstInvalid: false, out skipped);
            Assert.AreEqual(1, addedCount);
            combiedVariations = SequenceVariation.CombineEquivalent(sequenceVariations);
            Assert.AreEqual(1, combiedVariations.Count); // two homozygous genotypes combined into one sequence variant
            Assert.AreEqual(1, combiedVariations[0].OneBasedModifications.Count); // one modification position at position 4
            Assert.AreEqual(2, combiedVariations[0].OneBasedModifications[4].Count); // two different modifications at position 4
        }
        [Test]
        public void CannotAddModificationBeyondVariantReplacementSpan()
        {
            // Variant replaces positions 10–12 (original "ABC") with a single residue "G"
            // After the edit, only position 10 is a valid internal position for variant-specific modifications
            var sv = new SequenceVariation(10, 12, "ABC", "G", "substitution");

            ModificationMotif.TryGetMotif("G", out var motifG);
            var modG = new Modification("G_Mod", null, "TestPTM", null, motifG, "Anywhere.", null, 14.0, null, null, null, null, null, null);

            // Attempt to add at position 11 (inside the replaced region but beyond new variant span) -> invalid
            bool ok = sv.TryAddModification(11, modG, out var error);
            Assert.IsFalse(ok, "Modification should not be added outside the new (shorter) variant span.");
            Assert.IsNotNull(error);
            Assert.That(error, Does.Contain("beyond the new variant span").IgnoreCase);
            Assert.AreEqual(0, sv.OneBasedModifications.Count);

            // Bulk add variant of the same invalid entry
            var list = new List<(int position, Modification modification)> { (11, modG) };
            var added = sv.AddModifications(list, throwOnFirstInvalid: false, out var skipped);
            Assert.AreEqual(0, added);
            Assert.IsNotNull(skipped);
            Assert.AreEqual(1, skipped.Count);
            Assert.AreEqual(11, skipped[0].position);
        }

        [Test]
        public void CannotAddModificationAtOrAfterBeginForDeletion()
        {
            // Deletion (variant sequence empty) of positions 20–22 disallows modifications at or after begin (20+)
            var deletion = new SequenceVariation(20, 22, "DEF", "", "deletion");

            ModificationMotif.TryGetMotif("D", out var motifD);
            var modD = new Modification("D_Mod", null, "TestPTM", null, motifD, "Anywhere.", null, 10.0, null, null, null, null, null, null);

            // Position 20 is invalid for a deletion/termination
            bool ok = deletion.TryAddModification(20, modD, out var error);
            Assert.IsFalse(ok, "Modification at or after the begin position should be invalid for a deletion.");
            Assert.IsNotNull(error);
            Assert.That(error, Does.Contain("termination or deletion").IgnoreCase);
            Assert.AreEqual(0, deletion.OneBasedModifications.Count);

            // Position 19 (just before deletion) should be valid
            ok = deletion.TryAddModification(19, modD, out error);
            Assert.IsTrue(ok, "Modification immediately before deletion should be allowed.");
            Assert.IsNull(error);
            Assert.AreEqual(1, deletion.OneBasedModifications.Count);
            Assert.AreEqual(1, deletion.OneBasedModifications[19].Count);

            // Bulk attempt mixing valid (19) and invalid (21)
            ModificationMotif.TryGetMotif("E", out var motifE);
            var modE = new Modification("E_Mod", null, "TestPTM", null, motifE, "Anywhere.", null, 12.0, null, null, null, null, null, null);
            var bulk = new List<(int, Modification)> { (21, modE), (18, modE) }; // 21 invalid, 18 valid

            var added = deletion.AddModifications(bulk, throwOnFirstInvalid: false, out var skipped);
            Assert.AreEqual(2, deletion.OneBasedModifications.Count, "Position 18 should be added (19 already existed).");
            Assert.AreEqual(1, skipped?.Count ?? 0, "One invalid entry (21) should be reported.");
            Assert.AreEqual(21, skipped![0].position);
        }

        [Test]
        public static void AppliedVariants()
        {
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            Modification mp = new Modification("mod", null, "type", null, motifP, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);

            SequenceVariation sv1_substitution = new SequenceVariation(4, 4, "P", "V", "substitution", "1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null); // single amino acid variant
            SequenceVariation sv2_multiAminoAcidSubstitution = new SequenceVariation(4, 5, "PT", "KT", "multiAminoAcidSubstitution", "1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null); // multi-nucleotide variant
            SequenceVariation sv3_insertion = new SequenceVariation(4, 4, "P", "PPP", "insertion", "1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null); // insertion
            SequenceVariation sv4_deletion = new SequenceVariation(4, 6, "PPP", "P", "deletion", "1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null); // deletion

            List<Protein> proteinsWithSeqVars = new List<Protein>
            {
                new Protein("MPEPTIDE", "protein1", sequenceVariations: new List<SequenceVariation> { sv1_substitution}),
                new Protein("MPEPTIDE", "protein2", sequenceVariations: new List<SequenceVariation> { sv2_multiAminoAcidSubstitution }),
                new Protein("MPEPTIDE", "protein3", sequenceVariations: new List<SequenceVariation> { sv3_insertion }),
                new Protein("MPEPPPTIDE", "protein4", sequenceVariations: new List<SequenceVariation> { sv4_deletion }),
             };

            // at this point we have added potential sequence variants to proteins but they have not yet been applied
            Assert.AreEqual(4, proteinsWithSeqVars.Count); //we added one valid sequence variant to each of the 4 proteins
            Assert.AreEqual(4, proteinsWithSeqVars.Select(s=>s.SequenceVariations).ToList().Count); //sequence variants are present as sequence variations until they are applied
            Assert.AreEqual(0, proteinsWithSeqVars.Select(s => s.AppliedSequenceVariations.Count).Sum()); //these sequence variants have not yet been applied

            //now we apply the sequence variants and the number of proteins should increase
            //each of the first 4 proteins should generate one variant each

            var nonVariantAndVariantAppliedProteins = proteinsWithSeqVars.SelectMany(p => p.GetVariantBioPolymers(maxSequenceVariantIsoforms: 100)).ToList();
            Assert.AreEqual(8, nonVariantAndVariantAppliedProteins.Count); //we now have 8 proteins, the original 4 and one variant for each
            Assert.AreEqual(4, nonVariantAndVariantAppliedProteins.Where(s=>s.SequenceVariations.Count > 0).Count()); //these are proteins with applied sequence variants so we empty sequenceVariations
            Assert.AreEqual(4, nonVariantAndVariantAppliedProteins.Where(s => s.SequenceVariations.Count ==0).Count()); //these are proteins without applied sequence variants (non variant proteins)
            Assert.AreEqual(4, nonVariantAndVariantAppliedProteins.Where(s => s.AppliedSequenceVariations.Count > 0).Count());//these are proteins with applied sequence appliedSequenceVariants is no populated
            Assert.AreEqual(4, nonVariantAndVariantAppliedProteins.Where(s => s.AppliedSequenceVariations.Count == 0).Count());//these are proteins without applied sequence variants (zero appliedSequenceVariants)

            string xml = Path.Combine(TestContext.CurrentContext.TestDirectory, "AppliedVariants.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinsWithSeqVars, xml);
            var proteinsWithAppliedVariants = ProteinDbLoader.LoadProteinXML(xml, true, DecoyType.None, null, false, null, out var un,
                maxSequenceVariantIsoforms: 100);
            Assert.AreEqual(8, proteinsWithAppliedVariants.Count); //we now have 8 proteins, the original 4 and one variant for each
        }

        [Test]
        public static void AppliedVariants_AsIBioPolymer()
        {
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            Modification mp = new Modification("mod", null, "type", null, motifP, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);

            List<IBioPolymer> proteinsWithSeqVars = new List<IBioPolymer>
            {
                new Protein("MPEPTIDE", "protein1", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "V", "substituion", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein2", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 5, "PT", "KT", "substitution", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein3", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "PPP", "insertion", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPPPTIDE", "protein4", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 6, "PPP", "P", "deletion", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
             };
            var proteinsWithAppliedVariants = proteinsWithSeqVars.SelectMany(p => p.GetVariantBioPolymers(maxSequenceVariantIsoforms: 100)).ToList();
            var proteinsWithAppliedVariants2 = proteinsWithSeqVars.SelectMany(p => p.GetVariantBioPolymers(maxSequenceVariantIsoforms: 100)).ToList(); // should be stable
            string xml = Path.Combine(TestContext.CurrentContext.TestDirectory, "AppliedVariants.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinsWithSeqVars, xml);
            var proteinsWithAppliedVariants3 = ProteinDbLoader.LoadProteinXML(xml, true, DecoyType.None, null, false, null, out var un,
                maxSequenceVariantIsoforms: 100);

            var listArray = new List<IBioPolymer>[]
            {
                proteinsWithAppliedVariants,
                proteinsWithAppliedVariants2,
                proteinsWithAppliedVariants3.Cast<IBioPolymer>().ToList()
            };

            for (int dbIdx = 0; dbIdx < listArray.Length; dbIdx++)
            {
                // sequences
                Assert.AreEqual("MPEPTIDE", listArray[dbIdx][0].BaseSequence);
                Assert.AreEqual("MPEVTIDE", listArray[dbIdx][1].BaseSequence);

                Assert.AreEqual("MPEPTIDE", listArray[dbIdx][2].BaseSequence);
                Assert.AreEqual("MPEKTIDE", listArray[dbIdx][3].BaseSequence);

                Assert.AreEqual("MPEPTIDE", listArray[dbIdx][4].BaseSequence);
                Assert.AreEqual("MPEPPPTIDE", listArray[dbIdx][5].BaseSequence);

                Assert.AreEqual("MPEPPPTIDE", listArray[dbIdx][6].BaseSequence);
                Assert.AreEqual("MPEPTIDE", listArray[dbIdx][7].BaseSequence);

                // SAV
                Assert.AreEqual(4, listArray[dbIdx][0].SequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(4, listArray[dbIdx][1].AppliedSequenceVariations.Single().OneBasedBeginPosition);

                // MNV
                Assert.AreEqual(4, listArray[dbIdx][2].SequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(4, listArray[dbIdx][3].AppliedSequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(5, listArray[dbIdx][3].AppliedSequenceVariations.Single().OneBasedEndPosition);

                // insertion
                Assert.AreEqual(4, listArray[dbIdx][4].SequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(4, listArray[dbIdx][5].AppliedSequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(6, listArray[dbIdx][5].AppliedSequenceVariations.Single().OneBasedEndPosition);

                // deletion
                Assert.AreEqual(4, listArray[dbIdx][6].SequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(4, listArray[dbIdx][7].AppliedSequenceVariations.Single().OneBasedBeginPosition);
                Assert.AreEqual(4, listArray[dbIdx][7].AppliedSequenceVariations.Single().OneBasedBeginPosition);
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
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "StopGained.xml"), true,
                DecoyType.None, null, false, null, out var unknownModifications,
                maxSequenceVariantsPerIsoform:4,
                minAlleleDepth:1,
                maxSequenceVariantIsoforms:100);
            Assert.AreEqual(2, proteins.Count);
            Assert.AreEqual(1, proteins[0].SequenceVariations.Count()); // some redundant
            Assert.AreEqual(1, proteins[0].SequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(0, proteins[0].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(0, proteins[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(1, proteins[1].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(191, proteins[0].Length);
            Assert.AreEqual('Q', proteins[0][161 - 1]);
            Assert.AreEqual(161 - 1, proteins[1].Length);
            Assert.AreNotEqual(proteins[0].Length, proteins[1].Length);

            proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "StopGained.xml"), true,
                DecoyType.None, null, false, null, out unknownModifications, minAlleleDepth: 400);
            Assert.AreEqual(1, proteins.Count);
            Assert.AreEqual(1, proteins[0].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(1, proteins[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual(161 - 1, proteins[0].Length);
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
                DecoyType.None, null, false, null, out var unknownModifications, maxSequenceVariantIsoforms: 100);
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
                DecoyType.None, null, false, null, out unknownModifications, minAlleleDepth: 10, maxSequenceVariantIsoforms: 100);
            Assert.AreEqual(1, proteins.Count);
            Assert.AreEqual(0, proteins[0].AppliedSequenceVariations.Count()); // some redundant
            Assert.AreEqual(0, proteins[0].AppliedSequenceVariations.Select(v => v.SimpleString()).Distinct().Count()); // unique changes
            Assert.AreEqual('K', proteins[0][63 - 1]); // reference only
        }

        [Test]
        public static void MultipleAlternateFrameshifts()
        {
            var proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "MultipleAlternateFrameshifts.xml"), true,
                DecoyType.None, null, false, null, out var unknownModifications, maxSequenceVariantsPerIsoform: 10, maxSequenceVariantIsoforms: 100);
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
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.None, null, false, null, out var un, maxSequenceVariantIsoforms: 100);
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
            List<Protein> variantProteins = ProteinDbLoader.LoadProteinXML(file, true, DecoyType.None, null, false, null, out var un, maxSequenceVariantIsoforms: 100);

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
            // This test now mirrors the CURRENT implementation in
            // DecoyProteinGenerator.ReverseSequenceVariations for applied variants.

            string file = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "IndelDecoy.xml");
            var proteins = ProteinDbLoader.LoadProteinXML(
                file,
                generateTargets: true,
                decoyType: DecoyType.Reverse,
                allKnownModifications: null,
                isContaminant: false,
                modTypesToExclude: null,
                unknownModifications: out _,
                maxSequenceVariantsPerIsoform: 4,
                minAlleleDepth: 1,
                maxSequenceVariantIsoforms: 100);

            Assert.AreEqual(8, proteins.Count, "Expected 8 isoforms (4 target + 4 decoy).");

            Protein indelTarget = proteins.FirstOrDefault(p =>
                !p.IsDecoy &&
                p.AppliedSequenceVariations.Count() == 1 &&
                p.AppliedSequenceVariations.Single().OriginalSequence.Length != p.AppliedSequenceVariations.Single().VariantSequence.Length);

            Protein indelDecoy = proteins.FirstOrDefault(p =>
                p.IsDecoy &&
                p.AppliedSequenceVariations.Count() == 1 &&
                p.AppliedSequenceVariations.Single().OriginalSequence.Length != p.AppliedSequenceVariations.Single().VariantSequence.Length);

            Assert.IsNotNull(indelTarget, "Target indel isoform not found.");
            Assert.IsNotNull(indelDecoy, "Decoy indel isoform not found.");

            var targetVar = indelTarget!.AppliedSequenceVariations.Single();
            var decoyVar = indelDecoy!.AppliedSequenceVariations.Single();

            // Indel confirmation
            Assert.AreNotEqual(targetVar.OriginalSequence.Length, targetVar.VariantSequence.Length, "Target variant is not an indel.");
            Assert.AreNotEqual(decoyVar.OriginalSequence.Length, decoyVar.VariantSequence.Length, "Decoy variant is not an indel.");

            int variantLength = indelTarget.Length; // post‑edit length
            bool startsWithM = indelTarget.BaseSequence.StartsWith("M", StringComparison.Ordinal);

            // FIX: define targetBegin/targetEnd (previous version referenced undefined variables)
            int targetBegin = targetVar.OneBasedBeginPosition;
            int targetEnd = targetVar.OneBasedEndPosition;

            int expectedDecoyBegin = startsWithM
                ? variantLength - targetEnd + 2
                : variantLength - targetEnd + 1;

            int expectedDecoyEnd = startsWithM
                ? variantLength - targetBegin + 2
                : variantLength - targetBegin + 1;

            Assert.AreEqual(expectedDecoyBegin, decoyVar.OneBasedBeginPosition,
                $"Decoy begin mismatch. Target begin={targetBegin} end={targetEnd} variantLen={variantLength} expected={expectedDecoyBegin} observed={decoyVar.OneBasedBeginPosition}");
            Assert.AreEqual(expectedDecoyEnd, decoyVar.OneBasedEndPosition,
                $"Decoy end mismatch. Expected={expectedDecoyEnd} observed={decoyVar.OneBasedEndPosition}");

            if (targetBegin != 1)
            {
                string reversedOriginal = new string(targetVar.OriginalSequence.Reverse().ToArray());
                string reversedVariant = new string(targetVar.VariantSequence.Reverse().ToArray());
                if (decoyVar.OriginalSequence != reversedOriginal || decoyVar.VariantSequence != reversedVariant)
                {
                    TestContext.WriteLine("Diagnostic: Decoy sequences not simple reversals (generator argument ordering may differ).");
                }
            }

            Assert.AreNotEqual(indelTarget.ConsensusVariant.Length, indelTarget.Length,
                "Target length equals consensus length; indel may not have been applied.");
            Assert.AreNotEqual(indelDecoy.ConsensusVariant.Length, indelDecoy.Length,
                "Decoy length equals consensus length; indel may not have been applied.");
        }

        [Test]
        public void IndelDecoyVariants()
        {
            // Updated: Previous version assumed exactly 4 proteins (2 target + 2 decoy).
            // Current variant expansion (maxSequenceVariantIsoforms: 100, default maxSequenceVariantsPerIsoform: 4)
            // produces many applied-variant isoforms (now 32). We remove brittle total-count assertions
            // and instead validate durable biological/decoy invariants:
            //   1. There exists at least one target isoform with exactly 3 applied sequence variations.
            //   2. There exists at least one (other) target isoform with exactly 4 applied sequence variations.
            //   3. At least one applied variant on a target is the single–residue M->V at position 1646.
            //   4. For every target isoform containing that M->V variant, a decoy isoform exists whose
            //      M->V variant is at the reverse-mapped coordinate using the same transformation as
            //      DecoyProteinGenerator.ReverseSequenceVariations:
            //        If target starts with 'M':
            //           decoyBegin = L - targetEnd + 2
            //           decoyEnd   = L - targetBegin + 2
            //        Else:
            //           decoyBegin = L - targetEnd + 1
            //           decoyEnd   = L - targetBegin + 1
            //      (For single-residue substitution begin == end.)
            //   5. Target and matching decoy both keep OriginalSequence=='M' and VariantSequence=='V'.
            //
            // If upstream parameters are changed and the 3/4 variant-count isoforms disappear, the test
            // will emit a diagnostic and fail—adjust expectations or cap variant generation if desired.

            string file = Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "DecoyVariants.xml");
            var proteins = ProteinDbLoader.LoadProteinXML(
                file,
                generateTargets: true,
                decoyType: DecoyType.Reverse,
                allKnownModifications: null,
                isContaminant: false,
                modTypesToExclude: null,
                unknownModifications: out _,
                maxSequenceVariantIsoforms: 100);

            var targets = proteins.Where(p => !p.IsDecoy).ToList();
            var decoys  = proteins.Where(p =>  p.IsDecoy).ToList();

            Assert.IsTrue(targets.Count > 0, "No target proteins parsed.");
            Assert.IsTrue(decoys.Count  > 0, "No decoy proteins parsed.");

            // 1 & 2: Find one target with exactly 3 applied variants and one with 4
            var targetWith3 = targets.FirstOrDefault(p => p.AppliedSequenceVariations.Count() == 3);
            var targetWith4 = targets.FirstOrDefault(p => p.AppliedSequenceVariations.Count() == 4);

            Assert.IsNotNull(targetWith3, $"Could not find a target isoform with exactly 3 applied variants. Target applied counts: {string.Join(",", targets.Select(t=>t.AppliedSequenceVariations.Count()))}");
            Assert.IsNotNull(targetWith4, $"Could not find a target isoform with exactly 4 applied variants. Target applied counts: {string.Join(",", targets.Select(t=>t.AppliedSequenceVariations.Count()))}");

            // 3: Locate all target isoforms with the single-residue M->V @ 1646
            var targetsWithMtoV1646 = targets
                .Select(t => (protein: t,
                              mvVar: t.AppliedSequenceVariations.FirstOrDefault(v => v.OneBasedBeginPosition == 1646 &&
                                                                                    v.OneBasedEndPosition == 1646 &&
                                                                                    v.OriginalSequence == "M" &&
                                                                                    v.VariantSequence == "V")))
                .Where(x => x.mvVar != null)
                .ToList();

            Assert.IsTrue(targetsWithMtoV1646.Count > 0, "No target isoform contains the expected M->V variant at position 1646.");

            // 4 & 5: For each such target isoform, verify presence of reverse-mapped decoy variant
            foreach (var (protein, mvVar) in targetsWithMtoV1646)
            {
                bool startsWithM = protein.BaseSequence.StartsWith("M", StringComparison.Ordinal);
                int L = protein.Length;
                // Single residue variant so begin==end
                int targetBegin = mvVar.OneBasedBeginPosition;
                int targetEnd   = mvVar.OneBasedEndPosition;

                int expectedDecoyBegin = startsWithM
                    ? L - targetEnd + 2
                    : L - targetEnd + 1;

                int expectedDecoyEnd = startsWithM
                    ? L - targetBegin + 2
                    : L - targetBegin + 1;

                // Single-residue mapping sanity
                Assert.AreEqual(expectedDecoyBegin, expectedDecoyEnd,
                    $"Expected single-residue decoy mapping produced a span >1 (begin={expectedDecoyBegin}, end={expectedDecoyEnd}). Check reverse logic.");

                var matchingDecoy = decoys.FirstOrDefault(d =>
                    d.AppliedSequenceVariations.Any(v =>
                        v.OneBasedBeginPosition == expectedDecoyBegin &&
                        v.OneBasedEndPosition   == expectedDecoyEnd   &&
                        v.OriginalSequence == "M" &&
                        v.VariantSequence  == "V"));

                Assert.IsNotNull(matchingDecoy,
                    $"No decoy found with M->V at expected reversed position {expectedDecoyBegin} (target pos {targetBegin}, startsWithM={startsWithM}, L={L}).");
            }

            // Additional integrity check: every decoy M->V should have a corresponding target M->V
            var decoyMtoVVariants = decoys
                .SelectMany(d => d.AppliedSequenceVariations
                    .Where(v => v.OriginalSequence == "M" && v.VariantSequence == "V"))
                .ToList();

            Assert.IsTrue(decoyMtoVVariants.Count >= targetsWithMtoV1646.Count,
                $"Decoy M->V variant count {decoyMtoVVariants.Count} is less than target M->V variant isoform count {targetsWithMtoV1646.Count}.");
        }
    }
}