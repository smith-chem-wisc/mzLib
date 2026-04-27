using Chemistry;
using NUnit.Framework;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using CollectionAssert = NUnit.Framework.Legacy.CollectionAssert;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MzLibUtil;
using Omics;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SequenceConversion;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;
using Stopwatch = System.Diagnostics.Stopwatch;
using Transcriptomics;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class PeptideWithSetModsFullSequenceTests
    {
        private static Stopwatch Stopwatch { get; set; }

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

        /// <summary>
        /// CRITICAL: Tests parsing of complex modification strings with nested brackets.
        /// Modifications like "Cation:Fe[III]" contain brackets that could break parsing.
        /// Essential for correctly reading user-specified and database modifications.
        /// </summary>
        [Test]
        public static void TestHardToParseModifiedSequence()
        {
            string fullSequence = "PE[Metal:Cation:Fe[III] on X]PTIDE";

            ModificationMotif.TryGetMotif("X", out var motif);

            Modification mod = new Modification(_originalId: "Cation:Fe[III]", _modificationType: "Metal",
                _monoisotopicMass: 1, _locationRestriction: "Anywhere.", _target: motif);

            Dictionary<string, Modification> mods = new Dictionary<string, Modification> { { "Cation:Fe[III] on X", mod } };

            PeptideWithSetModifications pep = new PeptideWithSetModifications(fullSequence, mods);

            Assert.That(pep.AllModsOneIsNterminus.Count == 1);
            var annotatedMod = pep.AllModsOneIsNterminus.First();
            Assert.That(annotatedMod.Key == 3);
            Assert.That(annotatedMod.Value.IdWithMotif == "Cation:Fe[III] on X");
            Assert.That(annotatedMod.Value.OriginalId == "Cation:Fe[III]");
            Assert.That(annotatedMod.Value.ModificationType == "Metal");

            fullSequence = "[Metal:Cation:Fe[III] on X]PE[Metal:Cation:Fe[III] on X]PTIDE[Metal:Cation:Fe[III] on X]";
            pep = new PeptideWithSetModifications(fullSequence, mods);
            Assert.That(pep.AllModsOneIsNterminus.Count == 3);
            Assert.That(pep.AllModsOneIsNterminus.Keys.ToList().SequenceEqual(new int[] { 1, 3, 8 }));
        }

        /// <summary>
        /// CRITICAL: Tests correct parsing of C-terminal vs side-chain modifications.
        /// A modification on the last residue's side chain vs the peptide C-terminus are
        /// chemically distinct. Incorrect parsing would cause wrong mass calculations.
        /// </summary>
        [Test]
        public static void TestCTermAndLastSideChainModParsing()
        {
            string fullSequenceBothMods = "PEPTIDE[Mod:MyMod on E]-[PeptideCTermMod:MyCTermMod on E]";
            string fullSequenceCTermOnly = "PEPTIDE-[PeptideCTermMod:MyCTermMod on E]";
            string fullSequenceSideChainOnly = "PEPTIDE[Mod:MyMod on E]";

            ModificationMotif.TryGetMotif("E", out var motif);

            Modification mod = new Modification(_originalId: "MyMod", _modificationType: "Mod",
                _monoisotopicMass: 1, _locationRestriction: "Anywhere.", _target: motif);

            Modification cTermMod = new Modification(_originalId: "MyCTermMod", _modificationType: "PeptideCTermMod",
                _monoisotopicMass: 1, _locationRestriction: "Peptide C-terminal.", _target: motif);

            Dictionary<string, Modification> mods = new Dictionary<string, Modification>
            {
                { "MyMod on E", mod },
                { "MyCTermMod on E", cTermMod }
            };

            PeptideWithSetModifications pepBothMods = new PeptideWithSetModifications(fullSequenceBothMods, mods);
            PeptideWithSetModifications pepCterm = new PeptideWithSetModifications(fullSequenceCTermOnly, mods);
            PeptideWithSetModifications pepSideChain = new PeptideWithSetModifications(fullSequenceSideChainOnly, mods);

            Assert.That(pepBothMods.AllModsOneIsNterminus.Count == 2);
            Assert.That(pepBothMods.AllModsOneIsNterminus.Keys.SequenceEqual(new int[] { 8, 9 }));
            Assert.That(pepBothMods.AllModsOneIsNterminus[8].IdWithMotif == "MyMod on E");
            Assert.That(pepBothMods.AllModsOneIsNterminus[9].IdWithMotif == "MyCTermMod on E");
            Assert.That(pepCterm.AllModsOneIsNterminus.Count == 1);
            Assert.That(pepCterm.AllModsOneIsNterminus.Keys.SequenceEqual(new int[] { 9 }));
            Assert.That(pepCterm.AllModsOneIsNterminus[9].IdWithMotif == "MyCTermMod on E");
            Assert.That(pepSideChain.AllModsOneIsNterminus.Count == 1);
            Assert.That(pepSideChain.AllModsOneIsNterminus.Keys.SequenceEqual(new int[] { 8 }));
            Assert.That(pepSideChain.AllModsOneIsNterminus[8].IdWithMotif == "MyMod on E");
        }

        /// <summary>
        /// DEFENSIVE: Tests error handling for malformed peptide sequences.
        /// Validates that ambiguous sequences, bad modifications, and nonexistent mods
        /// throw MzLibException rather than crashing or returning invalid peptides.
        /// </summary>
        [Test]
        public static void BreakDeserializationMethod()
        {
            Assert.Throws<MzLibUtil.MzLibException>(() => new PeptideWithSetModifications("|", new Dictionary<string, Modification>())); // ambiguous
            Assert.Throws<MzLibUtil.MzLibException>(() => new PeptideWithSetModifications("[]", new Dictionary<string, Modification>())); // bad mod
            Assert.Throws<MzLibUtil.MzLibException>(() => new PeptideWithSetModifications("A[:mod]", new Dictionary<string, Modification>())); // nonexistent mod
        }

        /// <summary>
        /// CRITICAL: Tests MostAbundantMonoisotopicMass calculation accuracy.
        /// Compares against Protein Prospector values. Accurate mass calculation
        /// is essential for precursor mass matching in database searches.
        /// </summary>
        [Test]
        public static void CheckMostAbundantMonoisotopicMass()
        {
            PeptideWithSetModifications small_pep = new PeptideWithSetModifications(new Protein("PEPTIDE", "ACCESSION"), new DigestionParams(protease: "trypsin"), 1, 7, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            double small_pep_most_abundant_mass_prospector = 800.36724 - 1.0079;
            Assert.That(small_pep.MostAbundantMonoisotopicMass, Is.EqualTo(small_pep_most_abundant_mass_prospector).Within(0.01));

            PeptideWithSetModifications large_pep = new PeptideWithSetModifications(new Protein("PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE", "ACCESSION"), new DigestionParams(protease: "trypsin"), 1, 42, CleavageSpecificity.Full, null, 0, new Dictionary<int, Modification>(), 0, null);
            double large_pep_most_abundant_mass_prospector = 4709.12020 - 1.0079;
            Assert.That(large_pep.MostAbundantMonoisotopicMass, Is.EqualTo(large_pep_most_abundant_mass_prospector).Within(0.01));
        }

        /// <summary>
        /// CRITICAL: Tests sequence variant string generation for reporting.
        /// The SequenceVariantString method must correctly format variant annotations
        /// including N-terminal mods, middle mods, and truncated variants for output.
        /// Essential for human-readable variant peptide reporting.
        /// </summary>
        [Test]
        public static void TestSeqVarString()
        {
            Protein protein = new Protein("MACDEFGHIK", "test");

            // mod on N-terminus
            PeptideWithSetModifications pepe = new PeptideWithSetModifications(protein, new DigestionParams(), 1, 10, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification> { { 1, new Modification("mod on M", "mod", "mod", "mod") } }, 0);
            SequenceVariation sv1Before = new SequenceVariation(1, 1, "A", "M", ""); // n-terminal mod goes before the sequence
            Assert.AreEqual("A1[mod:mod on M]M", pepe.SequenceVariantString(sv1Before, true));

            // mod in middle
            PeptideWithSetModifications pepe2 = new PeptideWithSetModifications(protein, new DigestionParams(), 2, 10, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification> { { 2, new Modification("mod on A", "mod", "mod", "mod") } }, 0);
            SequenceVariation sv4MissenseBeginning = new SequenceVariation(2, 2, "V", "A", ""); // missense at beginning
            Assert.AreEqual("V2A[mod:mod on A]", pepe2.SequenceVariantString(sv4MissenseBeginning, true));

            // truncated seqvar doesn't truncate in string report (using applied variation correctly)
            PeptideWithSetModifications pepe3 = new PeptideWithSetModifications(protein, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification>(), 0);
            SequenceVariation svvvv = new SequenceVariation(7, 10, "GHM", "GHIK", ""); // insertion
            Assert.AreEqual("GHM7GHIK", pepe3.SequenceVariantString(svvvv, true));

            Protein protein2 = new Protein("WACDEFGHIK", "test");

            //variant starts at protein start but peptide does not
            PeptideWithSetModifications pepe4 = new PeptideWithSetModifications(protein2, new DigestionParams(), 4, 8, CleavageSpecificity.Unknown, "", 0, new Dictionary<int, Modification>(), 0);
            SequenceVariation variant = new SequenceVariation(1, 10, "MABCDEFGHIJKLMNOP", "WACDEFGHIK", ""); // frameshift
            Assert.AreEqual("MABCDEFGHIJKLMNOP1WACDEFGHIK", pepe4.SequenceVariantString(variant, true));
        }

        /// <summary>
        /// CRITICAL: Comprehensive test for variant identification and string methods.
        /// Tests all variant types (SAV, MNV, insertion, deletion, stop-gain, stop-loss)
        /// for correct intersection detection, identification logic, and string formatting.
        /// Essential for proteogenomic variant peptide analysis.
        /// </summary>
        [Test]
        public static void TestIdentifyandStringMethods()
        {
            ModificationMotif.TryGetMotif("V", out ModificationMotif motifV);
            Modification mv = new Modification("mod", null, "type", null, motifV, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            Modification mp = new Modification("mod", null, "type", null, motifP, "Anywhere.", null, 42.01, new Dictionary<string, IList<string>>(), null, null, null, null, null);
            Dictionary<int, Modification> modV = new Dictionary<int, Modification>();
            modV.Add(4, mv);
            Dictionary<int, Modification> modP = new Dictionary<int, Modification>();
            modP.Add(5, mp);

            Dictionary<int, List<Modification>> proteinPMods = new Dictionary<int, List<Modification>>();
            proteinPMods.Add(4, new List<Modification>() { mp });

            List<Protein> proteins = new List<Protein>
            {
                new Protein("MPEPTIDE", "protein0", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "V", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein1", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 5, "PT", "KT", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein2", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "PPP", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPPPTIDE", "protein3", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 6, "PPP", "P", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPKPKTIDE", "protein4", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 7, "PKPK", "PK", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTAIDE", "protein5",sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 6, "PTA", "KT", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEKKAIDE", "protein6", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 6, "KKA", "K", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein7", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "V", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", new Dictionary<int, List<Modification>> {{ 4, new[] { mv }.ToList() } }) }),
                new Protein("MPEPTIDE", "protein8",sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "PPP", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", new Dictionary<int, List<Modification>> {{ 5, new[] { mp }.ToList() } }) }),
                new Protein("MPEPTIDEPEPTIDE", "protein9", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 15, "PTIDEPEPTIDE", "PPP", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein10", oneBasedModifications: proteinPMods ,sequenceVariations: new List<SequenceVariation> { new SequenceVariation(4, 4, "P", "V", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein11", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(5, 5, "T", "*", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }), //stop-gain (can identify)
                new Protein("MPEKTIDE", "protein12", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(5, 5, "T", "*", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }), //stop-gain (can't identify)
                new Protein("MPEPTIPEPEPTIPE", "protein13", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(7, 7, "P", "D", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }),
                new Protein("MPEPTIDE", "protein14", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(8, 9, "E", "EK", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }), //peptide becomes longer, and cleavage site is created but cannot be identified
                new Protein("MPEPTIDE", "protein15", sequenceVariations: new List<SequenceVariation> { new SequenceVariation(9, 13, "*", "KMPEP", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=G||||||||||||||||\tGT:AD:DP\t1/1:30,30:30", null) }), // stop loss at end of original protein that cannot be identified
            };

            DigestionParams dp = new DigestionParams(minPeptideLength: 2);
            DigestionParams dp2 = new DigestionParams(protease: "Asp-N", minPeptideLength: 2);
            DigestionParams dp3 = new DigestionParams(protease: "Lys-N", minPeptideLength: 2);

            var protein0_variant = proteins.ElementAt(0).GetVariantBioPolymers().ElementAt(0);
            var protein1_variant = proteins.ElementAt(1).GetVariantBioPolymers().ElementAt(0);
            var protein2_variant = proteins.ElementAt(2).GetVariantBioPolymers().ElementAt(0);
            var protein3_variant = proteins.ElementAt(3).GetVariantBioPolymers().ElementAt(0);
            var protein4_variant = proteins.ElementAt(4).GetVariantBioPolymers().ElementAt(0);
            var protein5_variant = proteins.ElementAt(5).GetVariantBioPolymers().ElementAt(0);
            var protein6_variant = proteins.ElementAt(6).GetVariantBioPolymers().ElementAt(0);
            var protein7_variant = proteins.ElementAt(7).GetVariantBioPolymers().ElementAt(0);
            var protein8_variant = proteins.ElementAt(8).GetVariantBioPolymers().ElementAt(0);
            var protein9_variant = proteins.ElementAt(9).GetVariantBioPolymers().ElementAt(0);
            var protein10_variant = proteins.ElementAt(10).GetVariantBioPolymers().ElementAt(0);
            var protein11_variant = proteins.ElementAt(11).GetVariantBioPolymers().ElementAt(0);
            var protein12_variant = proteins.ElementAt(12).GetVariantBioPolymers().ElementAt(0);
            var protein13_variant = proteins.ElementAt(13).GetVariantBioPolymers().ElementAt(0);
            var protein14_variant = proteins.ElementAt(14).GetVariantBioPolymers().ElementAt(0);
            var protein15_variant = proteins.ElementAt(15).GetVariantBioPolymers().ElementAt(0);

            List<Modification> digestMods = new List<Modification>();

            var protein0_peptide = protein0_variant.Digest(dp, digestMods, digestMods).ElementAt(0);
            var protein0_peptide2 = protein0_variant.Digest(dp2, digestMods, digestMods).ElementAt(0);
            var protein1_peptide = protein1_variant.Digest(dp, digestMods, digestMods).ElementAt(2);
            var protein2_peptide = protein2_variant.Digest(dp, digestMods, digestMods).ElementAt(0);
            var protein3_peptide = protein3_variant.Digest(dp, digestMods, digestMods).ElementAt(0);
            var protein4_peptide = protein4_variant.Digest(dp, digestMods, digestMods).ElementAt(2);
            var protein5_peptide = protein5_variant.Digest(dp, digestMods, digestMods).ElementAt(2);
            var protein6_peptide = protein6_variant.Digest(dp, digestMods, digestMods).ElementAt(2);
            var protein7_peptide = protein7_variant.Digest(dp, digestMods, digestMods).ElementAt(1);
            var protein8_peptide = protein8_variant.Digest(dp, digestMods, digestMods).ElementAt(1);
            var protein9_peptide = protein9_variant.Digest(dp, digestMods, digestMods).ElementAt(0);
            var protein10_peptide = protein10_variant.Digest(dp, digestMods, digestMods).ElementAt(0);
            var protein11_peptide = protein11_variant.Digest(dp2, digestMods, digestMods).ElementAt(0);
            var protein11_peptide2 = protein11_variant.Digest(dp, digestMods, digestMods).ElementAt(0);
            var protein12_peptide = protein12_variant.Digest(dp, digestMods, digestMods).ElementAt(0);
            var protein13_peptide = protein13_variant.Digest(dp2, digestMods, digestMods).ElementAt(0);
            var protein14_peptide = protein14_variant.Digest(dp3, digestMods, digestMods).ElementAt(0);
            var protein15_peptide = protein15_variant.Digest(dp3, digestMods, digestMods).ElementAt(0);

            Assert.AreEqual((true, true), protein0_peptide.IntersectsAndIdentifiesVariation(protein0_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein0_peptide2.IntersectsAndIdentifiesVariation(protein0_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein1_peptide.IntersectsAndIdentifiesVariation(protein1_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein2_peptide.IntersectsAndIdentifiesVariation(protein2_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein3_peptide.IntersectsAndIdentifiesVariation(protein3_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((false, false), protein4_peptide.IntersectsAndIdentifiesVariation(protein4_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein5_peptide.IntersectsAndIdentifiesVariation(protein5_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((false, true), protein6_peptide.IntersectsAndIdentifiesVariation(protein6_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein7_peptide.IntersectsAndIdentifiesVariation(protein7_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein8_peptide.IntersectsAndIdentifiesVariation(protein8_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein9_peptide.IntersectsAndIdentifiesVariation(protein9_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, true), protein10_peptide.IntersectsAndIdentifiesVariation(protein10_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((false, true), protein11_peptide.IntersectsAndIdentifiesVariation(protein11_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((false, true), protein11_peptide2.IntersectsAndIdentifiesVariation(protein11_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((false, false), protein12_peptide.IntersectsAndIdentifiesVariation(protein12_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((false, true), protein13_peptide.IntersectsAndIdentifiesVariation(protein13_variant.AppliedSequenceVariations.ElementAt(0)));
            Assert.AreEqual((true, false), protein14_peptide.IntersectsAndIdentifiesVariation(protein14_variant.AppliedSequenceVariations.ElementAt(0)));// the peptide crosses the variant but the newly genrated cleavage site makes the same peptide as without the variant
            Assert.AreEqual((false, false), protein15_peptide.IntersectsAndIdentifiesVariation(protein15_variant.AppliedSequenceVariations.ElementAt(0)));// the peptide does not cross the variant, and the stop loss adds addition amino acids, but it creates the same peptide as without the variant

            Assert.AreEqual("P4V", protein0_peptide.SequenceVariantString(protein0_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("P4V", protein0_peptide2.SequenceVariantString(protein0_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("PT4KT", protein1_peptide.SequenceVariantString(protein1_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("P4PPP", protein2_peptide.SequenceVariantString(protein2_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("PPP4P", protein3_peptide.SequenceVariantString(protein3_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("PTA4KT", protein5_peptide.SequenceVariantString(protein5_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("KKA4K", protein6_peptide.SequenceVariantString(protein6_variant.AppliedSequenceVariations.ElementAt(0), false));
            Assert.AreEqual("P4V[type:mod on V]", protein7_peptide.SequenceVariantString(protein7_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("P4PP[type:mod on P]P", protein8_peptide.SequenceVariantString(protein8_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("PTIDEPEPTIDE4PPP", protein9_peptide.SequenceVariantString(protein9_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("P4V", protein10_peptide.SequenceVariantString(protein10_variant.AppliedSequenceVariations.ElementAt(0), true));
            Assert.AreEqual("T5*", protein11_peptide.SequenceVariantString(protein11_variant.AppliedSequenceVariations.ElementAt(0), false));
            Assert.AreEqual("T5*", protein11_peptide2.SequenceVariantString(protein11_variant.AppliedSequenceVariations.ElementAt(0), false));
            Assert.AreEqual("P7D", protein13_peptide.SequenceVariantString(protein13_variant.AppliedSequenceVariations.ElementAt(0), false));
        }

        /// <summary>
        /// CRITICAL: Tests EssentialSequence generation for peptide grouping.
        /// Essential sequences include only user-specified modification types and are
        /// used for peptide-level grouping and quantification. Regression test against
        /// expected output file ensures consistency across versions.
        /// </summary>
        [Test]
        public static void TestPeptideWithSetModsEssentialSequence()
        {
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            List<Modification> UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();

            Dictionary<string, int> modsToWrite = new Dictionary<string, int>();
            modsToWrite.Add("UniProt",0);

            var proteinXml = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "humanGAPDH.xml"), true, DecoyType.None, UniProtPtms, false, null, out var unknownMod);
            var gapdh = proteinXml[0];

            var gapdhPeptides = gapdh.Digest(new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Variable), UniProtPtms, new List<Modification>());

            List<string> allSequences = new List<string>();
            foreach (var peptide in gapdhPeptides)
            {
                allSequences.Add(peptide.EssentialSequence(modsToWrite));
            }

            var expectedFullStrings = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "essentialSequences.txt"));

            CollectionAssert.AreEquivalent(expectedFullStrings, allSequences.ToArray());
        }

        [Test]
        public static void TestPeptideWithSetModsEssentialSequence_SequenceSerializerMatchesExtension()
        {
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            List<Modification> UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList();

            Dictionary<string, int> modsToWrite = new Dictionary<string, int>
            {
                { "UniProt", 0 }
            };

            var proteinXml = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "humanGAPDH.xml"), true, DecoyType.None, UniProtPtms, false, null, out var unknownMod);
            var gapdh = proteinXml[0];

            var gapdhPeptides = gapdh.Digest(new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Variable), UniProtPtms, new List<Modification>());

            var serializer = new EssentialSequenceSerializer(modsToWrite);
            var service = SequenceConversionService.Default;

            List<string> extensionResults = new List<string>();
            List<string> serializerResults = new List<string>();
            List<string> serviceFullResults = new List<string>();
            List<string> serviceStepResults = new List<string>();
            List<string> serviceAutoFullResults = new List<string>();
            List<string> serviceAutoStepResults = new List<string>();
            foreach (var peptide in gapdhPeptides)
            {
                extensionResults.Add(peptide.EssentialSequence(modsToWrite));

                var canonical = peptide.ToCanonicalSequenceBuilder().Build();
                var serialized = serializer.Serialize(canonical);
                serializerResults.Add(serialized ?? string.Empty);

                var serviceFullConversion = service.Convert(peptide.FullSequence, "mzLib", "essential");
                serviceFullResults.Add(serviceFullConversion ?? string.Empty);

                var canon = service.Parse(peptide.FullSequence, "mzLib");
                var serviceStepConversion = service.Serialize(canon.Value, "essential");
                serviceStepResults.Add(serviceStepConversion ?? string.Empty);

                var serviceAutoFull = service.ConvertAutoDetect(peptide.FullSequence, "essential");
                serviceAutoFullResults.Add(serviceAutoFull ?? string.Empty);

                var serviceAutoStep = service.ParseAutoDetect(peptide.FullSequence);
                var serviceAutoStepConversion = service.Serialize(serviceAutoStep.Value, "essential");
                serviceAutoStepResults.Add(serviceAutoStepConversion ?? string.Empty);
            }

            CollectionAssert.AreEquivalent(extensionResults, serializerResults);
            CollectionAssert.AreEquivalent(extensionResults, serviceFullResults);
            CollectionAssert.AreEquivalent(extensionResults, serviceStepResults);
            CollectionAssert.AreEquivalent(extensionResults, serviceAutoFullResults);
            CollectionAssert.AreEquivalent(extensionResults, serviceAutoStepResults);


            var expectedFullStrings = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "essentialSequences.txt"));
            CollectionAssert.AreEquivalent(expectedFullStrings, serializerResults.ToArray());
            CollectionAssert.AreEquivalent(expectedFullStrings, serviceFullResults.ToArray());
            CollectionAssert.AreEquivalent(expectedFullStrings, serviceStepResults.ToArray());
            CollectionAssert.AreEquivalent(expectedFullStrings, serviceAutoFullResults.ToArray());
            CollectionAssert.AreEquivalent(expectedFullStrings, serviceAutoStepResults.ToArray());
        }
        /// <summary>
        /// CRITICAL: Tests FullSequence and FullSequenceWithMassShift generation.
        /// These string representations are used for peptide identification, reporting,
        /// and spectral library matching. Regression test against expected files
        /// ensures output format stability across versions.
        /// </summary>
        [Test]
        public static void TestPeptideWithSetModsFullSequence()
        {
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            List<Modification> UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"), formalChargesDictionary).ToList(); 
            var proteinXml = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "humanGAPDH.xml"), true, DecoyType.None, UniProtPtms, false, null, out var unknownMod);
            var gapdh = proteinXml[0];

            var gapdhPeptides = gapdh.Digest(new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, minPeptideLength:1, initiatorMethionineBehavior:InitiatorMethionineBehavior.Variable),UniProtPtms,new List<Modification>());
            
            List<string> allSequences = new List<string>();
            foreach (var peptide in gapdhPeptides)
            {
                allSequences.Add(peptide.FullSequence);
            }

            var expectedFullStrings = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "fullSequences.txt"));
            CollectionAssert.AreEquivalent(expectedFullStrings,allSequences.ToArray());

            allSequences.Clear();
            foreach (var peptide in gapdhPeptides)
            {
                allSequences.Add(peptide.FullSequenceWithMassShift());
            }

            var expectedFullStringsWithMassShifts = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "fullSequencesWithMassShift.txt"));
            CollectionAssert.AreEquivalent(expectedFullStringsWithMassShifts, allSequences.ToArray());
        }

        /// <summary>
        /// CRITICAL: Tests modification parsing from FullSequence strings.
        /// Validates that GetModificationDictionaryFromFullSequence and
        /// GetModificationsFromFullSequence correctly parse modifications from
        /// peptide strings. Essential for reading results files and spectral libraries.
        /// </summary>
        [Test]
        public static void TestIBioPolymerWithSetModsModificationFromFullSequence()
        {
            Dictionary<string, Modification> un = new Dictionary<string, Modification>();
            var psiModDeserialized = Loaders.LoadPsiMod(Path.Combine(TestContext.CurrentContext.TestDirectory, "PSI-MOD.obo2.xml"));
            Dictionary<string, int> formalChargesDictionary = Loaders.GetFormalChargesDictionary(psiModDeserialized);
            List<Modification> UniProtPtms = Loaders.LoadUniprot(Path.Combine(TestContext.CurrentContext.TestDirectory, "ptmlist2.txt"),
                    formalChargesDictionary).ToList();
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, "DatabaseTests", "cRAP_databaseGPTMD.xml"),
                true, DecoyType.None, UniProtPtms, false, new string[] { "exclude_me" }, out un);
            var allKnownModDict = UniProtPtms.ToDictionary(p => p.IdWithMotif, p => p);
            var digestionParameters = new DigestionParams(maxModsForPeptides: 3);

            foreach (Protein p in proteins)
            {
                List<PeptideWithSetModifications> digestedPeptides =
                    p.Digest(digestionParameters, [], [], null, null).ToList();
                // take the most modified peptide by base sequence and ensure all methods function properly
                foreach (var targetPeptide in digestedPeptides
                             .Where(pep => pep.FullSequence.Contains('['))
                             .GroupBy(pep => pep.BaseSequence)
                             .Select(pepGroup => pepGroup.MaxBy(pep => pep.AllModsOneIsNterminus.Count)))
                {
                    var startResidue = targetPeptide.OneBasedStartResidue;
                    var endResidue = targetPeptide.OneBasedEndResidue;

                    // Pull our expected modifications based upon parent protein object with a maximum value of DigestionParameters.MaxMods
                    // A bunch of logic to count the number of expected modifications based upon the xml database entries
                    int expectedModCount = 0;
                    foreach (var modDictEntry in p.OneBasedPossibleLocalizedModifications
                                 .Where(mod => mod.Key >= startResidue && mod.Key <= endResidue))
                    {
                        if (modDictEntry.Value.Count > 1)
                        {
                            var locRestrictions = modDictEntry.Value.Select(mod => mod.LocationRestriction).ToList();

                            if (locRestrictions.AllSame())
                            {
                                if (locRestrictions.First() == "Anywhere.")
                                    expectedModCount++;
                                else if (locRestrictions.First() == "N-terminal." && modDictEntry.Key == startResidue)
                                    expectedModCount++;
                            }
                            else if (modDictEntry.Value.Select(mod => mod.LocationRestriction).Contains("Anywhere.")
                                     && modDictEntry.Value.Select(mod => mod.LocationRestriction)
                                         .Contains("N-terminal."))
                            {
                                expectedModCount++;
                                if (modDictEntry.Key == startResidue)
                                    expectedModCount++;
                            }
                        }
                        else
                        {
                            switch (modDictEntry.Value.First().LocationRestriction)
                            {
                                case "Anywhere.":
                                case "N-terminal." when modDictEntry.Key == startResidue:
                                    expectedModCount++;
                                    break;
                            }
                        }
                    }

                    expectedModCount = Math.Min(expectedModCount, digestionParameters.MaxMods);

                    var expectedModifications = p.OneBasedPossibleLocalizedModifications.Where(mod =>
                        mod.Key >= startResidue &&
                        mod.Key <= endResidue).SelectMany(mod => mod.Value).ToList();

                    // Parse modifications from PWSM and two IBioPolymerWithSetMods methods
                    var pwsmModDict = targetPeptide.AllModsOneIsNterminus;
                    var bpwsmModDict = IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence(targetPeptide.FullSequence, allKnownModDict);
                    var bpwsmModList = IBioPolymerWithSetMods.GetModificationsFromFullSequence(targetPeptide.FullSequence, allKnownModDict);

                    // Ensure all methods are in agreement by modification count
                    Assert.AreEqual(pwsmModDict.Count, expectedModCount);
                    Assert.AreEqual(bpwsmModDict.Count, expectedModCount);
                    Assert.AreEqual(bpwsmModList.Count, expectedModCount);

                    // Ensure all methods are in agreement by modification identify
                    foreach (var pwsmModification in pwsmModDict.Values)
                        Assert.Contains(pwsmModification, expectedModifications);
                    foreach (var pwsmModification in bpwsmModDict.Values)
                        Assert.Contains(pwsmModification, expectedModifications);
                    foreach (var pwsmModification in bpwsmModList)
                        Assert.Contains(pwsmModification, expectedModifications);
                }
            }
        }

        /// <summary>
        /// CRITICAL: Tests parsing of nucleotide substitution annotations in full sequences.
        /// The ParseSubstitutedFullSequence method must correctly apply amino acid
        /// substitutions from GPTMD-discovered variants while preserving other modifications.
        /// Essential for variant peptide sequence reconstruction.
        /// </summary>
        [Test]
        public static void TestGetSubstitutedFullSequence()
        {
            //It should take care of multiple substitutions
            string test1 = "F[1 nucleotide substitution:F->Y on F]SIMGGGLA[1 nucleotide substitution:A->S on A]DR";
            string expected1 = "YSIMGGGLSDR";
            var actual1 = IBioPolymerWithSetMods.ParseSubstitutedFullSequence(test1);
            Assert.That(actual1, Is.EqualTo(expected1));

            //It should not change other modifications
            string test2 = "SANH[1 nucleotide substitution:H->L on H]M[Common Variable:Oxidation on M]AGHWVAISGAAGGLGSLAVQYAK";
            string expected2 = "SANLM[Common Variable:Oxidation on M]AGHWVAISGAAGGLGSLAVQYAK";
            var actual2 = IBioPolymerWithSetMods.ParseSubstitutedFullSequence(test2);
            Assert.That(actual2, Is.EqualTo(expected2));

            //It should work on 2 nucleotide substitutions
            string test3 = "S[2+ nucleotide substitution:S->E on S]AAADRLNLTSGHLNAGR";
            string expected3 = "EAAADRLNLTSGHLNAGR";
            var actual3 = IBioPolymerWithSetMods.ParseSubstitutedFullSequence(test3);
            Assert.That(actual3, Is.EqualTo(expected3));
        }
    }
}
