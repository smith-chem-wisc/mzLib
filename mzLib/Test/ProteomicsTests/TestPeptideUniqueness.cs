using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Proteomics;

namespace Test.ProteomicsTests
{
    /// <summary>
    /// Classifying a peptide as unique to one sequence, shared among isoforms of one gene, or shared
    /// across genes, with I and L treated as one residue.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestPeptideUniqueness
    {
        private const string Human = "Homo sapiens";

        private static Protein Entry(string sequence, string accession, string gene = null, string organism = Human,
            bool isDecoy = false, bool isContaminant = false, params string[] ensemblGenes) =>
            new Protein(sequence, accession, organism,
                gene == null ? null : new List<Tuple<string, string>> { new("primary", gene) },
                isDecoy: isDecoy, isContaminant: isContaminant,
                databaseReferences: ensemblGenes.Select((g, i) => new DatabaseReference("Ensembl", $"ENST0000000000{i}",
                    new List<Tuple<string, string>> { new("gene ID", g + ".1") })).ToList());

        private static PeptideUniqueness One(string peptide, params Protein[] proteins) =>
            PeptideUniquenessClassifier.Classify(new[] { peptide }, proteins).Single();

        [Test]
        public void PeptideInOneProtein_IsUnique_AndNamesItsGene()
        {
            var result = One("SAMPLER", Entry("MKSAMPLERK", "P00001", "GENEA"), Entry("MKOTHERSEQK", "P00002", "GENEB"));

            Assert.That(result.Sharing, Is.EqualTo(PeptideSharing.Unique));
            Assert.That(result.Accessions, Is.EqualTo(new[] { "P00001" }));
            Assert.That(result.SharedGeneKeys, Is.EqualTo(new[] { "entry:P00001", "gene:Homo sapiens:GENEA" }));
        }

        [Test]
        public void IdenticalSequencesUnderTwoAccessions_AreOneSequence_SoTheirPeptideIsUnique()
        {
            var result = One("SAMPLER", Entry("MKSAMPLERK", "P00001", "GENEA"), Entry("MKSAMPLERK", "P00002", "GENEB"));

            Assert.That(result.Sharing, Is.EqualTo(PeptideSharing.Unique));
            Assert.That(result.Accessions, Is.EqualTo(new[] { "P00001", "P00002" }));
        }

        [Test]
        public void TwoIsoformsOfOneEntry_ShareWithinGene_ThroughTheEntryAccession()
        {
            // A FASTA isoform carries no gene name here; the entry accession alone joins it to its canonical.
            var result = One("SAMPLER", Entry("MKSAMPLERK", "P00001", "GENEA"), Entry("MKSAMPLERAAK", "P00001-2"));

            Assert.That(result.Sharing, Is.EqualTo(PeptideSharing.SharedWithinGene));
            Assert.That(result.SharedGeneKeys, Is.EqualTo(new[] { "entry:P00001" }));
        }

        [Test]
        public void TwoEntriesOfOneEnsemblGene_ShareWithinGene()
        {
            var result = One("SAMPLER",
                Entry("MKSAMPLERK", "P00001", ensemblGenes: "ENSG00000000001"),
                Entry("MKSAMPLERAAK", "Q00002", ensemblGenes: "ENSG00000000001"));

            Assert.That(result.Sharing, Is.EqualTo(PeptideSharing.SharedWithinGene));
            Assert.That(result.SharedGeneKeys, Is.EqualTo(new[] { "ensembl:ENSG00000000001" }));
        }

        [Test]
        public void DifferentGenes_AreSharedAcrossGenes_WithNoSharedKeys()
        {
            var result = One("SAMPLER", Entry("MKSAMPLERK", "P00001", "GENEA"), Entry("MKSAMPLERAAK", "P00002", "GENEB"));

            Assert.That(result.Sharing, Is.EqualTo(PeptideSharing.SharedAcrossGenes));
            Assert.That(result.Accessions, Is.EqualTo(new[] { "P00001", "P00002" }));
            Assert.That(result.SharedGeneKeys, Is.Empty);
        }

        [Test]
        public void OneGeneNameInTwoOrganisms_IsTwoGenes()
        {
            // Bovine serum albumin is a contaminant in human samples; its gene symbol must not merge with ours.
            var result = One("SAMPLER",
                Entry("MKSAMPLERK", "P02768", "ALB"),
                Entry("MKSAMPLERAAK", "P02769", "ALB", organism: "Bos taurus", isContaminant: true));

            Assert.That(result.Sharing, Is.EqualTo(PeptideSharing.SharedAcrossGenes));
        }

        [Test]
        public void IsoleucineAndLeucine_AreOneResidue()
        {
            // The peptide says I, one protein says L, the other says I: both contain it.
            var result = One("SAMPIER", Entry("MKSAMPLERK", "P00001", "GENEA"), Entry("MKSAMPIERAAK", "P00002", "GENEB"));

            Assert.That(result.Sharing, Is.EqualTo(PeptideSharing.SharedAcrossGenes));
            Assert.That(result.Peptide, Is.EqualTo("SAMPIER"));
        }

        [Test]
        public void SequencesDifferingOnlyByIAndL_AreOneSequence()
        {
            var result = One("SAMPLER", Entry("MKSAMPLERIK", "P00001", "GENEA"), Entry("MKSAMPLERLK", "P00002", "GENEB"));

            Assert.That(result.Sharing, Is.EqualTo(PeptideSharing.Unique));
        }

        [Test]
        public void ContainmentIgnoresProteaseRules()
        {
            // SAMPLER sits after a P here, where trypsin would not cut; it is still in the protein.
            var result = One("SAMPLER", Entry("MKSAMPLERK", "P00001", "GENEA"), Entry("MPSAMPLERP", "P00002", "GENEB"));

            Assert.That(result.Sharing, Is.EqualTo(PeptideSharing.SharedAcrossGenes));
        }

        [Test]
        public void DecoysAreIgnored_ContaminantsAreNot()
        {
            var decoy = One("SAMPLER", Entry("MKSAMPLERK", "P00001", "GENEA"),
                Entry("MKSAMPLERK", "DECOY_P00002", "GENEB", isDecoy: true));
            var contaminant = One("SAMPLER", Entry("MKSAMPLERK", "P00001", "GENEA"),
                Entry("MKSAMPLERAAK", "CON__P00002", "GENEB", isContaminant: true));

            Assert.That(decoy.Sharing, Is.EqualTo(PeptideSharing.Unique));
            Assert.That(decoy.Accessions, Is.EqualTo(new[] { "P00001" }));
            Assert.That(contaminant.Sharing, Is.EqualTo(PeptideSharing.SharedAcrossGenes));
        }

        [Test]
        public void AContaminantIdenticalToATarget_IsOneSequenceWithIt_SoThePeptideIsUniqueAndListsBoth()
        {
            var result = One("SAMPLER", Entry("MKSAMPLERK", "P00001", "GENEA"),
                Entry("MKSAMPLERK", "CON__P00002", "GENEB", isContaminant: true));

            Assert.That(result.Sharing, Is.EqualTo(PeptideSharing.Unique));
            Assert.That(result.Accessions, Is.EqualTo(new[] { "CON__P00002", "P00001" }));
        }

        [Test]
        public void AOneResiduePeptide_DoesNotHideLongerOnes_AtAnyKeyLength()
        {
            // One index per key length: the 1-mer, a 5-mer, a 12-mer and a 14-mer are all found.
            var proteins = new[]
            {
                Entry("MKWPEPTIDEKLONGERPEPTIDEKR", "P00001", "GENEA"),
                Entry("MKPEPTLDEKGGR", "P00002", "GENEB")
            };

            var results = PeptideUniquenessClassifier.Classify(
                new[] { "W", "GGR", "LONGERPEPTLD", "LONGERPEPTIDEK", "ABSENTPEPTIDEK" }, proteins);

            Assert.That(results.Select(r => r.Sharing), Is.EqualTo(new[]
            {
                PeptideSharing.Unique, PeptideSharing.Unique, PeptideSharing.Unique, PeptideSharing.Unique,
                PeptideSharing.NotInDatabase
            }));
            Assert.That(results.Select(r => r.Accessions.Single()).Take(2), Is.EqualTo(new[] { "P00001", "P00002" }));
        }

        [Test]
        public void AbsentPeptide_IsNotInDatabase()
        {
            var result = One("NOWHERE", Entry("MKSAMPLERK", "P00001", "GENEA"));

            Assert.That(result.Sharing, Is.EqualTo(PeptideSharing.NotInDatabase));
            Assert.That(result.Accessions, Is.Empty);
            Assert.That(result.SharedGeneKeys, Is.Empty);
        }

        [Test]
        public void OneResultPerInput_InOrder_DuplicatesAndIAndLVariantsKeptVerbatim()
        {
            var proteins = new[] { Entry("MKSAMPLERKPEPTLDEK", "P00001", "GENEA") };

            var results = PeptideUniquenessClassifier.Classify(
                new[] { "PEPTIDEK", "SAMPLER", "PEPTLDEK", "PEPTIDEK", "ABSENT" }, proteins);

            Assert.That(results.Select(r => r.Peptide),
                Is.EqualTo(new[] { "PEPTIDEK", "SAMPLER", "PEPTLDEK", "PEPTIDEK", "ABSENT" }));
            Assert.That(results.Select(r => r.Sharing), Is.EqualTo(new[]
            {
                PeptideSharing.Unique, PeptideSharing.Unique, PeptideSharing.Unique, PeptideSharing.Unique,
                PeptideSharing.NotInDatabase
            }));
        }

        [Test]
        public void PeptidesOfMixedLengths_AreAllFound_IncludingOneLongerThanTheIndexKey()
        {
            // Each peptide is keyed on its own prefix, capped at twelve residues; longer peptides are verified in full.
            var proteins = new[]
            {
                Entry("MKAKLONGPEPTIDESEQUENCEKGGK", "P00001", "GENEA"),
                Entry("MKLONGPEPTIDESEQUENCEXXK", "P00002", "GENEB")
            };

            var results = PeptideUniquenessClassifier.Classify(
                new[] { "AK", "LONGPEPTIDESEQUENCEK", "LONGPEPTIDESEQUENCE", "GGK" }, proteins);

            Assert.That(results.Select(r => r.Sharing), Is.EqualTo(new[]
            {
                PeptideSharing.Unique, PeptideSharing.Unique, PeptideSharing.SharedAcrossGenes, PeptideSharing.Unique
            }));
        }

        [Test]
        public void APeptideSeenTwiceInOneProtein_ListsTheProteinOnce()
        {
            var result = One("SAMPLER", Entry("SAMPLERSAMPLER", "P00001", "GENEA"));

            Assert.That(result.Accessions, Is.EqualTo(new[] { "P00001" }));
        }

        [Test]
        public void ANonLetterInTheProtein_BreaksTheWindow()
        {
            var result = One("SAMPLER", Entry("MKSAM*PLERK", "P00001", "GENEA"));

            Assert.That(result.Sharing, Is.EqualTo(PeptideSharing.NotInDatabase));
        }

        [Test]
        public void ACustomGeneKey_Replaces_TheDefault()
        {
            var proteins = new[] { Entry("MKSAMPLERK", "P00001", "GENEA"), Entry("MKSAMPLERAAK", "P00002", "GENEB") };

            var result = PeptideUniquenessClassifier.Classify(new[] { "SAMPLER" }, proteins, _ => new[] { "family:X" })
                .Single();

            Assert.That(result.Sharing, Is.EqualTo(PeptideSharing.SharedWithinGene));
            Assert.That(result.SharedGeneKeys, Is.EqualTo(new[] { "family:X" }));
        }

        [Test]
        public void DefaultGeneKeys_CarryEveryEnsemblGene_TheOrganismQualifiedName_AndTheEntry()
        {
            var histone = Entry("MSGRGK", "P62805-3", "H4C1", ensemblGenes: new[] { "ENSG00000000002", "ENSG00000000001" });

            Assert.That(PeptideUniquenessClassifier.DefaultGeneKeys(histone), Is.EqualTo(new[]
            {
                "ensembl:ENSG00000000001", "ensembl:ENSG00000000002", "gene:Homo sapiens:H4C1", "entry:P62805"
            }));
        }

        [Test]
        public void DefaultGeneKeys_UseAnUnrecognizedAccessionVerbatim()
        {
            Assert.That(PeptideUniquenessClassifier.DefaultGeneKeys(Entry("MSGRGK", "CON__P02769")),
                Is.EqualTo(new[] { "entry:CON__P02769" }));
        }

        [Test]
        public void NoPeptides_GivesNoResults()
        {
            Assert.That(PeptideUniquenessClassifier.Classify(Array.Empty<string>(), new[] { Entry("MKSAMPLERK", "P00001") }),
                Is.Empty);
        }

        [Test]
        public void InvalidInput_Throws()
        {
            var proteins = new[] { Entry("MKSAMPLERK", "P00001") };

            Assert.Throws<ArgumentNullException>(() => PeptideUniquenessClassifier.Classify(null, proteins));
            Assert.Throws<ArgumentNullException>(() => PeptideUniquenessClassifier.Classify(new[] { "SAMPLER" }, null));
            Assert.Throws<ArgumentNullException>(() => PeptideUniquenessClassifier.Classify(new string[] { null }, proteins));
            Assert.Throws<ArgumentNullException>(() => PeptideUniquenessClassifier.Classify(new[] { "SAMPLER" }, new Protein[] { null }));
            Assert.Throws<ArgumentNullException>(() => PeptideUniquenessClassifier.DefaultGeneKeys(null));
            foreach (var bad in new[] { "", "sampler", "SAM[Oxidation]PLER", "SAMP LER" })
            {
                Assert.Throws<ArgumentException>(() => PeptideUniquenessClassifier.Classify(new[] { bad }, proteins), bad);
            }
        }
    }
}
