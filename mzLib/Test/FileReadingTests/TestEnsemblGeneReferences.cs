using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Proteomics;
using UsefulProteomicsDatabases;

namespace Test.FileReadingTests
{
    /// <summary>
    /// A UniProt entry's Ensembl gene ids are already parsed and stored -- ProteinXmlEntry keeps every
    /// &lt;dbReference&gt; generically -- and nothing reads them. Protein.EnsemblGeneReferences is a
    /// derived view over DatabaseReferences, following Protein.GoTerms.
    ///
    /// The shape differs from NcbiTaxonomyId, and the tests pin why:
    ///  - UniProt writes one Ensembl dbReference per TRANSCRIPT, with the gene id in a
    ///    &lt;property type="gene ID"&gt;, so the gene is not the reference's Id and a
    ///    FirstOrDefault(...)?.Id projection returns a transcript;
    ///  - one accession can map to several genes (the core histones span many loci), so the view
    ///    returns every gene and never picks one;
    ///  - UniProt's gene ids are versioned ("ENSG00000111640.15"); the stable id and the version are
    ///    both kept, because joins want the first and audits want the second.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestEnsemblGeneReferences
    {
        private static string Data(params string[] parts) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, Path.Combine(parts));

        private static Protein LoadGapdh() => ProteinDbLoader
            .LoadProteinXML(Data("DatabaseTests", "humanGAPDH.xml"), true, DecoyType.None, null, false, null, out _)
            .First(p => p.Accession == "P04406");

        private static Protein WithReferences(params DatabaseReference[] references) =>
            new Protein("PEPTIDEK", "P00001", databaseReferences: references.ToList());

        private static DatabaseReference Ensembl(string transcriptId, params (string Type, string Value)[] properties) =>
            new DatabaseReference("Ensembl", transcriptId,
                properties.Select(p => new Tuple<string, string>(p.Type, p.Value)).ToList());

        [Test]
        public void EnsemblGeneReferences_ParsedFromUniProtXml()
        {
            var protein = LoadGapdh();

            Assert.That(protein.EnsemblGeneReferences.Count, Is.EqualTo(5), "humanGAPDH.xml carries 5 Ensembl transcripts");
            Assert.That(protein.EnsemblGeneIds, Is.EqualTo(new[] { "ENSG00000111640" }), "five transcripts, one gene");

            var first = protein.EnsemblGeneReferences.Single(r => r.TranscriptId == "ENST00000229239.10");
            Assert.That(first.ProteinId, Is.EqualTo("ENSP00000229239.5"));
            Assert.That(first.VersionedGeneId, Is.EqualTo("ENSG00000111640.15"));
            Assert.That(first.GeneId, Is.EqualTo("ENSG00000111640"));
            Assert.That(first.GeneVersion, Is.EqualTo(15));
        }

        [Test]
        public void EnsemblGeneReferences_ExcludeNonEnsemblDatabaseReferences()
        {
            var protein = WithReferences(
                new DatabaseReference("GO", "GO:0005737", new List<Tuple<string, string>>()),
                new DatabaseReference("RefSeq", "NP_002037.2", new List<Tuple<string, string>>()),
                Ensembl("ENST00000000001.1", ("gene ID", "ENSG00000000001.1")));

            Assert.That(protein.EnsemblGeneReferences.Select(r => r.TranscriptId), Is.EqualTo(new[] { "ENST00000000001.1" }));
        }

        [Test]
        public void EnsemblGeneIds_EveryGeneOfAMultiGeneAccession_NeverAPick()
        {
            var protein = WithReferences(
                Ensembl("ENST00000000003.1", ("gene ID", "ENSG00000000009.2")),
                Ensembl("ENST00000000001.1", ("gene ID", "ENSG00000000002.4")),
                Ensembl("ENST00000000002.1", ("gene ID", "ENSG00000000002.4")));

            Assert.That(protein.EnsemblGeneIds, Is.EqualTo(new[] { "ENSG00000000002", "ENSG00000000009" }),
                "distinct stable ids, ordinal order, so output does not depend on XML order");
        }

        [Test]
        public void EnsemblGeneReferences_PropertiesMatchedByType_NeverByPosition()
        {
            var geneFirst = WithReferences(Ensembl("ENST00000000001.1",
                ("gene ID", "ENSG00000000001.3"), ("protein sequence ID", "ENSP00000000001.2")));
            var proteinFirst = WithReferences(Ensembl("ENST00000000001.1",
                ("protein sequence ID", "ENSP00000000001.2"), ("gene ID", "ENSG00000000001.3")));

            foreach (var protein in new[] { geneFirst, proteinFirst })
            {
                var r = protein.EnsemblGeneReferences.Single();
                Assert.That(r.VersionedGeneId, Is.EqualTo("ENSG00000000001.3"));
                Assert.That(r.ProteinId, Is.EqualTo("ENSP00000000001.2"));
            }
        }

        [Test]
        public void EnsemblGeneReferences_ReferenceWithoutAGeneId_IsNotAGeneLink()
        {
            var protein = WithReferences(Ensembl("ENST00000000001.1", ("protein sequence ID", "ENSP00000000001.2")));

            Assert.That(protein.EnsemblGeneReferences, Is.Empty);
            Assert.That(protein.EnsemblGeneIds, Is.Empty);
        }

        [Test]
        public void EnsemblGeneReferences_VersionSplitOnlyOnANumericSuffix()
        {
            var protein = WithReferences(
                Ensembl("ENST00000000001.1", ("gene ID", "ENSG00000000001")),
                Ensembl("ENST00000000002.1", ("gene ID", "ENSG00000000002.x")));

            var unversioned = protein.EnsemblGeneReferences.Single(r => r.TranscriptId == "ENST00000000001.1");
            Assert.That(unversioned.GeneId, Is.EqualTo("ENSG00000000001"));
            Assert.That(unversioned.GeneVersion, Is.Null, "absent is not version 0");

            var odd = protein.EnsemblGeneReferences.Single(r => r.TranscriptId == "ENST00000000002.1");
            Assert.That(odd.GeneId, Is.EqualTo("ENSG00000000002.x"), "not a version, so not split");
            Assert.That(odd.GeneVersion, Is.Null);
        }

        [Test]
        public void EnsemblGeneReferences_RepeatedTranscriptCountedOnce()
        {
            // ProteinDbLoader's duplicate-entry merge unions references, so one transcript can arrive twice.
            var protein = WithReferences(
                Ensembl("ENST00000000001.1", ("gene ID", "ENSG00000000001.1")),
                Ensembl("ENST00000000001.1", ("gene ID", "ENSG00000000001.1")));

            Assert.That(protein.EnsemblGeneReferences.Count, Is.EqualTo(1));
        }

        [Test]
        public void EnsemblGeneReferences_ReferencesWithoutATranscriptId_AreNotRepeatsOfEachOther()
        {
            var protein = WithReferences(
                Ensembl(null, ("gene ID", "ENSG00000000001.1")),
                Ensembl(null, ("gene ID", "ENSG00000000002.1")));

            Assert.That(protein.EnsemblGeneIds, Is.EqualTo(new[] { "ENSG00000000001", "ENSG00000000002" }),
                "no transcript id is not a shared transcript id, so the second gene link is not dropped");
        }

        [Test]
        public void Decoys_CarryNoEnsemblGeneReferences()
        {
            var decoy = ProteinDbLoader
                .LoadProteinXML(Data("DatabaseTests", "humanGAPDH.xml"), true, DecoyType.Reverse, null, false, null, out _)
                .First(p => p.IsDecoy);

            Assert.That(decoy.EnsemblGeneReferences, Is.Empty, "a decoy must never resolve to a real gene");
        }

        [Test]
        public void EnsemblGeneReferences_EmptyRatherThanNullWhenNothingIsAnnotated()
        {
            var protein = new Protein("PEPTIDEK", "P00001");

            Assert.That(protein.EnsemblGeneReferences, Is.Not.Null.And.Empty);
            Assert.That(protein.EnsemblGeneIds, Is.Not.Null.And.Empty);
        }

        [Test]
        public void EnsemblDatabaseReferenceType_IsReExportedFromTheLoader()
        {
            Assert.That(ProteinDbLoader.EnsemblDatabaseReferenceType, Is.EqualTo(Protein.EnsemblDatabaseReferenceType));
            Assert.That(Protein.EnsemblDatabaseReferenceType, Is.EqualTo("Ensembl"));
        }
    }
}
