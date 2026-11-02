using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.IO.Compression;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using UsefulProteomicsDatabases;

namespace Test.DatabaseTests
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestFastaHeaderPresets
    {
        // Deflines from NCBI efetch and rest.uniprot.org (2026-11-01). The gi| and ref| forms are the legacy layout
        // built from each record's gi number (esummary) and title. Each constant here is the thrL leader peptide.
        private const string ThrL = "MKRISTTITTTITITTGNGAG";
        private const string RefSeqDefline = ">NP_414542.1 thr operon leader peptide [Escherichia coli str. K-12 substr. MG1655]";
        private const string RefSeqPipedDefline = ">ref|NP_414542.1| thr operon leader peptide [Escherichia coli str. K-12 substr. MG1655]";
        private const string RefSeqGiDefline = ">gi|16127995|ref|NP_414542.1| thr operon leader peptide [Escherichia coli str. K-12 substr. MG1655]";
        private const string GenBankDefline = ">AAA24672.1 thrL [Escherichia coli]";
        private const string GenBankGiDefline = ">gi|290477|gb|AAA24672.1| thrL [Escherichia coli]";
        private const string GenBankSevenDigitDefline = ">CAN8157071.1 thr operon leader peptide [Escherichia coli]";
        private const string MultispeciesDefline = ">WP_001386572.1 MULTISPECIES: thr operon leader peptide [Bacteria]";
        private const string UniProtDefline = ">sp|P0AD86|LPT_ECOLI thr operon leader peptide OS=Escherichia coli (strain K12) OX=83333 GN=thrL PE=3 SV=1";

        private string _tempDir;

        [SetUp]
        public void SetUp()
        {
            _tempDir = Path.Combine(TestContext.CurrentContext.WorkDirectory, "TestFastaHeaderPresets_" + Guid.NewGuid().ToString("N"));
            Directory.CreateDirectory(_tempDir);
        }

        [TearDown]
        public void TearDown()
        {
            if (Directory.Exists(_tempDir))
                Directory.Delete(_tempDir, true);
        }

        private string WriteFasta(params string[] deflines)
        {
            string path = Path.Combine(_tempDir, Guid.NewGuid().ToString("N") + ".fasta");
            File.WriteAllLines(path, deflines.SelectMany(d => new[] { d, ThrL }));
            return path;
        }

        #region Public presets

        [Test]
        public void EnumValuesDoNotShift()
        {
            Assert.That((int)FastaHeaderType.UniProt, Is.EqualTo(0));
            Assert.That((int)FastaHeaderType.Ensembl, Is.EqualTo(1));
            Assert.That((int)FastaHeaderType.Gencode, Is.EqualTo(2));
            Assert.That((int)FastaHeaderType.Unknown, Is.EqualTo(3));
            Assert.That((int)FastaHeaderType.Ncbi, Is.EqualTo(4));
        }

        [Test]
        public void GetFieldRegexesReturnsTheNamedRegexesForEachFormat()
        {
            var uniprot = ProteinDbLoader.GetFieldRegexes(FastaHeaderType.UniProt);
            Assert.That(uniprot, Is.SameAs(ProteinDbLoader.UniprotFieldRegexes));
            Assert.That(uniprot.Accession, Is.SameAs(ProteinDbLoader.UniprotAccessionRegex));
            Assert.That(uniprot.FullName, Is.SameAs(ProteinDbLoader.UniprotFullNameRegex));
            Assert.That(uniprot.Name, Is.SameAs(ProteinDbLoader.UniprotNameRegex));
            Assert.That(uniprot.GeneName, Is.SameAs(ProteinDbLoader.UniprotGeneNameRegex));
            Assert.That(uniprot.Organism, Is.SameAs(ProteinDbLoader.UniprotOrganismRegex));
            Assert.That(uniprot.OrganismId, Is.SameAs(ProteinDbLoader.UniprotOrganismIdRegex));

            var ensembl = ProteinDbLoader.GetFieldRegexes(FastaHeaderType.Ensembl);
            Assert.That(ensembl, Is.SameAs(ProteinDbLoader.EnsemblFieldRegexes));
            Assert.That(ensembl, Is.EqualTo(new FastaHeaderFieldRegexBundle(ProteinDbLoader.EnsemblAccessionRegex,
                ProteinDbLoader.EnsemblFullNameRegex, null, ProteinDbLoader.EnsemblGeneNameRegex, null, null)));

            var gencode = ProteinDbLoader.GetFieldRegexes(FastaHeaderType.Gencode);
            Assert.That(gencode, Is.SameAs(ProteinDbLoader.GencodeFieldRegexes));
            Assert.That(gencode, Is.EqualTo(new FastaHeaderFieldRegexBundle(ProteinDbLoader.GencodeAccessionRegex,
                ProteinDbLoader.GencodeFullNameRegex, null, ProteinDbLoader.GencodeGeneNameRegex, null, null)));

            var ncbi = ProteinDbLoader.GetFieldRegexes(FastaHeaderType.Ncbi);
            Assert.That(ncbi, Is.SameAs(ProteinDbLoader.NcbiFieldRegexes));
            Assert.That(ncbi, Is.EqualTo(new FastaHeaderFieldRegexBundle(ProteinDbLoader.NcbiAccessionRegex,
                ProteinDbLoader.NcbiFullNameRegex, null, null, ProteinDbLoader.NcbiOrganismRegex, null)));
        }

        [TestCase(FastaHeaderType.Unknown)]
        [TestCase((FastaHeaderType)99)]
        public void GetFieldRegexesRefusesATypeWithNoPreset(FastaHeaderType headerType)
        {
            Assert.That(() => ProteinDbLoader.GetFieldRegexes(headerType), Throws.TypeOf<ArgumentOutOfRangeException>());
        }

        #endregion

        #region Differential: the public presets reproduce the switch they replaced, on every FASTA fixture

        // git ls-files Test | grep -iE '\.(fasta|fa|faa)(\.gz)?$' at 564e7286.
        private static readonly string[] FastaFixtures =
        {
            "DataFiles/humanInsulin.fasta",
            "DatabaseTests/Isoform.fasta",
            "DatabaseTests/Isoform.fasta.gz",
            "DatabaseTests/TwoHumanHistone_mimic_Retain3.fasta",
            "DatabaseTests/bad.fasta",
            "DatabaseTests/bad3.fasta",
            "DatabaseTests/bad4.fasta",
            "DatabaseTests/blank.fasta",
            "DatabaseTests/ensembl_prrc2a.fa",
            "DatabaseTests/fasta.fasta",
            "DatabaseTests/gencode_mmp20.fa",
            "DatabaseTests/isoformTest.fasta",
            "DatabaseTests/ncbi_gi_legacy.fasta",
            "DatabaseTests/test_ensembl.pep.all.fasta",
            "DatabaseTests/uniprot_aifm1.fasta",
            "Transcriptomics/TestData/20mer1.fasta",
            "Transcriptomics/TestData/20mer1.fasta.gz",
            "Transcriptomics/TestData/Burkmers_ModomicsSequences.fasta",
            "Transcriptomics/TestData/Burkmers_PrimarySequence.fasta",
            "Transcriptomics/TestData/ModomicsUnmodifiedTrimmed.fasta",
            "Transcriptomics/TestData/ModomicsUnmodifiedTrimmed.fasta.gz",
            "Transcriptomics/TestData/TRNAdbUnmodifiedTrimmed.fasta",
            "Transcriptomics/TestData/TestDatabase_Ensembl.GRCh38.ncrna.fa",
        };

        // Caller-supplied regexes for the non-accession fields, which the presets must keep or replace exactly as before.
        private static readonly FastaHeaderFieldRegex CallerFullName = new("fullName", @"^>(\S{1,5})", 0, 1);
        private static readonly FastaHeaderFieldRegex CallerName = new("name", @"^>(\S{1,4})", 0, 1);
        private static readonly FastaHeaderFieldRegex CallerGeneName = new("geneName", @"^>(\S{1,3})", 0, 1);
        private static readonly FastaHeaderFieldRegex CallerOrganism = new("organism", @"(\S+)$", 0, 1);
        private static readonly FastaHeaderFieldRegex CallerOrganismId = new("organismId", @"(\d{3})", 0, 1);

        private static IEnumerable<TestCaseData> FixtureCases() =>
            FastaFixtures.SelectMany(f => new[]
            {
                new TestCaseData(f, "none").SetArgDisplayNames(f, "none"),
                new TestCaseData(f, "nameOrganismOrganismId").SetArgDisplayNames(f, "nameOrganismOrganismId"),
                new TestCaseData(f, "allButAccession").SetArgDisplayNames(f, "allButAccession"),
            });

        private static string FixturePath(string relative) =>
            Path.Combine(new[] { TestContext.CurrentContext.TestDirectory }.Concat(relative.Split('/')).ToArray());

        private static List<string> HeaderLines(string path)
        {
            using Stream file = File.OpenRead(path);
            using Stream stream = path.EndsWith(".gz") ? new GZipStream(file, CompressionMode.Decompress) : file;
            using var reader = new StreamReader(stream);
            var headers = new List<string>();
            for (string line = reader.ReadLine(); line != null; line = reader.ReadLine())
                if (line.StartsWith(">"))
                    headers.Add(line);
            return headers;
        }

        // Verbatim copy of DetectFastaHeaderFormat before Ncbi was added (564e7286).
        private static FastaHeaderType LegacyDetectFastaHeaderFormat(string line)
        {
            if (!line.StartsWith(">"))
            {
                return FastaHeaderType.Unknown;
            }
            if (line.Length >= 4 && char.IsLower(line[1]) && char.IsLower(line[2]) && line[3] == '|')
            {
                return FastaHeaderType.UniProt;
            }
            if (line.Contains("gene_symbol:", StringComparison.InvariantCultureIgnoreCase)
                || line.Contains("gene:", StringComparison.InvariantCultureIgnoreCase)
                || line.Contains("Acc:", StringComparison.InvariantCultureIgnoreCase)
                || line.Contains("description:", StringComparison.InvariantCultureIgnoreCase))
            {
                return FastaHeaderType.Ensembl;
            }
            if (line.Count(p => p == '|') == 7)
            {
                return FastaHeaderType.Gencode;
            }
            return FastaHeaderType.Unknown;
        }

        // Verbatim copy of the switch LoadProteinFasta used before GetFieldRegexes (564e7286, ProteinDbLoader.cs:311-333).
        private static void LegacySwitch(FastaHeaderType? HeaderType, ref FastaHeaderFieldRegex accessionRegex,
            ref FastaHeaderFieldRegex fullNameRegex, ref FastaHeaderFieldRegex nameRegex, ref FastaHeaderFieldRegex geneNameRegex,
            ref FastaHeaderFieldRegex organismRegex, ref FastaHeaderFieldRegex organismIdRegex)
        {
            switch (HeaderType)
            {
                case FastaHeaderType.UniProt:
                    accessionRegex = ProteinDbLoader.UniprotAccessionRegex;
                    fullNameRegex = ProteinDbLoader.UniprotFullNameRegex;
                    nameRegex = ProteinDbLoader.UniprotNameRegex;
                    organismRegex = ProteinDbLoader.UniprotOrganismRegex;
                    organismIdRegex ??= ProteinDbLoader.UniprotOrganismIdRegex;
                    geneNameRegex = ProteinDbLoader.UniprotGeneNameRegex;
                    break;

                case FastaHeaderType.Ensembl:
                    accessionRegex = ProteinDbLoader.EnsemblAccessionRegex;
                    fullNameRegex = ProteinDbLoader.EnsemblFullNameRegex;
                    geneNameRegex = ProteinDbLoader.EnsemblGeneNameRegex;
                    break;

                case FastaHeaderType.Gencode:
                    accessionRegex = ProteinDbLoader.GencodeAccessionRegex;
                    fullNameRegex = ProteinDbLoader.GencodeFullNameRegex;
                    geneNameRegex = ProteinDbLoader.GencodeGeneNameRegex;
                    break;
            }
        }

        private static string Describe(Protein p) => string.Join("\t",
            p.Accession ?? "<null>", p.Name ?? "<null>", p.FullName ?? "<null>", p.Organism ?? "<null>",
            string.Join(";", p.GeneNames.Select(g => g.Item1 + "=" + g.Item2)),
            string.Join(";", (p.DatabaseReferences ?? new List<DatabaseReference>()).Select(d => d.Type + "=" + d.Id)),
            p.IsDecoy, p.IsContaminant, p.IsEntrapment, p.BaseSequence);

        [TestCaseSource(nameof(FixtureCases))]
        public void PresetsParseEveryFixtureAsTheReplacedSwitchDid(string fixture, string callerFields)
        {
            string path = FixturePath(fixture);
            Assert.That(File.Exists(path), Is.True, "fixture not copied to the test output: " + fixture);

            FastaHeaderFieldRegex fullName = null, name = null, geneName = null, organism = null, organismId = null;
            if (callerFields != "none")
            {
                name = CallerName;
                organism = CallerOrganism;
                organismId = CallerOrganismId;
            }
            if (callerFields == "allButAccession")
            {
                fullName = CallerFullName;
                geneName = CallerGeneName;
            }

            List<string> headers = HeaderLines(path);
            if (headers.Count > 0 && LegacyDetectFastaHeaderFormat(headers[0]) == FastaHeaderType.Unknown)
            {
                Assert.That(() => ProteinDbLoader.LoadProteinFasta(path, true, DecoyType.Reverse, false, out _,
                        null, fullName, name, geneName, organism, organismIdRegex: organismId),
                    Throws.TypeOf<MzLibException>().With.Message.EqualTo("Unknown fasta header format: " + headers[0]));
                return;
            }

            // Old arm: resolve the regexes with the legacy switch and pass them explicitly, which bypasses detection.
            FastaHeaderFieldRegex legacyAccession = null, legacyFullName = fullName, legacyName = name,
                legacyGeneName = geneName, legacyOrganism = organism, legacyOrganismId = organismId;
            FastaHeaderType? detected = headers.Count > 0 ? LegacyDetectFastaHeaderFormat(headers[0]) : null;
            LegacySwitch(detected, ref legacyAccession, ref legacyFullName, ref legacyName, ref legacyGeneName,
                ref legacyOrganism, ref legacyOrganismId);

            var expected = ProteinDbLoader.LoadProteinFasta(path, true, DecoyType.Reverse, false, out var expectedErrors,
                legacyAccession, legacyFullName, legacyName, legacyGeneName, legacyOrganism, organismIdRegex: legacyOrganismId);
            var actual = ProteinDbLoader.LoadProteinFasta(path, true, DecoyType.Reverse, false, out var actualErrors,
                null, fullName, name, geneName, organism, organismIdRegex: organismId);

            Assert.That(actual.Select(Describe), Is.EqualTo(expected.Select(Describe)));
            Assert.That(actualErrors, Is.EqualTo(expectedErrors));
        }

        [TestCaseSource(nameof(FastaFixtures))]
        public void NcbiDetectionTakesNoFixtureHeaderFromAnotherFormat(string fixture)
        {
            foreach (string header in HeaderLines(FixturePath(fixture)))
            {
                FastaHeaderType before = LegacyDetectFastaHeaderFormat(header);
                FastaHeaderType after = ProteinDbLoader.DetectFastaHeaderFormat(header);
                if (before == FastaHeaderType.Unknown && after == FastaHeaderType.Ncbi)
                    continue;
                Assert.That(after, Is.EqualTo(before), header);
            }
        }

        [Test]
        public void TheOnlyFixtureHeaderNowDetectedAsNcbiIsTheRefOnlyLine()
        {
            var ncbi = FastaFixtures.SelectMany(f => HeaderLines(FixturePath(f)))
                .Where(h => ProteinDbLoader.DetectFastaHeaderFormat(h) == FastaHeaderType.Ncbi);
            Assert.That(ncbi, Is.EqualTo(new[] { ">ref|NP_414557.1| homoserine kinase [Escherichia coli str. K-12 substr. MG1655]" }));
        }

        #endregion

        #region NCBI

        [TestCase(RefSeqDefline)]
        [TestCase(RefSeqPipedDefline)]
        [TestCase(GenBankDefline)]
        [TestCase(GenBankSevenDigitDefline)]
        [TestCase(MultispeciesDefline)]
        [TestCase(">XP_011541469.1 protein Shroom1 isoform X1 [Homo sapiens]")]
        [TestCase(">CAA23472.1 unnamed protein product (chloroplast) [Nicotiana tabacum]")]
        [TestCase(">NP_414542.1")]
        public void NcbiDeflinesAreDetected(string defline)
        {
            Assert.That(ProteinDbLoader.DetectFastaHeaderFormat(defline), Is.EqualTo(FastaHeaderType.Ncbi));
        }

        // gi| and gb| begin with two lowercase letters and a pipe, which the UniProt check claims first.
        [TestCase(RefSeqGiDefline, FastaHeaderType.UniProt)]
        [TestCase(GenBankGiDefline, FastaHeaderType.UniProt)]
        [TestCase(UniProtDefline, FastaHeaderType.UniProt)]
        [TestCase(">NP_414542 thr operon leader peptide [Escherichia coli str. K-12 substr. MG1655]", FastaHeaderType.Unknown)]
        [TestCase(">lcl|U00096.3_prot_AAC73112.1_1 [gene=thrL] [locus_tag=b0001] [db_xref=UniProtKB/Swiss-Prot:P0AD86] [protein=thr operon leader peptide] [protein_id=AAC73112.1] [location=190..255] [gbkey=CDS]", FastaHeaderType.Unknown)]
        [TestCase(">tdbR00000016", FastaHeaderType.Unknown)]
        public void OtherDeflinesAreNotDetectedAsNcbi(string defline, FastaHeaderType expected)
        {
            Assert.That(ProteinDbLoader.DetectFastaHeaderFormat(defline), Is.EqualTo(expected));
        }

        // Synthetic: an NCBI accession followed by Ensembl keys, or in the first of eight Gencode fields, keeps the earlier format.
        [TestCase(">NP_414542.1 pep:known gene:thrL", FastaHeaderType.Ensembl)]
        [TestCase(">NP_414542.1|t|g|-|-|thrL-201|thrL|21", FastaHeaderType.Gencode)]
        public void NcbiDetectionIsTriedAfterTheExistingFormats(string defline, FastaHeaderType expected)
        {
            Assert.That(ProteinDbLoader.DetectFastaHeaderFormat(defline), Is.EqualTo(expected));
        }

        [TestCase(RefSeqDefline, "NP_414542.1", "thr operon leader peptide", "Escherichia coli str. K-12 substr. MG1655")]
        [TestCase(RefSeqPipedDefline, "NP_414542.1", "thr operon leader peptide", "Escherichia coli str. K-12 substr. MG1655")]
        [TestCase(RefSeqGiDefline, "NP_414542.1", "thr operon leader peptide", "Escherichia coli str. K-12 substr. MG1655")]
        [TestCase(GenBankDefline, "AAA24672.1", "thrL", "Escherichia coli")]
        [TestCase(GenBankGiDefline, "AAA24672.1", "thrL", "Escherichia coli")]
        [TestCase(GenBankSevenDigitDefline, "CAN8157071.1", "thr operon leader peptide", "Escherichia coli")]
        [TestCase(MultispeciesDefline, "WP_001386572.1", "MULTISPECIES: thr operon leader peptide", "Bacteria")]
        [TestCase(">CAA23472.1 unnamed protein product (chloroplast) [Nicotiana tabacum]", "CAA23472.1", "unnamed protein product (chloroplast)", "Nicotiana tabacum")]
        [TestCase(">NP_414542.1", "NP_414542.1", null, null)]
        public void NcbiPresetParsesEachDeflineShape(string defline, string accession, string fullName, string organism)
        {
            var ncbi = ProteinDbLoader.NcbiFieldRegexes;
            Assert.That(ncbi.Accession.ApplyRegex(defline), Is.EqualTo(accession));
            Assert.That(ncbi.FullName.ApplyRegex(defline), Is.EqualTo(fullName));
            Assert.That(ncbi.Organism.ApplyRegex(defline), Is.EqualTo(organism));
        }

        [Test]
        public void AnNcbiFastaThatUsedToThrowNowLoadsByAutoDetection()
        {
            string path = WriteFasta(RefSeqDefline, GenBankDefline, GenBankSevenDigitDefline, MultispeciesDefline);

            var proteins = ProteinDbLoader.LoadProteinFasta(path, true, DecoyType.None, false, out var errors);

            Assert.That(proteins.Select(p => p.Accession),
                Is.EqualTo(new[] { "NP_414542.1", "AAA24672.1", "CAN8157071.1", "WP_001386572.1" }));
            Assert.That(proteins.Select(p => p.FullName),
                Is.EqualTo(new[] { "thr operon leader peptide", "thrL", "thr operon leader peptide", "MULTISPECIES: thr operon leader peptide" }));
            Assert.That(proteins.Select(p => p.Organism),
                Is.EqualTo(new[] { "Escherichia coli str. K-12 substr. MG1655", "Escherichia coli", "Escherichia coli", "Bacteria" }));
            Assert.That(proteins.Select(p => p.Name), Is.All.Null);
            Assert.That(proteins.SelectMany(p => p.GeneNames), Is.Empty);
            Assert.That(proteins.Select(p => p.BaseSequence), Is.All.EqualTo(ThrL));
            Assert.That(errors, Is.Empty);
        }

        [Test]
        public void ARefPipedFirstHeaderAutoDetectsAsNcbi()
        {
            string path = WriteFasta(RefSeqPipedDefline);

            var protein = ProteinDbLoader.LoadProteinFasta(path, true, DecoyType.None, false, out var errors).Single();

            Assert.That(protein.Accession, Is.EqualTo("NP_414542.1"));
            Assert.That(protein.FullName, Is.EqualTo("thr operon leader peptide"));
            Assert.That(protein.Organism, Is.EqualTo("Escherichia coli str. K-12 substr. MG1655"));
            Assert.That(errors, Is.Empty);
        }

        [Test]
        public void TheNcbiPresetKeepsTheCallersRegexesForFieldsItDoesNotParse()
        {
            string path = WriteFasta(RefSeqDefline);

            var protein = ProteinDbLoader.LoadProteinFasta(path, true, DecoyType.None, false, out _,
                nameRegex: CallerName, geneNameRegex: CallerGeneName, organismIdRegex: CallerOrganismId).Single();

            Assert.That(protein.Name, Is.EqualTo("NP_4"));
            Assert.That(protein.GeneNames.Single().Item2, Is.EqualTo("NP_"));
            Assert.That(protein.DatabaseReferences.Single().Id, Is.EqualTo("414"));
        }

        [Test]
        public void TheNcbiPresetParsesLegacyGiDeflinesWhenChosenExplicitly()
        {
            string path = WriteFasta(RefSeqGiDefline, GenBankGiDefline);
            var ncbi = ProteinDbLoader.NcbiFieldRegexes;

            var proteins = ProteinDbLoader.LoadProteinFasta(path, true, DecoyType.None, false, out var errors,
                ncbi.Accession, ncbi.FullName, ncbi.Name, ncbi.GeneName, ncbi.Organism, organismIdRegex: ncbi.OrganismId);

            Assert.That(proteins.Select(p => p.Accession), Is.EqualTo(new[] { "NP_414542.1", "AAA24672.1" }));
            Assert.That(proteins.Select(p => p.FullName), Is.EqualTo(new[] { "thr operon leader peptide", "thrL" }));
            Assert.That(proteins.Select(p => p.Organism), Is.EqualTo(new[] { "Escherichia coli str. K-12 substr. MG1655", "Escherichia coli" }));
            Assert.That(errors, Is.Empty);
        }

        #endregion
    }
}
