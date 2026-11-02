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
        #region Public presets

        [Test]
        public void EnumValuesDoNotShift()
        {
            Assert.That((int)FastaHeaderType.UniProt, Is.EqualTo(0));
            Assert.That((int)FastaHeaderType.Ensembl, Is.EqualTo(1));
            Assert.That((int)FastaHeaderType.Gencode, Is.EqualTo(2));
            Assert.That((int)FastaHeaderType.Unknown, Is.EqualTo(3));
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

        // Verbatim copy of DetectFastaHeaderFormat at 564e7286.
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

        #endregion
    }
}
