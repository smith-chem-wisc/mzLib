using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests.ExternalFileReading
{
    /// <summary>
    /// The test file is a 37-row excerpt of a DIA-NN 1.9 long-format report.tsv from
    /// MSV000095360 (PBMC DIA on an Orbitrap Astral). Rows were chosen to cover the cases the
    /// reader has to get right: three runs, charges 1 through 3, carbamidomethylated peptides,
    /// protein groups holding one and four accessions, and groups with no gene name.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestDiaNnResultFiles
    {
        private const string DiaNnReportPath = @"FileReadingTests\ExternalFileTypes\DiaNn_LongFormat_report.tsv";

        private static string TestFilePath =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, DiaNnReportPath);

        private string _outputDirectory;

        [OneTimeSetUp]
        public void SetUp()
        {
            _outputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\DiaNnReadWriteTests");
            Directory.CreateDirectory(_outputDirectory);
        }

        [OneTimeTearDown]
        public void TearDown()
        {
            Directory.Delete(_outputDirectory, true);
        }

        [Test]
        public void TestDiaNnReportLoadsAndCountCorrect()
        {
            DiaNnReportFile file = new DiaNnReportFile(TestFilePath);

            Assert.That(file.Count(), Is.EqualTo(37));
            Assert.That(file.CanRead(TestFilePath));
            Assert.That(file.FileType, Is.EqualTo(SupportedFileType.DiaNnReport));
            Assert.That(file.Software, Is.EqualTo(Software.DiaNn));
        }

        /// <summary>
        /// The report is recognized by its header rather than its name, so a file whose name gives
        /// away nothing must still be identified as a DIA-NN report.
        /// </summary>
        [Test]
        public void TestDiaNnReportIsDetectedByHeaderRegardlessOfFileName()
        {
            string renamedPath = Path.Combine(_outputDirectory, "some_unrelated_name.tsv");
            File.Copy(TestFilePath, renamedPath, true);

            Assert.That(renamedPath.ParseFileType(), Is.EqualTo(SupportedFileType.DiaNnReport));

            IQuantifiableResultFile quantifiable = FileReader.ReadQuantifiableResultFile(renamedPath);
            Assert.That(quantifiable, Is.TypeOf<DiaNnReportFile>());
            Assert.That(quantifiable.GetQuantifiableResults().Count(), Is.EqualTo(37));
        }

        [Test]
        public void TestDiaNnReportFirstRecordCorrect()
        {
            DiaNnPrecursor first = new DiaNnReportFile(TestFilePath).First();

            Assert.That(first.SpectraFilePath, Is.EqualTo(
                "/common/hendricksn/iphronesis/platform_workspace/jobs/SJOB12004/thermo_converter/20240316_PBMCs_24min_3Th-1_centroid.mzML"));
            Assert.That(first.Run, Is.EqualTo("20240316_PBMCs_24min_3Th-1_centroid"));
            Assert.That(first.ProteinGroup, Is.EqualTo("Q8NFD5"));
            Assert.That(first.ProteinIds, Is.EqualTo("Q8NFD5;Q9Y651"));
            Assert.That(first.ProteinNames, Is.EqualTo("ARI1B_HUMAN"));
            Assert.That(first.Genes, Is.EqualTo("ARID1B"));
            Assert.That(first.ProteinGroupQuantity, Is.EqualTo(20851.74609));
            Assert.That(first.ProteinGroupNormalized, Is.EqualTo(22475.39062));
            Assert.That(first.ProteinGroupMaxLfq, Is.EqualTo(24653.61523));
            Assert.That(first.GenesQuantity, Is.EqualTo(20851.74609));
            Assert.That(first.ModifiedSequence, Is.EqualTo("AAAAAAAAAAR"));
            Assert.That(first.BaseSequence, Is.EqualTo("AAAAAAAAAAR"));
            Assert.That(first.PrecursorId, Is.EqualTo("AAAAAAAAAAR2"));
            Assert.That(first.PrecursorCharge, Is.EqualTo(2));
            Assert.That(first.QValue, Is.EqualTo(0.0001727389172));
            Assert.That(first.PosteriorErrorProbability, Is.EqualTo(0.003596759867));
            Assert.That(first.GlobalQValue, Is.EqualTo(0.0717228204));
            Assert.That(first.ProteinQValue, Is.EqualTo(0.0001707067277));
            Assert.That(first.ProteinGroupQValue, Is.EqualTo(0.0001623640273));
            Assert.That(first.Proteotypic, Is.False);
            Assert.That(first.PrecursorQuantity, Is.EqualTo(5359.077148));
            Assert.That(first.PrecursorNormalized, Is.EqualTo(5263.249023));
            Assert.That(first.RetentionTime, Is.EqualTo(2.962105513));
            Assert.That(first.RetentionTimeStart, Is.EqualTo(2.865492821));
            Assert.That(first.RetentionTimeStop, Is.EqualTo(3.059067249));
            Assert.That(first.IndexedRetentionTime, Is.EqualTo(-34.94962311));
            Assert.That(first.PredictedRetentionTime, Is.EqualTo(2.874304056));
            Assert.That(first.Ms1Area, Is.EqualTo(173254.0625));
            Assert.That(first.CScore, Is.EqualTo(0.997564137));
            Assert.That(first.DecoyCScore, Is.EqualTo(0.05423271656));
            Assert.That(first.Ms2ScanNumber, Is.EqualTo(18291));
            Assert.That(first.PrecursorMz, Is.EqualTo(443.2495117));
            Assert.That(first.LibraryIndex, Is.EqualTo(26));
            Assert.That(first.IonMobility, Is.EqualTo(0));
        }

        /// <summary>
        /// The per-fragment columns are 40% of the file's bytes and nothing in the quantification
        /// path reads them, so they are skipped unless asked for.
        /// </summary>
        [Test]
        public void TestFragmentColumnsAreSkippedUnlessRequested()
        {
            DiaNnPrecursor withoutFragments = new DiaNnReportFile(TestFilePath).First();
            Assert.That(withoutFragments.FragmentInfo, Is.Null);
            Assert.That(withoutFragments.FragmentQuantRaw, Is.Null);
            Assert.That(withoutFragments.FragmentQuantCorrected, Is.Null);
            Assert.That(withoutFragments.FragmentCorrelations, Is.Null);

            // Everything else is unaffected
            Assert.That(withoutFragments.PrecursorId, Is.EqualTo("AAAAAAAAAAR2"));
            Assert.That(withoutFragments.Ms2ScanNumber, Is.EqualTo(18291));

            DiaNnPrecursor withFragments = new DiaNnReportFile(TestFilePath, readFragmentColumns: true).First();
            Assert.That(withFragments.FragmentInfo, Does.StartWith("y7^1/601.3427124;y6^1/530.305603;"));
            Assert.That(withFragments.FragmentQuantRaw, Does.StartWith("709.9885254;1161.098877;"));
            Assert.That(withFragments.FragmentCorrelations, Does.StartWith("0.4798017144;0.6663590074;"));
            Assert.That(withFragments.PrecursorId, Is.EqualTo("AAAAAAAAAAR2"));
        }

        /// <summary>
        /// The strings that repeat on every row -- run name, spectra path, protein columns, and the
        /// sequences -- are pooled so the same text is stored once rather than once per row.
        /// </summary>
        [Test]
        public void TestRepeatedStringsAreShared()
        {
            var records = new DiaNnReportFile(TestFilePath)
                .Where(precursor => precursor.Run == "20240316_PBMCs_24min_3Th-1_centroid")
                .ToList();

            Assert.That(records.Count, Is.GreaterThan(1));
            Assert.That(records.All(r => ReferenceEquals(r.Run, records[0].Run)));
            Assert.That(records.All(r => ReferenceEquals(r.SpectraFilePath, records[0].SpectraFilePath)));
        }

        [Test]
        public void TestDiaNnReportLastRecordCorrect()
        {
            DiaNnPrecursor last = new DiaNnReportFile(TestFilePath).Last();

            Assert.That(last.Run, Is.EqualTo("20240316_PBMCs_24min_3Th-1_centroid"));
            Assert.That(last.ProteinGroup, Is.EqualTo("Q6ZSR9"));
            Assert.That(last.ProteinIds, Is.EqualTo("Q6ZSR9"));
            Assert.That(last.ProteinNames, Is.EqualTo("YJ005_HUMAN"));
            Assert.That(last.Genes, Is.Empty);
            Assert.That(last.ModifiedSequence, Is.EqualTo("EAVSGPMAGKPFRPQSLSK"));
            Assert.That(last.PrecursorId, Is.EqualTo("EAVSGPMAGKPFRPQSLSK3"));
            Assert.That(last.PrecursorCharge, Is.EqualTo(3));
            Assert.That(last.QValue, Is.EqualTo(0.008654053323));
            Assert.That(last.Proteotypic, Is.True);
            Assert.That(last.RetentionTime, Is.EqualTo(9.183587074));
            Assert.That(last.Ms2ScanNumber, Is.EqualTo(56935));
            Assert.That(last.PrecursorMz, Is.EqualTo(663.0198975));
            Assert.That(last.DecoyCScore, Is.EqualTo(-10000000));
        }

        /// <summary>
        /// The Genes column and its three aggregate columns are blank whenever DIA-NN's database
        /// has no gene for the group. Those must read back as null rather than as a quantity of zero,
        /// since a blank means "not reported" and zero would be a real measurement.
        /// </summary>
        [Test]
        public void TestBlankNumericColumnsReadAsNull()
        {
            DiaNnPrecursor last = new DiaNnReportFile(TestFilePath).Last();

            Assert.That(last.GenesQuantity, Is.Null);
            Assert.That(last.GenesNormalized, Is.Null);
            Assert.That(last.GenesMaxLfq, Is.Null);
            Assert.That(last.GenesMaxLfqUnique, Is.Null);

            // Translated.Quality is blank on every row of this report
            Assert.That(new DiaNnReportFile(TestFilePath).All(p => p.TranslatedQuality is null));
        }

        /// <summary>
        /// DIA-NN writes modifications as UniMod accessions. Those are passed through verbatim,
        /// which keeps the record lossless; the stripped sequence supplies the base sequence.
        /// </summary>
        [Test]
        public void TestModifiedSequenceIsPassedThroughVerbatim()
        {
            DiaNnPrecursor modified = new DiaNnReportFile(TestFilePath)
                .First(p => p.BaseSequence == "AAAACLDK");

            Assert.That(modified.ModifiedSequence, Is.EqualTo("AAAAC(UniMod:4)LDK"));
            Assert.That(modified.FullSequence, Is.EqualTo("AAAAC(UniMod:4)LDK"));
            Assert.That(modified.BaseSequence, Is.EqualTo("AAAACLDK"));
        }

        [Test]
        public void TestIQuantifiableRecordImplementation()
        {
            DiaNnPrecursor first = new DiaNnReportFile(TestFilePath).First();

            // The Run column is used as the file name, not the File.Name path, which points at
            // wherever the search was run rather than at the caller's copy of the spectra file.
            Assert.That(first.FileName, Is.EqualTo("20240316_PBMCs_24min_3Th-1_centroid"));
            Assert.That(first.FullSequence, Is.EqualTo("AAAAAAAAAAR"));
            Assert.That(first.ChargeState, Is.EqualTo(2));
            Assert.That(first.IsDecoy, Is.False);

            // DIA-NN reports RT in minutes, so IQuantifiableRecord.RetentionTime needs no conversion
            Assert.That(first.RetentionTime, Is.EqualTo(2.962105513));

            // DIA-NN reports no peptide mass, so it comes from the precursor m/z and charge:
            // 2 * 443.2495117 - 2 * 1.007276466879
            Assert.That(first.MonoisotopicMass, Is.EqualTo(884.4844705).Within(1e-6));
        }

        /// <summary>
        /// Protein.Ids, Genes, and Protein.Names are independently sorted sets rather than parallel
        /// lists, so a gene is only claimed when the group names exactly one.
        /// </summary>
        [Test]
        public void TestProteinGroupInfosParsing()
        {
            var allResults = new DiaNnReportFile(TestFilePath).ToList();

            // One accession, one gene: unambiguous, so everything is populated
            DiaNnPrecursor oneProtein = allResults.First(p => p.BaseSequence == "AAAAAAAAAVSR");
            Assert.That(oneProtein.ProteinGroupInfos.Count, Is.EqualTo(1));
            Assert.That(oneProtein.ProteinGroupInfos.Single(), Is.EqualTo(("Q96JP5", "ZFP91", "HUMAN")));

            // Two accessions but one reported gene: that gene belongs to one of them and the report
            // does not say which (Q8NFD5 is ARID1B, Q9Y651 is not), so neither claims it. The
            // organism is the same either way and survives.
            DiaNnPrecursor twoAccessionsOneGene = allResults.First(p => p.BaseSequence == "AAAAAAAAAAR");
            Assert.That(twoAccessionsOneGene.Genes, Is.EqualTo("ARID1B"));
            Assert.That(twoAccessionsOneGene.ProteinGroupInfos, Is.EqualTo(new List<(string, string, string)>
            {
                ("Q8NFD5", "", "HUMAN"),
                ("Q9Y651", "", "HUMAN"),
            }));

            // No gene reported at all; accession and organism still come through
            DiaNnPrecursor noGene = allResults.First(p => p.BaseSequence == "EAVSGPMAGKPFRPQSLSK");
            Assert.That(noGene.ProteinGroupInfos.Single(), Is.EqualTo(("Q6ZSR9", "", "HUMAN")));

            // Every record yields at least one protein group, which FlashLFQ requires
            Assert.That(allResults.All(p => p.ProteinGroupInfos.Count > 0));
        }

        /// <summary>
        /// The three protein columns are independently sorted, de-duplicated sets, so zipping them
        /// by index invents pairings that are simply wrong. In the report this snippet came from,
        /// group "P06753;P07951;P09493;P67936" carries genes "TPM1;TPM2;TPM3;TPM4", but that
        /// report's own single-accession rows prove P06753 is TPM3 and P09493 is TPM1 -- so a
        /// positional zip gets both backwards. Rather than guess, no gene is claimed.
        /// </summary>
        [Test]
        public void TestNoGeneIsClaimedWhenTheGroupNamesSeveral()
        {
            var allResults = new DiaNnReportFile(TestFilePath).ToList();

            DiaNnPrecursor fourProteins = allResults.First(p => p.BaseSequence == "AADESER");
            Assert.That(fourProteins.ProteinIds, Is.EqualTo("P06753;P07951;P09493;P67936"));
            Assert.That(fourProteins.Genes, Is.EqualTo("TPM1;TPM2;TPM3;TPM4"));

            Assert.That(fourProteins.ProteinGroupInfos, Is.EqualTo(new List<(string, string, string)>
            {
                ("P06753", "", "HUMAN"),
                ("P07951", "", "HUMAN"),
                ("P09493", "", "HUMAN"),
                ("P67936", "", "HUMAN"),
            }));

            // Specifically, the two that a positional zip would transpose
            Assert.That(fourProteins.ProteinGroupInfos.Single(p => p.proteinAccessions == "P06753").geneName,
                Is.Not.EqualTo("TPM1"));
            Assert.That(fourProteins.ProteinGroupInfos.Single(p => p.proteinAccessions == "P09493").geneName,
                Is.Not.EqualTo("TPM3"));

            // A two-accession group with two different genes is the same situation
            DiaNnPrecursor twoGenes = allResults.First(p => p.BaseSequence == "AAAVADK");
            Assert.That(twoGenes.Genes, Is.EqualTo("ITPRID2;TADA3"));
            Assert.That(twoGenes.ProteinGroupInfos.All(p => p.geneName == ""));
            Assert.That(twoGenes.ProteinGroupInfos.Select(p => p.proteinAccessions),
                Is.EqualTo(new[] { "O75528", "P28290" }));
        }

        [Test]
        public void TestFileNameToFilePathExactMatch()
        {
            DiaNnReportFile file = new DiaNnReportFile(TestFilePath);

            List<string> fullFilePaths = new()
            {
                @"D:\Data\20240315_PBMCs_24min_3Th-1_centroid.mzML",
                @"D:\Data\20240316_PBMCs_24min_3Th-1_centroid.mzML",
                @"D:\Data\20240317_PBMCs_24min_3Th-1_centroid.mzML",
                @"D:\Data\SomethingUnrelated.mzML",
            };

            Dictionary<string, string> matched = file.FileNameToFilePath(fullFilePaths);

            Assert.That(matched.Count, Is.EqualTo(3));
            Assert.That(matched["20240315_PBMCs_24min_3Th-1_centroid"],
                Is.EqualTo(@"D:\Data\20240315_PBMCs_24min_3Th-1_centroid.mzML"));
            Assert.That(matched["20240316_PBMCs_24min_3Th-1_centroid"],
                Is.EqualTo(@"D:\Data\20240316_PBMCs_24min_3Th-1_centroid.mzML"));
            Assert.That(matched["20240317_PBMCs_24min_3Th-1_centroid"],
                Is.EqualTo(@"D:\Data\20240317_PBMCs_24min_3Th-1_centroid.mzML"));
        }

        /// <summary>
        /// DIA-NN names a run after whatever file it was handed, which is routinely a converted
        /// copy of the file the caller still has. Here the caller holds the original .raw files
        /// while DIA-NN searched "_centroid" mzMLs, so matching has to tolerate the suffix.
        /// </summary>
        [Test]
        public void TestFileNameToFilePathMatchesConvertedFileNames()
        {
            DiaNnReportFile file = new DiaNnReportFile(TestFilePath);

            List<string> fullFilePaths = new()
            {
                @"D:\Data\20240315_PBMCs_24min_3Th-1.raw",
                @"D:\Data\20240316_PBMCs_24min_3Th-1.raw",
                @"D:\Data\20240317_PBMCs_24min_3Th-1.raw",
            };

            Dictionary<string, string> matched = file.FileNameToFilePath(fullFilePaths);

            Assert.That(matched.Count, Is.EqualTo(3));
            Assert.That(matched["20240316_PBMCs_24min_3Th-1_centroid"],
                Is.EqualTo(@"D:\Data\20240316_PBMCs_24min_3Th-1.raw"));
        }

        /// <summary>
        /// A run with no corresponding spectra file is left out of the dictionary rather than
        /// mapped to an arbitrary file.
        /// </summary>
        [Test]
        public void TestFileNameToFilePathOmitsUnmatchedRuns()
        {
            DiaNnReportFile file = new DiaNnReportFile(TestFilePath);

            Dictionary<string, string> matched = file.FileNameToFilePath(new List<string>
            {
                @"D:\Data\20240316_PBMCs_24min_3Th-1_centroid.mzML",
            });

            Assert.That(matched.Count, Is.EqualTo(1));
            Assert.That(matched.ContainsKey("20240315_PBMCs_24min_3Th-1_centroid"), Is.False);

            Assert.That(file.FileNameToFilePath(new List<string>()), Is.Empty);
        }

        /// <summary>
        /// A path that strips to an empty name is a prefix of every run name, so it must not be
        /// allowed to match everything.
        /// </summary>
        [Test]
        public void TestFileNameToFilePathIgnoresEmptyFileNames()
        {
            DiaNnReportFile file = new DiaNnReportFile(TestFilePath);

            Dictionary<string, string> matched = file.FileNameToFilePath(new List<string> { @"D:\Data\" });

            Assert.That(matched, Is.Empty);
        }

        /// <summary>
        /// A fuzzy match has to stop at a delimiter. Plain containment lets the file "A1" match the
        /// run "A10", and 96-well plate naming produces exactly that collision -- which would
        /// quietly attribute two runs' peptides to one spectra file.
        /// </summary>
        [Test]
        public void TestFileNameToFilePathDoesNotMatchAcrossANameBoundary()
        {
            DiaNnReportFile plateFile = new DiaNnReportFile(TestFilePath)
            {
                Results = new List<DiaNnPrecursor>
                {
                    new() { Run = "A1", ProteinIds = "P1", PrecursorCharge = 2 },
                    new() { Run = "A10", ProteinIds = "P1", PrecursorCharge = 2 },
                }
            };

            Dictionary<string, string> matched = plateFile.FileNameToFilePath(new List<string> { @"D:\d\A10.raw" });

            Assert.That(matched.Count, Is.EqualTo(1));
            Assert.That(matched["A10"], Is.EqualTo(@"D:\d\A10.raw"));
            Assert.That(matched.ContainsKey("A1"), Is.False);
        }

        /// <summary>
        /// Two runs must never share one spectra file: FlashLFQ would merge their identifications
        /// and quantify both as a single file. Dropping the contested pair makes the problem surface
        /// as a named error instead of a silently wrong result.
        /// </summary>
        [Test]
        public void TestFileNameToFilePathDropsRunsThatContestTheSameFile()
        {
            DiaNnReportFile contested = new DiaNnReportFile(TestFilePath)
            {
                Results = new List<DiaNnPrecursor>
                {
                    new() { Run = "sample_a", ProteinIds = "P1", PrecursorCharge = 2 },
                    new() { Run = "sample_b", ProteinIds = "P1", PrecursorCharge = 2 },
                }
            };

            Assert.That(contested.FileNameToFilePath(new List<string> { @"D:\d\sample.raw" }), Is.Empty);
        }

        /// <summary>
        /// Two equally good candidates give no reason to prefer either, and picking by list order
        /// would make the result depend on how the caller happened to sort its spectra files.
        /// </summary>
        [Test]
        public void TestFileNameToFilePathLeavesAmbiguousRunsUnmatched()
        {
            DiaNnReportFile ambiguous = new DiaNnReportFile(TestFilePath)
            {
                Results = new List<DiaNnPrecursor>
                {
                    new() { Run = "sample", ProteinIds = "P1", PrecursorCharge = 2 },
                }
            };

            var forward = ambiguous.FileNameToFilePath(new List<string> { @"D:\d\sample_rep1.raw", @"D:\d\sample_rep2.raw" });
            var reversed = ambiguous.FileNameToFilePath(new List<string> { @"D:\d\sample_rep2.raw", @"D:\d\sample_rep1.raw" });

            Assert.That(forward, Is.Empty);
            Assert.That(forward, Is.EqualTo(reversed));
        }

        /// <summary>
        /// DIA-NN writes report.pr_matrix.tsv into the same directory as the main report, and it
        /// carries Precursor.Id and Stripped.Sequence too -- but one column per run instead of
        /// File.Name, so it is not readable as a long-format report.
        /// </summary>
        [Test]
        public void TestMatrixReportsAreNotMistakenForTheMainReport()
        {
            string prMatrixPath = Path.Combine(_outputDirectory, "sample.pr_matrix.tsv");
            File.WriteAllLines(prMatrixPath, new[]
            {
                "Protein.Group\tProtein.Ids\tProtein.Names\tGenes\tFirst.Protein.Description\tProteotypic\t" +
                "Stripped.Sequence\tModified.Sequence\tPrecursor.Charge\tPrecursor.Id\t/data/run1.mzML",
                "Q8NFD5\tQ8NFD5\tARI1B_HUMAN\tARID1B\tdescription\t1\tAAAAAAAAAAR\tAAAAAAAAAAR\t2\tAAAAAAAAAAR2\t123.4",
            });

            Assert.Throws<MzLibException>(() => prMatrixPath.ParseFileType());
        }

        /// <summary>
        /// DIA-NN commonly runs on Linux and writes '/'-delimited paths, which .NET's path helpers
        /// would not split when the tests run on Windows.
        /// </summary>
        [Test]
        [TestCase("/common/jobs/SJOB12004/20240316_PBMCs_centroid.mzML", "20240316_PBMCs_centroid")]
        [TestCase(@"D:\Data\20240316_PBMCs_centroid.mzML", "20240316_PBMCs_centroid")]
        [TestCase("20240316_PBMCs_centroid.mzML", "20240316_PBMCs_centroid")]
        [TestCase("20240316_PBMCs_centroid", "20240316_PBMCs_centroid")]
        [TestCase("", "")]
        public void TestStripPathAndExtension(string path, string expected)
        {
            Assert.That(DiaNnPrecursor.StripPathAndExtension(path), Is.EqualTo(expected));
        }

        /// <summary>
        /// The full path through FlashLFQ's adapter: a DIA-NN report plus the spectra files it
        /// refers to should produce one Identification per row.
        /// </summary>
        [Test]
        public void TestMakeIdentifications()
        {
            DiaNnReportFile file = new DiaNnReportFile(TestFilePath);

            List<SpectraFileInfo> spectraFiles = new()
            {
                new SpectraFileInfo(@"D:\Data\20240315_PBMCs_24min_3Th-1_centroid.mzML", "PBMC", 0, 0, 0),
                new SpectraFileInfo(@"D:\Data\20240316_PBMCs_24min_3Th-1_centroid.mzML", "PBMC", 1, 0, 0),
                new SpectraFileInfo(@"D:\Data\20240317_PBMCs_24min_3Th-1_centroid.mzML", "PBMC", 2, 0, 0),
            };

            List<Identification> identifications = MzLibExtensions.MakeIdentifications(file, spectraFiles);

            Assert.That(identifications.Count, Is.EqualTo(37));

            Identification first = identifications[0];
            Assert.That(first.BaseSequence, Is.EqualTo("AAAAAAAAAAR"));
            Assert.That(first.ModifiedSequence, Is.EqualTo("AAAAAAAAAAR"));
            Assert.That(first.Ms2RetentionTimeInMinutes, Is.EqualTo(2.962105513));
            Assert.That(first.MonoisotopicMass, Is.EqualTo(884.4844705).Within(1e-6));
            Assert.That(first.PrecursorChargeState, Is.EqualTo(2));
            Assert.That(first.FileInfo.FullFilePathWithExtension,
                Is.EqualTo(@"D:\Data\20240316_PBMCs_24min_3Th-1_centroid.mzML"));
            Assert.That(first.ProteinGroups.Count, Is.EqualTo(2));

            ProteinGroup proteinGroup = first.ProteinGroups.First(pg => pg.ProteinGroupName == "Q8NFD5");
            Assert.That(proteinGroup.Organism, Is.EqualTo("HUMAN"));

            // A single-accession group does carry its gene through to FlashLFQ
            Identification unambiguous = identifications.First(id => id.BaseSequence == "AAAAAAAAAVSR");
            ProteinGroup zfp91 = unambiguous.ProteinGroups.Single();
            Assert.That(zfp91.ProteinGroupName, Is.EqualTo("Q96JP5"));
            Assert.That(zfp91.GeneName, Is.EqualTo("ZFP91"));

            // The four-accession tropomyosin group carries all four proteins through
            Identification fourProteins = identifications.First(id => id.BaseSequence == "AADESER");
            Assert.That(fourProteins.ProteinGroups.Count, Is.EqualTo(4));

            // Decoys are not reported as rows in DIA-NN's main report, so everything is usable
            // for protein quantification
            Assert.That(identifications.All(id => !id.IsDecoy));
            Assert.That(identifications.All(id => id.UseForProteinQuant));
        }

        /// <summary>
        /// FlashLFQ gates MBR donor selection on Identification.QValue, so leaving it at the default
        /// zero would make every DIA-NN precursor a top-confidence donor. The value has to be
        /// Global.Q.Value rather than Q.Value: DIA-NN scopes Q.Value to a single run and has already
        /// filtered it below 1%, so it carries no information across runs.
        /// </summary>
        [Test]
        public void TestQValuesReachIdentifications()
        {
            DiaNnReportFile file = new DiaNnReportFile(TestFilePath);

            List<SpectraFileInfo> spectraFiles = new()
            {
                new SpectraFileInfo(@"D:\Data\20240315_PBMCs_24min_3Th-1_centroid.mzML", "PBMC", 0, 0, 0),
                new SpectraFileInfo(@"D:\Data\20240316_PBMCs_24min_3Th-1_centroid.mzML", "PBMC", 1, 0, 0),
                new SpectraFileInfo(@"D:\Data\20240317_PBMCs_24min_3Th-1_centroid.mzML", "PBMC", 2, 0, 0),
            };

            DiaNnPrecursor firstRecord = file.First();
            Assert.That(firstRecord.GlobalQValue, Is.EqualTo(0.0717228204));
            Assert.That(firstRecord.PosteriorErrorProbability, Is.EqualTo(0.003596759867));
            Assert.That(firstRecord.CScore, Is.EqualTo(0.997564137));

            Identification first = MzLibExtensions.MakeIdentifications(file, spectraFiles)[0];
            Assert.That(first.QValue, Is.EqualTo(0.0717228204));
            Assert.That(first.PsmScore, Is.EqualTo(0.997564137));

            // usePepQValue swaps in the posterior error probability instead
            Identification firstUsingPep = MzLibExtensions.MakeIdentifications(file, spectraFiles, usePepQValue: true)[0];
            Assert.That(firstUsingPep.QValue, Is.EqualTo(0.003596759867));
            Assert.That(firstUsingPep.PsmScore, Is.EqualTo(0.997564137));

            // Not everything collapses to zero, and this row would have passed a 1% donor gate it
            // does not deserve: its run-specific Q.Value is 0.00017 but experiment-wide it is 7.2%
            Assert.That(firstRecord.QValue, Is.LessThan(0.01));
            Assert.That(first.QValue, Is.GreaterThan(0.01));

            var allIdentifications = MzLibExtensions.MakeIdentifications(file, spectraFiles);
            Assert.That(allIdentifications.Any(id => id.QValue > 0));
            Assert.That(allIdentifications.All(id => id.PsmScore > 0));
        }

        /// <summary>
        /// Round-tripping the fragment columns requires asking for them; that is the case this
        /// covers, since a file loaded without them writes them back empty by design.
        /// </summary>
        [Test]
        public void TestDiaNnReportReadWriteRoundTrip()
        {
            DiaNnReportFile original = new DiaNnReportFile(TestFilePath, readFragmentColumns: true);
            string outputPath = Path.Combine(_outputDirectory, "roundTrip_report.tsv");

            original.WriteResults(outputPath);
            Assert.That(File.Exists(outputPath));

            DiaNnReportFile rewritten = new DiaNnReportFile(outputPath, readFragmentColumns: true);
            Assert.That(rewritten.Count(), Is.EqualTo(original.Count()));

            foreach (var (before, after) in original.Zip(rewritten))
            {
                Assert.That(after.Run, Is.EqualTo(before.Run));
                Assert.That(after.PrecursorId, Is.EqualTo(before.PrecursorId));
                Assert.That(after.ModifiedSequence, Is.EqualTo(before.ModifiedSequence));
                Assert.That(after.BaseSequence, Is.EqualTo(before.BaseSequence));
                Assert.That(after.ProteinIds, Is.EqualTo(before.ProteinIds));
                Assert.That(after.Genes, Is.EqualTo(before.Genes));
                Assert.That(after.PrecursorCharge, Is.EqualTo(before.PrecursorCharge));
                Assert.That(after.PrecursorMz, Is.EqualTo(before.PrecursorMz));
                Assert.That(after.RetentionTime, Is.EqualTo(before.RetentionTime));
                Assert.That(after.QValue, Is.EqualTo(before.QValue));
                Assert.That(after.PrecursorQuantity, Is.EqualTo(before.PrecursorQuantity));
                Assert.That(after.GenesQuantity, Is.EqualTo(before.GenesQuantity));
                Assert.That(after.Ms2ScanNumber, Is.EqualTo(before.Ms2ScanNumber));
                Assert.That(after.FragmentInfo, Is.EqualTo(before.FragmentInfo));
                Assert.That(after.Proteotypic, Is.EqualTo(before.Proteotypic));
            }

            // The file we wrote is still a DIA-NN report as far as the reader is concerned
            Assert.That(outputPath.ParseFileType(), Is.EqualTo(SupportedFileType.DiaNnReport));
        }

        /// <summary>
        /// Proteotypic is written by DIA-NN as 0 or 1. CsvHelper's default boolean converter reads
        /// that correctly but writes "False"/"True", which reads back fine here while quietly
        /// breaking every other tool that consumes the file.
        /// </summary>
        [Test]
        public void TestProteotypicRoundTripsAsZeroOrOne()
        {
            DiaNnReportFile original = new DiaNnReportFile(TestFilePath, readFragmentColumns: true);
            string outputPath = Path.Combine(_outputDirectory, "proteotypic_report.tsv");
            original.WriteResults(outputPath);

            string[] originalLines = File.ReadAllLines(TestFilePath);
            string[] writtenLines = File.ReadAllLines(outputPath);
            int column = Array.IndexOf(originalLines[0].Split('\t'), "Proteotypic");

            Assert.That(column, Is.GreaterThanOrEqualTo(0));
            for (int row = 1; row < originalLines.Length; row++)
            {
                Assert.That(writtenLines[row].Split('\t')[column],
                    Is.EqualTo(originalLines[row].Split('\t')[column]));
            }
        }

        /// <summary>
        /// WriteResults appends the extension when the destination does not already carry it, but
        /// must not do so for a path that is already a .tsv. This type's "extension" is the whole
        /// conventional name "report.tsv", so the usual CanRead check would treat a perfectly good
        /// "myResults.tsv" as unnamed and turn it into "myResults.tsvreport.tsv".
        /// </summary>
        [Test]
        public void TestWriteResultsAppendsExtensionOnlyWhenMissing()
        {
            DiaNnReportFile original = new DiaNnReportFile(TestFilePath);

            string noExtension = Path.Combine(_outputDirectory, "noExtension");
            original.WriteResults(noExtension);
            Assert.That(File.Exists(noExtension + SupportedFileType.DiaNnReport.GetFileExtension()));

            string alreadyTsv = Path.Combine(_outputDirectory, "myResults.tsv");
            original.WriteResults(alreadyTsv);
            Assert.That(File.Exists(alreadyTsv));
            Assert.That(File.Exists(alreadyTsv + SupportedFileType.DiaNnReport.GetFileExtension()), Is.False);
        }
    }
}
