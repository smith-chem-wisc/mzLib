using System;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.IO.Compression;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using Readers.ExternalResults.IndividualResultRecords;
using Readers.ExternalResults.ResultFiles;

namespace Test.FileReadingTests.ExternalFileReading
{
    /// <summary>
    /// MzIdentMLResultFile over the PRIDE captures in DataFiles (see PRIDE_MZID_PROVENANCE.md) and the published
    /// HUPO-PSI examples. Expected values were read out of each file's XML, not off the reader.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class TestMzIdentMLResultFile
    {
        private static string DataFile(string name) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", name);

        private static MzIdentMLResultFile Read(string name) => new(DataFile(name));

        private string _outputDirectory;

        [OneTimeSetUp]
        public void OneTimeSetUp()
        {
            _outputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "MzIdentMLResultFileTests");
            Directory.CreateDirectory(_outputDirectory);
        }

        [OneTimeTearDown]
        public void OneTimeTearDown()
        {
            if (Directory.Exists(_outputDirectory))
            {
                Directory.Delete(_outputDirectory, true);
            }
        }

        /// <summary>
        /// Every SpectrumIdentificationItem in the file is either a record or a skipped match, and the only items
        /// skipped here are crosslinks (MS:1002511), counted in the XML. multiple_spectra_per_id has three
        /// SpectrumIdentificationLists, so its count also pins that every list is read.
        /// </summary>
        [TestCase("SmallCalibratible_Yeast.mzID", 73, 0)]
        [TestCase("PXD078927_msgf_1_1_0.mzid", 12, 0)]
        [TestCase("PXD078927_msgf_1_1_0.mzid.gz", 12, 0)]
        [TestCase("PXD000783_scaffold_1_1_0.mzid", 3, 0)]
        [TestCase("PXD019591_mascotparser_1_1_0.mzid", 7, 0)]
        [TestCase("PXD019733_proteomediscoverer_1_1_0.mzid", 3, 0)]
        [TestCase("PXD000783_scaffold_mods_1_1_0.mzid", 6, 0)]
        [TestCase("PXD019591_mascotparser_mods_1_1_0.mzid", 15, 0)]
        [TestCase("PXD070193_xifdr_1_2_0.mzid", 0, 4)]
        [TestCase("OpenxQuest_example_1_2_0.mzid", 0, 16)]
        [TestCase("mzidLib_xtandem_fdr_1_2_0.mzid", 5, 0)]
        [TestCase("multiple_spectra_per_id_1_3_0.mzid", 4, 4)]
        [TestCase("noncovalently_assoc_1_3_0.mzid", 2, 0)]
        public void EveryItemIsARecordOrASkippedCrosslink(string fileName, int records, int crosslinks)
        {
            var file = Read(fileName);

            Assert.Multiple(() =>
            {
                Assert.That(file.Results, Has.Count.EqualTo(records));
                Assert.That(file.SkippedMatches, Has.Count.EqualTo(crosslinks));
                Assert.That(file.SkippedMatches.Select(s => s.Reason), Is.All.EqualTo("crosslink identification"));
            });
        }

        [Test]
        public void LoadsAndReportsTheFileTypeAndSoftware()
        {
            var plain = Read("PXD078927_msgf_1_1_0.mzid");
            var compressed = Read("PXD078927_msgf_1_1_0.mzid.gz");

            Assert.Multiple(() =>
            {
                Assert.That(plain.FileType, Is.EqualTo(SupportedFileType.MzIdentML));
                Assert.That(compressed.FileType, Is.EqualTo(SupportedFileType.MzIdentMLGz));
                Assert.That(plain.Software, Is.EqualTo(Software.MzIdentML));
                Assert.That(new MzIdentMLResultFile().Software, Is.EqualTo(Software.MzIdentML));
                Assert.That(plain.CanRead(plain.FilePath), Is.True);
                Assert.That(compressed.CanRead(compressed.FilePath), Is.True);
            });
        }

        /// <summary>
        /// The gzip file is the plain fixture compressed, so the two must read identically; the compressed one is
        /// streamed through GZipStream rather than written out first.
        /// </summary>
        [Test]
        public void CompressedFileReadsTheSameRecordsAsThePlainFile()
        {
            var plain = Read("PXD078927_msgf_1_1_0.mzid").Results;
            var compressed = Read("PXD078927_msgf_1_1_0.mzid.gz").Results;

            Assert.That(compressed.Select(r => (r.SpectrumId, r.FullSequence, r.Rank, r.QValue, r.Accession)),
                Is.EqualTo(plain.Select(r => (r.SpectrumId, r.FullSequence, r.Rank, r.QValue, r.Accession))));
        }

        /// <summary>
        /// MS-GF+ on a Thermo mzML: the scan number comes from the Thermo nativeID, the q-value from
        /// MS:1002054 MS-GF:QValue (a child of MS:1002354), and a peptide shared between two antigens reports both.
        /// </summary>
        [Test]
        public void MsgfRecord_ThermoNativeIdChildQValueAndSharedPeptide()
        {
            var records = Read("PXD078927_msgf_1_1_0.mzid").Results;
            var first = records[0];
            var second = records[1];

            Assert.Multiple(() =>
            {
                Assert.That(first.SpectrumId, Is.EqualTo("controllerType=0 controllerNumber=1 scan=14316"));
                Assert.That(first.OneBasedScanNumber, Is.EqualTo(14316));
                Assert.That(first.FileNameWithoutExtension, Is.EqualTo("MHC1_Pep"));
                Assert.That(first.SpectraFileLocation, Is.EqualTo("/mnt/c/Users/USER/Documents/LCMS_data/mzML_file/3rd_repeat/MHC1_Pep.mzML"));
                Assert.That(first.BaseSequence, Is.EqualTo("HSNLNDATYQRT"));
                Assert.That(first.FullSequence, Is.EqualTo("HSNLNDATYQRT"));
                Assert.That(first.AllModsOneIsNterminus, Is.Empty);
                Assert.That(first.Accession, Is.EqualTo("sp|P2_H1N1|Antigen_P2_H1N1|sp|P5_H3N2|Antigen_P5_H3N2"));
                Assert.That(first.ChargeState, Is.EqualTo(2));
                Assert.That(first.Rank, Is.EqualTo(1));
                Assert.That(first.PassThreshold, Is.True);
                Assert.That(first.QValue, Is.EqualTo(0.0).Within(1e-12));
                Assert.That(first.Scores["MS-GF:RawScore"], Is.EqualTo(115));
                Assert.That(first.Scores["MS-GF:SpecEValue"], Is.EqualTo(3.041116E-14).Within(1e-20));

                // same spectrum, second-ranked peptide
                Assert.That(second.OneBasedScanNumber, Is.EqualTo(14316));
                Assert.That(second.Rank, Is.EqualTo(2));
                Assert.That(second.QValue, Is.EqualTo(0.6588785).Within(1e-9));

                Assert.That(records.Count(r => r.IsDecoy), Is.EqualTo(6));
                Assert.That(records.All(r => r.QValue.HasValue), Is.True);
            });
        }

        /// <summary>
        /// Proteome Discoverer's spectrumID "scan=N file=K" names the RAW file's SpectraData by id, although the
        /// result itself references the merged mzML. The record reports the RAW file, 1K_4D.raw (id 454).
        /// </summary>
        [Test]
        public void ProteomeDiscovererRecord_FileIdInTheNativeIdNamesTheRawFile()
        {
            var records = Read("PXD019733_proteomediscoverer_1_1_0.mzid").Results;

            Assert.Multiple(() =>
            {
                Assert.That(records[0].SpectrumId, Is.EqualTo("scan=81409 file=454"));
                Assert.That(records[0].OneBasedScanNumber, Is.EqualTo(81409));
                Assert.That(records[0].FileNameWithoutExtension, Is.EqualTo("1K_4D"));
                Assert.That(records[0].SpectraFileLocation, Does.EndWith(@"\1K_4D.raw"));
                Assert.That(records[2].SpectrumId, Is.EqualTo("scan=81596 file=455"));
                Assert.That(records[2].FileNameWithoutExtension, Is.EqualTo("1K_4E"));
                Assert.That(records[0].QValue, Is.Null);
                Assert.That(records[0].Scores["Proteome Discoverer Delta Score"], Is.EqualTo(0.0106).Within(1e-12));
            });
        }

        /// <summary>
        /// Mascot Parser on an MGF: "index=16332" is the zero-based position in the peak list, so the record
        /// reports 16333, which is also the query number the title starts with. The instrument scan (22014) is
        /// only in the title, and is not parsed out of it. A Windows UNC location still yields a file name.
        /// </summary>
        [Test]
        public void MascotParserRecord_IndexNativeIdIsAOneBasedPeakListPosition()
        {
            var first = Read("PXD019591_mascotparser_1_1_0.mzid").Results[0];

            Assert.Multiple(() =>
            {
                Assert.That(first.SpectrumId, Is.EqualTo("index=16332"));
                Assert.That(first.OneBasedScanNumber, Is.EqualTo(16333));
                Assert.That(first.SpectrumTitle, Does.StartWith("16333: Scan 22014 "));
                Assert.That(first.FileNameWithoutExtension, Is.EqualTo("QE2_PharmacoDB_TP_Sample_16B.raw.-1"));
                Assert.That(first.PassThreshold, Is.False);
                Assert.That(first.Scores["Mascot:score"], Is.EqualTo(11.94).Within(1e-9));
                Assert.That(first.QValue, Is.Null);
            });
        }

        /// <summary>
        /// Scaffold names the spectra file with backslashes. The location must split on them on any OS.
        /// </summary>
        [Test]
        public void ScaffoldRecord_WindowsLocationSplitsOnAnyOperatingSystem()
        {
            var first = Read("PXD000783_scaffold_1_1_0.mzid").Results[0];

            Assert.Multiple(() =>
            {
                Assert.That(first.OneBasedScanNumber, Is.EqualTo(2));
                Assert.That(first.FileNameWithoutExtension, Is.EqualTo("mascot_daemon_merge_F008897.mzid_mascot_daemon_merge_F008897"));
                Assert.That(first.Scores["Scaffold:Peptide Probability"], Is.EqualTo(0.92933995).Within(1e-9));
            });
        }

        /// <summary>
        /// The Scaffold/Mascot capture keeps the first result carrying each modification kind. Locations are
        /// mzIdentML's (0 the N-terminus, 1..length a residue), keys are AllModsOneIsNterminus's (1 the
        /// N-terminus, residue r at r + 2 counting from zero). The N-terminal acetyl has no residues attribute
        /// and must resolve to Unimod's N-terminal entry, not a side-chain acetylation.
        /// </summary>
        [Test]
        public void ScaffoldModifications_ResolveToUnimodAtTheRightPosition()
        {
            var records = Read("PXD000783_scaffold_mods_1_1_0.mzid").Results;
            string Mods(MzIdentMLRecord r) =>
                string.Join("; ", r.AllModsOneIsNterminus.OrderBy(m => m.Key).Select(m => $"{m.Key}:{m.Value.IdWithMotif}"));

            Assert.Multiple(() =>
            {
                Assert.That(Mods(records[0]), Is.EqualTo("4:Phospho on S"));           // SLS(3)FGGR
                Assert.That(Mods(records[1]), Is.EqualTo("2:Carbamidomethyl on C"));   // C(1)SHLGLPIQGK
                Assert.That(Mods(records[2]), Is.EqualTo("7:Phospho on T; 11:Oxidation on M")); // AAAAAT(6)PSPM(10)VK
                Assert.That(Mods(records[3]), Is.EqualTo("2:Phospho on S; 4:Deamidated on Q"));  // S(1)IQ(3)LIR
                Assert.That(Mods(records[5]), Is.EqualTo("1:Acetyl on X; 2:Phospho on T"));     // N-term, T(1)AILER

                Assert.That(records[5].AllModsOneIsNterminus[1].LocationRestriction, Does.Contain("N-terminal"));
                Assert.That(records[5].BaseSequence, Is.EqualTo("TAILER"));
                Assert.That(records[5].FullSequence, Is.EqualTo("[Unimod:Acetyl on X]T[Unimod:Phospho on T]AILER"));
                Assert.That(records[2].FullSequence, Is.EqualTo("AAAAAT[Unimod:Phospho on T]PSPM[Unimod:Oxidation on M]VK"));

                // two peptides for one spectrum, both kept
                Assert.That(records[3].SpectrumId, Is.EqualTo(records[4].SpectrumId));
            });
        }

        /// <summary>
        /// Mascot Parser writes pyro-glutamate at location 0 on a peptide starting with Q; Unimod has it only as a
        /// Peptide N-terminal modification on Q, which is what it must resolve to.
        /// </summary>
        [Test]
        public void MascotParserModifications_NTerminalPyroGluResolvesToTheTerminalEntry()
        {
            var records = Read("PXD019591_mascotparser_mods_1_1_0.mzid").Results;
            var pyroGlu = records.Single(r => r.BaseSequence == "QGITKSAPLR");
            var carbamidomethyl = records.Single(r => r.BaseSequence == "SSHAVELACR");

            Assert.Multiple(() =>
            {
                Assert.That(pyroGlu.AllModsOneIsNterminus.Keys, Is.EquivalentTo(new[] { 1, 5 }));
                Assert.That(pyroGlu.AllModsOneIsNterminus[1].IdWithMotif, Is.EqualTo("Gln->pyro-Glu on Q"));
                Assert.That(pyroGlu.AllModsOneIsNterminus[1].LocationRestriction, Is.EqualTo("Peptide N-terminal."));
                Assert.That(pyroGlu.AllModsOneIsNterminus[5].IdWithMotif, Is.EqualTo("Phospho on T"));
                Assert.That(pyroGlu.Rank, Is.EqualTo(2));

                Assert.That(carbamidomethyl.AllModsOneIsNterminus.Keys, Is.EquivalentTo(new[] { 10 }));
                Assert.That(carbamidomethyl.PassThreshold, Is.True);

                // ranks 1..5 of one spectrum are all records
                Assert.That(records.Where(r => r.SpectrumId == "index=12903").Select(r => r.Rank), Is.EqualTo(new[] { 1, 2, 3, 4, 5 }));
            });
        }

        [Test]
        public void SkippedMatchesAreAvailableWithoutTouchingResults()
        {
            var file = Read("PXD070193_xifdr_1_2_0.mzid");

            Assert.That(file.SkippedMatches.Select(s => s.SpectrumIdentificationItemId),
                Is.EqualTo(new[] { "SII_10_1", "SII_10_2", "SII_56_1", "SII_56_2" }));
            Assert.That(file.SkippedMatches[0].SpectrumId, Is.EqualTo("index=5896"));
        }

        /// <summary>
        /// ResultFile reloads whenever Results is empty, so a file whose every item is skipped is read once per
        /// path rather than on every access. Replacing the file with garbage after the first read proves it.
        /// </summary>
        [Test]
        public void AFileWithNoRecordsIsNotReadAgainOnEveryAccess()
        {
            string copy = Path.Combine(_outputDirectory, "all_crosslinks.mzid");
            File.Copy(DataFile("PXD070193_xifdr_1_2_0.mzid"), copy, true);
            var file = new MzIdentMLResultFile(copy);
            Assert.That(file.Results, Is.Empty);

            File.WriteAllText(copy, "not xml");

            Assert.That(() => file.Results, Throws.Nothing);
            Assert.That(file.SkippedMatches, Has.Count.EqualTo(4));
        }

        [Test]
        public void TheFactoryReadsBothExtensions()
        {
            var plain = FileReader.ReadResultFile(DataFile("PXD078927_msgf_1_1_0.mzid"));
            var compressed = FileReader.ReadResultFile(DataFile("PXD078927_msgf_1_1_0.mzid.gz"));

            Assert.That(plain, Is.TypeOf<MzIdentMLResultFile>());
            Assert.That(compressed, Is.TypeOf<MzIdentMLResultFile>());
            Assert.That(((MzIdentMLResultFile)compressed).Results, Has.Count.EqualTo(12));
        }

        [Test]
        public void WritingIsNotSupported()
        {
            var file = Read("PXD078927_msgf_1_1_0.mzid");

            Assert.That(() => file.WriteResults(Path.Combine(_outputDirectory, "out.mzid")),
                Throws.TypeOf<NotSupportedException>());
        }

        [TestCase("not xml at all", "not_xml.mzid")]
        [TestCase("<?xml version=\"1.0\"?><Unrelated xmlns=\"urn:x\"/>", "wrong_root.mzid")]
        public void AnUnreadableDocumentThrowsMzLibExceptionNamingTheFile(string contents, string fileName)
        {
            string path = Path.Combine(_outputDirectory, fileName);
            File.WriteAllText(path, contents);

            var e = Assert.Throws<MzLibException>(() => _ = new MzIdentMLResultFile(path).Results);
            Assert.That(e!.Message, Does.StartWith($"Could not read mzIdentML file '{path}'"));
        }

        [Test]
        public void ACorruptGzipThrowsMzLibException()
        {
            string path = Path.Combine(_outputDirectory, "corrupt.mzid.gz");
            var bytes = File.ReadAllBytes(DataFile("PXD078927_msgf_1_1_0.mzid.gz"));
            File.WriteAllBytes(path, bytes.Take(bytes.Length / 2).Concat(new byte[64]).ToArray());

            Assert.Throws<MzLibException>(() => _ = new MzIdentMLResultFile(path).Results);
        }

        /// <summary>
        /// Results returns nothing for a missing file (ResultFile only loads a file that exists), so only a
        /// direct LoadResults sees it, and sees it as a missing file rather than an unreadable one.
        /// </summary>
        [TestCase("missing.mzid")]
        [TestCase("missing.mzid.gz")]
        public void LoadingAMissingFileThrowsFileNotFoundException(string fileName)
        {
            string path = Path.Combine(_outputDirectory, fileName);

            Assert.Throws<FileNotFoundException>(() => new MzIdentMLResultFile(path).LoadResults());
        }

        /// <summary>
        /// An I/O fault on a file that exists is a read failure like any other, so it names the file. A file
        /// another process holds open without sharing stands in for a truncated or yanked one.
        /// </summary>
        [TestCase("PXD078927_msgf_1_1_0.mzid")]
        [TestCase("PXD078927_msgf_1_1_0.mzid.gz")]
        public void AnIOFaultThrowsMzLibExceptionNamingTheFile(string fileName)
        {
            string path = Path.Combine(_outputDirectory, "locked_" + fileName);
            File.Copy(DataFile(fileName), path, overwrite: true);

            using (new FileStream(path, FileMode.Open, FileAccess.ReadWrite, FileShare.None))
            {
                var e = Assert.Throws<MzLibException>(() => _ = new MzIdentMLResultFile(path).Results);
                Assert.That(e!.Message, Does.StartWith($"Could not read mzIdentML file '{path}'"));
                Assert.That(e.InnerException, Is.InstanceOf<IOException>());
            }
        }

        /// <summary>
        /// A plain .mzid is read from disk, not copied into memory first as a .gz must be. The document is
        /// padded with comments the deserializer skips, so reading it allocates far less than its size unless
        /// the whole file is buffered.
        /// </summary>
        [Test]
        public void APlainFileIsNotBufferedWholeIntoMemory()
        {
            string source = File.ReadAllText(DataFile("PXD078927_msgf_1_1_0.mzid"));
            int body = source.IndexOf("<MzIdentML", StringComparison.Ordinal);
            string padding = string.Concat(Enumerable.Repeat("<!--" + new string('x', 1000) + "-->\n", 40_000));
            string path = Path.Combine(_outputDirectory, "padded.mzid");
            File.WriteAllText(path, source[..body] + padding + source[body..]);
            long size = new FileInfo(path).Length;

            _ = Read("PXD078927_msgf_1_1_0.mzid").Results; // serializer and Unimod warm-up
            long before = GC.GetAllocatedBytesForCurrentThread();
            var file = new MzIdentMLResultFile(path);
            Assert.That(file.Results, Has.Count.EqualTo(12));
            long allocated = GC.GetAllocatedBytesForCurrentThread() - before;

            Assert.That(allocated, Is.LessThan(size / 2), $"allocated {allocated:N0} bytes reading a {size:N0}-byte file");
        }

        [TestCase(null, -1)]
        [TestCase("", -1)]
        [TestCase("scan=42", 42)]
        [TestCase("controllerType=0 controllerNumber=1 scan=7", 7)]
        [TestCase("scan=81409 file=454", 81409)]
        [TestCase("index=0", 1)]
        [TestCase("index=16332", 16333)]
        [TestCase("mzMLid=controllerType=0", -1)]
        [TestCase("index=12x", -1)]
        [TestCase("scan=1,scan=2", 1)]
        public void OneBasedScanNumber_FromTheNativeId(string spectrumId, int expected)
        {
            Assert.That(MzIdentMLResultFile.OneBasedScanNumberOf(spectrumId), Is.EqualTo(expected));
        }

        #region Synthetic documents

        /// <summary>
        /// SYNTHETIC, not a capture: one result, one item, and a peptide carrying <paramref name="modifications"/>
        /// verbatim, for the skip reasons no published file exercises.
        /// </summary>
        private string WriteSynthetic(string fileName, string sequence, string modifications,
            string itemCvParams = "", string peptideRef = "PEP_1", string substitution = "")
        {
            string document = $@"<?xml version=""1.0"" encoding=""utf-8""?>
<MzIdentML xmlns=""http://psidev.info/psi/pi/mzIdentML/1.1"" id=""synthetic"" version=""1.1.0"">
  <SequenceCollection>
    <DBSequence id=""DB_1"" accession=""P00001"" searchDatabase_ref=""SDB_1"" />
    <Peptide id=""PEP_1"">
      <PeptideSequence>{sequence}</PeptideSequence>{modifications}{substitution}
    </Peptide>
    <PeptideEvidence id=""PE_1"" dBSequence_ref=""DB_1"" peptide_ref=""PEP_1"" isDecoy=""false"" />
  </SequenceCollection>
  <DataCollection>
    <Inputs>
      <SpectraData id=""SD_1"" location=""C:\data\run.mgf"">
        <FileFormat><cvParam cvRef=""PSI-MS"" accession=""MS:1001062"" name=""Mascot MGF format"" /></FileFormat>
        <SpectrumIDFormat><cvParam cvRef=""PSI-MS"" accession=""MS:1000774"" name=""multiple peak list nativeID format"" /></SpectrumIDFormat>
      </SpectraData>
    </Inputs>
    <AnalysisData>
      <SpectrumIdentificationList id=""SIL_1"">
        <SpectrumIdentificationResult id=""SIR_1"" spectrumID=""index=4"" spectraData_ref=""SD_1"">
          <SpectrumIdentificationItem id=""SII_1"" chargeState=""2"" experimentalMassToCharge=""500.5"" peptide_ref=""{peptideRef}"" rank=""1"" passThreshold=""true"">
            <PeptideEvidenceRef peptideEvidence_ref=""PE_1"" />{itemCvParams}
          </SpectrumIdentificationItem>
        </SpectrumIdentificationResult>
      </SpectrumIdentificationList>
    </AnalysisData>
  </DataCollection>
</MzIdentML>";
            string path = Path.Combine(_outputDirectory, fileName);
            File.WriteAllText(path, document);
            return path;
        }

        private static string Modification(string location, string accession, string residues = "", string mass = "") =>
            $@"
      <Modification{(location == null ? "" : $@" location=""{location}""")}{(residues.Length > 0 ? $@" residues=""{residues}""" : "")}{(mass.Length > 0 ? $@" monoisotopicMassDelta=""{mass}""" : "")}>
        <cvParam cvRef=""UNIMOD"" accession=""{accession}"" name=""mod"" />
      </Modification>";

        [TestCase("no_unimod.mzid", "PEPTIDE", "2", "MS:1001460", "modification MS:1001460 at location 2 has no UNIMOD accession")]
        [TestCase("unknown_unimod.mzid", "PEPTIDE", "2", "UNIMOD:999999", "UNIMOD:999999 at location 2 does not resolve to a Unimod modification")]
        [TestCase("bad_unimod.mzid", "PEPTIDE", "2", "UNIMOD:abc", "modification UNIMOD:abc at location 2 has no UNIMOD accession")]
        [TestCase("past_the_end.mzid", "PEPTIDE", "9", "UNIMOD:35", "modification location 9 is not on the peptide")]
        [TestCase("negative.mzid", "PEPTIDE", "-1", "UNIMOD:35", "modification location -1 is not on the peptide")]
        [TestCase("no_location.mzid", "PEPTIDE", null, "UNIMOD:35", "modification has no location")]
        public void AnItemWithAModificationThatCannotBeResolvedIsSkipped(
            string fileName, string sequence, string location, string accession, string reason)
        {
            var file = new MzIdentMLResultFile(WriteSynthetic(fileName, sequence, Modification(location, accession)));

            Assert.That(file.Results, Is.Empty);
            Assert.That(file.SkippedMatches.Single(), Is.EqualTo(new MzIdentMLSkippedMatch("SII_1", "index=4", reason)));
        }

        [Test]
        public void AModificationWhoseOnlyTermHasNoAccessionIsSkipped()
        {
            const string mod = @"
      <Modification location=""2"">
        <cvParam cvRef=""PSI-MS"" name=""a modification with no accession"" />
      </Modification>";
            var file = new MzIdentMLResultFile(WriteSynthetic("no_accession.mzid", "PEPTIDE", mod));

            Assert.That(file.SkippedMatches.Single().Reason,
                Is.EqualTo("modification (no accession) at location 2 has no UNIMOD accession"));
        }

        /// <summary>
        /// spectrumID is required by the schema but not by the deserializer; without it there is no scan number.
        /// </summary>
        [Test]
        public void AResultWithNoSpectrumIdHasNoScanNumber()
        {
            string path = WriteSynthetic("no_spectrum_id.mzid", "PEPTIDE", "");
            File.WriteAllText(path, File.ReadAllText(path).Replace(@"spectrumID=""index=4"" ", ""));

            var file = new MzIdentMLResultFile(path);
            var record = file.Results.Single();

            Assert.That(record.SpectrumId, Is.Empty);
            Assert.That(record.OneBasedScanNumber, Is.EqualTo(-1));
            Assert.That(record.FileNameWithoutExtension, Is.EqualTo("run"));
        }

        [Test]
        public void TwoModificationsAtOnePositionAreSkipped()
        {
            string mods = Modification("4", "UNIMOD:21", "T") + Modification("4", "UNIMOD:1", "T");
            var file = new MzIdentMLResultFile(WriteSynthetic("two_at_one.mzid", "PEPTIDE", mods));

            Assert.That(file.SkippedMatches.Single().Reason, Is.EqualTo("more than one modification at location 4"));
        }

        [Test]
        public void ASubstitutionIsSkipped()
        {
            const string substitution = @"
      <SubstitutionModification originalResidue=""E"" replacementResidue=""D"" location=""2"" />";
            var file = new MzIdentMLResultFile(WriteSynthetic("substitution.mzid", "PEPTIDE", "", substitution: substitution));

            Assert.That(file.SkippedMatches.Single().Reason, Is.EqualTo("peptide carries a substitution modification"));
        }

        [Test]
        public void AnUnresolvedPeptideReferenceIsSkipped()
        {
            var file = new MzIdentMLResultFile(WriteSynthetic("dangling.mzid", "PEPTIDE", "", peptideRef: "PEP_MISSING"));

            Assert.That(file.SkippedMatches.Single().Reason, Is.EqualTo("peptide reference does not resolve"));
        }

        /// <summary>
        /// Location length + 1 is the C-terminus, keyed length + 2; the residue a terminal modification sits on
        /// comes from the sequence, since residues is optional.
        /// </summary>
        [Test]
        public void TerminalModificationsAreKeyedAtBothEnds()
        {
            string mods = Modification("0", "UNIMOD:1") + Modification("8", "UNIMOD:2") + Modification("4", "UNIMOD:21", "T", "79.966331");
            var record = new MzIdentMLResultFile(WriteSynthetic("termini.mzid", "PEPTIDE", mods)).Results.Single();

            Assert.Multiple(() =>
            {
                Assert.That(record.AllModsOneIsNterminus.Keys, Is.EquivalentTo(new[] { 1, 5, 9 }));
                Assert.That(record.AllModsOneIsNterminus[1].IdWithMotif, Is.EqualTo("Acetyl on X"));
                Assert.That(record.AllModsOneIsNterminus[9].IdWithMotif, Is.EqualTo("Amidated on X"));
                Assert.That(record.AllModsOneIsNterminus[9].LocationRestriction, Does.Contain("C-terminal"));
                Assert.That(record.FullSequence, Is.EqualTo("[Unimod:Acetyl on X]PEPT[Unimod:Phospho on T]IDE-[Unimod:Amidated on X]"));
                Assert.That(record.FileNameWithoutExtension, Is.EqualTo("run"));
                Assert.That(record.OneBasedScanNumber, Is.EqualTo(5));
                Assert.That(record.CalculatedMassToCharge, Is.Null);
                Assert.That(record.SpectrumTitle, Is.Null);
                Assert.That(record.Match.SpectrumIdentificationItemId, Is.EqualTo("SII_1"));
            });
        }

        /// <summary>
        /// Unimod has methylation both on lysine side chains and at the peptide C-terminus. A modification at
        /// location length + 1 is the C-terminal one even though the peptide ends in K.
        /// </summary>
        [Test]
        public void ACTerminalModificationOnALysineIsTheTerminalEntry()
        {
            var record = new MzIdentMLResultFile(WriteSynthetic("cterm_methyl.mzid", "PEPTIDEK", Modification("9", "UNIMOD:34"))).Results.Single();

            Assert.That(record.AllModsOneIsNterminus.Keys, Is.EquivalentTo(new[] { 10 }));
            Assert.That(record.AllModsOneIsNterminus[10].IdWithMotif, Is.EqualTo("Methyl on X"));
            Assert.That(record.AllModsOneIsNterminus[10].LocationRestriction, Does.Contain("C-terminal"));
        }

        /// <summary>
        /// A peptide found twice in one protein has two evidence entries naming the same accession; it is
        /// reported once.
        /// </summary>
        [Test]
        public void AnAccessionNamedByTwoEvidenceEntriesIsReportedOnce()
        {
            string path = WriteSynthetic("repeated_protein.mzid", "PEPTIDE", "");
            File.WriteAllText(path, File.ReadAllText(path)
                .Replace(@"peptide_ref=""PEP_1"" isDecoy=""false"" />",
                    @"peptide_ref=""PEP_1"" isDecoy=""false"" />
    <PeptideEvidence id=""PE_2"" dBSequence_ref=""DB_1"" peptide_ref=""PEP_1"" isDecoy=""false"" />")
                .Replace(@"<PeptideEvidenceRef peptideEvidence_ref=""PE_1"" />",
                    @"<PeptideEvidenceRef peptideEvidence_ref=""PE_1"" /><PeptideEvidenceRef peptideEvidence_ref=""PE_2"" />"));

            var record = new MzIdentMLResultFile(path).Results.Single();

            Assert.That(record.Match.PeptideEvidence, Has.Count.EqualTo(2));
            Assert.That(record.Accession, Is.EqualTo("P00001"));
        }

        /// <summary>
        /// The spectra file comes from the SpectraData location, then its name; a result whose SpectraData does
        /// not resolve has no file, and a nativeID naming an unknown file keeps the referenced one.
        /// </summary>
        [TestCase(@"location=""C:\data\run.mgf""", "index=4", "run", @"C:\data\run.mgf")]
        [TestCase(@"location=""/data/deep/run.two.mgf""", "index=4", "run.two", "/data/deep/run.two.mgf")]
        [TestCase(@"name=""named.mzML""", "index=4", "named", "named.mzML")]
        [TestCase(@"location=""C:\data\run.mgf""", "scan=4 file=UNKNOWN", "run", @"C:\data\run.mgf")]
        public void TheSpectraFileComesFromTheSpectraDataEntry(string spectraDataAttributes, string spectrumId,
            string expectedName, string expectedLocation)
        {
            string path = WriteSynthetic("spectra_" + expectedName + spectrumId.Length + ".mzid", "PEPTIDE", "");
            File.WriteAllText(path, File.ReadAllText(path)
                .Replace(@"location=""C:\data\run.mgf""", spectraDataAttributes)
                .Replace(@"spectrumID=""index=4""", $@"spectrumID=""{spectrumId}"""));

            var record = new MzIdentMLResultFile(path).Results.Single();

            Assert.That(record.FileNameWithoutExtension, Is.EqualTo(expectedName));
            Assert.That(record.SpectraFileLocation, Is.EqualTo(expectedLocation));
        }

        [TestCase("false", false)]
        [TestCase("true", true)]
        public void AnItemIsADecoyOnlyWhenEveryProteinItMapsToIsADecoy(string firstEvidenceIsDecoy, bool expected)
        {
            string path = WriteSynthetic($"decoy_{firstEvidenceIsDecoy}.mzid", "PEPTIDE", "");
            File.WriteAllText(path, File.ReadAllText(path)
                .Replace(@"peptide_ref=""PEP_1"" isDecoy=""false"" />",
                    $@"peptide_ref=""PEP_1"" isDecoy=""{firstEvidenceIsDecoy}"" />
    <DBSequence id=""DB_2"" accession=""DECOY_P00001"" searchDatabase_ref=""SDB_1"" />
    <PeptideEvidence id=""PE_2"" dBSequence_ref=""DB_2"" peptide_ref=""PEP_1"" isDecoy=""true"" />")
                .Replace(@"<PeptideEvidenceRef peptideEvidence_ref=""PE_1"" />",
                    @"<PeptideEvidenceRef peptideEvidence_ref=""PE_1"" /><PeptideEvidenceRef peptideEvidence_ref=""PE_2"" />"));

            var record = new MzIdentMLResultFile(path).Results.Single();

            Assert.That(record.IsDecoy, Is.EqualTo(expected));
            Assert.That(record.Accession, Is.EqualTo("P00001|DECOY_P00001"));
        }

        /// <summary>
        /// An item that references no peptide evidence names no protein and cannot be called a decoy.
        /// </summary>
        [Test]
        public void AnItemWithNoPeptideEvidenceHasNoAccessionAndIsNotADecoy()
        {
            string path = WriteSynthetic("no_evidence.mzid", "PEPTIDE", "");
            File.WriteAllText(path, File.ReadAllText(path).Replace(@"<PeptideEvidenceRef peptideEvidence_ref=""PE_1"" />", ""));

            var record = new MzIdentMLResultFile(path).Results.Single();

            Assert.That(record.Accession, Is.Empty);
            Assert.That(record.IsDecoy, Is.False);
            Assert.That(record.BaseSequence, Is.EqualTo("PEPTIDE"));
        }

        [Test]
        public void AResultWhoseSpectraDataDoesNotResolveHasNoFile()
        {
            string path = WriteSynthetic("no_spectradata.mzid", "PEPTIDE", "");
            File.WriteAllText(path, File.ReadAllText(path).Replace(@"spectraData_ref=""SD_1""", @"spectraData_ref=""SD_MISSING"""));

            var record = new MzIdentMLResultFile(path).Results.Single();

            Assert.That(record.FileNameWithoutExtension, Is.Empty);
            Assert.That(record.SpectraFileLocation, Is.Empty);
        }

        /// <summary>
        /// A decoy only when every protein the item maps to is a decoy; the q-value is read only from a numeric
        /// q-value term, and a score with no numeric value is left out of Scores.
        /// </summary>
        [Test]
        public void ScoresAndQValueIgnoreNonNumericValues()
        {
            const string cvParams = @"
            <cvParam cvRef=""PSI-MS"" accession=""MS:1002354"" name=""PSM-level q-value"" value=""n/a"" />
            <cvParam cvRef=""PSI-MS"" accession=""MS:1001363"" name=""peptide unique to one protein"" />
            <cvParam cvRef=""PSI-MS"" accession=""MS:1001171"" name=""Mascot:score"" value=""40.5"" />
            <cvParam cvRef=""PSI-MS"" accession=""MS:1001171"" name=""Mascot:score"" value=""1"" />
            <userParam name=""engine rank"" value=""3"" />";
            var record = new MzIdentMLResultFile(WriteSynthetic("scores.mzid", "PEPTIDE", "", cvParams)).Results.Single();

            Assert.Multiple(() =>
            {
                Assert.That(record.QValue, Is.Null);
                Assert.That(record.Scores.Keys, Is.EquivalentTo(new[] { "Mascot:score", "engine rank" }));
                Assert.That(record.Scores["Mascot:score"], Is.EqualTo(40.5));
                Assert.That(record.IsDecoy, Is.False);
                Assert.That(record.Accession, Is.EqualTo("P00001"));
            });
        }

        #endregion
    }
}
