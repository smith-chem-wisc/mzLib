using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using MzIdentML;
using MzLibUtil;
using NUnit.Framework;

namespace Test.MzIdentML
{
    [TestFixture]
    public class MzidIdentificationsTests
    {

        private string FilePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "SmallCalibratible_Yeast.mzID");
        private MzidIdentifications mzid;

        [SetUp]
        public void Setup()
        {
            mzid = new MzidIdentifications(FilePath);
        }

        private static MzidIdentifications Read(string fileName) =>
            new MzidIdentifications(Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", fileName));

        /// <summary>
        /// Calls every accessor for every PSM and every modification in a file.
        ///
        /// The assertions are limited to the ones that can actually fail. Each accessor starts from a
        /// sentinel (null, or -1 for the two numeric modification accessors) and only overwrites it
        /// once a PeptideEvidenceRef resolves against SequenceCollection, so "not the sentinel" is the
        /// real check -- it fails whenever a lookup falls through. ProteinFullName and
        /// CalculatedMassToCharge are read but not constrained: both are legitimately absent from some
        /// of these files, so any bound on them would pass no matter what the reader did.
        ///
        /// Sweeping rather than reading one PSM matters because the lookups are linear scans, so a
        /// single PSM only ever takes the matching side of the id comparison inside them.
        /// </summary>
        private static (int Psms, int Decoys, int WithQValue) ReadEveryPsm(MzidIdentifications ids)
        {
            int psmCount = 0, decoys = 0, withQValue = 0;

            // Count is itself a version cascade, so it is read once rather than per iteration
            int results = ids.Count;

            for (int sir = 0; sir < results; sir++)
            {
                Assert.That(ids.Ms2SpectrumID(sir), Is.Not.Null.And.Not.Empty);

                int psms = ids.NumPSMsFromScan(sir);
                Assert.That(psms, Is.GreaterThan(0));

                for (int sii = 0; sii < psms; sii++)
                {
                    psmCount++;

                    Assert.That(ids.PeptideSequenceWithoutModifications(sir, sii), Is.Not.Null.And.Not.Empty);
                    Assert.That(ids.ProteinAccession(sir, sii), Is.Not.Null.And.Not.Empty);
                    Assert.That(ids.StartResidueInProtein(sir, sii), Is.Not.Null.And.Not.Empty);
                    Assert.That(ids.EndResidueInProtein(sir, sii), Is.Not.Null.And.Not.Empty);
                    Assert.That(ids.ChargeState(sir, sii), Is.GreaterThan(0));
                    Assert.That(ids.ExperimentalMassToCharge(sir, sii), Is.GreaterThan(0));

                    _ = ids.ProteinFullName(sir, sii);
                    _ = ids.CalculatedMassToCharge(sir, sii);

                    if (ids.IsDecoy(sir, sii)) decoys++;
                    if (ids.QValue(sir, sii) != -1) withQValue++;

                    for (int mod = 0; mod < ids.NumModifications(sir, sii); mod++)
                    {
                        Assert.That(ids.ModificationAcession(sir, sii, mod), Is.Not.Null.And.Not.Empty);
                        Assert.That(ids.ModificationDictionary(sir, sii, mod), Is.Not.Null.And.Not.Empty);

                        // -1 is the not-found sentinel for both. 0 is a real location: N-terminal.
                        Assert.That(ids.ModificationLocation(sir, sii, mod), Is.Not.EqualTo(-1));
                        Assert.That(ids.ModificationMass(sir, sii, mod), Is.Not.EqualTo(-1));
                    }
                }
            }

            return (psmCount, decoys, withQValue);
        }

        /// <summary>
        /// PSM, decoy and q-value counts were read out of each file's XML rather than off the reader,
        /// so a regression that silently drops PSMs or stops resolving PeptideEvidence fails here.
        ///
        /// The counts are per SpectrumIdentificationList[0], which is all MzidIdentifications ever
        /// looks at -- see MzIdentML130_OnlyTheFirstSpectrumIdentificationListIsRead.
        /// </summary>
        [TestCase("SmallCalibratible_Yeast.mzID", 73, 0, 73)]
        [TestCase("OpenxQuest_example_1_2_0.mzid", 16, 10, 0)]
        [TestCase("mzidLib_xtandem_fdr_1_2_0.mzid", 5, 0, 0)]
        [TestCase("multiple_spectra_per_id_1_3_0.mzid", 2, 0, 0)]
        [TestCase("noncovalently_assoc_1_3_0.mzid", 2, 0, 0)]
        // MS-GF+ does report a q-value on all 12 PSMs, as MS:1002054 "MS-GF:QValue". That term is_a
        // MS:1002354, but QValue matches only the parent accession, so the reader sees none of them.
        // 0 pins the current behaviour; it does not mean the file has no q-values.
        [TestCase("PXD078927_msgf_1_1_0.mzid", 12, 6, 0)]
        [TestCase("PXD000783_scaffold_1_1_0.mzid", 3, 0, 0)]
        [TestCase("PXD019591_mascotparser_1_1_0.mzid", 7, 0, 0)]
        [TestCase("PXD019733_proteomediscoverer_1_1_0.mzid", 3, 0, 0)]
        [TestCase("PXD070193_xifdr_1_2_0.mzid", 4, 0, 0)]
        public void EveryPsmIsReadable(string fileName, int psms, int decoys, int withQValue)
        {
            var swept = ReadEveryPsm(Read(fileName));

            Assert.Multiple(() =>
            {
                Assert.That(swept.Psms, Is.EqualTo(psms));
                Assert.That(swept.Decoys, Is.EqualTo(decoys));
                Assert.That(swept.WithQValue, Is.EqualTo(withQValue));
            });
        }

        /// <summary>
        /// The PXD fixtures are real files downloaded from PRIDE Archive, each cut down to its first
        /// few SpectrumIdentificationResults plus the DBSequence, Peptide and PeptideEvidence entries
        /// those results reference. Every kept line is byte-for-byte from the published file; the only
        /// text not in the original is the closing tags after the last kept result. Sources, checksums
        /// and the full trim rule are in DataFiles/PRIDE_MZID_PROVENANCE.md. They pin what the published,
        /// synthetic examples above cannot: the CV names real writers use.
        ///
        /// Ms2SpectrumID decides between spectrumID and the spectrum title from SpectraData's FileFormat.
        /// It matched that term by display name only, and two of the four writers spell it differently
        /// from the PSI-MS name while carrying the same accession:
        ///
        /// - MS-GF+ (PXD078927) writes "mzML file" for MS:1000584, so the mzML branch was skipped and
        ///   the method returned null for every result.
        /// - Scaffold (PXD000783) writes "Mascot MGF file" for MS:1001062, with the same null outcome.
        ///
        /// And the MGF branch read the result's FIRST cvParam as the title. Scaffold happens to write the
        /// title first; Mascot Parser (PXD019591) writes "Mascot:identity threshold" first, so the method
        /// returned the identity threshold ("17") instead of the title, with no error.
        ///
        /// Proteome Discoverer (PXD019733) names the mzML format exactly and was already right; it is
        /// here to pin that the match widening does not disturb an exact-name file.
        ///
        /// xiFDR (PXD070193, mzIdentML 1.2.0, cut to two results) writes no spectrum title at all, only
        /// "peak list scans", so it takes the fallback to the first cvParam -- the value that was returned
        /// before -- rather than changing what a title-less file reports.
        /// </summary>
        [TestCase("PXD078927_msgf_1_1_0.mzid", 0, "controllerType=0 controllerNumber=1 scan=14316")]
        [TestCase("PXD078927_msgf_1_1_0.mzid", 2, "controllerType=0 controllerNumber=1 scan=14164")]
        [TestCase("PXD000783_scaffold_1_1_0.mzid", 0, "Locus:1.1.1.3846.2 File:\"20130114-cw15-HILIC-polyMAC-biorep2-F6.wiff\"")]
        [TestCase("PXD019591_mascotparser_1_1_0.mzid", 0, @"16333: Scan 22014 (rt=62.2092) [\\qmcr.qmul.ac.uk\rfs\HAEMONC\Pedro's Lab\Mass Spec\PharmacoDB\RAW\TP\QE2_PharmacoDB_TP_Sample_16B.raw]")]
        [TestCase("PXD019591_mascotparser_1_1_0.mzid", 2, @"3735: Scan 7629 (rt=31.44) [\\qmcr.qmul.ac.uk\rfs\HAEMONC\Pedro's Lab\Mass Spec\PharmacoDB\RAW\TP\QE2_PharmacoDB_TP_Sample_16B.raw]")]
        [TestCase("PXD019733_proteomediscoverer_1_1_0.mzid", 0, "scan=81409 file=454")]
        [TestCase("PXD070193_xifdr_1_2_0.mzid", 1, "8426")]
        public void Ms2SpectrumID_RealWriters_RecogniseTheFormatByAccessionAndReadTheTitleByAccession(
            string fileName, int sirIndex, string expected)
        {
            Assert.That(Read(fileName).Ms2SpectrumID(sirIndex), Is.EqualTo(expected));
        }

        /// <summary>
        /// Proteome Discoverer (PXD019733) scores each SpectrumIdentificationItem with a userParam and no
        /// cvParam at all, so the item's cvParam array deserializes as null. QValue filtered that array
        /// for MS:1002354 without a null check and threw ArgumentNullException for every PSM in the file
        /// -- 861,250 of them in the published original. An item with no cvParams has no q-value, which
        /// is the -1 the method already returns when the term is absent.
        /// </summary>
        [Test]
        public void QValue_ItemWithNoCvParams_IsAbsentRatherThanAThrow()
        {
            var pd = Read("PXD019733_proteomediscoverer_1_1_0.mzid");

            Assert.That(pd.QValue(0, 0), Is.EqualTo(-1).Within(1e-9));
        }

        private const string Mzid110 = "http://psidev.info/psi/pi/mzIdentML/1.1";
        private const string Mzid111 = "http://psidev.info/psi/pi/mzIdentML/1.1.1";
        private const string Mzid120 = "http://psidev.info/psi/pi/mzIdentML/1.2";
        private const string Mzid130 = "http://psidev.info/psi/pi/mzIdentML/1.3";

        /// <summary>
        /// SYNTHETIC, not a capture: the smallest document that reaches Ms2SpectrumID and QValue. Every
        /// accessor fix in this file is repeated once per schema version, and the real PRIDE fixtures are
        /// all 1.1.0 or 1.2.0, so the other version arms need a document no writer we have produced.
        /// The one result carries <paramref name="resultCvParams"/> verbatim (empty for none) and its one
        /// item carries <paramref name="itemCvParams"/>; by default none at all, which is the Proteome
        /// Discoverer shape.
        /// </summary>
        private static string SyntheticDocument(string ns, string formatAccession, string formatName, string resultCvParams,
            string itemCvParams = "") =>
            $@"<?xml version=""1.0"" encoding=""utf-8""?>
<MzIdentML xmlns=""{ns}"" id=""synthetic"" version=""1.1.0"">
  <DataCollection>
    <Inputs>
      <SpectraData id=""SD_1"" location=""spectra"">
        <FileFormat>
          <cvParam cvRef=""PSI-MS"" accession=""{formatAccession}"" name=""{formatName}"" />
        </FileFormat>
      </SpectraData>
    </Inputs>
    <AnalysisData>
      <SpectrumIdentificationList id=""SIL_1"">
        <SpectrumIdentificationResult id=""SIR_1"" spectrumID=""index=7"" spectraData_ref=""SD_1"">
          <SpectrumIdentificationItem id=""SII_1"" chargeState=""2"" experimentalMassToCharge=""500.5"" rank=""1"" passThreshold=""true"">
{itemCvParams}
            <userParam name=""Percolator q-Value"" value=""0.001"" />
          </SpectrumIdentificationItem>{resultCvParams}
        </SpectrumIdentificationResult>
      </SpectrumIdentificationList>
    </AnalysisData>
  </DataCollection>
</MzIdentML>";

        private static MzidIdentifications ReadSynthetic(string contents, string fileName)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, fileName);
            File.WriteAllText(path, contents);
            try
            {
                return new MzidIdentifications(path);
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// The null guard on the item's cvParam array was added to every version arm, but the Proteome
        /// Discoverer capture only reaches the 1.1.0 one. Without the guard each of these throws
        /// ArgumentNullException instead of reporting the q-value as absent.
        /// </summary>
        [TestCase(Mzid111)]
        [TestCase(Mzid120)]
        [TestCase(Mzid130)]
        public void QValue_ItemWithNoCvParams_IsAbsentInEveryVersion(string ns)
        {
            var ids = ReadSynthetic(SyntheticDocument(ns, "MS:1000584", "mzML format", ""), "synthetic_qvalue.mzid");

            Assert.That(ids.QValue(0, 0), Is.EqualTo(-1).Within(1e-9));
        }

        /// <summary>
        /// An MGF-backed result with no cvParam has no title to return. Reading position 0 of a null
        /// array threw; the method now reports the ID as missing, in every version arm.
        /// </summary>
        [TestCase(Mzid110)]
        [TestCase(Mzid111)]
        [TestCase(Mzid120)]
        [TestCase(Mzid130)]
        public void Ms2SpectrumID_MgfResultWithNoCvParams_IsNullRatherThanAThrow(string ns)
        {
            var ids = ReadSynthetic(SyntheticDocument(ns, "MS:1001062", "Mascot MGF format", ""), "synthetic_notitle.mzid");

            string id = "not read";
            Assert.That(() => id = ids.Ms2SpectrumID(0), Throws.Nothing);
            Assert.That(id, Is.Null);
        }

        private static readonly string[] EveryNamespace = { Mzid110, Mzid111, Mzid120, Mzid130 };

        private static readonly (string Accession, string Name, string Expected)[] FormatsUnderWritersNames =
        {
            ("MS:1000563", "Thermo RAW file", "index=7"),
            ("MS:1000584", "mzML file", "index=7"),
            ("MS:1001062", "Mascot MGF file", "title of spectrum 7"),
        };

        private static IEnumerable<TestCaseData> FormatsUnderWritersNamesInEveryVersion() =>
            from ns in EveryNamespace
            from format in FormatsUnderWritersNames
            select new TestCaseData(ns, format.Accession, format.Name, format.Expected);

        /// <summary>
        /// Every version arm carries its own copy of both Ms2SpectrumID fixes, and IsFileFormat is one
        /// shared static, so line coverage reads 100% whether or not each call site is right. Only the
        /// 1.1.0 arm is reached by a real writer's display name, so each arm is driven here: a format
        /// recognised by accession under a writer's name, and a title found by accession when another
        /// term precedes it. The spectrumID is "index=7" and the leading term is the identity threshold,
        /// so reading the wrong branch or the first cvParam gives a different answer.
        /// </summary>
        [TestCaseSource(nameof(FormatsUnderWritersNamesInEveryVersion))]
        public void Ms2SpectrumID_EveryVersion_RecognisesTheFormatAndTitleByAccession(
            string ns, string formatAccession, string formatName, string expected)
        {
            const string cvParams = @"
          <cvParam cvRef=""PSI-MS"" accession=""MS:1001371"" name=""Mascot:identity threshold"" value=""17"" />
          <cvParam cvRef=""PSI-MS"" accession=""MS:1000796"" name=""spectrum title"" value=""title of spectrum 7"" />";
            var ids = ReadSynthetic(SyntheticDocument(ns, formatAccession, formatName, cvParams), "synthetic_formats.mzid");

            Assert.That(ids.Ms2SpectrumID(0), Is.EqualTo(expected));
        }

        /// <summary>
        /// A title term that is present but unusable is a missing title. The fallback to the first cvParam
        /// is only for results that carry no title term at all; taking it here returns the identity
        /// threshold, which is the defect the title lookup exists to fix. The value attribute is optional
        /// on every cvParam, and PXD019591 already carries a valueless one on an item.
        /// </summary>
        [TestCase(@"<cvParam cvRef=""PSI-MS"" accession=""MS:1000796"" name=""spectrum title"" />", null)]
        [TestCase(@"<cvParam cvRef=""PSI-MS"" accession=""MS:1000796"" name=""spectrum title"" value="""" />", null)]
        [TestCase(@"<cvParam cvRef=""PSI-MS"" accession=""MS:1001416"" name=""spectrum title"" value=""title of spectrum 7"" />", "title of spectrum 7")]
        public void Ms2SpectrumID_TitleTermIsNotConfusedWithTheFirstCvParam(string titleTerm, string expected)
        {
            string cvParams = @"
          <cvParam cvRef=""PSI-MS"" accession=""MS:1001371"" name=""Mascot:identity threshold"" value=""17"" />
          " + titleTerm;
            var ids = ReadSynthetic(SyntheticDocument(Mzid110, "MS:1001062", "Mascot MGF file", cvParams), "synthetic_title.mzid");

            Assert.That(ids.Ms2SpectrumID(0), Is.EqualTo(expected));
        }

        /// <summary>
        /// A q-value term with no value is absent, not 0: Convert.ToDouble reads a null string as 0, which
        /// is indistinguishable from the most confident q-value there is. The value parse is shared, but
        /// each version arm passes its own term into it, so all four are driven.
        /// </summary>
        [TestCaseSource(nameof(EveryNamespace))]
        public void QValue_TermWithNoValue_IsAbsentRatherThanZero(string ns)
        {
            const string itemCvParams = @"
            <cvParam cvRef=""PSI-MS"" accession=""MS:1002354"" name=""PSM-level q-value"" />";
            var ids = ReadSynthetic(SyntheticDocument(ns, "MS:1000584", "mzML format", "", itemCvParams), "synthetic_novalue.mzid");

            Assert.That(ids.QValue(0, 0), Is.EqualTo(-1).Within(1e-9));
        }

        /// <summary>
        /// Control for the test above, so it cannot pass by never reading the term.
        /// </summary>
        [Test]
        public void QValue_TermWithValue_IsRead()
        {
            const string itemCvParams = @"
            <cvParam cvRef=""PSI-MS"" accession=""MS:1002354"" name=""PSM-level q-value"" value=""0.0125"" />";
            var ids = ReadSynthetic(SyntheticDocument(Mzid130, "MS:1000584", "mzML format", "", itemCvParams), "synthetic_value.mzid");

            Assert.That(ids.QValue(0, 0), Is.EqualTo(0.0125).Within(1e-9));
        }

        /// <summary>
        /// There was no 1.2.0 fixture, so the 1.2.0 level of every cascade was only ever reached as the
        /// null dereference that falls through to 1.3.0. This is
        /// examples/1_2examples/crosslinking/OpenxQuest_example.mzid from the HUPO-PSI repository.
        ///
        /// Its SpectraData is mzML, so Ms2SpectrumID returns spectrumID verbatim. As with the 1.3.0
        /// test below, MatchedIons and MatchedIonCounts are not covered: the file carries no
        /// fragmentation table.
        /// </summary>
        [Test]
        public void MzIdentML120_IsReadThroughEveryVersionCascade()
        {
            var mzid120 = Read("OpenxQuest_example_1_2_0.mzid");

            Assert.Multiple(() =>
            {
                Assert.That(mzid120.ParentTolerance, Is.TypeOf<PpmTolerance>());
                Assert.That(mzid120.ParentTolerance.Value, Is.EqualTo(10.0).Within(1e-9));
                Assert.That(mzid120.FragmentTolerance.Value, Is.EqualTo(0.2).Within(1e-9));

                Assert.That(mzid120.Count, Is.EqualTo(1));
                Assert.That(mzid120.NumPSMsFromScan(0), Is.EqualTo(16));
                Assert.That(mzid120.Ms2SpectrumID(0), Is.EqualTo("scan=1,scan=2"));

                Assert.That(mzid120.CalculatedMassToCharge(0, 0), Is.EqualTo(718.396192605738).Within(1e-9));
                Assert.That(mzid120.ExperimentalMassToCharge(0, 0), Is.EqualTo(672.374450683594).Within(1e-9));
                Assert.That(mzid120.ChargeState(0, 0), Is.EqualTo(3));

                // this file writes isDecoy as "0"/"1" rather than "false"/"true"; both are valid
                // xs:boolean, and 5 of its 8 PeptideEvidence entries are decoys
                Assert.That(mzid120.IsDecoy(0, 0), Is.False);
                Assert.That(mzid120.IsDecoy(0, 4), Is.True);

                // no MS:1002354, reported as -1
                Assert.That(mzid120.QValue(0, 0), Is.EqualTo(-1).Within(1e-9));

                Assert.That(mzid120.PeptideSequenceWithoutModifications(0, 0), Is.EqualTo("SPAIIFIDELDAIGTKR"));
                Assert.That(mzid120.ProteinAccession(0, 0), Is.EqualTo("sp|O14126|PRS6A_SCHPO"));
                Assert.That(mzid120.StartResidueInProtein(0, 0), Is.EqualTo("277"));
                Assert.That(mzid120.EndResidueInProtein(0, 0), Is.EqualTo("293"));

                // its DBSequence entries carry no name attribute
                Assert.That(mzid120.ProteinFullName(0, 0), Is.Empty);

                Assert.That(mzid120.NumModifications(0, 0), Is.EqualTo(1));
                Assert.That(mzid120.ModificationAcession(0, 0, 0), Is.EqualTo("UNIMOD:1020"));
                Assert.That(mzid120.ModificationDictionary(0, 0, 0), Is.EqualTo("UNIMOD"));
                Assert.That(mzid120.ModificationLocation(0, 0, 0), Is.EqualTo(16));

                // the cvParam has no value attribute and the Modification no monoisotopicMassDelta,
                // which the accessors report as null and 0 rather than as absent
                Assert.That(mzid120.ModificationValue(0, 0, 0), Is.Null);
                Assert.That(mzid120.ModificationMass(0, 0, 0), Is.Zero);
            });
        }

        /// <summary>
        /// examples/1_2examples/protein_inference/mzidLib_peaklist2a_plus_ecoli_versus_unimod_full_xtandem_fdr_threshold_groups.mzid
        /// from the HUPO-PSI repository, shortened on the way in. It reaches three 1.2.0 branches
        /// OpenxQuest does not:
        ///
        /// - both tolerances are in daltons, so they read back as AbsoluteTolerance, not PpmTolerance
        /// - its SpectraData is Mascot MGF, so Ms2SpectrumID reads the SpectrumIdentificationResult's
        ///   spectrum title (MS:1000796) instead of spectrumID
        /// - its peptides map to up to 29 PeptideEvidence entries, so the accessors that concatenate
        ///   across shared proteins actually concatenate
        /// </summary>
        [Test]
        public void MzIdentML120_DaltonTolerancesMgfSpectraAndSharedPeptides()
        {
            var mzid120 = Read("mzidLib_xtandem_fdr_1_2_0.mzid");

            Assert.Multiple(() =>
            {
                Assert.That(mzid120.ParentTolerance, Is.TypeOf<AbsoluteTolerance>());
                Assert.That(mzid120.ParentTolerance.Value, Is.EqualTo(0.2).Within(1e-9));
                Assert.That(mzid120.FragmentTolerance, Is.TypeOf<AbsoluteTolerance>());
                Assert.That(mzid120.FragmentTolerance.Value, Is.EqualTo(0.3).Within(1e-9));

                Assert.That(mzid120.Count, Is.EqualTo(5));
                Assert.That(mzid120.NumPSMsFromScan(0), Is.EqualTo(1));

                // the MS:1000796 spectrum title, found by accession; spectrumID is "index=12"
                Assert.That(mzid120.Ms2SpectrumID(0), Is.EqualTo("Locus:11.1.1.4652.4 File:\"R1 p450 iTRAQ QS CEX11.wiff\""));

                Assert.That(mzid120.PeptideSequenceWithoutModifications(0, 0), Is.EqualTo("MPYTNAVIHEVQR"));
                Assert.That(mzid120.ChargeState(0, 0), Is.EqualTo(3));
                Assert.That(mzid120.ExperimentalMassToCharge(0, 0), Is.EqualTo(567.9671).Within(1e-9));
                Assert.That(mzid120.CalculatedMassToCharge(0, 0), Is.EqualTo(567.966917).Within(1e-9));

                // this PSM has ten PeptideEvidence entries across ten proteins. ProteinAccession
                // reports only the first; the other three report all ten, joined with " or "
                Assert.That(mzid120.ProteinAccession(0, 0), Does.StartWith("sp|P24457|CP2DB_MOUSE"));
                Assert.That(mzid120.StartResidueInProtein(0, 0), Is.EqualTo("356 or 356 or 356 or 356 or 356 or 356 or 356 or 194 or 194 or 194"));
                Assert.That(mzid120.EndResidueInProtein(0, 0), Is.EqualTo("368 or 368 or 368 or 368 or 368 or 368 or 368 or 206 or 206 or 206"));
                Assert.That(mzid120.ProteinFullName(0, 0).Split(" or "), Has.Length.EqualTo(10));
                Assert.That(mzid120.ProteinFullName(0, 0), Does.StartWith("sp|P24457|CP2DB_MOUSE Cytochrome P450 2D11"));

                // it reports MS:1001868 rather than the MS:1002354 QValue looks for
                Assert.That(mzid120.QValue(0, 0), Is.EqualTo(-1).Within(1e-9));

                Assert.That(mzid120.NumModifications(0, 0), Is.EqualTo(1));
                Assert.That(mzid120.ModificationAcession(0, 0, 0), Is.EqualTo("UNIMOD:214"));
                Assert.That(mzid120.ModificationDictionary(0, 0, 0), Is.EqualTo("UNIMOD"));
                Assert.That(mzid120.ModificationMass(0, 0, 0), Is.EqualTo(144.10201).Within(1e-9));

                // N-terminal iTRAQ, and 0 rather than -1 is what distinguishes it from a failed lookup
                Assert.That(mzid120.ModificationLocation(0, 0, 0), Is.Zero);
            });
        }

        /// <summary>
        /// examples/1_3examples/crosslinking/multiple_spectra_per_id_1_3_0_draft.mzid from the HUPO-PSI
        /// repository, renamed on the way in. It is the only published 1.3.0 example whose SpectraData
        /// is mzML rather than Mascot MGF, so it is the only one that takes Ms2SpectrumID's first
        /// branch, and its DBSequence entries carry names so ProteinFullName returns something.
        /// </summary>
        [Test]
        public void MzIdentML130_MzmlSpectraDataTakesTheSpectrumIdBranch()
        {
            var mzid130 = Read("multiple_spectra_per_id_1_3_0.mzid");

            Assert.Multiple(() =>
            {
                Assert.That(mzid130.ParentTolerance.Value, Is.EqualTo(3.0).Within(1e-9));
                Assert.That(mzid130.FragmentTolerance.Value, Is.EqualTo(10.0).Within(1e-9));

                Assert.That(mzid130.NumPSMsFromScan(0), Is.EqualTo(2));

                // spectrumID verbatim, not a cvParam
                Assert.That(mzid130.Ms2SpectrumID(0), Is.EqualTo("index=1"));

                Assert.That(mzid130.ExperimentalMassToCharge(0, 0), Is.EqualTo(210.093).Within(1e-9));
                Assert.That(mzid130.ChargeState(0, 0), Is.EqualTo(3));
                Assert.That(mzid130.PeptideSequenceWithoutModifications(0, 0), Is.EqualTo("PEPK"));
                Assert.That(mzid130.ProteinAccession(0, 0), Is.EqualTo("PA"));
                Assert.That(mzid130.ProteinFullName(0, 0), Is.EqualTo("Protein A"));
                Assert.That(mzid130.StartResidueInProtein(0, 0), Is.EqualTo("11"));

                Assert.That(mzid130.NumModifications(0, 0), Is.EqualTo(1));
                Assert.That(mzid130.ModificationAcession(0, 0, 0), Is.EqualTo("MS:1003393"));
                Assert.That(mzid130.ModificationValue(0, 0, 0), Is.EqualTo("DSSO_crosslink_donor"));
                Assert.That(mzid130.ModificationDictionary(0, 0, 0), Is.EqualTo("PSI-MS"));
                Assert.That(mzid130.ModificationLocation(0, 0, 0), Is.EqualTo(4));
                Assert.That(mzid130.ModificationMass(0, 0, 0), Is.EqualTo(158.003765).Within(1e-9));

                Assert.That(mzid130.PeptideSequenceWithoutModifications(0, 1), Is.EqualTo("TIDEK"));
                Assert.That(mzid130.ProteinAccession(0, 1), Is.EqualTo("PB"));
                Assert.That(mzid130.ProteinFullName(0, 1), Is.EqualTo("Protein B"));
            });
        }

        /// <summary>
        /// Every accessor indexes SpectrumIdentificationList[0], on all four version branches, so a
        /// file with more than one list is read only in part. This fixture has three -- sil_HCD,
        /// sil_ETD and sil_MS3, holding 1, 1 and 4 SpectrumIdentificationResults -- and Count reports
        /// the 1 in sil_HCD.
        ///
        /// Asserting it so the limitation is visible and so a future fix has to change this test
        /// rather than quietly change what callers see. Pre-existing and not specific to 1.3.0.
        /// </summary>
        [Test]
        public void MzIdentML130_OnlyTheFirstSpectrumIdentificationListIsRead()
        {
            Assert.That(Read("multiple_spectra_per_id_1_3_0.mzid").Count, Is.EqualTo(1));
        }

        /// <summary>
        /// The Stream constructor runs the same version cascade as the path constructor over a buffered copy, so
        /// every fixture, one per schema version and every writer, must read identically either way.
        /// </summary>
        [TestCase("SmallCalibratible_Yeast.mzID")]
        [TestCase("OpenxQuest_example_1_2_0.mzid")]
        [TestCase("multiple_spectra_per_id_1_3_0.mzid")]
        [TestCase("PXD078927_msgf_1_1_0.mzid")]
        [TestCase("PXD019733_proteomediscoverer_1_1_0.mzid")]
        public void StreamConstructor_ReadsLikeThePathConstructor(string fileName)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", fileName);
            var fromPath = new MzidIdentifications(path);
            using var stream = File.OpenRead(path);
            var fromStream = new MzidIdentifications(stream);

            Assert.Multiple(() =>
            {
                Assert.That(fromStream.Count, Is.EqualTo(fromPath.Count));
                Assert.That(fromStream.Ms2SpectrumID(0), Is.EqualTo(fromPath.Ms2SpectrumID(0)));
                Assert.That(fromStream.GetSpectrumMatches().Select(m => m.SpectrumIdentificationItemId),
                    Is.EqualTo(fromPath.GetSpectrumMatches().Select(m => m.SpectrumIdentificationItemId)));
            });
        }

        /// <summary>
        /// A GZipStream cannot seek, and the cascade opens the document once per version it tries, so the
        /// constructor must buffer it. The caller's stream is read but left open.
        /// </summary>
        [Test]
        public void StreamConstructor_ReadsAForwardOnlyGzipStreamAndLeavesItOpen()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "PXD078927_msgf_1_1_0.mzid.gz");
            using var file = File.OpenRead(path);
            using var gzip = new GZipStream(file, CompressionMode.Decompress);

            var ids = new MzidIdentifications(gzip);

            Assert.That(ids.GetSpectrumMatches().Count(), Is.EqualTo(12));
            Assert.That(gzip.CanRead, Is.True);
            Assert.That(file.CanRead, Is.True);
        }

        [Test]
        public void StreamConstructor_RejectsANullStream()
        {
            Assert.That(() => new MzidIdentifications((Stream)null), Throws.ArgumentNullException);
        }

        /// <summary>
        /// The last arm of the cascade, the legacy ".../1.1.0" namespace mzLib used to write, is reached through
        /// the Stream constructor too.
        /// </summary>
        [Test]
        public void StreamConstructor_ReadsTheLegacyNamespace()
        {
            string document = SyntheticDocument("http://psidev.info/psi/pi/mzIdentML/1.1.0", "MS:1000584", "mzML format", "");
            using var stream = new MemoryStream(System.Text.Encoding.UTF8.GetBytes(document));

            var ids = new MzidIdentifications(stream);

            Assert.That(ids.Ms2SpectrumID(0), Is.EqualTo("index=7"));
        }

        /// <summary>
        /// SYNTHETIC, not a capture. Two SpectrumIdentificationLists, and in the first result two items: one
        /// with every optional attribute and param a record reads, one with neither its own peptide_ref nor
        /// calculatedMassToCharge, whose evidence names its peptide. The second list's result references a
        /// SpectraData that does not exist, and its item references evidence that does not exist.
        /// </summary>
        private static string RichSyntheticDocument(string ns) =>
            $@"<?xml version=""1.0"" encoding=""utf-8""?>
<MzIdentML xmlns=""{ns}"" id=""synthetic"" version=""1.1.0"">
  <SequenceCollection>
    <DBSequence id=""DB_1"" accession=""P12345"" searchDatabase_ref=""SDB_1"" />
    <DBSequence id=""DB_2"" accession=""DECOY_P12345"" searchDatabase_ref=""SDB_1"" />
    <Peptide id=""PEP_1"">
      <PeptideSequence>PEPTIDE</PeptideSequence>
      <Modification location=""0"" monoisotopicMassDelta=""42.010565"">
        <cvParam cvRef=""UNIMOD"" accession=""UNIMOD:1"" name=""Acetyl"" />
      </Modification>
      <Modification location=""4"" residues=""T"">
        <cvParam cvRef=""UNIMOD"" accession=""UNIMOD:21"" name=""Phospho"" />
        <cvParam cvRef=""PSI-MS"" accession=""MS:1001524"" name=""fragment neutral loss"" value=""0"" unitCvRef=""UO"" unitAccession=""UO:0000221"" unitName=""dalton"" />
      </Modification>
      <SubstitutionModification originalResidue=""E"" replacementResidue=""D"" location=""7"" />
    </Peptide>
    <Peptide id=""PEP_2"">
      <PeptideSequence>SAMPLER</PeptideSequence>
    </Peptide>
    <PeptideEvidence id=""PE_1"" dBSequence_ref=""DB_1"" peptide_ref=""PEP_1"" start=""10"" end=""16"" pre=""K"" post=""A"" isDecoy=""false"" />
    <PeptideEvidence id=""PE_2"" dBSequence_ref=""DB_2"" peptide_ref=""PEP_2"" isDecoy=""true"" />
    <PeptideEvidence id=""PE_3"" dBSequence_ref=""DB_MISSING"" peptide_ref=""PEP_2"" isDecoy=""true"" />
  </SequenceCollection>
  <DataCollection>
    <Inputs>
      <SpectraData id=""SD_1"" name=""run one"" location=""C:\data\run1.raw"">
        <FileFormat><cvParam cvRef=""PSI-MS"" accession=""MS:1000563"" name=""Thermo RAW format"" /></FileFormat>
        <SpectrumIDFormat><cvParam cvRef=""PSI-MS"" accession=""MS:1000768"" name=""Thermo nativeID format"" /></SpectrumIDFormat>
      </SpectraData>
      <SpectraData id=""SD_2"" location=""run2.mgf"" />
    </Inputs>
    <AnalysisData>
      <SpectrumIdentificationList id=""SIL_1"">
        <SpectrumIdentificationResult id=""SIR_1"" spectrumID=""scan=7"" spectraData_ref=""SD_1"">
          <SpectrumIdentificationItem id=""SII_1"" chargeState=""2"" experimentalMassToCharge=""500.5"" calculatedMassToCharge=""500.25"" peptide_ref=""PEP_1"" rank=""1"" passThreshold=""true"">
            <PeptideEvidenceRef peptideEvidence_ref=""PE_1"" />
            <cvParam cvRef=""PSI-MS"" accession=""MS:1002354"" name=""PSM-level q-value"" value=""0.01"" />
            <userParam name=""engine score"" value=""12.5"" type=""xsd:double"" />
          </SpectrumIdentificationItem>
          <SpectrumIdentificationItem id=""SII_2"" chargeState=""3"" experimentalMassToCharge=""333.5"" rank=""2"" passThreshold=""false"">
            <PeptideEvidenceRef peptideEvidence_ref=""PE_2"" />
            <PeptideEvidenceRef peptideEvidence_ref=""PE_3"" />
          </SpectrumIdentificationItem>
          <cvParam cvRef=""PSI-MS"" accession=""MS:1000796"" name=""spectrum title"" value=""the title"" />
          <userParam name=""result note"" value=""kept"" />
        </SpectrumIdentificationResult>
      </SpectrumIdentificationList>
      <SpectrumIdentificationList id=""SIL_2"">
        <SpectrumIdentificationResult id=""SIR_2"" spectrumID=""index=3"" spectraData_ref=""SD_MISSING"">
          <SpectrumIdentificationItem id=""SII_3"" chargeState=""1"" experimentalMassToCharge=""700"" peptide_ref=""PEP_MISSING"" rank=""1"" passThreshold=""true"">
            <PeptideEvidenceRef peptideEvidence_ref=""PE_MISSING"" />
          </SpectrumIdentificationItem>
        </SpectrumIdentificationResult>
      </SpectrumIdentificationList>
    </AnalysisData>
  </DataCollection>
</MzIdentML>";

        /// <summary>
        /// GetSpectrumMatches is written once per schema version, so each version's copy is driven through
        /// every field it projects, every list, and every reference that does or does not resolve.
        /// </summary>
        [TestCaseSource(nameof(EveryNamespace))]
        public void GetSpectrumMatches_EveryVersion_ResolvesEveryReferenceInEveryList(string ns)
        {
            var matches = ReadSynthetic(RichSyntheticDocument(ns), "synthetic_matches.mzid").GetSpectrumMatches().ToList();

            Assert.That(matches.Select(m => m.SpectrumIdentificationItemId), Is.EqualTo(new[] { "SII_1", "SII_2", "SII_3" }));
            var full = matches[0];
            var viaEvidence = matches[1];
            var unresolved = matches[2];

            Assert.Multiple(() =>
            {
                Assert.That(full.SpectrumIdentificationListId, Is.EqualTo("SIL_1"));
                Assert.That(full.SpectrumIdentificationResultId, Is.EqualTo("SIR_1"));
                Assert.That(full.SpectrumId, Is.EqualTo("scan=7"));
                Assert.That(full.SpectraData, Is.EqualTo(new MzidSpectraData("SD_1", "run one", @"C:\data\run1.raw",
                    new CvParam("PSI-MS", "MS:1000563", "Thermo RAW format", ""),
                    new CvParam("PSI-MS", "MS:1000768", "Thermo nativeID format", ""))));
                Assert.That(full.SpectrumTitle, Is.EqualTo("the title"));
                Assert.That(full.Rank, Is.EqualTo(1));
                Assert.That(full.PassThreshold, Is.True);
                Assert.That(full.ChargeState, Is.EqualTo(2));
                Assert.That(full.ExperimentalMassToCharge, Is.EqualTo(500.5));
                Assert.That(full.CalculatedMassToCharge, Is.EqualTo(500.25));
                Assert.That(full.PeptideSequence, Is.EqualTo("PEPTIDE"));
                Assert.That(full.HasSubstitutionModifications, Is.True);
                Assert.That(full.Modifications, Has.Count.EqualTo(2));
                Assert.That(full.Modifications[0].Location, Is.EqualTo(0));
                Assert.That(full.Modifications[0].Residues, Is.Empty);
                Assert.That(full.Modifications[0].MonoisotopicMassDelta, Is.EqualTo(42.010565));
                Assert.That(full.Modifications[0].CvParams.Single().Accession, Is.EqualTo("UNIMOD:1"));
                Assert.That(full.Modifications[1].Residues, Is.EqualTo(new[] { "T" }));
                Assert.That(full.Modifications[1].MonoisotopicMassDelta, Is.Null);
                Assert.That(full.Modifications[1].CvParams[1], Is.EqualTo(
                    new CvParam("PSI-MS", "MS:1001524", "fragment neutral loss", "0", "UO", "UO:0000221", "dalton")));
                Assert.That(full.PeptideEvidence.Single(), Is.EqualTo(new MzidPeptideEvidence("PE_1", false, "P12345", 10, 16, "K", "A")));
                Assert.That(full.ItemCvParams.Single().Value, Is.EqualTo("0.01"));
                Assert.That(full.ItemUserParams.Single(), Is.EqualTo(new MzidUserParam("engine score", "12.5", "xsd:double")));
                Assert.That(full.ResultCvParams.Single().Accession, Is.EqualTo("MS:1000796"));
                Assert.That(full.ResultUserParams.Single(), Is.EqualTo(new MzidUserParam("result note", "kept", null)));

                // no peptide_ref of its own: the first evidence names the peptide; an unknown DBSequence gives no accession
                Assert.That(viaEvidence.PeptideSequence, Is.EqualTo("SAMPLER"));
                Assert.That(viaEvidence.CalculatedMassToCharge, Is.Null);
                Assert.That(viaEvidence.Modifications, Is.Empty);
                Assert.That(viaEvidence.HasSubstitutionModifications, Is.False);
                Assert.That(viaEvidence.PeptideEvidence.Select(e => (e.DBSequenceAccession, e.IsDecoy, e.Start)),
                    Is.EqualTo(new (string, bool, int?)[] { ("DECOY_P12345", true, null), (null, true, null) }));

                // the second list is read, and references that resolve to nothing are reported as absent
                Assert.That(unresolved.SpectrumIdentificationListId, Is.EqualTo("SIL_2"));
                Assert.That(unresolved.SpectraData, Is.Null);
                Assert.That(unresolved.SpectrumTitle, Is.Null);
                Assert.That(unresolved.PeptideSequence, Is.Null);
                Assert.That(unresolved.PeptideEvidence, Is.Empty);
                Assert.That(unresolved.ItemCvParams, Is.Empty);
            });
        }

        /// <summary>
        /// The title is the value of MS:1000796, or of the obsolete MS:1001416 it replaced; an empty value is no
        /// title. Unlike Ms2SpectrumID, nothing falls back to another cvParam.
        /// </summary>
        [TestCase(@"accession=""MS:1001416"" name=""spectrum title"" value=""old title""", "old title")]
        [TestCase(@"accession=""MS:1000796"" name=""spectrum title"" value=""""", null)]
        [TestCase(@"accession=""MS:1001371"" name=""Mascot:identity threshold"" value=""17""", null)]
        public void GetSpectrumMatches_SpectrumTitle(string titleTerm, string expected)
        {
            string document = RichSyntheticDocument(Mzid110).Replace(
                @"accession=""MS:1000796"" name=""spectrum title"" value=""the title""", titleTerm);

            var match = ReadSynthetic(document, "synthetic_title_terms.mzid").GetSpectrumMatches().First();

            Assert.That(match.SpectrumTitle, Is.EqualTo(expected));
        }

        [TestCaseSource(nameof(EveryNamespace))]
        public void GetSpectraData_EveryVersion_ListsEveryEntry(string ns)
        {
            var spectraData = ReadSynthetic(RichSyntheticDocument(ns), "synthetic_spectradata.mzid").GetSpectraData();

            Assert.That(spectraData.Select(d => d.Id), Is.EqualTo(new[] { "SD_1", "SD_2" }));
            Assert.That(spectraData[1].FileFormat, Is.Null);
            Assert.That(spectraData[1].Location, Is.EqualTo("run2.mgf"));
        }

        /// <summary>
        /// A document with no SequenceCollection at all still enumerates its items, with nothing resolved.
        /// </summary>
        [TestCaseSource(nameof(EveryNamespace))]
        public void GetSpectrumMatches_EveryVersion_ToleratesAMissingSequenceCollection(string ns)
        {
            var match = ReadSynthetic(SyntheticDocument(ns, "MS:1000584", "mzML format", ""), "synthetic_noseq.mzid")
                .GetSpectrumMatches().Single();

            Assert.Multiple(() =>
            {
                Assert.That(match.SpectrumId, Is.EqualTo("index=7"));
                Assert.That(match.SpectraData.Id, Is.EqualTo("SD_1"));
                Assert.That(match.PeptideSequence, Is.Null);
                Assert.That(match.ItemUserParams.Single().Name, Is.EqualTo("Percolator q-Value"));
            });
        }

        private static IEnumerable<TestCaseData> MissingSectionsInEveryVersion() =>
            from ns in EveryNamespace
            from body in new[]
            {
                // no DataCollection at all
                "",
                // DataCollection with neither Inputs nor AnalysisData
                "<DataCollection />",
                // Inputs with no SpectraData; a list with no results; a result with no items and no spectraData_ref
                @"<DataCollection>
    <Inputs />
    <AnalysisData>
      <SpectrumIdentificationList id=""SIL_EMPTY"" />
      <SpectrumIdentificationList id=""SIL_1"">
        <SpectrumIdentificationResult id=""SIR_EMPTY"" spectrumID=""index=0"" />
      </SpectrumIdentificationList>
    </AnalysisData>
  </DataCollection>",
            }
            select new TestCaseData(ns, body);

        /// <summary>
        /// Every section GetSpectrumMatches and GetSpectraData walk is optional to the deserializer, and each
        /// version's copy treats an absent one as empty rather than throwing.
        /// </summary>
        [TestCaseSource(nameof(MissingSectionsInEveryVersion))]
        public void GetSpectrumMatches_EveryVersion_TreatsMissingSectionsAsEmpty(string ns, string body)
        {
            string document = $@"<?xml version=""1.0"" encoding=""utf-8""?>
<MzIdentML xmlns=""{ns}"" id=""synthetic"" version=""1.1.0"">
  {body}
</MzIdentML>";
            var ids = ReadSynthetic(document, "synthetic_sections.mzid");

            Assert.That(ids.GetSpectrumMatches(), Is.Empty);
            Assert.That(ids.GetSpectraData(), Is.Empty);
        }

        /// <summary>
        /// Every attribute of a cvParam is optional to the deserializer; an absent one reads as empty, the
        /// CvParam default.
        /// </summary>
        [Test]
        public void GetSpectrumMatches_ACvParamWithOnlyAnAccessionReadsTheRestAsEmpty()
        {
            string document = RichSyntheticDocument(Mzid110).Replace(
                @"<cvParam cvRef=""PSI-MS"" accession=""MS:1002354"" name=""PSM-level q-value"" value=""0.01"" />",
                @"<cvParam accession=""MS:1002354"" />");

            var cv = ReadSynthetic(document, "synthetic_bare_cv.mzid").GetSpectrumMatches().First().ItemCvParams.Single();

            Assert.That(cv, Is.EqualTo(new CvParam("", "MS:1002354", "", "")));
        }

        /// <summary>
        /// A repeated id keeps its first entry, and an entry with no id is not indexed.
        /// </summary>
        [Test]
        public void GetSpectrumMatches_ARepeatedIdResolvesToItsFirstEntry()
        {
            string document = RichSyntheticDocument(Mzid110)
                .Replace(@"<DBSequence id=""DB_2"" accession=""DECOY_P12345""",
                    @"<DBSequence id=""DB_1"" accession=""SHADOWED"" searchDatabase_ref=""SDB_1"" /><DBSequence accession=""NO_ID"" /><DBSequence id=""DB_2"" accession=""DECOY_P12345""");

            var match = ReadSynthetic(document, "synthetic_repeated.mzid").GetSpectrumMatches().First();

            Assert.That(match.PeptideEvidence.Single().DBSequenceAccession, Is.EqualTo("P12345"));
        }

        [Test]
        public void Count_IsGreaterThanZero()
        {
            Assert.That(mzid.Count, Is.GreaterThan(0));
            Assert.That(mzid.Count, Is.EqualTo(65));
        }

        /// <summary>
        /// A 1.3.0 file declares the http://psidev.info/psi/pi/mzIdentML/1.3 namespace, so deserialization
        /// against the 1.1.0, 1.1.1 and 1.2.0 types all fail before the 1.3.0 branch is reached.
        ///
        /// Every accessor carries its own independent version cascade, so each one has to be exercised
        /// separately -- a cascade missing its 1.3.0 branch ends at a null 1.2.0 field and throws
        /// NullReferenceException rather than returning a wrong answer. MatchedIons and MatchedIonCounts
        /// are the two not covered here: no published 1.3.0 example carries a fragmentation table.
        ///
        /// The fixture is examples/1_3examples/crosslinking/noncovalently_assoc_1_3_0_draft.mzid from the
        /// HUPO-PSI repository, which validates against the released 1.3.0 schema.
        /// </summary>
        [Test]
        public void MzIdentML130_IsReadThroughEveryVersionCascade()
        {
            var mzid130 = new MzidIdentifications(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "noncovalently_assoc_1_3_0.mzid"));

            Assert.Multiple(() =>
            {
                Assert.That(mzid130.ParentTolerance.Value, Is.EqualTo(3.0).Within(1e-9));
                Assert.That(mzid130.FragmentTolerance.Value, Is.EqualTo(20.0).Within(1e-9));

                Assert.That(mzid130.Count, Is.EqualTo(1));
                Assert.That(mzid130.NumPSMsFromScan(0), Is.EqualTo(2));
                Assert.That(mzid130.Ms2SpectrumID(0), Is.EqualTo("13773"));

                Assert.That(mzid130.CalculatedMassToCharge(0, 0), Is.EqualTo(1392.567094980103).Within(1e-9));
                Assert.That(mzid130.ExperimentalMassToCharge(0, 0), Is.EqualTo(1392.897440641436).Within(1e-9));
                Assert.That(mzid130.ChargeState(0, 0), Is.EqualTo(3));
                Assert.That(mzid130.IsDecoy(0, 0), Is.False);

                // the fixture declares no q-value, and the accessor reports that as -1
                Assert.That(mzid130.QValue(0, 0), Is.EqualTo(-1).Within(1e-9));

                Assert.That(mzid130.PeptideSequenceWithoutModifications(0, 0), Is.EqualTo("AYALMTDIHWDDCFCR"));
                Assert.That(mzid130.ProteinAccession(0, 0), Is.EqualTo("P15640"));
                Assert.That(mzid130.ProteinFullName(0, 0), Does.StartWith("PUR2_ECOLI"));
                Assert.That(mzid130.StartResidueInProtein(0, 0), Is.EqualTo("401"));
                Assert.That(mzid130.EndResidueInProtein(0, 0), Is.EqualTo("416"));

                Assert.That(mzid130.NumModifications(0, 0), Is.EqualTo(3));
                Assert.That(mzid130.ModificationAcession(0, 0, 0), Is.EqualTo("MS:1003393"));
                Assert.That(mzid130.ModificationValue(0, 0, 0), Is.EqualTo("ox"));
                Assert.That(mzid130.ModificationDictionary(0, 0, 0), Is.EqualTo("PSI-MS"));
                Assert.That(mzid130.ModificationLocation(0, 0, 0), Is.EqualTo(5));
                Assert.That(mzid130.ModificationMass(0, 0, 0), Is.EqualTo(15.99491).Within(1e-9));
            });
        }

        [Test]
        public void ParentTolerance_IsNotNull()
        {
            Assert.That(mzid.ParentTolerance, Is.Not.Null);
            Assert.That(mzid.ParentTolerance.Within(mzid.ParentTolerance.Value,5.0));
        }

        [Test]
        public void FragmentTolerance_IsNotNull()
        {
            Assert.That(mzid.FragmentTolerance, Is.Not.Null);
            Assert.That(mzid.FragmentTolerance.Within(mzid.FragmentTolerance.Value, 20.0));
        }

        [Test]
        public void PeptideSequenceWithoutModifications_ReturnsExpectedValue()
        {
            // Example indices, adjust as needed for your file
            var expectedSequence = "KAPAGGAADAAAK";
            Assert.That(mzid.PeptideSequenceWithoutModifications(0, 0), Is.EqualTo(expectedSequence));
        }

        [Test]
        public void ProteinAccession_ReturnsExpectedValue()
        {
            // Example indices, adjust as needed for your file
            var expectedAccession = "P46672";
            Assert.That(mzid.ProteinAccession(0, 0), Is.EqualTo(expectedAccession));
        }

        [Test]
        public void ProteinFullName_ReturnsExpectedValue()
        {
            // Example indices, adjust as needed for your file
            var expectedFullName = "tRNA-aminoacylation cofactor ARC1";
            Assert.That(mzid.ProteinFullName(0, 0), Is.EqualTo(expectedFullName));
        }

    }
}
