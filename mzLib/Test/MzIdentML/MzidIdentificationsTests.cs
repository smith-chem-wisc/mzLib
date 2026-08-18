using System.IO;
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

            for (int sir = 0; sir < ids.Count; sir++)
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
        /// - its SpectraData is Mascot MGF, so Ms2SpectrumID falls through to the SpectrumIdentification-
        ///   Result's first cvParam instead of reading spectrumID
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

                // the first cvParam happens to be the spectrum title; spectrumID is "index=12"
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
