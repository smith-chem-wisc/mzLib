using System.IO;
using MzIdentML;
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
