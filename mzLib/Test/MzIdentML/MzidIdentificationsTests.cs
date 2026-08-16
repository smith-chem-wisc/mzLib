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
        /// against the 1.1.0, 1.1.1 and 1.2.0 types all fail before the 1.3.0 branch is reached. Reading a
        /// value through several different accessors is the point: each one has its own version cascade, so
        /// a missing 1.3.0 branch in any of them surfaces as a NullReferenceException rather than a wrong
        /// answer. The fixture is the HUPO-PSI non-covalent association example.
        /// </summary>
        [Test]
        public void MzIdentML130_IsReadThroughTheVersionCascade()
        {
            var mzid130 = new MzidIdentifications(
                Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "noncovalently_assoc_1_3_0.mzid"));

            Assert.That(mzid130.Count, Is.EqualTo(1));
            Assert.That(mzid130.PeptideSequenceWithoutModifications(0, 0), Is.EqualTo("AYALMTDIHWDDCFCR"));
            Assert.That(mzid130.ChargeState(0, 0), Is.EqualTo(3));
            Assert.That(mzid130.ExperimentalMassToCharge(0, 0), Is.EqualTo(1392.897440641436).Within(1e-9));
            Assert.That(mzid130.CalculatedMassToCharge(0, 0), Is.EqualTo(1392.567094980103).Within(1e-9));
            Assert.That(mzid130.ProteinAccession(0, 0), Is.EqualTo("P15640"));
            Assert.That(mzid130.IsDecoy(0, 0), Is.False);
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
