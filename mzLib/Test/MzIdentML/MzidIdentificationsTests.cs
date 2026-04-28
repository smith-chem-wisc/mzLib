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
