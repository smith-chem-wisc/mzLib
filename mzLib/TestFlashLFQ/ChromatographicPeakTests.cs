using FlashLFQ;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TestFlashLFQ
{
    public class ChromatographicPeakTests
    {
        private ChromatographicPeak CreateChromatographicPeak()
        {
            // Create a sample SpectraFileInfo
            SpectraFileInfo spectraFileInfo = new SpectraFileInfo("sampleFile", "A", 1, 1, 1);

            // Create a sample Identification
            Identification identification = new Identification(spectraFileInfo, "MPEPTIDE", "M[Oxidation]PEPTIDE", 100, 10, 2, new List<ProteinGroup>());

            // Create a ChromatographicPeak instance
            ChromatographicPeak chromatographicPeak = new ChromatographicPeak(identification, false, spectraFileInfo);

            IndexedMassSpectralPeak peak1 = new IndexedMassSpectralPeak(100, 300, 1, 9.5);
            IndexedMassSpectralPeak peak2 = new IndexedMassSpectralPeak(100, 300, 1, 10.5);

            // Add sample IsotopicEnvelopes
            chromatographicPeak.IsotopicEnvelopes = new List<IsotopicEnvelope>()
            {
                new IsotopicEnvelope(peak1, 2, 300, 1),
                new IsotopicEnvelope(peak2, 2, 300, 1)
            };

            return chromatographicPeak;
        }


        [Test]
        public void TestResolveIdentifications()
        {
            // Arrange
            ChromatographicPeak chromatographicPeak = CreateChromatographicPeak();

            // Act
            chromatographicPeak.ResolveIdentifications();

            // Assert
            Assert.AreEqual(1, chromatographicPeak.NumIdentificationsByBaseSeq);
            Assert.AreEqual(1, chromatographicPeak.NumIdentificationsByFullSeq);
        }

        [Test]
        public void TestToString()
        {
            // Arrange
            ChromatographicPeak chromatographicPeak = CreateChromatographicPeak();

            // Act
            string result = chromatographicPeak.ToString();

            // Assert
            string expected = "sampleFile\tMPEPTIDE\tM[Oxidation]PEPTIDE\t\t\t100\t10\t2\t51.007276466879\t0\t-\t-\t-\t-\t-\t0\tMSMS\t\t\t1\t1\t1\t0\tNaN\tFalse\tFalse";
            Assert.AreEqual(expected, result);
        }
    }
}
