using FlashLFQ;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.Linq;
using MassSpectrometry;
using IsotopicEnvelope = FlashLFQ.IsotopicEnvelope;

namespace Test.FlashLFQ
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
            ChromatographicPeak chromatographicPeak = new ChromatographicPeak(identification, spectraFileInfo);

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

        [Test]
        public void TestDetectionType()
        {
            // create the sample ChromatographicPeaks with different DetectionTypes
            SpectraFileInfo spectraFileInfo = new SpectraFileInfo("sampleFile", "A", 1, 1, 1);
            Identification identification = new Identification(spectraFileInfo, "MPEPTIDE", "Mz[Oxidation]PEPTIDE", 100, 10, 2, new List<ProteinGroup>());
            IndexedMassSpectralPeak peak1 = new IndexedMassSpectralPeak(100, 300, 1, 9.5);
            IndexedMassSpectralPeak peak2 = new IndexedMassSpectralPeak(100, 300, 1, 10.5);

            ChromatographicPeak peak_MSMS = new ChromatographicPeak(identification, spectraFileInfo);
            MbrChromatographicPeak peak_MBR = new MbrChromatographicPeak(identification, spectraFileInfo, 1, false);
            ChromatographicPeak peak_IsoTecker_MBR = new ChromatographicPeak(identification, spectraFileInfo, DetectionType.IsoTrack_MBR);
            ChromatographicPeak peak_IsoTecker_Amb = new ChromatographicPeak(identification, spectraFileInfo, DetectionType.IsoTrack_Ambiguous);

            peak_MSMS.IsotopicEnvelopes = new List<IsotopicEnvelope>()
            {
                new IsotopicEnvelope(peak1, 2, 300, 1),
                new IsotopicEnvelope(peak2, 2, 300, 1)
            };
            peak_MBR.IsotopicEnvelopes = new List<IsotopicEnvelope>()
            {
                new IsotopicEnvelope(peak1, 2, 300, 1),
                new IsotopicEnvelope(peak2, 2, 300, 1)
            };
            peak_IsoTecker_MBR.IsotopicEnvelopes = new List<IsotopicEnvelope>()
            {
                new IsotopicEnvelope(peak1, 2, 300, 1),
                new IsotopicEnvelope(peak2, 2, 300, 1)
            };
            peak_IsoTecker_Amb.IsotopicEnvelopes = new List<IsotopicEnvelope>()
            {
                new IsotopicEnvelope(peak1, 2, 300, 1),
                new IsotopicEnvelope(peak2, 2, 300, 1)
            };

            // Test the DetectionType of the peaks
            Assert.AreEqual(peak_MSMS.DetectionType, DetectionType.MSMS);
            Assert.AreEqual(peak_MBR.DetectionType, DetectionType.MBR);
            Assert.AreEqual(peak_IsoTecker_MBR.DetectionType, DetectionType.IsoTrack_MBR);
            Assert.AreEqual(peak_IsoTecker_Amb.DetectionType, DetectionType.IsoTrack_Ambiguous);
            
            List<string> str_MSMS = peak_MSMS.ToString().Split("\t").ToList();
            List<string> str_MBR = peak_MBR.ToString().Split("\t").ToList();
            List<string> str_IsoTecker_MBR = peak_IsoTecker_MBR.ToString().Split("\t").ToList();
            List<string> str_IsoTecker_Amb = peak_IsoTecker_Amb.ToString().Split("\t").ToList();

            // Test the DetectionType output
            Assert.AreEqual(str_MSMS[16], "MSMS");
            Assert.AreEqual(str_MBR[16], "MBR");
            Assert.AreEqual(str_IsoTecker_MBR[16], "MBR_IsoTrack");
            Assert.AreEqual(str_IsoTecker_Amb[16], "IsoTrack_Ambiguous");

            // The MS2Retention time will only be printed for MSMS peaks
            Assert.AreEqual(str_MSMS[6], "10");
            Assert.AreEqual(str_MBR[6], "");
            Assert.AreEqual(str_IsoTecker_MBR[6], "");
            Assert.AreEqual(str_IsoTecker_Amb[6], "");

        }
    }
}
