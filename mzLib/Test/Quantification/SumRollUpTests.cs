using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.ExperimentalDesign;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Quantification;
using Quantification.Strategies;

namespace Test.Quantification
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class SumRollUpTests
    {
        private class TestExperimentalDesign : IExperimentalDesign
        {
            public Dictionary<string, ISampleInfo[]> FileNameSampleInfoDictionary { get; }

            public TestExperimentalDesign(Dictionary<string, ISampleInfo[]> dict)
            {
                FileNameSampleInfoDictionary = dict;
            }
        }

        [Test]
        public void RollUpSpectralMatches_SumsIntensitiesCorrectly()
        {
            // Setup: Create 2 files with 3 TMT channels each
            string file1 = "file1.raw";
            string file2 = "file2.raw";

            // The files represent two fractions, the channels are the same in both files
            var file1Samples = new ISampleInfo[]
            {
                // Showing one example with all the arguments labeled for clarity
                new IsobaricQuantSampleInfo(file1,
                    condition: "Reference", 
                    biologicalReplicate: 0, 
                    technicalReplicate: 0, 
                    fraction: 0, 
                    plexId: 0, 
                    channelLabel: "126", 
                    reporterIonMz: 126.0, 
                    isReferenceChannel: true),
                new IsobaricQuantSampleInfo(file1, "Control", 0, 0, 0, 0, "127N", 127.1, false),
                new IsobaricQuantSampleInfo(file1, "Treatment", 0, 0, 0, 0, "127C", 127.2, false)
            };

            var file2Samples = new ISampleInfo[]
            {
                new IsobaricQuantSampleInfo(file2, "Reference", 0, 0, fraction: 1, 0, "126", 126.0, true),
                new IsobaricQuantSampleInfo(file2, "Control", 0, 0, fraction: 1, 0, "127N", 127.1, false),
                new IsobaricQuantSampleInfo(file2, "Treatment", 0, 0, fraction: 1, 0, "127C", 127.2, false)
            };

            var expDesign = new TestExperimentalDesign(new Dictionary<string, ISampleInfo[]>
            {
                { file1, file1Samples },
                { file2, file2Samples }
            });

            // Create 2 peptides
            var protein = new Protein("SAMPLERPEPTIDEK", "P1");
            var sampler = new PeptideWithSetModifications(protein, null, 1, 7, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            var peptidek = new PeptideWithSetModifications(protein, null, 8, 15, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);

            var peptides = new List<IBioPolymerWithSetMods> { sampler, peptidek };

            // Create 8 spectral matches:
            // - 4 for peptide1 (two from each file)
            // - 4 for peptide2 (two from each file)
            var sm1 = new BaseSpectralMatch(file1, 1, 100.0, sampler.FullSequence, sampler.BaseSequence, new[] { sampler }) // SAMPLER - File 1
            {
                QuantValues = new double[] { 1000.0, 2000.0, 3000.0 } // file1 channels
            };

            var sm2 = new BaseSpectralMatch(file2, 1, 95.0, sampler.FullSequence, sampler.BaseSequence, new[] { sampler }) // SAMPLER - File 2
            {
                QuantValues = new double[] { 1500.0, 2500.0, 3500.0 } // file2 channels
            };

            var sm3 = new BaseSpectralMatch(file1, 2, 90.0, peptidek.FullSequence, peptidek.BaseSequence, new[] { peptidek }) //PEPTIDEK - File 1
            {
                QuantValues = new double[] { 500.0, 1000.0, 1500.0 } // file1 channels
            };

            var sm4 = new BaseSpectralMatch(file2, 2, 85.0, peptidek.FullSequence, peptidek.BaseSequence, new[] { peptidek }) //PEPTIDEK - File 2
            {
                QuantValues = new double[] { 750.0, 1250.0, 1750.0 } // file2 channels
            };

            var sm5 = new BaseSpectralMatch(file1, 3, 100.0, sampler.FullSequence, sampler.BaseSequence, new[] { sampler }) // SAMPLER - File 1
            {
                QuantValues = new double[] { 1000.0, 2000.0, 3000.0 } // file1 channels
            };

            var sm6 = new BaseSpectralMatch(file2, 3, 95.0, sampler.FullSequence, sampler.BaseSequence, new[] { sampler }) // SAMPLER - File 2
            {
                QuantValues = new double[] { 1500.0, 2500.0, 3500.0 } // file2 channels
            };

            var sm7 = new BaseSpectralMatch(file1, 4, 90.0, peptidek.FullSequence, peptidek.BaseSequence, new[] { peptidek }) //PEPTIDEK - File 1
            {
                QuantValues = new double[] { 500.0, 1000.0, 1500.0 } // file1 channels
            };

            var sm8 = new BaseSpectralMatch(file2, 4, 85.0, peptidek.FullSequence, peptidek.BaseSequence, new[] { peptidek }) //PEPTIDEK - File 2
            {
                QuantValues = new double[] { 750.0, 1250.0, 1750.0 } // file2 channels
            };

            var spectralMatches = new List<ISpectralMatch> { sm1, sm2, sm3, sm4, sm5, sm6, sm7, sm8 };

            // Execute
            var rollUp = new SumRollUp();
            var smMatrix = QuantificationEngine.Pivot(spectralMatches, expDesign);
            var map = QuantificationEngine.GetPsmToPeptideMap(smMatrix, peptides);

            var result = rollUp.RollUp(smMatrix, map);

            // Assert that the SampleInfoArray was assembled correctly
            Assert.That(result.ColumnKeys.Count, Is.EqualTo(6)); // 2 files x 3 channels each = 6 columns
            Assert.That(result.ColumnKeys[0].FullFilePathWithExtension, Is.EqualTo(file1));

            // Check that the entire list is correct
            var expectedSampleInfos = new List<ISampleInfo>
            {
                file1Samples[0], // file1, channel 126
                file1Samples[1], // file1, channel 127N
                file1Samples[2], // file1, channel 127C
                file2Samples[0], // file2, channel 126
                file2Samples[1], // file2, channel 127N
                file2Samples[2]  // file2, channel 127C
            };
            CollectionAssert.AreEqual(expectedSampleInfos, result.ColumnKeys.ToList());

            // Check that the intensities were summed correctly       

            // Expected intensities for SAMPLER:
            //     file1: (1000+1000), (2000+2000), (3000+3000) = 2000, 4000, 6000
            //     file2: (1500+1500), (2500+2500), (3500+3500) = 3000, 5000, 7000
            var samplerRow = result.GetRow(sampler);
            Assert.That(samplerRow[0], Is.EqualTo(2000.0)); // file1, channel 126
            Assert.That(samplerRow[1], Is.EqualTo(4000.0)); // file1, channel 127N
            Assert.That(samplerRow[2], Is.EqualTo(6000.0)); // file1, channel 127C
            Assert.That(samplerRow[3], Is.EqualTo(3000.0)); // file2, channel 126
            Assert.That(samplerRow[4], Is.EqualTo(5000.0)); // file2, channel 127N
            Assert.That(samplerRow[5], Is.EqualTo(7000.0)); // file2, channel 127C

            // Expected intensities for PEPTIDEK:
            //     file1: (500+500), (1000+1000), (1500+1500) = 1000, 2000, 3000
            //     file2: (750+750), (1250+1250), (1750+1750) = 1500, 2500, 3500
            var peptide2Row = result.GetRow(peptidek);
            Assert.That(peptide2Row[0], Is.EqualTo(1000.0));  // file1, channel 126
            Assert.That(peptide2Row[1], Is.EqualTo(2000.0)); // file1, channel 127N
            Assert.That(peptide2Row[2], Is.EqualTo(3000.0)); // file1, channel 127C
            Assert.That(peptide2Row[3], Is.EqualTo(1500.0));  // file2, channel 126
            Assert.That(peptide2Row[4], Is.EqualTo(2500.0)); // file2, channel 127N
            Assert.That(peptide2Row[5], Is.EqualTo(3500.0)); // file2, channel 127C
        }
    }
}
