using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
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
using Test.Omics;
using Omics.SpectralMatch;

namespace Test.Quantification
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class QuantificationTests
    {
        #region Helper Methods

        /// <summary>
        /// Creates a test experimental design with specified number of files and channels
        /// </summary>
        private static TestExperimentalDesign CreateTestExperimentalDesign(int numFiles, int numChannels, bool multipleFractions = false)
        {
            var dict = new Dictionary<string, ISampleInfo[]>();

            for (int i = 0; i < numFiles; i++)
            {
                string fileName = $"file{i + 1}.raw";
                int fraction = multipleFractions ? i : 0;
                var samples = new ISampleInfo[numChannels];

                for (int j = 0; j < numChannels; j++)
                {
                    string condition = j == 0 ? "Reference" : (j == 1 ? "Control" : "Treatment");
                    bool isReference = j == 0;
                    string channelLabel = $"Channel_{j}";
                    double reporterIonMz = 126.0 + j * 0.1;

                    samples[j] = new IsobaricQuantSampleInfo(
                        fileName, condition, 0, 0, fraction, 0,
                        channelLabel, reporterIonMz, isReference);
                }

                dict[fileName] = samples;
            }

            return new TestExperimentalDesign(dict);
        }

        /// <summary>
        /// Creates test proteins with digested peptides
        /// </summary>
        private static List<IBioPolymerGroup> CreateTestProteins(out List<IBioPolymerWithSetMods> peptides)
        {
            var proteins = new List<Protein>
            {
                new Protein("PEPTIDEK", "P1"),
                new Protein("SAMPLERTIDE", "P2"),
                new Protein("SEQUENCEK", "P3")
            };

            peptides = new List<IBioPolymerWithSetMods>();
            var bioPolymerGroups = new List<IBioPolymerGroup>();

            foreach (var protein in proteins)
            {
                var digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5);
                var digestedPeptides = protein.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();

                // Create a protein group for each protein using BioPolymerGroup constructor
                var proteinGroup = new BioPolymerGroup(
                    new HashSet<IBioPolymer> { protein },
                    new HashSet<IBioPolymerWithSetMods>(digestedPeptides),
                    new HashSet<IBioPolymerWithSetMods>(digestedPeptides));

                bioPolymerGroups.Add(proteinGroup);
                peptides.AddRange(digestedPeptides);
            }

            return bioPolymerGroups;
        }

        /// <summary>
        /// Creates test spectral matches with quantification values
        /// </summary>
        private static List<ISpectralMatch> CreateTestSpectralMatches(List<IBioPolymerWithSetMods> peptides, 
            int numFiles, int numChannels)
        {
            var matches = new List<ISpectralMatch>();
            int scanNumber = 1;

            foreach (var peptide in peptides)
            {
                for (int fileIndex = 0; fileIndex < numFiles; fileIndex++)
                {
                    string fileName = $"file{fileIndex + 1}.raw";

                    // Create 2 matches per peptide per file for better coverage
                    for (int matchNum = 0; matchNum < 2; matchNum++)
                    {
                        var quantValues = new double[numChannels];
                        for (int j = 0; j < numChannels; j++)
                        {
                            // Create varying intensities: base intensity * channel multiplier * file multiplier
                            quantValues[j] = (1000.0 + scanNumber * 100) * (j + 1) * (fileIndex + 1);
                        }

                        var match = new MockSpectralMatch(
                            fileName,
                            peptide.FullSequence,
                            peptide.BaseSequence,
                            100.0 - matchNum * 5,
                            scanNumber++,
                            new[] { peptide })
                        {
                            Intensities = quantValues
                        };

                        matches.Add(match);
                    }
                }
            }

            return matches;
        }

        #endregion

        #region Existing Test

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
            var sm1 = new MockSpectralMatch(file1, sampler.FullSequence, sampler.BaseSequence, 100.0, 1, new[] { sampler }) // SAMPLER - File 1
            {
                Intensities = new double[] { 1000.0, 2000.0, 3000.0 } // file1 channels
            };

            var sm2 = new MockSpectralMatch(file2, sampler.FullSequence, sampler.BaseSequence, 95.0, 1, new[] { sampler }) // SAMPLER - File 2
            {
                Intensities = new double[] { 1500.0, 2500.0, 3500.0 } // file2 channels
            };

            var sm3 = new MockSpectralMatch(file1, peptidek.FullSequence, peptidek.BaseSequence, 90.0, 2, new[] { peptidek }) //PEPTIDEK - File 1
            {
                Intensities = new double[] { 500.0, 1000.0, 1500.0 } // file1 channels
            };

            var sm4 = new MockSpectralMatch(file2, peptidek.FullSequence, peptidek.BaseSequence, 85.0, 2, new[] { peptidek }) //PEPTIDEK - File 2
            {
                Intensities = new double[] { 750.0, 1250.0, 1750.0 } // file2 channels
            };

            var sm5 = new MockSpectralMatch(file1, sampler.FullSequence, sampler.BaseSequence, 100.0, 3, new[] { sampler }) // SAMPLER - File 1
            {
                Intensities = new double[] { 1000.0, 2000.0, 3000.0 } // file1 channels
            };

            var sm6 = new MockSpectralMatch(file2, sampler.FullSequence, sampler.BaseSequence, 95.0, 3, new[] { sampler }) // SAMPLER - File 2
            {
                Intensities = new double[] { 1500.0, 2500.0, 3500.0 } // file2 channels
            };

            var sm7 = new MockSpectralMatch(file1, peptidek.FullSequence, peptidek.BaseSequence, 90.0, 4, new[] { peptidek }) //PEPTIDEK - File 1
            {
                Intensities = new double[] { 500.0, 1000.0, 1500.0 } // file1 channels
            };

            var sm8 = new MockSpectralMatch(file2, peptidek.FullSequence, peptidek.BaseSequence, 85.0, 4, new[] { peptidek }) //PEPTIDEK - File 2
            {
                Intensities = new double[] { 750.0, 1250.0, 1750.0 } // file2 channels
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

        #endregion

        #region New QuantificationEngine Tests

        [Test]
        public void QuantificationEngine_CompleteRun_SuccessfullyExecutes()
        {
            // Arrange
            int numFiles = 2;
            int numChannels = 3;

            var experimentalDesign = CreateTestExperimentalDesign(numFiles, numChannels);
            var proteinGroups = CreateTestProteins(out var peptides);
            var spectralMatches = CreateTestSpectralMatches(peptides, numFiles, numChannels);

            var parameters = QuantificationParameters.GetSimpleParameters();
            parameters.WriteRawInformation = false;
            parameters.WritePeptideInformation = false;
            parameters.WriteProteinInformation = false;

            var engine = new QuantificationEngine(
                parameters,
                experimentalDesign,
                spectralMatches,
                peptides.ToList(),
                proteinGroups);

            // Act
            var result = engine.Run();

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Summary, Is.EqualTo("Quantification completed successfully."));

            // Verify the engine processed the data through all steps
            // by checking that we can access the pivot functionality
            var smMatrix = QuantificationEngine.Pivot(spectralMatches, experimentalDesign);
            Assert.That(smMatrix, Is.Not.Null);
            Assert.That(smMatrix.RowKeys.Count, Is.GreaterThan(0));
            Assert.That(smMatrix.ColumnKeys.Count, Is.EqualTo(numFiles * numChannels));

            // Verify PSM to peptide mapping works
            var peptideMap = QuantificationEngine.GetPsmToPeptideMap(smMatrix, peptides.ToList());
            Assert.That(peptideMap, Is.Not.Null);
            Assert.That(peptideMap.Keys.Count, Is.EqualTo(peptides.Count));

            // Verify each peptide has associated PSMs
            foreach (var peptide in peptides)
            {
                Assert.That(peptideMap.ContainsKey(peptide), Is.True);
                Assert.That(peptideMap[peptide].Count, Is.GreaterThan(0));
            }
        }

        [Test]
        public void QuantificationEngine_ValidationTests_ProperlyHandlesErrors()
        {
            // Test 1: Null experimental design
            var proteinGroups = CreateTestProteins(out var peptides);
            var spectralMatches = CreateTestSpectralMatches(peptides, 1, 2);
            var parameters = QuantificationParameters.GetSimpleParameters();

            var engine1 = new QuantificationEngine(
                parameters,
                null,  // null experimental design
                spectralMatches,
                peptides.ToList(),
                proteinGroups);

            var result1 = engine1.Run();
            Assert.That(result1.Summary, Does.Contain("Experimental design is null"));

            // Test 2: Empty spectral matches
            var experimentalDesign = CreateTestExperimentalDesign(1, 2);
            var engine2 = new QuantificationEngine(
                parameters,
                experimentalDesign,
                new List<ISpectralMatch>(),  // empty spectral matches
                peptides.ToList(),
                proteinGroups);

            var result2 = engine2.Run();
            Assert.That(result2.Summary, Does.Contain("No spectral matches"));

            // Test 3: Empty peptides
            var engine3 = new QuantificationEngine(
                parameters,
                experimentalDesign,
                spectralMatches,
                new List<IBioPolymerWithSetMods>(),  // empty peptides
                proteinGroups);

            var result3 = engine3.Run();
            Assert.That(result3.Summary, Does.Contain("No modified biopolymers"));

            // Test 4: Empty protein groups
            var engine4 = new QuantificationEngine(
                parameters,
                experimentalDesign,
                spectralMatches,
                peptides.ToList(),
                new List<IBioPolymerGroup>());  // empty protein groups

            var result4 = engine4.Run();
            Assert.That(result4.Summary, Does.Contain("No biopolymer groups"));
        }

        [Test]
        public void QuantificationEngine_ComplexScenario_MultipleProteinsAndFiles()
        {
            // Arrange - Create a more complex scenario with multiple proteins and fractions
            int numFiles = 3;
            int numChannels = 4;

            var experimentalDesign = CreateTestExperimentalDesign(numFiles, numChannels, multipleFractions: true);

            // Create multiple proteins with overlapping peptides
            var protein1 = new Protein("PEPTIDEKSAMPLER", "P1");
            var protein2 = new Protein("SAMPLERSEQUENCEK", "P2");
            var protein3 = new Protein("SEQUENCEKRTIDE", "P3");

            var allPeptides = new List<IBioPolymerWithSetMods>();
            var proteinGroups = new List<IBioPolymerGroup>();

            var digestionParams = new DigestionParams(maxMissedCleavages: 1, minPeptideLength: 5);

            foreach (var protein in new[] { protein1, protein2, protein3 })
            {
                var peptides = protein.Digest(digestionParams, new List<Modification>(), new List<Modification>()).ToList();
                var proteinGroup = new BioPolymerGroup(
                    new HashSet<IBioPolymer> { protein },
                    new HashSet<IBioPolymerWithSetMods>(peptides),
                    new HashSet<IBioPolymerWithSetMods>(peptides));
                proteinGroups.Add(proteinGroup);
                allPeptides.AddRange(peptides);
            }

            // Create spectral matches with varying intensities across files and channels
            var spectralMatches = CreateTestSpectralMatches(allPeptides, numFiles, numChannels);

            var parameters = QuantificationParameters.GetSimpleParameters();
            parameters.WriteRawInformation = false;
            parameters.WritePeptideInformation = false;
            parameters.WriteProteinInformation = false;

            var engine = new QuantificationEngine(
                parameters,
                experimentalDesign,
                spectralMatches,
                allPeptides,
                proteinGroups);

            // Act
            var result = engine.Run();

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Summary, Is.EqualTo("Quantification completed successfully."));

            // Verify the pivot operation produced the correct dimensions
            var smMatrix = QuantificationEngine.Pivot(spectralMatches, experimentalDesign);
            Assert.That(smMatrix.ColumnKeys.Count, Is.EqualTo(numFiles * numChannels));
            Assert.That(smMatrix.RowKeys.Count, Is.GreaterThan(0));

            // Verify peptide mapping captured all peptides
            var peptideMap = QuantificationEngine.GetPsmToPeptideMap(smMatrix, allPeptides);
            Assert.That(peptideMap.Keys.Count, Is.EqualTo(allPeptides.Distinct().Count()));

            // Verify protein mapping works
            var rollUpStrategy = new SumRollUp();
            var peptideMatrix = rollUpStrategy.RollUp(smMatrix, peptideMap);
            Assert.That(peptideMatrix, Is.Not.Null);
            Assert.That(peptideMatrix.RowKeys.Count, Is.GreaterThan(0));

            // Verify unique peptide to protein mapping
            var proteinMap = QuantificationEngine.GetUniquePeptideToProteinMap(peptideMatrix, proteinGroups);
            Assert.That(proteinMap, Is.Not.Null);
            Assert.That(proteinMap.Keys.Count, Is.EqualTo(proteinGroups.Count));

            // Verify each protein group has at least one mapped peptide
            foreach (var proteinGroup in proteinGroups)
            {
                Assert.That(proteinMap.ContainsKey(proteinGroup), Is.True);
                // Note: Some protein groups might have 0 peptides if all their peptides are shared
            }
        }

        [Test]
        public void RunAndReturnProteinMatrix_VariablePeptidesAndPSMs_QuantifiesCorrectly()
        {
            // Arrange - Create 1 file with 2 channels for simpler intensity tracking
            string file1 = "file1.raw";
            var file1Samples = new ISampleInfo[]
            {
                new IsobaricQuantSampleInfo(file1, "Control", 0, 0, 0, 0, "126", 126.0, false),
                new IsobaricQuantSampleInfo(file1, "Treatment", 0, 0, 0, 0, "127N", 127.1, false)
            };

            var expDesign = new TestExperimentalDesign(new Dictionary<string, ISampleInfo[]>
            {
                { file1, file1Samples }
            });

            // Create 4 proteins with variable peptide counts (0, 1, 2, 4 unique peptides)
            var protein1 = new Protein("ABCDEFGHIJK", "P1"); // Will have peptides (empty protein group)
            var protein2 = new Protein("LMNOPQR", "P2"); // Will have 1 unique peptide
            var protein3 = new Protein("STUVWXYZ", "P3"); // Will have 2 unique peptides
            var protein4 = new Protein("AAAABBBBCCCCDDDD", "P4"); // Will have 4 unique peptides

            // Manually create peptides to control the exact number
            var p2_peptide1 = new PeptideWithSetModifications(protein2, null, 1, 7, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0); // LMNOPQR

            var p3_peptide1 = new PeptideWithSetModifications(protein3, null, 1, 4, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0); // STUV
            var p3_peptide2 = new PeptideWithSetModifications(protein3, null, 5, 8, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0); // WXYZ

            var p4_peptide1 = new PeptideWithSetModifications(protein4, null, 1, 4, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0); // AAAA
            var p4_peptide2 = new PeptideWithSetModifications(protein4, null, 5, 8, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0); // BBBB
            var p4_peptide3 = new PeptideWithSetModifications(protein4, null, 9, 12, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0); // CCCC
            var p4_peptide4 = new PeptideWithSetModifications(protein4, null, 13, 16, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0); // DDDD

            // Create protein groups
            var proteinGroup1 = new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein1 },
                new HashSet<IBioPolymerWithSetMods>(), // 0 peptides
                new HashSet<IBioPolymerWithSetMods>());

            var proteinGroup2 = new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein2 },
                new HashSet<IBioPolymerWithSetMods> { p2_peptide1 }, // 1 peptide
                new HashSet<IBioPolymerWithSetMods> { p2_peptide1 });

            var proteinGroup3 = new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein3 },
                new HashSet<IBioPolymerWithSetMods> { p3_peptide1, p3_peptide2 }, // 2 peptides
                new HashSet<IBioPolymerWithSetMods> { p3_peptide1, p3_peptide2 });

            var proteinGroup4 = new BioPolymerGroup(
                new HashSet<IBioPolymer> { protein4 },
                new HashSet<IBioPolymerWithSetMods> { p4_peptide1, p4_peptide2, p4_peptide3, p4_peptide4 }, // 4 peptides (unique + shared)
                new HashSet<IBioPolymerWithSetMods> { p4_peptide2, p4_peptide3, p4_peptide4 }); // 3 unique peptides

            var proteinGroups = new List<IBioPolymerGroup> { proteinGroup1, proteinGroup2, proteinGroup3, proteinGroup4 };

            // All peptides except those from protein1
            var allPeptides = new List<IBioPolymerWithSetMods>
            {
                p2_peptide1, p3_peptide1, p3_peptide2, p4_peptide1, p4_peptide2, p4_peptide3, p4_peptide4
            };

            // Create spectral matches with variable PSM counts per peptide
            var spectralMatches = new List<ISpectralMatch>();

            // P2 peptide 1: 1 PSM with intensities [100, 200]
            spectralMatches.Add(new MockSpectralMatch(file1, p2_peptide1.FullSequence, p2_peptide1.BaseSequence, 100.0, 1, new[] { p2_peptide1 })
            {
                Intensities = new double[] { 100.0, 200.0 }
            });

            // P3 peptide 1: 2 PSMs with intensities [50, 100] and [50, 100] = sum [100, 200]
            spectralMatches.Add(new MockSpectralMatch(file1, p3_peptide1.FullSequence, p3_peptide1.BaseSequence, 95.0, 2, new[] { p3_peptide1 })
            {
                Intensities = new double[] { 50.0, 100.0 }
            });
            spectralMatches.Add(new MockSpectralMatch(file1, p3_peptide1.FullSequence, p3_peptide1.BaseSequence, 95.0, 3, new[] { p3_peptide1 })
            {
                Intensities = new double[] { 50.0, 100.0 }
            });

            // P3 peptide 2: 3 PSMs with intensities [30, 60], [30, 60], [40, 80] = sum [100, 200]
            spectralMatches.Add(new MockSpectralMatch(file1, p3_peptide2.FullSequence, p3_peptide2.BaseSequence, 90.0, 4, new[] { p3_peptide2 })
            {
                Intensities = new double[] { 30.0, 60.0 }
            });
            spectralMatches.Add(new MockSpectralMatch(file1, p3_peptide2.FullSequence, p3_peptide2.BaseSequence, 90.0, 5, new[] { p3_peptide2 })
            {
                Intensities = new double[] { 30.0, 60.0 }
            });
            spectralMatches.Add(new MockSpectralMatch(file1, p3_peptide2.FullSequence, p3_peptide2.BaseSequence, 90.0, 6, new[] { p3_peptide2 })
            {
                Intensities = new double[] { 40.0, 80.0 }
            });

            // P4 peptide 1: 1 PSM [25, 50] (This one is non-unique, and the intensities shouldn't be included in the final protein quant)
            spectralMatches.Add(new MockSpectralMatch(file1, p4_peptide1.FullSequence, p4_peptide1.BaseSequence, 85.0, 7, new[] { p4_peptide1 })
            {
                Intensities = new double[] { 25.0, 50.0 } // Should be ignored
            });

            // P4 peptide 2: 2 PSMs [25, 50] and [25, 50] = sum [50, 100]
            spectralMatches.Add(new MockSpectralMatch(file1, p4_peptide2.FullSequence, p4_peptide2.BaseSequence, 85.0, 8, new[] { p4_peptide2 })
            {
                Intensities = new double[] { 25.0, 50.0 }
            });
            spectralMatches.Add(new MockSpectralMatch(file1, p4_peptide2.FullSequence, p4_peptide2.BaseSequence, 85.0, 9, new[] { p4_peptide2 })
            {
                Intensities = new double[] { 25.0, 50.0 }
            });

            // P4 peptide 3: 3 PSMs [10, 20], [10, 20], [5, 10] = sum [25, 50]
            spectralMatches.Add(new MockSpectralMatch(file1, p4_peptide3.FullSequence, p4_peptide3.BaseSequence, 80.0, 10, new[] { p4_peptide3 })
            {
                Intensities = new double[] { 10.0, 20.0 }
            });
            spectralMatches.Add(new MockSpectralMatch(file1, p4_peptide3.FullSequence, p4_peptide3.BaseSequence, 80.0, 11, new[] { p4_peptide3 })
            {
                Intensities = new double[] { 10.0, 20.0 }
            });
            spectralMatches.Add(new MockSpectralMatch(file1, p4_peptide3.FullSequence, p4_peptide3.BaseSequence, 80.0, 12, new[] { p4_peptide3 })
            {
                Intensities = new double[] { 5.0, 10.0 }
            });

            // P4 peptide 4: 4 PSMs [10, 20] each = sum [40, 80]
            spectralMatches.Add(new MockSpectralMatch(file1, p4_peptide4.FullSequence, p4_peptide4.BaseSequence, 75.0, 13, new[] { p4_peptide4 })
            {
                Intensities = new double[] { 10.0, 20.0 }
            });
            spectralMatches.Add(new MockSpectralMatch(file1, p4_peptide4.FullSequence, p4_peptide4.BaseSequence, 75.0, 14, new[] { p4_peptide4 })
            {
                Intensities = new double[] { 10.0, 20.0 }
            });
            spectralMatches.Add(new MockSpectralMatch(file1, p4_peptide4.FullSequence, p4_peptide4.BaseSequence, 75.0, 15, new[] { p4_peptide4 })
            {
                Intensities = new double[] { 10.0, 20.0 }
            });
            spectralMatches.Add(new MockSpectralMatch(file1, p4_peptide4.FullSequence, p4_peptide4.BaseSequence, 75.0, 16, new[] { p4_peptide4 })
            {
                Intensities = new double[] { 10.0, 20.0 }
            });

            // Add an orphaned spectral match (not in allPeptides list)
            var orphanPeptide = new PeptideWithSetModifications(protein1, null, 1, 5, CleavageSpecificity.Full, "", 0, new Dictionary<int, Modification>(), 0);
            spectralMatches.Add(new MockSpectralMatch(file1, orphanPeptide.FullSequence, orphanPeptide.BaseSequence, 70.0, 17, new[] { orphanPeptide })
            {
                Intensities = new double[] { 999.0, 999.0 } // Should be ignored
            });

            var parameters = QuantificationParameters.GetSimpleParameters();
            parameters.WriteRawInformation = false;
            parameters.WritePeptideInformation = false;
            parameters.WriteProteinInformation = false;

            var engine = new QuantificationEngine(
                parameters,
                expDesign,
                spectralMatches,
                allPeptides,
                proteinGroups);

            // Act
            var proteinMatrix = engine.RunAndReturnProteinMatrix();

            // Assert
            Assert.That(proteinMatrix, Is.Not.Null);
            Assert.That(proteinMatrix.RowKeys.Count, Is.EqualTo(4));
            Assert.That(proteinMatrix.ColumnKeys.Count, Is.EqualTo(2)); // 2 channels

            // Verify protein 1 (0 peptides) has 0 intensity
            var p1Row = proteinMatrix.GetRow(proteinGroup1);
            Assert.That(p1Row[0], Is.EqualTo(0.0)); // Channel 1
            Assert.That(p1Row[1], Is.EqualTo(0.0)); // Channel 2

            // Verify protein 2 (1 peptide, 1 PSM): [100, 200]
            var p2Row = proteinMatrix.GetRow(proteinGroup2);
            Assert.That(p2Row[0], Is.EqualTo(100.0)); // Channel 1
            Assert.That(p2Row[1], Is.EqualTo(200.0)); // Channel 2

            // Verify protein 3 (2 peptides): peptide1[100, 200] + peptide2[100, 200] = [200, 400]
            var p3Row = proteinMatrix.GetRow(proteinGroup3);
            Assert.That(p3Row[0], Is.EqualTo(200.0)); // Channel 1
            Assert.That(p3Row[1], Is.EqualTo(400.0)); // Channel 2

            // Verify protein 4 (3 unique peptides): [25, 50]*0 + [50, 100] + [25, 50] + [40, 80] = [115, 230]
            var p4Row = proteinMatrix.GetRow(proteinGroup4);
            Assert.That(p4Row[0], Is.EqualTo(115.0)); // Channel 1
            Assert.That(p4Row[1], Is.EqualTo(230.0)); // Channel 2

            // Verify the orphaned PSM was not included in any protein quantification
            // (already verified by checking specific protein intensities match expected values)
        }

        #endregion

        private class TestExperimentalDesign : IExperimentalDesign
        {
            public Dictionary<string, ISampleInfo[]> FileNameSampleInfoDictionary { get; }

            public TestExperimentalDesign(Dictionary<string, ISampleInfo[]> dict)
            {
                FileNameSampleInfoDictionary = dict;
            }
        }
    }
}
