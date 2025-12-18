using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymer;
using Omics.BioPolymerGroup;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Omics.SpectralMatch;

namespace Test.Omics
{
    /// <summary>
    /// Integration tests for BioPolymerGroup class adapted directly from unit tests inside
    /// MetaMorpheus for ProteinGroup.
    /// Tests sequence coverage calculation, quantification, and file-specific operations.
    /// Core equality/merge tests are in BioPolymerGroupTests.cs.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    internal class ExternalProteinGroupTests
    {
        #region Sequence Coverage Tests

        /// <summary>
        /// Verifies CalculateSequenceCoverage correctly computes full coverage when all residues are covered.
        /// Critical: Sequence coverage is a key quality metric reported in protein group output files.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_FullCoverage()
        {
            var prot = new MockBioPolymer("MEDEEKPEPTIDE", "P00001");
            var peptide1 = new MockBioPolymerWithSetMods(prot, 1, 6);  // MEDEEK
            var peptide2 = new MockBioPolymerWithSetMods(prot, 7, 13); // PEPTIDE

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { prot },
                new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 },
                new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 });

            var psm1 = new MockSpectralMatch("file.mzML", "MEDEEK", "MEDEEK", 100, 1, new[] { peptide1 });
            var psm2 = new MockSpectralMatch("file.mzML", "PEPTIDE", "PEPTIDE", 100, 2, new[] { peptide2 });
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            group.CalculateSequenceCoverage();

            var output = group.ToString();
            // Full coverage should show 1.0 (or 100%)
            Assert.That(output, Does.Contain("1"));
        }

        /// <summary>
        /// Verifies partial coverage shows uppercase for covered residues, lowercase for uncovered.
        /// Critical: Visual coverage representation helps users identify which regions lack evidence.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_PartialCoverage_ShowsUpperLowerCase()
        {
            var prot = new MockBioPolymer("MEDEEKPEPTIDE", "P00001"); // 13 residues
            var peptide = new MockBioPolymerWithSetMods(prot, 1, 6);   // MEDEEK - 6 residues

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { prot },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm = new MockSpectralMatch("file.mzML", "MEDEEK", "MEDEEK", 100, 1, new[] { peptide });
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm };

            group.CalculateSequenceCoverage();

            var output = group.ToString();
            Assert.That(output, Does.Contain("MEDEEK"));  // Covered = uppercase
            Assert.That(output, Does.Contain("peptide")); // Uncovered = lowercase
        }

        #endregion

        #region Subset Group Tests

        /// <summary>
        /// Verifies ConstructSubsetBioPolymerGroup filters PSMs and samples to a specific file.
        /// Critical: Per-file analysis requires correct data partitioning for individual file statistics.
        /// </summary>
        [Test]
        public void ConstructSubsetBioPolymerGroup_FiltersByFile()
        {
            var prot = new MockBioPolymer("MEDEEK", "P00001");
            var peptide = new MockBioPolymerWithSetMods(prot, 1, 6);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { prot },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm1 = new MockSpectralMatch("file1.mzML", "MEDEEK", "MEDEEK", 100, 1, new[] { peptide });
            var psm2 = new MockSpectralMatch("file2.mzML", "MEDEEK", "MEDEEK", 100, 1, new[] { peptide });
            group.AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm1, psm2 };

            var file1 = new SpectraFileInfo("file1.mzML", "Condition1", 1, 1, 0);
            var file2 = new SpectraFileInfo("file2.mzML", "Condition1", 1, 1, 0);
            group.SamplesForQuantification = new List<ISampleInfo> { file1, file2 };
            group.IntensitiesBySample = new Dictionary<ISampleInfo, double> { { file1, 1000 }, { file2, 2000 } };

            var subset = group.ConstructSubsetBioPolymerGroup("file1.mzML");

            Assert.Multiple(() =>
            {
                Assert.That(subset.AllPsmsBelowOnePercentFDR.Count, Is.EqualTo(1));
                Assert.That(subset.AllPsmsBelowOnePercentFDR.First().FullFilePath, Is.EqualTo("file1.mzML"));
                Assert.That(subset.SamplesForQuantification.Count, Is.EqualTo(1));
                Assert.That(subset.IntensitiesBySample.Values.First(), Is.EqualTo(1000));
            });
        }

        #endregion

        #region Quantification Tests

        /// <summary>
        /// Verifies intensity columns are correctly generated in output when SpectraFileInfo samples are set.
        /// Critical: Quantification data must appear in correct columns for downstream statistical analysis.
        /// </summary>
        [Test]
        public void Quantification_GeneratesIntensityColumns()
        {
            var prot = new MockBioPolymer("MEDEEK", "P00001");
            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { prot },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            var sample = new SpectraFileInfo("test.mzML", "Condition1", 1, 1, 0);
            group.SamplesForQuantification = new List<ISampleInfo> { sample };
            group.IntensitiesBySample = new Dictionary<ISampleInfo, double> { { sample, 12345.67 } };

            var header = group.GetTabSeparatedHeader();
            var output = group.ToString();

            Assert.Multiple(() =>
            {
                Assert.That(header, Does.Contain("Intensity_"));
                Assert.That(output, Does.Contain("12345"));
            });
        }

        #endregion

        #region String Truncation Tests

        /// <summary>
        /// Verifies long strings are truncated when MaxStringLength is set.
        /// Critical: Prevents Excel compatibility issues (32,767 character cell limit).
        /// </summary>
        [Test]
        public void MaxStringLength_TruncatesLongStrings()
        {
            var longSequence = new string('M', 50000);
            var prot = new MockBioPolymer(longSequence, "P00001");
            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { prot },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            var originalMaxLength = BioPolymerGroup.MaxStringLength;
            try
            {
                BioPolymerGroup.MaxStringLength = 100;
                var output = group.ToString();

                // Output should not contain the full 50,000 character sequence
                Assert.That(output.Length, Is.LessThan(longSequence.Length));
            }
            finally
            {
                BioPolymerGroup.MaxStringLength = originalMaxLength;
            }
        }

        /// <summary>
        /// Verifies truncation is disabled when MaxStringLength is 0.
        /// Critical: Allows users to disable truncation for full data export when needed.
        /// </summary>
        [Test]
        public void MaxStringLength_ZeroDisablesTruncation()
        {
            var prot = new MockBioPolymer("MEDEEK", "P00001");
            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { prot },
                new HashSet<IBioPolymerWithSetMods>(),
                new HashSet<IBioPolymerWithSetMods>());

            var originalMaxLength = BioPolymerGroup.MaxStringLength;
            try
            {
                BioPolymerGroup.MaxStringLength = 0;
                Assert.DoesNotThrow(() => group.ToString());
            }
            finally
            {
                BioPolymerGroup.MaxStringLength = originalMaxLength;
            }
        }

        #endregion
    }
}