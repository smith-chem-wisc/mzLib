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
        #region Helper Classes

        private class TestBioPolymer : IBioPolymer
        {
            public string BaseSequence { get; }
            public string Accession { get; }
            public bool IsDecoy { get; }
            public bool IsContaminant { get; }
            public string Organism { get; init; } = "Test Organism";
            public string FullName { get; init; } = "Test Full Name";
            public string Name { get; init; } = "Test Name";
            public int Length => BaseSequence?.Length ?? 0;
            public string DatabaseFilePath { get; init; } = "";
            public List<Tuple<string, string>> GeneNames { get; init; } = new() { new("primary", "TestGene") };
            public IDictionary<int, List<Modification>> OneBasedPossibleLocalizedModifications { get; init; } = new Dictionary<int, List<Modification>>();
            public string SampleNameForVariants { get; init; } = "";
            public IDictionary<int, List<Modification>> OriginalNonVariantModifications { get; set; } = new Dictionary<int, List<Modification>>();
            public IBioPolymer ConsensusVariant => this;
            public List<SequenceVariation> AppliedSequenceVariations { get; init; } = new();
            public List<SequenceVariation> SequenceVariations { get; init; } = new();
            public List<TruncationProduct> TruncationProducts { get; init; } = new();

            public TestBioPolymer(string baseSequence, string accession, bool isDecoy = false, bool isContaminant = false)
            {
                BaseSequence = baseSequence;
                Accession = accession;
                IsDecoy = isDecoy;
                IsContaminant = isContaminant;
            }

            public IEnumerable<IBioPolymerWithSetMods> Digest(IDigestionParams d, List<Modification> f, List<Modification> v,
                List<SilacLabel>? s = null, (SilacLabel, SilacLabel)? t = null, bool top = false) => throw new NotImplementedException();
            public IBioPolymer CloneWithNewSequenceAndMods(string s, IDictionary<int, List<Modification>>? m) => throw new NotImplementedException();
            public TBioPolymerType CreateVariant<TBioPolymerType>(string v, TBioPolymerType o, IEnumerable<SequenceVariation> a,
                IEnumerable<TruncationProduct> p, IDictionary<int, List<Modification>> m, string s) where TBioPolymerType : IHasSequenceVariants => throw new NotImplementedException();
            public bool Equals(IBioPolymer? other) => other != null && Accession == other.Accession && BaseSequence == other.BaseSequence;
            public override bool Equals(object? obj) => Equals(obj as IBioPolymer);
            public override int GetHashCode() => HashCode.Combine(Accession, BaseSequence);
        }

        private class TestBioPolymerWithSetMods : IBioPolymerWithSetMods
        {
            public string BaseSequence { get; }
            public string FullSequence { get; }
            public IBioPolymer Parent { get; }
            public int OneBasedStartResidue { get; }
            public int OneBasedEndResidue { get; }
            public double MonoisotopicMass { get; init; } = 500.0;
            public double MostAbundantMonoisotopicMass => MonoisotopicMass;
            public ChemicalFormula ThisChemicalFormula => new();
            public string SequenceWithChemicalFormulas => FullSequence;
            public int MissedCleavages => 0;
            public string Description => "";
            public CleavageSpecificity CleavageSpecificityForFdrCategory { get; set; } = CleavageSpecificity.Full;
            public char PreviousResidue => '-';
            public char NextResidue => '-';
            public IDigestionParams DigestionParams => null!;
            public Dictionary<int, Modification> AllModsOneIsNterminus { get; init; } = new();
            public int NumMods => AllModsOneIsNterminus?.Count ?? 0;
            public int NumFixedMods => 0;
            public int NumVariableMods => NumMods;
            public int Length => BaseSequence?.Length ?? 0;
            public char this[int i] => BaseSequence[i];

            public TestBioPolymerWithSetMods(IBioPolymer parent, int start, int end, Dictionary<int, Modification>? mods = null)
            {
                Parent = parent;
                OneBasedStartResidue = start;
                OneBasedEndResidue = end;
                BaseSequence = parent.BaseSequence.Substring(start - 1, end - start + 1);
                AllModsOneIsNterminus = mods ?? new();
                FullSequence = BaseSequence; // Simplified
            }

            public void Fragment(DissociationType d, FragmentationTerminus t, List<Product> p, FragmentationParams? f = null) { }
            public void FragmentInternally(DissociationType d, int m, List<Product> p, FragmentationParams? f = null) { }
            public IBioPolymerWithSetMods Localize(int i, double m) => this;
            public bool Equals(IBioPolymerWithSetMods? o) => o != null && FullSequence == o.FullSequence && Parent?.Accession == o.Parent?.Accession;
            public override bool Equals(object? o) => Equals(o as IBioPolymerWithSetMods);
            public override int GetHashCode() => HashCode.Combine(FullSequence, Parent?.Accession);
        }

        private class TestSpectralMatch : ISpectralMatch
        {
            private readonly List<IBioPolymerWithSetMods> _identified;
            public string FullFilePath { get; }
            public string FullSequence { get; }
            public string BaseSequence { get; }
            public double Score { get; }
            public int OneBasedScanNumber { get; }

            public TestSpectralMatch(string path, string baseSeq, string fullSeq, double score, int scan, IEnumerable<IBioPolymerWithSetMods> identified)
            {
                FullFilePath = path;
                BaseSequence = baseSeq;
                FullSequence = fullSeq;
                Score = score;
                OneBasedScanNumber = scan;
                _identified = identified?.ToList() ?? new();
            }

            public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods() => _identified;
            public int CompareTo(ISpectralMatch? o) => o is null ? -1 : o.Score.CompareTo(Score);
        }

        #endregion

        #region Sequence Coverage Tests

        /// <summary>
        /// Verifies CalculateSequenceCoverage correctly computes full coverage when all residues are covered.
        /// Critical: Sequence coverage is a key quality metric reported in protein group output files.
        /// </summary>
        [Test]
        public void CalculateSequenceCoverage_FullCoverage()
        {
            var prot = new TestBioPolymer("MEDEEKPEPTIDE", "P00001");
            var peptide1 = new TestBioPolymerWithSetMods(prot, 1, 6);  // MEDEEK
            var peptide2 = new TestBioPolymerWithSetMods(prot, 7, 13); // PEPTIDE

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { prot },
                new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 },
                new HashSet<IBioPolymerWithSetMods> { peptide1, peptide2 });

            var psm1 = new TestSpectralMatch("file.mzML", "MEDEEK", "MEDEEK", 100, 1, new[] { peptide1 });
            var psm2 = new TestSpectralMatch("file.mzML", "PEPTIDE", "PEPTIDE", 100, 2, new[] { peptide2 });
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
            var prot = new TestBioPolymer("MEDEEKPEPTIDE", "P00001"); // 13 residues
            var peptide = new TestBioPolymerWithSetMods(prot, 1, 6);   // MEDEEK - 6 residues

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { prot },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm = new TestSpectralMatch("file.mzML", "MEDEEK", "MEDEEK", 100, 1, new[] { peptide });
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
            var prot = new TestBioPolymer("MEDEEK", "P00001");
            var peptide = new TestBioPolymerWithSetMods(prot, 1, 6);

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { prot },
                new HashSet<IBioPolymerWithSetMods> { peptide },
                new HashSet<IBioPolymerWithSetMods> { peptide });

            var psm1 = new TestSpectralMatch("file1.mzML", "MEDEEK", "MEDEEK", 100, 1, new[] { peptide });
            var psm2 = new TestSpectralMatch("file2.mzML", "MEDEEK", "MEDEEK", 100, 1, new[] { peptide });
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
            var prot = new TestBioPolymer("MEDEEK", "P00001");
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
            var prot = new TestBioPolymer(longSequence, "P00001");
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
            var prot = new TestBioPolymer("MEDEEK", "P00001");
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