using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Digestion;
using Omics.Modifications;
using Omics.SpectralMatch;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Quantification;
using Test.Omics;

namespace Test.Quantification
{
    /// <summary>
    /// Tests the unambiguous filter in GetPsmToPeptideMap. Before this, the map took
    /// GetIdentifiedBioPolymersWithSetMods().First() with the comment "Assumes unambiguous mapping",
    /// so an ambiguous match was attributed to whichever biopolymer came first and a match with no
    /// identification threw InvalidOperationException.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class AmbiguousMatchFilterTests
    {
        private const string File = "file1.raw";

        private class TestExperimentalDesign : IExperimentalDesign
        {
            public Dictionary<string, ISampleInfo[]> FileNameSampleInfoDictionary { get; }

            public TestExperimentalDesign(Dictionary<string, ISampleInfo[]> dict)
            {
                FileNameSampleInfoDictionary = dict;
            }
        }

        private static IBioPolymerWithSetMods Peptide(string sequence, string accession)
        {
            var protein = new Protein(sequence, accession);
            var digestionParams = new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5);
            return protein.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();
        }

        private static List<ISampleInfo> Columns() => new List<ISampleInfo>
        {
            new IsobaricQuantSampleInfo(File, "Control", 0, 0, 0, 0, "126", 126.0, false),
            new IsobaricQuantSampleInfo(File, "Treatment", 0, 0, 0, 0, "127N", 127.1, false)
        };

        private static TestExperimentalDesign Design(List<ISampleInfo> columns) =>
            new TestExperimentalDesign(new Dictionary<string, ISampleInfo[]> { [File] = columns.ToArray() });

        private static MockSpectralMatch Match(int scan, params IBioPolymerWithSetMods[] identified) =>
            new MockSpectralMatch(File, $"SEQ{scan}", $"SEQ{scan}", 100.0, scan, identified)
            {
                Intensities = new[] { 1000.0, 2000.0 }
            };

        private static QuantMatrix<ISpectralMatch> MatrixOf(params ISpectralMatch[] matches)
        {
            var columns = Columns();
            return new QuantMatrix<ISpectralMatch>(matches.ToList(), columns, Design(columns));
        }

        [Test]
        public void AmbiguousMatch_IsExcluded_NotAttributedToTheFirstCandidate()
        {
            var pepA = Peptide("PEPTIDEK", "P1");
            var pepB = Peptide("SEQUENCEK", "P2");

            var unambiguous = Match(1, pepA);
            var ambiguous = Match(2, pepA, pepB); // could be either -- so it is neither

            var matrix = MatrixOf(unambiguous, ambiguous);

            var map = QuantificationEngine.GetPsmToPeptideMap(
                matrix, new List<IBioPolymerWithSetMods> { pepA, pepB });

            Assert.Multiple(() =>
            {
                // Row 1 is the ambiguous match. Previously it landed on pepA, because pepA is first.
                Assert.That(map[pepA], Is.EqualTo(new List<int> { 0 }));
                Assert.That(map[pepB], Is.Empty);
            });
        }

        [Test]
        public void MatchWithNoIdentification_IsExcludedRatherThanThrowing()
        {
            var pepA = Peptide("PEPTIDEK", "P1");

            var matrix = MatrixOf(Match(1, pepA), Match(2));

            // .First() on an empty sequence threw InvalidOperationException here.
            Dictionary<IBioPolymerWithSetMods, List<int>> map = null;
            Assert.DoesNotThrow(() => map = QuantificationEngine.GetPsmToPeptideMap(
                matrix, new List<IBioPolymerWithSetMods> { pepA }));

            Assert.That(map[pepA], Is.EqualTo(new List<int> { 0 }));
        }

        [Test]
        public void SamePeptideNamedTwice_IsNotAmbiguous()
        {
            var pepA = Peptide("PEPTIDEK", "P1");

            // A match can name one peptide once per protein it maps to. That identifies a single
            // peptide in substance, so counting raw entries would throw away a good measurement.
            var matrix = MatrixOf(Match(1, pepA, pepA));

            var map = QuantificationEngine.GetPsmToPeptideMap(
                matrix, new List<IBioPolymerWithSetMods> { pepA });

            Assert.That(map[pepA], Is.EqualTo(new List<int> { 0 }));
        }

        [Test]
        public void NullAlongsideOneIdentification_IsNotAmbiguous()
        {
            var pepA = Peptide("PEPTIDEK", "P1");

            var matrix = MatrixOf(Match(1, pepA, null));

            var map = QuantificationEngine.GetPsmToPeptideMap(
                matrix, new List<IBioPolymerWithSetMods> { pepA });

            Assert.That(map[pepA], Is.EqualTo(new List<int> { 0 }));
        }

        [Test]
        public void TheReportedCountMatchesWhatWasActuallyDropped()
        {
            var pepA = Peptide("PEPTIDEK", "P1");
            var pepB = Peptide("SEQUENCEK", "P2");
            var peptides = new List<IBioPolymerWithSetMods> { pepA, pepB };

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { new Protein("MPEPTIDEK", "G1") },
                new HashSet<IBioPolymerWithSetMods>(peptides),
                new HashSet<IBioPolymerWithSetMods>(peptides));

            var matches = new List<ISpectralMatch>
            {
                Match(1, pepA),
                Match(2, pepA, pepA),   // one peptide, named twice -- kept
                Match(3, pepA, null),   // one peptide plus a null -- kept
                Match(4, pepA, pepB),   // genuinely ambiguous -- dropped
                Match(5)                // no identification -- dropped
            };

            var columns = Columns();
            var engine = new QuantificationEngine(
                QuantificationParameters.GetSimpleParameters(), Design(columns), matches, peptides,
                new List<IBioPolymerGroup> { group });

            var result = engine.Run();

            // The count and the filter run off one predicate, so this pins them together.
            var matrix = MatrixOf(matches.ToArray());
            var kept = QuantificationEngine.GetPsmToPeptideMap(matrix, peptides).Values.Sum(v => v.Count);

            Assert.Multiple(() =>
            {
                Assert.That(result.AmbiguousSpectralMatchesExcluded, Is.EqualTo(2));
                Assert.That(kept, Is.EqualTo(matches.Count - result.AmbiguousSpectralMatchesExcluded));
            });
        }

        [Test]
        public void UnambiguousMatches_AreUnaffected()
        {
            var pepA = Peptide("PEPTIDEK", "P1");
            var pepB = Peptide("SEQUENCEK", "P2");

            var matrix = MatrixOf(Match(1, pepA), Match(2, pepB), Match(3, pepA));

            var map = QuantificationEngine.GetPsmToPeptideMap(
                matrix, new List<IBioPolymerWithSetMods> { pepA, pepB });

            Assert.Multiple(() =>
            {
                Assert.That(map[pepA], Is.EqualTo(new List<int> { 0, 2 }));
                Assert.That(map[pepB], Is.EqualTo(new List<int> { 1 }));
            });
        }

        [Test]
        public void Engine_ReportsHowManyMatchesItExcluded()
        {
            var pepA = Peptide("PEPTIDEK", "P1");
            var pepB = Peptide("SEQUENCEK", "P2");
            var peptides = new List<IBioPolymerWithSetMods> { pepA, pepB };

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { new Protein("MPEPTIDEK", "G1") },
                new HashSet<IBioPolymerWithSetMods>(peptides),
                new HashSet<IBioPolymerWithSetMods>(peptides));

            var matches = new List<ISpectralMatch>
            {
                Match(1, pepA),
                Match(2, pepB),
                Match(3, pepA, pepB), // ambiguous
                Match(4)              // no identification
            };

            var columns = Columns();
            var parameters = QuantificationParameters.GetSimpleParameters();
            var engine = new QuantificationEngine(
                parameters, Design(columns), matches, peptides,
                new List<IBioPolymerGroup> { group });

            var result = engine.Run();

            Assert.Multiple(() =>
            {
                Assert.That(result.Success, Is.True, result.Summary);
                // Dropping these is correct; dropping them silently is not.
                Assert.That(result.AmbiguousSpectralMatchesExcluded, Is.EqualTo(2));
                Assert.That(result.Summary, Does.Contain("2 spectral match(es) were excluded"));
            });
        }

        [Test]
        public void Engine_SaysNothingExtraWhenNothingWasExcluded()
        {
            var pepA = Peptide("PEPTIDEK", "P1");
            var peptides = new List<IBioPolymerWithSetMods> { pepA };

            var group = new BioPolymerGroup(
                new HashSet<IBioPolymer> { new Protein("MPEPTIDEK", "G1") },
                new HashSet<IBioPolymerWithSetMods>(peptides),
                new HashSet<IBioPolymerWithSetMods>(peptides));

            var columns = Columns();
            var engine = new QuantificationEngine(
                QuantificationParameters.GetSimpleParameters(),
                Design(columns),
                new List<ISpectralMatch> { Match(1, pepA) },
                peptides,
                new List<IBioPolymerGroup> { group });

            var result = engine.Run();

            Assert.Multiple(() =>
            {
                Assert.That(result.AmbiguousSpectralMatchesExcluded, Is.Zero);
                Assert.That(result.Summary, Is.EqualTo("Quantification completed successfully."));
            });
        }
    }
}
