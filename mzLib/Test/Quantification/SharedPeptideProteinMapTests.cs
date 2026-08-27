using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Digestion;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Quantification;
using Omics.SpectralMatch;
using Test.Omics;

namespace Test.Quantification
{
    /// <summary>
    /// Tests for QuantificationEngine.GetAllPeptideToProteinMap -- the shared-peptide half of
    /// protein roll-up, which threw NotImplementedException until now, so
    /// UseSharedPeptidesForProteinQuant = true could not run at all.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class SharedPeptideProteinMapTests
    {
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

        private static BioPolymerGroup Group(
            string accession,
            IEnumerable<IBioPolymerWithSetMods> all,
            IEnumerable<IBioPolymerWithSetMods> unique)
        {
            return new BioPolymerGroup(
                new HashSet<IBioPolymer> { new Protein("MPEPTIDEK", accession) },
                new HashSet<IBioPolymerWithSetMods>(all),
                new HashSet<IBioPolymerWithSetMods>(unique));
        }

        /// <summary>The group's delivered intensities, keyed by channel label rather than sample object.</summary>
        private static Dictionary<string, double> ByChannel(IBioPolymerGroup group)
        {
            return group.IntensitiesBySample.ToDictionary(
                kvp => ((IsobaricQuantSampleInfo)kvp.Key).ChannelLabel,
                kvp => kvp.Value);
        }

        private static QuantMatrix<IBioPolymerWithSetMods> MatrixOf(params IBioPolymerWithSetMods[] peptides)
        {
            const string file = "file1.raw";
            var columns = new List<ISampleInfo>
            {
                new IsobaricQuantSampleInfo(file, "Control", 0, 0, 0, 0, "126", 126.0, false)
            };
            var design = new TestExperimentalDesign(
                new Dictionary<string, ISampleInfo[]> { [file] = columns.ToArray() });

            return new QuantMatrix<IBioPolymerWithSetMods>(peptides.ToList(), columns, design);
        }

        [Test]
        public void SharedPeptide_CountsTowardEveryGroupItBelongsTo()
        {
            var pepA = Peptide("PEPTIDEK", "P1");
            var pepShared = Peptide("SEQUENCEK", "P2");
            var pepC = Peptide("ANOTHERK", "P3");

            // pepShared is in both groups; pepA is unique to g1 and pepC unique to g2.
            var g1 = Group("G1", new[] { pepA, pepShared }, new[] { pepA });
            var g2 = Group("G2", new[] { pepShared, pepC }, new[] { pepC });
            var groups = new List<IBioPolymerGroup> { g1, g2 };

            var matrix = MatrixOf(pepA, pepShared, pepC);

            var map = QuantificationEngine.GetAllPeptideToProteinMap(matrix, groups);

            Assert.Multiple(() =>
            {
                // Row 1 (the shared peptide) appears in both lists. That is the whole difference
                // from the unique map, where each row index appears at most once.
                Assert.That(map[g1], Is.EqualTo(new List<int> { 0, 1 }));
                Assert.That(map[g2], Is.EqualTo(new List<int> { 1, 2 }));
            });
        }

        [Test]
        public void UniqueMap_ExcludesTheSharedPeptide()
        {
            var pepA = Peptide("PEPTIDEK", "P1");
            var pepShared = Peptide("SEQUENCEK", "P2");
            var pepC = Peptide("ANOTHERK", "P3");

            var g1 = Group("G1", new[] { pepA, pepShared }, new[] { pepA });
            var g2 = Group("G2", new[] { pepShared, pepC }, new[] { pepC });
            var groups = new List<IBioPolymerGroup> { g1, g2 };

            var matrix = MatrixOf(pepA, pepShared, pepC);

            var unique = QuantificationEngine.GetUniquePeptideToProteinMap(matrix, groups);
            var all = QuantificationEngine.GetAllPeptideToProteinMap(matrix, groups);

            Assert.Multiple(() =>
            {
                // Contrast, on the same input: the shared row is absent from the unique map and
                // present in both lists of the all map.
                Assert.That(unique[g1], Is.EqualTo(new List<int> { 0 }));
                Assert.That(unique[g2], Is.EqualTo(new List<int> { 2 }));
                Assert.That(all[g1].Count + all[g2].Count, Is.EqualTo(4));
                Assert.That(unique[g1].Count + unique[g2].Count, Is.EqualTo(2));
            });
        }

        [Test]
        public void EveryGroupGetsAnEntry_EvenWithNoPeptidesInTheMatrix()
        {
            var pepA = Peptide("PEPTIDEK", "P1");
            var absent = Peptide("SEQUENCEK", "P2");

            var g1 = Group("G1", new[] { pepA }, new[] { pepA });
            var g2 = Group("G2", new[] { absent }, new[] { absent });
            var groups = new List<IBioPolymerGroup> { g1, g2 };

            // Only pepA is in the matrix -- g2's peptide was never passed to the engine.
            var matrix = MatrixOf(pepA);

            var map = QuantificationEngine.GetAllPeptideToProteinMap(matrix, groups);

            Assert.Multiple(() =>
            {
                Assert.That(map.Keys, Has.Count.EqualTo(2));
                Assert.That(map[g1], Is.EqualTo(new List<int> { 0 }));
                // An empty list, not a missing key -- the roll-up strategy iterates the map.
                Assert.That(map[g2], Is.Empty);
            });
        }

        [Test]
        public void RowIndicesAreSorted()
        {
            var pep1 = Peptide("PEPTIDEK", "P1");
            var pep2 = Peptide("SEQUENCEK", "P2");
            var pep3 = Peptide("ANOTHERK", "P3");
            var pep4 = Peptide("MYFAVORITEK", "P4");

            var group = Group("G1", new[] { pep4, pep1, pep3, pep2 }, new[] { pep1 });
            var groups = new List<IBioPolymerGroup> { group };

            var matrix = MatrixOf(pep1, pep2, pep3, pep4);

            var map = QuantificationEngine.GetAllPeptideToProteinMap(matrix, groups);

            // AllBioPolymersWithSetMods is a HashSet, whose enumeration order is not guaranteed
            // stable. Sorting keeps roll-up from depending on it.
            Assert.That(map[group], Is.EqualTo(new List<int> { 0, 1, 2, 3 }));
        }

        /// <summary>
        /// Builds the same two-group scenario each time -- fresh groups, because the engine writes its
        /// results onto them -- and returns each group's per-channel totals.
        /// </summary>
        private static (Dictionary<string, double> g1, Dictionary<string, double> g2, QuantificationResults result)
            RunSharedPeptideScenario(bool useSharedPeptides)
        {
            const string file = "file1.raw";
            var columns = new ISampleInfo[]
            {
                new IsobaricQuantSampleInfo(file, "Control", 0, 0, 0, 0, "126", 126.0, false),
                new IsobaricQuantSampleInfo(file, "Treatment", 0, 0, 0, 0, "127N", 127.1, false)
            };
            var design = new TestExperimentalDesign(
                new Dictionary<string, ISampleInfo[]> { [file] = columns });

            var pepA = Peptide("PEPTIDEK", "P1");
            var pepShared = Peptide("SEQUENCEK", "P2");
            var pepC = Peptide("ANOTHERK", "P3");
            var peptides = new List<IBioPolymerWithSetMods> { pepA, pepShared, pepC };

            var g1 = Group("G1", new[] { pepA, pepShared }, new[] { pepA });
            var g2 = Group("G2", new[] { pepShared, pepC }, new[] { pepC });
            var groups = new List<IBioPolymerGroup> { g1, g2 };

            // Distinct intensities per peptide, so the expected protein totals are only reachable if
            // the shared peptide really is counted in both groups. With equal values the totals match
            // even when the shared peptide is dropped from both.
            var intensityByPeptide = new Dictionary<IBioPolymerWithSetMods, double[]>
            {
                [pepA] = new[] { 100.0, 200.0 },
                [pepShared] = new[] { 30.0, 60.0 },
                [pepC] = new[] { 7.0, 14.0 },
            };

            int scan = 1;
            var matches = peptides
                .Select(p => (ISpectralMatch)new MockSpectralMatch(
                    file, p.FullSequence, p.BaseSequence, 100.0, scan++, new[] { p })
                {
                    Intensities = intensityByPeptide[p]
                })
                .ToList();

            var parameters = QuantificationParameters.GetSimpleParameters();
            parameters.UseSharedPeptidesForProteinQuant = useSharedPeptides;

            var engine = new QuantificationEngine(parameters, design, matches, peptides, groups);

            // Before this change, Run() threw NotImplementedException here with the flag on -- it was
            // settable and unusable.
            var result = engine.Run();

            return (ByChannel(g1), ByChannel(g2), result);
        }

        [Test]
        public void Engine_RunsWithSharedPeptidesEnabled()
        {
            var (g1ByChannel, g2ByChannel, result) = RunSharedPeptideScenario(useSharedPeptides: true);

            Assert.Multiple(() =>
            {
                Assert.That(result.Success, Is.True, result.Summary);

                // G1 = pepA + pepShared = 100+30 and 200+60. G2 = pepShared + pepC = 30+7 and 60+14.
                // Neither is reachable without the shared peptide -- the companion test below pins
                // what the same input produces when it is excluded.
                Assert.That(g1ByChannel["126"], Is.EqualTo(130.0));
                Assert.That(g1ByChannel["127N"], Is.EqualTo(260.0));
                Assert.That(g2ByChannel["126"], Is.EqualTo(37.0));
                Assert.That(g2ByChannel["127N"], Is.EqualTo(74.0));
            });
        }

        [Test]
        public void Engine_WithSharedPeptidesOff_LeavesTheSharedPeptideOut()
        {
            var (g1ByChannel, g2ByChannel, result) = RunSharedPeptideScenario(useSharedPeptides: false);

            Assert.Multiple(() =>
            {
                Assert.That(result.Success, Is.True, result.Summary);

                // Same input, unique peptides only: G1 is pepA alone and G2 is pepC alone. Every value
                // differs from the shared-peptide run, so the pair of tests distinguishes the two
                // behaviours rather than merely agreeing with each other.
                Assert.That(g1ByChannel["126"], Is.EqualTo(100.0));
                Assert.That(g1ByChannel["127N"], Is.EqualTo(200.0));
                Assert.That(g2ByChannel["126"], Is.EqualTo(7.0));
                Assert.That(g2ByChannel["127N"], Is.EqualTo(14.0));
            });
        }

        [Test]
        public void NoGroups_GivesAnEmptyMap()
        {
            var matrix = MatrixOf(Peptide("PEPTIDEK", "P1"));

            var map = QuantificationEngine.GetAllPeptideToProteinMap(matrix, new List<IBioPolymerGroup>());

            Assert.That(map, Is.Empty);
        }

        [Test]
        public void EmptyMatrix_GivesEveryGroupAnEmptyList()
        {
            var pepA = Peptide("PEPTIDEK", "P1");
            var group = Group("G1", new[] { pepA }, new[] { pepA });

            var matrix = MatrixOf();

            var map = QuantificationEngine.GetAllPeptideToProteinMap(
                matrix, new List<IBioPolymerGroup> { group });

            Assert.That(map[group], Is.Empty);
        }
    }
}
