using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Modifications;
using Omics.SpectralMatch;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Quantification;
using Test.Omics;

namespace Test.Quantification
{
    /// <summary>
    /// One peptide sequence found in several proteins, or twice in one protein (mzLib #1280). A search names one
    /// candidate per place the sequence occurs, and PeptideWithSetModifications equality includes the parent accession
    /// and the start residue, so those candidates are unequal objects. Before this, the engine compared them as objects:
    /// such a match counted as ambiguous and was dropped, every copy became its own peptide row, and a peptide counted
    /// toward a protein group only through the caller's UniqueBioPolymersWithSetMods -- which MetaMorpheus fills before
    /// parsimony, and so leaves empty for every group of indistinguishable proteins. Such a group got no value at all.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class SharedSequenceQuantificationTests
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

        private static readonly ISampleInfo[] Channels =
        {
            new IsobaricQuantSampleInfo(File, "Control", 0, 0, 0, 0, "126", 126.0, false),
            new IsobaricQuantSampleInfo(File, "Treatment", 0, 0, 0, 0, "127N", 127.1, false)
        };

        private static TestExperimentalDesign Design() =>
            new TestExperimentalDesign(new Dictionary<string, ISampleInfo[]> { [File] = Channels });

        /// <summary>The protein and every tryptic peptide of it, in sequence order.</summary>
        private static (Protein protein, List<IBioPolymerWithSetMods> peptides) Digest(string sequence, string accession)
        {
            var protein = new Protein(sequence, accession);
            var peptides = protein
                .Digest(new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 5),
                    new List<Modification>(), new List<Modification>())
                .Cast<IBioPolymerWithSetMods>()
                .OrderBy(p => p.OneBasedStartResidue)
                .ToList();
            return (protein, peptides);
        }

        private static IBioPolymerWithSetMods Of(List<IBioPolymerWithSetMods> peptides, string sequence) =>
            peptides.Single(p => p.FullSequence == sequence);

        private static MockSpectralMatch Match(int scan, double[] intensities, params IBioPolymerWithSetMods[] identified) =>
            new MockSpectralMatch(File, identified[0].FullSequence, identified[0].BaseSequence, 100.0, scan, identified)
            {
                Intensities = intensities
            };

        private static QuantMatrix<ISpectralMatch> MatrixOf(params ISpectralMatch[] matches) =>
            new QuantMatrix<ISpectralMatch>(matches.ToList(), Channels.ToList(), Design());

        private static BioPolymerGroup Group(IEnumerable<IBioPolymer> proteins,
            IEnumerable<IBioPolymerWithSetMods> all, IEnumerable<IBioPolymerWithSetMods> unique) =>
            new BioPolymerGroup(new HashSet<IBioPolymer>(proteins),
                new HashSet<IBioPolymerWithSetMods>(all), new HashSet<IBioPolymerWithSetMods>(unique));

        /// <summary>Intensities by channel label; empty when there are none.</summary>
        private static Dictionary<string, double> ByChannel(IDictionary<ISampleInfo, double> intensities) =>
            (intensities ?? new Dictionary<ISampleInfo, double>()).ToDictionary(
                kvp => ((IsobaricQuantSampleInfo)kvp.Key).ChannelLabel, kvp => kvp.Value);

        private static Dictionary<string, double> ByChannel(IBioPolymerGroup group) => ByChannel(group.IntensitiesBySample);

        private static Dictionary<string, double> PeptideByChannel(QuantificationResults result, string sequence) =>
            ByChannel(result.PeptideIntensities.Single(kvp => kvp.Key.FullSequence == sequence).Value);

        private static QuantificationResults Run(List<ISpectralMatch> matches, List<IBioPolymerWithSetMods> peptides,
            List<IBioPolymerGroup> groups, bool useSharedPeptides)
        {
            var parameters = QuantificationParameters.GetSimpleParameters();
            parameters.UseSharedPeptidesForProteinQuant = useSharedPeptides;
            return new QuantificationEngine(parameters, Design(), matches, peptides, groups).Run();
        }

        [Test]
        public void OneSequenceInTwoProteins_IsOneIdentification()
        {
            var inP1 = Of(Digest("PEPTIDEK", "P1").peptides, "PEPTIDEK");
            var inP2 = Of(Digest("PEPTIDEK", "P2").peptides, "PEPTIDEK");
            Assert.That(inP1, Is.Not.EqualTo(inP2), "the fixture only means anything if the two copies are unequal");

            var map = QuantificationEngine.GetPsmToPeptideMap(
                MatrixOf(Match(1, new[] { 1000.0, 2000.0 }, inP1, inP2)),
                new List<IBioPolymerWithSetMods> { inP1, inP2 });

            Assert.Multiple(() =>
            {
                // Previously two keys, both empty: the match was dropped as ambiguous.
                Assert.That(map.Keys.Select(p => p.FullSequence), Is.EqualTo(new[] { "PEPTIDEK" }));
                Assert.That(map.Values.Single(), Is.EqualTo(new List<int> { 0 }));
            });
        }

        [Test]
        public void OneSequenceTwiceInOneProtein_IsOneIdentification()
        {
            var (_, peptides) = Digest("PEPTIDEKPEPTIDEK", "P1");
            Assert.That(peptides.Select(p => p.OneBasedStartResidue), Is.EqualTo(new[] { 1, 9 }),
                "the fixture only means anything if the protein yields the sequence at two places");

            var map = QuantificationEngine.GetPsmToPeptideMap(
                MatrixOf(Match(1, new[] { 1000.0, 2000.0 }, peptides[0], peptides[1])), peptides);

            Assert.Multiple(() =>
            {
                Assert.That(map.Keys.Select(p => p.FullSequence), Is.EqualTo(new[] { "PEPTIDEK" }));
                Assert.That(map.Values.Single(), Is.EqualTo(new List<int> { 0 }));
            });
        }

        [Test]
        public void TwoDifferentSequences_AreStillAmbiguous_EvenInOneProtein()
        {
            var (_, peptides) = Digest("PEPTIDEKAGLLVEDK", "P1");

            var map = QuantificationEngine.GetPsmToPeptideMap(
                MatrixOf(Match(1, new[] { 1000.0, 2000.0 }, peptides[0], peptides[1])), peptides);

            // Two sequences for one spectrum is a real ambiguity, whichever proteins they come from.
            Assert.That(map.Values.SelectMany(v => v), Is.Empty);
        }

        [Test]
        public void TheRowIsTheSameCopy_WhicheverOrderTheCallerListsThem()
        {
            var inP1 = Of(Digest("PEPTIDEK", "P1").peptides, "PEPTIDEK");
            var inP2 = Of(Digest("PEPTIDEK", "P2").peptides, "PEPTIDEK");
            var matrix = MatrixOf(Match(1, new[] { 1000.0, 2000.0 }, inP1, inP2));

            var forward = QuantificationEngine.GetPsmToPeptideMap(matrix, new List<IBioPolymerWithSetMods> { inP1, inP2 });
            var backward = QuantificationEngine.GetPsmToPeptideMap(matrix, new List<IBioPolymerWithSetMods> { inP2, inP1 });

            // The lowest accession, so the key a caller finds in PeptideIntensities does not depend on list order.
            Assert.Multiple(() =>
            {
                Assert.That(forward.Keys.Single(), Is.SameAs(inP1));
                Assert.That(backward.Keys.Single(), Is.SameAs(inP1));
            });
        }

        /// <summary>
        /// Calmodulin's case: two proteins with the same sequence, one group, and an empty unique list, which is what
        /// MetaMorpheus passes for indistinguishable proteins. Every value differs, so each total is reachable only
        /// one way: match 3 has already been narrowed to a single copy, and with shared peptides on, a sequence
        /// counted once per copy would double it.
        /// </summary>
        [TestCase(false)]
        [TestCase(true)]
        public void IndistinguishableProteins_TheirGroupIsQuantified(bool useSharedPeptides)
        {
            var (protein1, copies1) = Digest("PEPTIDEKAGLLVEDK", "P1");
            var (protein2, copies2) = Digest("PEPTIDEKAGLLVEDK", "P2");
            var peptides = copies1.Concat(copies2).ToList();
            var group = Group(new[] { protein1, protein2 }, peptides, new IBioPolymerWithSetMods[0]);

            var matches = new List<ISpectralMatch>
            {
                Match(1, new[] { 1000.0, 2000.0 }, Of(copies1, "PEPTIDEK"), Of(copies2, "PEPTIDEK")),
                Match(2, new[] { 30.0, 60.0 }, Of(copies1, "AGLLVEDK"), Of(copies2, "AGLLVEDK")),
                Match(3, new[] { 5.0, 10.0 }, Of(copies2, "PEPTIDEK")),
            };

            var result = Run(matches, peptides, new List<IBioPolymerGroup> { group }, useSharedPeptides);

            Assert.Multiple(() =>
            {
                Assert.That(result.Success, Is.True, result.Summary);
                Assert.That(result.AmbiguousSpectralMatchesExcluded, Is.Zero, "previously 2: matches 1 and 2");
                Assert.That(result.Summary, Is.EqualTo("Quantification completed successfully."));

                // One row per sequence, not one per copy.
                Assert.That(result.PeptideIntensities.Keys.Select(p => p.FullSequence),
                    Is.EquivalentTo(new[] { "AGLLVEDK", "PEPTIDEK" }));
                Assert.That(PeptideByChannel(result, "PEPTIDEK"),
                    Is.EqualTo(new Dictionary<string, double> { ["126"] = 1005.0, ["127N"] = 2010.0 }));

                // Previously nothing with shared peptides off, and 5 and 10 with them on.
                Assert.That(ByChannel(group), Is.EqualTo(new Dictionary<string, double> { ["126"] = 1035.0, ["127N"] = 2070.0 }));
            });
        }

        /// <summary>
        /// The contrast: PEPTIDEK sits in two proteins that parsimony left in two different groups, so it is shared
        /// between groups. It is still one peptide, quantified once; whether a group counts it is the shared-peptide
        /// setting's call, as for any shared peptide.
        /// </summary>
        [TestCase(false, 30.0, 7.0)]
        [TestCase(true, 1030.0, 1007.0)]
        public void ASequenceInTwoGroups_IsSharedBetweenThem(bool useSharedPeptides, double group1, double group2)
        {
            var (protein1, copies1) = Digest("PEPTIDEKAGLLVEDK", "P1");
            var (protein2, copies2) = Digest("PEPTIDEKQLFEGASK", "P2");
            var peptides = copies1.Concat(copies2).ToList();
            var g1 = Group(new[] { protein1 }, copies1, new[] { Of(copies1, "AGLLVEDK") });
            var g2 = Group(new[] { protein2 }, copies2, new[] { Of(copies2, "QLFEGASK") });

            var matches = new List<ISpectralMatch>
            {
                Match(1, new[] { 1000.0, 2000.0 }, Of(copies1, "PEPTIDEK"), Of(copies2, "PEPTIDEK")),
                Match(2, new[] { 30.0, 60.0 }, Of(copies1, "AGLLVEDK")),
                Match(3, new[] { 7.0, 14.0 }, Of(copies2, "QLFEGASK")),
            };

            var result = Run(matches, peptides, new List<IBioPolymerGroup> { g1, g2 }, useSharedPeptides);

            Assert.Multiple(() =>
            {
                Assert.That(result.AmbiguousSpectralMatchesExcluded, Is.Zero, "previously 1: match 1");
                Assert.That(PeptideByChannel(result, "PEPTIDEK"),
                    Is.EqualTo(new Dictionary<string, double> { ["126"] = 1000.0, ["127N"] = 2000.0 }));
                Assert.That(ByChannel(g1)["126"], Is.EqualTo(group1));
                Assert.That(ByChannel(g1)["127N"], Is.EqualTo(2 * group1));
                Assert.That(ByChannel(g2)["126"], Is.EqualTo(group2));
                Assert.That(ByChannel(g2)["127N"], Is.EqualTo(2 * group2));
            });
        }

        /// <summary>
        /// A peptide that no other group lists counts as unique to its group, whatever the caller's unique list says.
        /// MetaMorpheus decides that list before parsimony, so it leaves out a peptide whose other protein parsimony
        /// later discarded; the engine judges among the groups it is given instead.
        /// </summary>
        [Test]
        public void UniquenessIsJudgedAmongTheGroups_NotTakenFromTheCallersList()
        {
            var (protein1, copies1) = Digest("PEPTIDEKAGLLVEDK", "P1");
            var group = Group(new[] { protein1 }, copies1, new[] { Of(copies1, "AGLLVEDK") });

            var matches = new List<ISpectralMatch>
            {
                Match(1, new[] { 1000.0, 2000.0 }, Of(copies1, "PEPTIDEK")),
                Match(2, new[] { 30.0, 60.0 }, Of(copies1, "AGLLVEDK")),
            };

            var result = Run(matches, copies1, new List<IBioPolymerGroup> { group }, useSharedPeptides: false);

            // Previously 30 and 60: PEPTIDEK was left out because the caller's unique list did not name it.
            Assert.That(ByChannel(group), Is.EqualTo(new Dictionary<string, double> { ["126"] = 1030.0, ["127N"] = 2060.0 }),
                result.Summary);
        }
    }
}
