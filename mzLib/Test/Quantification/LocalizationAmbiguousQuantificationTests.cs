using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics;
using Omics.BioPolymerGroup;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SpectralMatch;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Quantification;
using Test.Omics;

namespace Test.Quantification
{
    /// <summary>
    /// A spectral match ambiguous only in where its modifications sit: two or more full sequences of one base sequence,
    /// carrying the same modifications at different sites, explain the scan equally well. Such a match is quantified as
    /// one peptide row standing for all of its sites (<see cref="LocalizationAmbiguousPeptide"/>), and counts toward a
    /// protein group when that group is the only one listing its sites. Before, every match without one full sequence
    /// was left out, and a protein got nothing from the spectra that could not be localized.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class LocalizationAmbiguousQuantificationTests
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

        private static readonly DigestionParams Digestion = new(minPeptideLength: 1);

        private static Modification Mod(string id, string residue, string formula)
        {
            ModificationMotif.TryGetMotif(residue, out var motif);
            return new Modification(_originalId: id, _modificationType: "Test", _target: motif,
                _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula(formula));
        }

        private static readonly Modification Oxidation = Mod("Oxidation", "M", "O");
        private static readonly Modification Dioxidation = Mod("Dioxidation", "M", "O2");

        /// <summary>A peptide of <paramref name="protein"/>, with <paramref name="mods"/> keyed by 1-based residue in the peptide.</summary>
        private static IBioPolymerWithSetMods Peptide(Protein protein, int start, int end, params (int Residue, Modification Mod)[] mods) =>
            new PeptideWithSetModifications(protein, Digestion, start, end, CleavageSpecificity.Full, "", 0,
                mods.ToDictionary(m => m.Residue + 1, m => m.Mod), 0);

        private static MockSpectralMatch Match(int scan, double[] intensities, params IBioPolymerWithSetMods[] identified) =>
            new MockSpectralMatch(File, string.Join("|", identified.Select(p => p.FullSequence).Distinct()),
                identified[0].BaseSequence, 100.0, scan, identified)
            {
                Intensities = intensities
            };

        private static BioPolymerGroup Group(IEnumerable<IBioPolymer> proteins, IEnumerable<IBioPolymerWithSetMods> peptides) =>
            new BioPolymerGroup(new HashSet<IBioPolymer>(proteins), new HashSet<IBioPolymerWithSetMods>(peptides),
                new HashSet<IBioPolymerWithSetMods>());

        private static Dictionary<string, double> ByChannel(IDictionary<ISampleInfo, double> intensities) =>
            (intensities ?? new Dictionary<ISampleInfo, double>()).ToDictionary(
                kvp => ((IsobaricQuantSampleInfo)kvp.Key).ChannelLabel, kvp => kvp.Value);

        private static QuantificationResults Run(List<ISpectralMatch> matches, List<IBioPolymerWithSetMods> peptides,
            List<IBioPolymerGroup> groups, bool useSharedPeptides = false)
        {
            var parameters = QuantificationParameters.GetSimpleParameters();
            parameters.UseSharedPeptidesForProteinQuant = useSharedPeptides;
            var design = new TestExperimentalDesign(new Dictionary<string, ISampleInfo[]> { [File] = Channels });
            return new QuantificationEngine(parameters, design, matches, peptides, groups).Run();
        }

        /// <summary>
        /// The caller's peptide list: the engine refuses an empty one, and a localization-ambiguous row is built from the
        /// matches, not from this list, so these tests pass one unrelated peptide.
        /// </summary>
        private static List<IBioPolymerWithSetMods> Unrelated() =>
            new() { Peptide(new Protein("WWWWWK", "UNRELATED"), 1, 6) };

        private static string Name(params IBioPolymerWithSetMods[] sites) =>
            string.Join("|", sites.Select(s => s.FullSequence).OrderBy(s => s, StringComparer.Ordinal));

        /// <summary>
        /// The user's case: one protein, oxidation on M3 or on M5. The match is one peptide row named by both sites,
        /// and P1's group, the only group listing them, gets its reporter intensities.
        /// </summary>
        [Test]
        public void AMatchAmbiguousOnlyInSite_IsOneRowAndCountsForItsGroup()
        {
            var protein = new Protein("ACMDMEKGGGR", "P1");
            var at3 = Peptide(protein, 1, 7, (3, Oxidation));
            var at5 = Peptide(protein, 1, 7, (5, Oxidation));
            var group = Group(new[] { protein }, new[] { at3, at5 });

            var result = Run(new List<ISpectralMatch> { Match(1, new[] { 1000.0, 2000.0 }, at3, at5) },
                Unrelated(), new List<IBioPolymerGroup> { group });

            Assert.Multiple(() =>
            {
                Assert.That(result.AmbiguousSpectralMatchesExcluded, Is.Zero);
                var row = result.PeptideIntensities.Keys.OfType<LocalizationAmbiguousPeptide>().Single();
                Assert.That(row, Is.InstanceOf<LocalizationAmbiguousPeptide>());
                Assert.That(row.FullSequence, Is.EqualTo(Name(at3, at5)));
                Assert.That(row.BaseSequence, Is.EqualTo("ACMDMEK"));
                Assert.That(ByChannel(result.ProteinIntensities[group]), Is.EqualTo(new Dictionary<string, double> { ["126"] = 1000, ["127N"] = 2000 }));
            });
        }

        /// <summary>
        /// One scan is ambiguous between M1 and M3, another between M3 and M5: one base sequence and one oxidation, so one
        /// row, named by all three sites, holding both scans.
        /// </summary>
        [Test]
        public void MatchesAmbiguousOverDifferentSitesOfOnePeptide_AreOneRowNamedByEverySite()
        {
            var protein = new Protein("MDMEMKGGGR", "P1");
            var at1 = Peptide(protein, 1, 6, (1, Oxidation));
            var at3 = Peptide(protein, 1, 6, (3, Oxidation));
            var at5 = Peptide(protein, 1, 6, (5, Oxidation));

            var result = Run(new List<ISpectralMatch>
                {
                    Match(1, new[] { 1000.0, 2000.0 }, at1, at3),
                    Match(2, new[] { 10.0, 20.0 }, at5, at3),
                },
                Unrelated(), new List<IBioPolymerGroup> { Group(new[] { protein }, new[] { at1, at3, at5 }) });

            var row = result.PeptideIntensities.Single(kvp => kvp.Key is LocalizationAmbiguousPeptide);
            Assert.That(row.Key.FullSequence, Is.EqualTo(Name(at1, at3, at5)));
            Assert.That(ByChannel(row.Value), Is.EqualTo(new Dictionary<string, double> { ["126"] = 1010, ["127N"] = 2020 }));
        }

        /// <summary>
        /// A localized match keeps its own row; the ambiguous one beside it is a second row and does not add to the first.
        /// Each scan is its own reporter-ion measurement, so the two are not counted twice.
        /// </summary>
        [Test]
        public void ALocalizedMatchKeepsItsRowBesideTheAmbiguousOne()
        {
            var protein = new Protein("ACMDMEKGGGR", "P1");
            var at3 = Peptide(protein, 1, 7, (3, Oxidation));
            var at5 = Peptide(protein, 1, 7, (5, Oxidation));
            var group = Group(new[] { protein }, new[] { at3, at5 });

            var result = Run(new List<ISpectralMatch>
                {
                    Match(1, new[] { 1000.0, 2000.0 }, at3),
                    Match(2, new[] { 10.0, 20.0 }, at3, at5),
                },
                new List<IBioPolymerWithSetMods> { at3 }, new List<IBioPolymerGroup> { group });

            var rows = result.PeptideIntensities.ToDictionary(kvp => kvp.Key.FullSequence, kvp => ByChannel(kvp.Value));
            Assert.That(rows.Keys, Has.None.EqualTo(at5.FullSequence), "a site seen only in an ambiguous match gets no row of its own");
            Assert.That(rows.Keys, Is.EquivalentTo(new[] { at3.FullSequence, Name(at3, at5) }));
            Assert.That(rows[at3.FullSequence], Is.EqualTo(new Dictionary<string, double> { ["126"] = 1000, ["127N"] = 2000 }));
            Assert.That(rows[Name(at3, at5)], Is.EqualTo(new Dictionary<string, double> { ["126"] = 10, ["127N"] = 20 }));
            Assert.That(ByChannel(result.ProteinIntensities[group]), Is.EqualTo(new Dictionary<string, double> { ["126"] = 1010, ["127N"] = 2020 }));
        }

        /// <summary>
        /// ACMDMEK is in P1 and in P2, which are two groups. The ambiguous match could be either group's; crediting it
        /// would guess the group on top of the site, so it counts for neither -- with shared peptides off, and on.
        /// </summary>
        [TestCase(false)]
        [TestCase(true)]
        public void AMatchAmbiguousInSiteAcrossTwoGroups_CountsForNeither(bool useSharedPeptides)
        {
            var p1 = new Protein("ACMDMEKGGGR", "P1");
            var p2 = new Protein("ACMDMEKWWWR", "P2");
            var sites1 = new[] { Peptide(p1, 1, 7, (3, Oxidation)), Peptide(p1, 1, 7, (5, Oxidation)) };
            var sites2 = new[] { Peptide(p2, 1, 7, (3, Oxidation)), Peptide(p2, 1, 7, (5, Oxidation)) };
            var g1 = Group(new[] { p1 }, sites1);
            var g2 = Group(new[] { p2 }, sites2);

            var result = Run(new List<ISpectralMatch> { Match(1, new[] { 1000.0, 2000.0 }, sites1.Concat(sites2).ToArray()) },
                Unrelated(), new List<IBioPolymerGroup> { g1, g2 }, useSharedPeptides);

            Assert.That(result.PeptideIntensities.Keys.OfType<LocalizationAmbiguousPeptide>().Single(), Is.InstanceOf<LocalizationAmbiguousPeptide>());
            Assert.That(ByChannel(result.ProteinIntensities[g1]).Values, Has.None.GreaterThan(0));
            Assert.That(ByChannel(result.ProteinIntensities[g2]).Values, Has.None.GreaterThan(0));
        }

        /// <summary>
        /// The sites split across groups: the database annotates the oxidation on M3 only in P1 and on M5 only in P2, so
        /// one site belongs to each group. The row may be either group's, so it counts for neither, whichever site
        /// happens to sort first.
        /// </summary>
        [Test]
        public void AMatchWhoseSitesBelongToDifferentGroups_CountsForNeither()
        {
            var p1 = new Protein("ACMDMEKGGGR", "P1");
            var p2 = new Protein("ACMDMEKWWWR", "P2");
            var at3InP1 = Peptide(p1, 1, 7, (3, Oxidation));
            var at5InP2 = Peptide(p2, 1, 7, (5, Oxidation));
            var g1 = Group(new[] { p1 }, new[] { at3InP1 });
            var g2 = Group(new[] { p2 }, new[] { at5InP2 });

            var result = Run(new List<ISpectralMatch> { Match(1, new[] { 1000.0, 2000.0 }, at3InP1, at5InP2) },
                Unrelated(), new List<IBioPolymerGroup> { g1, g2 });

            Assert.That(result.PeptideIntensities.Keys.OfType<LocalizationAmbiguousPeptide>().Single().FullSequence, Is.EqualTo(Name(at3InP1, at5InP2)));
            Assert.That(ByChannel(result.ProteinIntensities[g1]).Values, Has.None.GreaterThan(0));
            Assert.That(ByChannel(result.ProteinIntensities[g2]).Values, Has.None.GreaterThan(0));
        }

        /// <summary>
        /// Several proteins, one group: P1 and P2 carry the same sequence and form the group P1|P2, so a match ambiguous in
        /// site and in protein still has every candidate in that one group, and counts for it.
        /// </summary>
        [Test]
        public void AMatchAmbiguousInSiteAcrossProteinsOfOneGroup_CountsForThatGroup()
        {
            var p1 = new Protein("ACMDMEKGGGR", "P1");
            var p2 = new Protein("ACMDMEKGGGR", "P2");
            var sites = new[]
            {
                Peptide(p1, 1, 7, (3, Oxidation)), Peptide(p1, 1, 7, (5, Oxidation)),
                Peptide(p2, 1, 7, (3, Oxidation)), Peptide(p2, 1, 7, (5, Oxidation)),
            };
            var group = Group(new[] { p1, p2 }, sites);

            var result = Run(new List<ISpectralMatch> { Match(1, new[] { 1000.0, 2000.0 }, sites) },
                Unrelated(), new List<IBioPolymerGroup> { group });

            Assert.That(result.PeptideIntensities.Keys.OfType<LocalizationAmbiguousPeptide>().Single().FullSequence.Split('|'), Has.Length.EqualTo(2), "two sites; two proteins add no names");
            Assert.That(ByChannel(result.ProteinIntensities[group]), Is.EqualTo(new Dictionary<string, double> { ["126"] = 1000, ["127N"] = 2000 }));
        }

        /// <summary>
        /// Two oxidations or one dioxidation: one total formula, but different modifications, so not a localization
        /// ambiguity. Still excluded, and counted.
        /// </summary>
        [Test]
        public void AMatchAmbiguousBetweenDifferentModifications_IsStillExcluded()
        {
            var protein = new Protein("ACMDMEKGGGR", "P1");
            var twoOxidations = Peptide(protein, 1, 7, (3, Oxidation), (5, Oxidation));
            var oneDioxidation = Peptide(protein, 1, 7, (3, Dioxidation));
            Assume.That(twoOxidations.MonoisotopicMass, Is.EqualTo(oneDioxidation.MonoisotopicMass).Within(1e-9));

            var result = Run(new List<ISpectralMatch> { Match(1, new[] { 1000.0, 2000.0 }, twoOxidations, oneDioxidation) },
                Unrelated(), new List<IBioPolymerGroup> { Group(new[] { protein }, new[] { twoOxidations, oneDioxidation }) });

            Assert.That(result.AmbiguousSpectralMatchesExcluded, Is.EqualTo(1));
            Assert.That(result.PeptideIntensities.Keys.OfType<LocalizationAmbiguousPeptide>(), Is.Empty);
        }

        /// <summary>
        /// A modification is what it is by its type, its name without the residue and its formula, not by its formula
        /// alone. Oxidation on M3 against Hydroxylation on D4 is one oxygen either way, but two modifications: the match
        /// cannot say which one it carries, so it is excluded as before. Phosphorylation on S or on T is one modification
        /// on two residues, so that match is one row.
        /// </summary>
        [Test]
        public void ModificationsAreComparedByIdentity_NotByFormulaAlone()
        {
            var protein = new Protein("ACMDMEKGGGR", "P1");
            var oxidationAtM3 = Peptide(protein, 1, 7, (3, Oxidation));
            var hydroxylationAtD4 = Peptide(protein, 1, 7, (4, Mod("Hydroxylation", "D", "O")));
            Assume.That(oxidationAtM3.MonoisotopicMass, Is.EqualTo(hydroxylationAtD4.MonoisotopicMass).Within(1e-9));

            var phosphoProtein = new Protein("ASDTEKGGGR", "P2");
            var phosphoAtS2 = Peptide(phosphoProtein, 1, 6, (2, Mod("Phosphorylation", "S", "HO3P")));
            var phosphoAtT4 = Peptide(phosphoProtein, 1, 6, (4, Mod("Phosphorylation", "T", "HO3P")));

            var result = Run(new List<ISpectralMatch>
                {
                    Match(1, new[] { 1000.0, 2000.0 }, oxidationAtM3, hydroxylationAtD4),
                    Match(2, new[] { 300.0, 400.0 }, phosphoAtS2, phosphoAtT4)
                },
                Unrelated(), new List<IBioPolymerGroup>
                {
                    Group(new[] { protein }, new[] { oxidationAtM3, hydroxylationAtD4 }),
                    Group(new[] { phosphoProtein }, new[] { phosphoAtS2, phosphoAtT4 })
                });

            Assert.Multiple(() =>
            {
                Assert.That(result.AmbiguousSpectralMatchesExcluded, Is.EqualTo(1), "only the oxidation/hydroxylation match");
                Assert.That(result.PeptideIntensities.Keys.OfType<LocalizationAmbiguousPeptide>().Select(p => p.FullSequence),
                    Is.EqualTo(new[] { Name(phosphoAtS2, phosphoAtT4) }));
                Assert.That(LocalizationAmbiguousPeptide.SiteFreeIdentity(oxidationAtM3),
                    Is.Not.EqualTo(LocalizationAmbiguousPeptide.SiteFreeIdentity(hydroxylationAtD4)));
                Assert.Throws<ArgumentException>(() => new LocalizationAmbiguousPeptide(new[] { oxidationAtM3, hydroxylationAtD4 }));
            });
        }

        [Test]
        public void ALocalizationAmbiguousPeptide_IsItsSitesAndClaimsNone()
        {
            var protein = new Protein("ACMDMEKGGGR", "P1");
            var at3 = Peptide(protein, 1, 7, (3, Oxidation));
            var at5 = Peptide(protein, 1, 7, (5, Oxidation));
            var peptide = new LocalizationAmbiguousPeptide(new[] { at5, at3 });

            Assert.Multiple(() =>
            {
                Assert.That(peptide.Sites, Is.EqualTo(new[] { at3, at5 }.OrderBy(s => s.FullSequence, StringComparer.Ordinal)));
                Assert.That(peptide.FullSequence, Is.EqualTo(Name(at3, at5)));
                Assert.That(peptide.MonoisotopicMass, Is.EqualTo(at3.MonoisotopicMass));
                Assert.That(peptide.AllModsOneIsNterminus, Is.Empty, "no site every candidate agrees on");
                Assert.That(peptide, Is.EqualTo(new LocalizationAmbiguousPeptide(new[] { at3, at5 })));
                Assert.That(peptide.GetHashCode(), Is.EqualTo(new LocalizationAmbiguousPeptide(new[] { at3, at5 }).GetHashCode()));
                Assert.That(peptide.Equals(at3), Is.False);
                Assert.Throws<NotSupportedException>(() => peptide.Fragment(DissociationType.HCD, FragmentationTerminus.Both, new List<Product>()));
                Assert.Throws<ArgumentException>(() => new LocalizationAmbiguousPeptide(new[] { at3 }));
                Assert.Throws<ArgumentException>(() => new LocalizationAmbiguousPeptide(new[] { at3, Peptide(protein, 1, 7) }));
            });
        }
    }
}
