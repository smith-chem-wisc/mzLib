using System.Linq;
using Chemistry;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using Omics.Modifications;
using Omics.Modifications.IO.Modomics;

namespace Test.Transcriptomics
{
    [TestFixture]
    public class ModomicsLoaderTests
    {
        [Test]
        public void LoadModomics_ParsesEntriesCorrectly()
        {
            var report = ModomicsLoader.LoadModomics();

            Assert.That(report.LoadedModifications, Is.Not.Empty);
            Assert.That(report.LoadedCount, Is.EqualTo(report.LoadedModifications.Count));
            Assert.That(report.LoadedModifications.All(m => !string.IsNullOrWhiteSpace(m.IdWithMotif)), Is.True);
            Assert.That(report.LoadedModifications.All(m => m.IdWithMotif.Contains(" on ")), Is.True);
        }

        [Test]
        public void LoadModomics_SkipsUnknownBases()
        {
            var report = ModomicsLoader.LoadModomics();

            Assert.That(report.LoadedModifications.Any(m => m.IdWithMotif.Contains("on X")), Is.False);
        }

        [Test]
        public void LoadModomics_CachesResults()
        {
            var firstLoad = ModomicsLoader.LoadModomics();
            var secondLoad = ModomicsLoader.LoadModomics();
            Assert.That(ReferenceEquals(firstLoad, secondLoad), Is.True);
        }


        [Test]
        public void LoadModomics_TracksNotYetRepresentableEntriesWithReasons()
        {
            var report = ModomicsLoader.LoadModomics();

            Assert.That(report.NotYetRepresentableEntries, Is.Not.Empty);

            // Every tracked entry remains identifiable (some source rows have no short name).
            Assert.That(report.NotYetRepresentableEntries.All(e => !string.IsNullOrWhiteSpace(e.Name) || !string.IsNullOrWhiteSpace(e.ShortName)), Is.True);

            // Unknown stand-in entries (e.g. "xG", "xX") carry no formula and must be tracked, not silently dropped.
            Assert.That(report.NotYetRepresentableEntries.Any(e => e.Reason == ModomicsRepresentationFailureReason.EmptyFormula), Is.True);

            // Generic "X" moieties (e.g. synthetic nucleotides such as M1Y) need per-nucleotide expansion.
            Assert.That(report.NotYetRepresentableEntries.Any(e => e.Reason == ModomicsRepresentationFailureReason.UnsupportedReferenceMoiety), Is.True);

            // Modified bases (e.g. preQ0base/preQ1base) require new residue definitions.
            Assert.That(report.NotYetRepresentableEntries.Any(e => e.Reason == ModomicsRepresentationFailureReason.BaseMoietyRequiresNewResidue), Is.True);

            // Real-base nucleotide formulas are anionic and protonation-inconsistent (e.g. pm1G).
            Assert.That(report.NotYetRepresentableEntries.Any(e => e.Reason == ModomicsRepresentationFailureReason.NucleotideProtonationAmbiguous), Is.True);

            // Terminus descriptors and generic-N cofactors (e.g. "5' diphosphate end", NAD) need the terminus model.
            Assert.That(report.NotYetRepresentableEntries.Any(e => e.Reason == ModomicsRepresentationFailureReason.TerminalStructureRequiresTerminusModel), Is.True);

            // Canonical nucleoside entries (e.g. "G on G") and pseudouridine impart no mass shift.
            Assert.That(report.NotYetRepresentableEntries.Any(e => e.Reason == ModomicsRepresentationFailureReason.NoMassShift), Is.True);
        }

        [Test]
        public void NucleotideEntriesAreTrackedRatherThanLoadedWithWrongMasses()
        {
            var report = ModomicsLoader.LoadModomics();

            // pm1G ("1N-methylguanosine-5'-monophosphate") is the anionic counterpart of m1G; it must not
            // load with the sugar-removal shift (~215 Da) that a naive nucleotide transform produces.
            Assert.That(report.LoadedModifications.Any(m => m.IdWithMotif == "pm1G on G"), Is.False);
            Assert.That(report.NotYetRepresentableEntries.Any(e => e.ShortName == "pm1G"
                && e.Reason == ModomicsRepresentationFailureReason.NucleotideProtonationAmbiguous), Is.True);

            // The neutral nucleoside twin carries the same chemistry.
            var m1G = report.LoadedModifications.SingleOrDefault(m => m.IdWithMotif == "m1G on G");
            Assert.That(m1G, Is.Not.Null);
            Assert.That(m1G!.ChemicalFormula.Equals(ChemicalFormula.ParseFormula("CH2")), Is.True);
        }

        [Test]
        public void LoadModomics_TerminalCapsUseExistingTerminalRestriction()
        {
            var report = ModomicsLoader.LoadModomics();

            Assert.That(report.TerminalModifications, Is.Not.Empty);
            Assert.That(report.TerminalModifications.All(m => m.LocationRestriction == "5'-terminal."), Is.True);
            Assert.That(report.TerminalModifications.All(m => m.ModificationType == "5' Terminal Cap"), Is.True);

            // Cap 0 (m7GpppN) is emitted once per reference moiety, mirroring the multi-target JSON data.
            var cap0Targets = report.TerminalModifications
                .Where(m => m.OriginalId == "m7GpppN")
                .Select(m => m.Target.ToString())
                .OrderBy(t => t)
                .ToList();
            CollectionAssert.AreEqual(new[] { "A", "C", "G", "U" }, cap0Targets);

            // The cap shift keeps the cap nucleoside and phosphate chain: m7GpppN (C16H23N5O17P3) minus the
            // generic ribose (C5H7O3).
            var cap0OnA = report.TerminalModifications.Single(m => m.IdWithMotif == "m7GpppN on A");
            Assert.That(cap0OnA.ChemicalFormula.Equals(ChemicalFormula.ParseFormula("C11H16N5O14P3")), Is.True);

            // The terminal restriction is enforced by the existing localization semantics: the cap fits the
            // first residue of the polymer, not an interior residue with a matching motif.
            Assert.That(ModificationLocalization.ModFits(cap0OnA, "AAAA", 1, 4, 1), Is.True);
            Assert.That(ModificationLocalization.ModFits(cap0OnA, "AAAA", 1, 4, 3), Is.False);
        }

        [Test]
        public void ModomicsModsAreRegisteredAsASeparateCategoryAlongsideCuratedMods()
        {
            // Same chemistry, two names: MODOMICS "Am" vs curated "2'-O-Methyladenosine" (both C1H2 on A).
            Assert.That(Mods.ModomicsRnaModifications.Any(m => m.IdWithMotif == "Am on A"), Is.True);
            Assert.That(Mods.MetaMorpheusRnaModifications.Any(m => m.IdWithMotif == "2'-O-Methyladenosine on A"), Is.True);
            Assert.That(Mods.AllRnaModsList.Any(m => m.IdWithMotif == "Am on A"), Is.True);
            Assert.That(Mods.AllRnaModsList.Any(m => m.IdWithMotif == "2'-O-Methyladenosine on A"), Is.True);
            Assert.That(Mods.AllKnownRnaModsDictionary.ContainsKey("Am on A"), Is.True);
            Assert.That(Mods.AllKnownRnaModsDictionary.ContainsKey("2'-O-Methyladenosine on A"), Is.True);
            Assert.That(Mods.ModsByConvention[ModificationNamingConvention.Modomics], Is.Not.Empty);

            // The chemistry-equivalent overlap is reported without removing either entry.
            var modomicsAm = Mods.ModomicsLoadReport.LoadedModifications.Single(m => m.IdWithMotif == "Am on A");
            var duplicate = Mods.ModomicsLoadReport.DuplicateModifications.SingleOrDefault(d => d.ModomicsModification == modomicsAm);
            Assert.That(duplicate, Is.Not.Null);
            Assert.That(duplicate!.ExistingModifications.Any(m => m.IdWithMotif == "2'-O-Methyladenosine on A"), Is.True);
        }

        [Test]
        public void AmbiguousRnaLookupDefaultsToMetaMorpheus()
        {
            // For any id present in both conventions, the combined RNA dictionary must resolve to the
            // curated (MetaMorpheus) entry, mirroring how protein conventions coexist.
            var collisions = Mods.ModomicsRnaModifications
                .Where(m => Mods.MetaMorpheusRnaModifications.Any(c => c.IdWithMotif == m.IdWithMotif))
                .ToList();

            foreach (var modomicsMod in collisions)
            {
                var resolved = Mods.GetModification(modomicsMod.IdWithMotif, false, true);
                Assert.That(resolved, Is.Not.Null, modomicsMod.IdWithMotif);
                Assert.That(resolved!.ModificationType, Is.Not.EqualTo("Modomics"), modomicsMod.IdWithMotif);
                Assert.That(resolved.ModificationType, Is.Not.EqualTo("5' Terminal Cap"), modomicsMod.IdWithMotif);
            }
        }

        [Test]
        public void ConventionLookupRoutesEachNamingStyleToItsOwnMods()
        {
            var modomicsMod = Mods.GetModification("m6A on A", ModificationNamingConvention.Modomics);
            Assert.That(modomicsMod, Is.Not.Null);
            Assert.That(modomicsMod!.ModificationType, Is.EqualTo("Modomics"));

            var curatedMod = Mods.GetModification("N6-methyladenosine on A", ModificationNamingConvention.MetaMorpheus_Rna);
            Assert.That(curatedMod, Is.Not.Null);
            Assert.That(curatedMod!.ModificationType, Is.Not.EqualTo("Modomics"));

            // "Closest to the input": each naming style resolves through the combined RNA dictionary too.
            var byModomicsName = Mods.GetModification("m6A on A", false, true);
            Assert.That(byModomicsName, Is.Not.Null);
            Assert.That(byModomicsName!.ModificationType, Is.EqualTo("Modomics"));

            var byCuratedName = Mods.GetModification("N6-methyladenosine on A", false, true);
            Assert.That(byCuratedName, Is.Not.Null);
            Assert.That(byCuratedName.ModificationType, Is.Not.EqualTo("Modomics"));
        }

        [Test]
        public void MethylsHaveCorrectFormula()
        {
            var mods = ModomicsLoader.LoadModomics().LoadedModifications
                .Where(p => !p.ModificationType.Contains("Cap"))
                .ToList();

            var singleMethylMods = mods.Where(m =>
                    System.Text.RegularExpressions.Regex.IsMatch(
                        m.IdWithMotif,
                        @"^m\d+[ACGU] on [ACGU]$"))
                    .ToList();

            // Exception to the normal CH2 rule: the 3-methylcytidine cation formula carries an additional proton.
            var m3C = singleMethylMods.FirstOrDefault(m => m.IdWithMotif.StartsWith("m3C on C"));
            Assert.That(m3C, Is.Not.Null, "Expected to find m3C modification");
            Assert.That(m3C!.ChemicalFormula.Equals(ChemicalFormula.ParseFormula("C1H3")), Is.True, "m3C should have formula C1H3");
            singleMethylMods.Remove(m3C);

            var expectedFormula = ChemicalFormula.ParseFormula("C1H2");
            CollectionAssert.AreEqual(
                Enumerable.Repeat(expectedFormula, singleMethylMods.Count),
                singleMethylMods.Select(p => p.ChemicalFormula)
            );

            var diMethylMods = mods.Where(m =>
                    System.Text.RegularExpressions.Regex.IsMatch(
                        m.IdWithMotif,
                        @"^m\d+[ACGU]m on [ACGU]$"))
                    .ToList();

            expectedFormula = ChemicalFormula.ParseFormula("C2H4");
            CollectionAssert.AreEqual(
                Enumerable.Repeat(expectedFormula, diMethylMods.Count).ToList(),
                diMethylMods.Select(p => p.ChemicalFormula).ToList()
            );
        }
    }
}
