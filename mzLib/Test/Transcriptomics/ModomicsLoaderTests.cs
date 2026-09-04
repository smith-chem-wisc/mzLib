using System.Linq;
using Chemistry;
using MassSpectrometry;
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
            Assert.That(report.LoadedModifications.Any(m => m.IdWithMotif == "1N-methylguanosine-5'-monophosphate on G"), Is.False);
            Assert.That(report.NotYetRepresentableEntries.Any(e => e.ShortName == "pm1G"
                && e.Reason == ModomicsRepresentationFailureReason.NucleotideProtonationAmbiguous), Is.True);

            // The neutral nucleoside twin carries the same chemistry.
            var m1G = report.LoadedModifications.SingleOrDefault(m => m.IdWithMotif == "1-methylguanosine on G");
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
                .Where(m => m.OriginalId == "N7-methyl-guanosine cap (cap 0)")
                .Select(m => m.Target.ToString())
                .OrderBy(t => t)
                .ToList();
            CollectionAssert.AreEqual(new[] { "A", "C", "G", "U" }, cap0Targets);

            // The cap shift keeps the cap nucleoside and phosphate chain: m7GpppN is published as an anion,
            // so it must be neutralized first and then have only the base-less ribose scaffold removed.
            var cap0OnA = report.TerminalModifications.Single(m => m.IdWithMotif == "N7-methyl-guanosine cap (cap 0) on A");
            Assert.That(cap0OnA.ChemicalFormula.Equals(ChemicalFormula.ParseFormula("C11H16N5O13P3")), Is.True);

            var methylTriphosphateCap = report.TerminalModifications.Single(m => m.IdWithMotif == "gamma-methyltriphosphate 5' cap on A");
            Assert.That(methylTriphosphateCap.ChemicalFormula.Equals(ChemicalFormula.ParseFormula("CH5O9P3")), Is.True);

            var methylMonophosphateCap = report.TerminalModifications.Single(m => m.IdWithMotif == "alpha-methylmonophosphate 5' cap on A");
            Assert.That(methylMonophosphateCap.ChemicalFormula.Equals(ChemicalFormula.ParseFormula("CH3O3P")), Is.True);

            var dimethylMonophosphateCap = report.TerminalModifications.Single(m => m.IdWithMotif == "alpha-dimethylmonophosphate 5' cap on A");
            Assert.That(dimethylMonophosphateCap.ChemicalFormula.Equals(ChemicalFormula.ParseFormula("C2H5O3P")), Is.True);

            // The terminal restriction is enforced by the existing localization semantics: the cap fits the
            // first residue of the polymer, not an interior residue with a matching motif.
            Assert.That(ModificationLocalization.ModFits(cap0OnA, "AAAA", 1, 4, 1), Is.True);
            Assert.That(ModificationLocalization.ModFits(cap0OnA, "AAAA", 1, 4, 3), Is.False);
        }

        [Test]
        public void ModomicsModsAreRegisteredAsASeparateCategoryAlongsideCuratedMods()
        {
            // Same chemistry, two names: MODOMICS "2'-O-methyladenosine" vs curated "2'-O-Methyladenosine"
            // (both C1H2 on A); the ids differ only in letter case, so both keys coexist.
            Assert.That(Mods.ModomicsRnaModifications.Any(m => m.IdWithMotif == "2'-O-methyladenosine on A"), Is.True);
            Assert.That(Mods.MetaMorpheusRnaModifications.Any(m => m.IdWithMotif == "2'-O-Methyladenosine on A"), Is.True);
            Assert.That(Mods.AllRnaModsList.Any(m => m.IdWithMotif == "2'-O-methyladenosine on A"), Is.True);
            Assert.That(Mods.AllRnaModsList.Any(m => m.IdWithMotif == "2'-O-Methyladenosine on A"), Is.True);
            Assert.That(Mods.AllKnownRnaModsDictionary.ContainsKey("2'-O-methyladenosine on A"), Is.True);
            Assert.That(Mods.AllKnownRnaModsDictionary.ContainsKey("2'-O-Methyladenosine on A"), Is.True);
            Assert.That(Mods.ModsByConvention[ModificationNamingConvention.Modomics], Is.Not.Empty);

            // The chemistry-equivalent overlap is reported without removing either entry.
            var modomicsAm = Mods.ModomicsLoadReport.LoadedModifications.Single(m => m.IdWithMotif == "2'-O-methyladenosine on A");
            var duplicate = Mods.ModomicsLoadReport.DuplicateModifications.SingleOrDefault(d => d.ModomicsModification == modomicsAm);
            Assert.That(duplicate, Is.Not.Null);
            Assert.That(duplicate!.ExistingModifications.Any(m => m.IdWithMotif == "2'-O-Methyladenosine on A"), Is.True);

            // 3-methylcytidine has no curated RnaMods.txt twin, so the Modomics registry itself is the second
            // surface that must preserve the neutral CH2 shift.
            var modomicsM3C = Mods.ModomicsRnaModifications.Single(m => m.IdWithMotif == "3-methylcytidine on C");
            Assert.That(modomicsM3C.ChemicalFormula.Equals(ChemicalFormula.ParseFormula("C1H2")), Is.True);
        }

        [Test]
        public void AmbiguousRnaLookupDefaultsToMetaMorpheus()
        {
            // MODOMICS "N6-methyladenosine on A" shares its id exactly with the curated mod, so at least
            // one real collision exists and the precedence rule is exercised, not vacuous.
            var collisions = Mods.ModomicsRnaModifications
                .Where(m => Mods.MetaMorpheusRnaModifications.Any(c => c.IdWithMotif == m.IdWithMotif))
                .ToList();
            Assert.That(collisions, Is.Not.Empty);

            // For any id present in both conventions, the combined RNA dictionary must resolve to the
            // curated (MetaMorpheus) entry, mirroring how protein conventions coexist.
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
            // "N6-methyladenosine on A" exists in both conventions; each convention bucket resolves it
            // to its own entry.
            var modomicsMod = Mods.GetModification("N6-methyladenosine on A", ModificationNamingConvention.Modomics);
            Assert.That(modomicsMod, Is.Not.Null);
            Assert.That(modomicsMod!.ModificationType, Is.EqualTo("Modomics"));

            var curatedMod = Mods.GetModification("N6-methyladenosine on A", ModificationNamingConvention.MetaMorpheus_Rna);
            Assert.That(curatedMod, Is.Not.Null);
            Assert.That(curatedMod!.ModificationType, Is.Not.EqualTo("Modomics"));

            // The shared id resolves to the curated (MetaMorpheus) entry in the combined RNA dictionary.
            var bySharedName = Mods.GetModification("N6-methyladenosine on A", false, true);
            Assert.That(bySharedName, Is.Not.Null);
            Assert.That(bySharedName!.ModificationType, Is.Not.EqualTo("Modomics"));

            // "Closest to the input": an unambiguous MODOMICS name resolves through the combined dictionary.
            var byModomicsName = Mods.GetModification("1-methyladenosine on A", false, true);
            Assert.That(byModomicsName, Is.Not.Null);
            Assert.That(byModomicsName!.ModificationType, Is.EqualTo("Modomics"));
        }

        [Test]
        public void ProductIonsBecomeAnyActivationTypeDiagnosticIons()
        {
            var report = ModomicsLoader.LoadModomics();

            // Diagnostic ions are attached only when the source lists them.
            Assert.That(report.LoadedModifications.Count(m => m.DiagnosticIons is not null), Is.GreaterThan(50));

            // The primary base ion is upgraded to the accurate neutral mass implied by the published nominal m/z:
            // 1-methyladenosine lists "150"; [1-methyladenine]+ = 150.078, so the stored neutral mass is 149.0707.
            var m1A = report.LoadedModifications.Single(m => m.IdWithMotif == "1-methyladenosine on A");
            Assert.That(m1A.DiagnosticIons, Is.Not.Null);
            Assert.That(m1A.DiagnosticIons![DissociationType.AnyActivationType], Has.Some.EqualTo(149.07).Within(0.001));

            // 5-methyl-2-thiouridine lists a single nominal ion, "143", stored as the corresponding neutral mass.
            var m5s2U = report.LoadedModifications.Single(m => m.IdWithMotif == "5-methyl-2-thiouridine on U");
            Assert.That(m5s2U.DiagnosticIons![DissociationType.AnyActivationType], Has.Some.EqualTo(142.021).Within(0.001));

            // Multi-ion entries upgrade only the primary ion; secondary fragments keep the published nominal
            // m/z values converted to neutral mass: 5-carboxymethyl-2-thiouridine lists "187/169/141".
            var cm5s2U = report.LoadedModifications.Single(m => m.IdWithMotif == "5-carboxymethyl-2-thiouridine on U");
            var cm5s2UIons = cm5s2U.DiagnosticIons![DissociationType.AnyActivationType];
            Assert.That(cm5s2UIons, Has.Some.EqualTo(186.01).Within(0.001));
            Assert.That(cm5s2UIons, Has.Some.EqualTo(167.993).Within(0.001));
            Assert.That(cm5s2UIons, Has.Some.EqualTo(139.993).Within(0.001));

            // The upgrade follows the measured topology rather than an all-on-base assumption:
            // 2'-O-methyladenosine's base ion is unmodified [adenine]+ (136 -> neutral 135.055), and
            // 3-methylcytidine's base carries only the CH2 of its neutral C1H2 shift (126 -> neutral 125.059).
            var am = report.LoadedModifications.Single(m => m.IdWithMotif == "2'-O-methyladenosine on A");
            Assert.That(am.DiagnosticIons![DissociationType.AnyActivationType], Has.Some.EqualTo(135.055).Within(0.001));
            Assert.That(am.DiagnosticIons[DissociationType.AnyActivationType], Has.None.EqualTo(149.07).Within(0.001));

            var m3C = report.LoadedModifications.Single(m => m.IdWithMotif == "3-methylcytidine on C");
            Assert.That(m3C.DiagnosticIons![DissociationType.AnyActivationType], Has.Some.EqualTo(125.059).Within(0.001));
        }

        [Test]
        public void BaseLossBehaviorIsDerivedFromProductIons()
        {
            var report = ModomicsLoader.LoadModomics();

            // Sugar-localized: 2'-O-methyladenosine fragments to unmodified [adenine]+ (136), so
            // nothing extra leaves with the base; the diagnostic ion is retained alongside the
            // derived behavior.
            var am = report.LoadedModifications.Single(m => m.IdWithMotif == "2'-O-methyladenosine on A");
            var amBaseMod = am as BaseModification;
            Assert.That(amBaseMod, Is.Not.Null);
            Assert.That(amBaseMod!.BaseLossType, Is.EqualTo(BaseLossBehavior.Default));
            Assert.That(amBaseMod.BaseLossModification, Is.Null);
            Assert.That(am.DiagnosticIons![DissociationType.AnyActivationType], Has.Some.EqualTo(135.055).Within(0.001));

            // Base-localized: N6-methyladenosine fragments to [N6-methyladenine]+ (150 = [adenine]+ + CH2),
            // so the modification departs with the base during base loss.
            var m6A = report.LoadedModifications.Single(m => m.IdWithMotif == "N6-methyladenosine on A");
            var m6ABaseMod = m6A as BaseModification;
            Assert.That(m6ABaseMod, Is.Not.Null);
            Assert.That(m6ABaseMod!.BaseLossType, Is.EqualTo(BaseLossBehavior.Modified));
            Assert.That(m6ABaseMod.BaseLossModification!.Equals(ChemicalFormula.ParseFormula("C1H2")), Is.True);
            Assert.That(m6A.DiagnosticIons![DissociationType.AnyActivationType], Has.Some.EqualTo(149.07).Within(0.001));

            // Base conversion: inosine targets adenine and fragments to [hypoxanthine]+
            // (137 = [adenine]+ + H-1N-1O1), so the full base difference leaves with the base.
            var inosine = report.LoadedModifications.Single(m => m.IdWithMotif == "inosine on A");
            var inosineBaseMod = inosine as BaseModification;
            Assert.That(inosineBaseMod, Is.Not.Null);
            Assert.That(inosineBaseMod!.BaseLossType, Is.EqualTo(BaseLossBehavior.Modified));
            Assert.That(inosineBaseMod.BaseLossModification!.Equals(ChemicalFormula.ParseFormula("H-1N-1O1")), Is.True);

            // Partially base-localized: N6,2'-O-dimethyladenosine's base ion is [N6-methyladenine]+
            // (150 = [adenine]+ + CH2), so only the N6 methyl departs with the base; C1H2 is the unique
            // sub-formula of the C2H4 shift matching that measurement.
            var m6Am = report.LoadedModifications.Single(m => m.IdWithMotif == "N6,2'-O-dimethyladenosine on A");
            var m6AmBaseMod = m6Am as BaseModification;
            Assert.That(m6AmBaseMod, Is.Not.Null);
            Assert.That(m6AmBaseMod!.BaseLossType, Is.EqualTo(BaseLossBehavior.Modified));
            Assert.That(m6AmBaseMod.BaseLossModification!.Equals(ChemicalFormula.ParseFormula("C1H2")), Is.True);
            Assert.That(m6Am.DiagnosticIons![DissociationType.AnyActivationType], Has.Some.EqualTo(149.07).Within(0.001));

            // Protonation-ambiguous cations resolve through measurement too: 3-methylcytidine's source
            // formula is cationic, but both its neutral shift and its base-localized portion are C1H2.
            var m3C = report.LoadedModifications.Single(m => m.IdWithMotif == "3-methylcytidine on C");
            var m3CBaseMod = m3C as BaseModification;
            Assert.That(m3CBaseMod, Is.Not.Null);
            Assert.That(m3CBaseMod!.BaseLossType, Is.EqualTo(BaseLossBehavior.Modified));
            Assert.That(m3CBaseMod.BaseLossModification!.Equals(ChemicalFormula.ParseFormula("C1H2")), Is.True);
            Assert.That(m3C.ChemicalFormula.Equals(ChemicalFormula.ParseFormula("C1H2")), Is.True);

            // Entries without published ions keep the plain representation.
            Assert.That(report.LoadedModifications.Any(m => m.DiagnosticIons is null), Is.True);
        }

        [Test]
        public void NucleosidePhosphateAnionsAreNeutralizedBeforeResidueShiftDerivation()
        {
            var report = ModomicsLoader.LoadModomics();

            var arP = report.LoadedModifications.Single(m => m.IdWithMotif == "2'-O-ribosyladenosine (phosphate) on A");
            Assert.That(arP.ChemicalFormula.Equals(ChemicalFormula.ParseFormula("C5H9O7P")), Is.True);

            var grP = report.LoadedModifications.Single(m => m.IdWithMotif == "2'-O-ribosylguanosine (phosphate) on G");
            Assert.That(grP.ChemicalFormula.Equals(ChemicalFormula.ParseFormula("C5H9O7P")), Is.True);
        }

        [Test]
        public void MethylsHaveCorrectFormula()
        {
            var mods = ModomicsLoader.LoadModomics().LoadedModifications
                .Where(p => !p.ModificationType.Contains("Cap"))
                .ToList();

            // Base methylations stated by full MODOMICS name (e.g. "1-methyladenosine", "N6-methyladenosine");
            // excludes ribose methyls ("2'-O-methyl..."), inosine/pseudouridine base conversions, and compound
            // names that merely contain "methyl" (methylthio, taurinomethyl, methyldihydrouridine, ...).
            var singleMethylMods = mods.Where(m =>
                    System.Text.RegularExpressions.Regex.IsMatch(
                        m.IdWithMotif,
                        @"^(N\d+-|\d+-)methyl(adenosine|cytidine|guanosine|uridine) on [ACGU]$"))
                    .ToList();
            Assert.That(singleMethylMods, Is.Not.Empty, "Expected to find base-methylated nucleosides");

            var expectedFormula = ChemicalFormula.ParseFormula("C1H2");
            CollectionAssert.AreEqual(
                Enumerable.Repeat(expectedFormula, singleMethylMods.Count),
                singleMethylMods.Select(p => p.ChemicalFormula)
            );

            // Base + 2'-O ribose dimethylations (e.g. "N6,2'-O-dimethyladenosine"), excluding
            // N,N-dimethyls and trimethyls.
            var diMethylMods = mods.Where(m =>
                    System.Text.RegularExpressions.Regex.IsMatch(
                        m.IdWithMotif,
                        @"^(\d+|N\d+),2'-O-dimethyl(adenosine|cytidine|guanosine|uridine) on [ACGU]$"))
                    .ToList();
            Assert.That(diMethylMods, Is.Not.Empty, "Expected to find base + ribose dimethylated nucleosides");

            expectedFormula = ChemicalFormula.ParseFormula("C2H4");
            CollectionAssert.AreEqual(
                Enumerable.Repeat(expectedFormula, diMethylMods.Count).ToList(),
                diMethylMods.Select(p => p.ChemicalFormula).ToList()
            );
        }
    }
}
