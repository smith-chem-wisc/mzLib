using System;
using System.Collections.Generic;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Omics.Modifications;
using Readers.ProForma;
using Tdp = TopDownProteomics.ProForma;

namespace Test.FileReadingTests.ProForma
{
    /// <summary>
    /// Layer-2: ProFormaTerm &lt;-&gt; mzLib (base sequence + AllModsOneIsNterminus). Covers
    /// per-residue name and accession (UNIMOD/MOD/RESID) modifications, motif-aware accession
    /// resolution, and the unsupported-feature guard. Terminal mods and mass/formula/glycan
    /// descriptors are later slices (see known-limitations.md).
    /// </summary>
    [TestFixture]
    internal class ProFormaConverterTests
    {
        private static Modification MakeMod(string name, char residue, double mass)
        {
            ModificationMotif.TryGetMotif(residue.ToString(), out var motif);
            return new Modification(_originalId: name, _modificationType: "testMods", _target: motif,
                _locationRestriction: "Anywhere.", _monoisotopicMass: mass);
        }

        private static Modification MakeMod(string name, char residue, double mass, string dbKey, string accession)
        {
            ModificationMotif.TryGetMotif(residue.ToString(), out var motif);
            return new Modification(_originalId: name, _modificationType: dbKey, _target: motif,
                _locationRestriction: "Anywhere.", _monoisotopicMass: mass,
                _databaseReference: new Dictionary<string, IList<string>> { [dbKey] = new List<string> { accession } });
        }

        // Terminal mods use the any-residue motif "X" and a terminal LocationRestriction, mirroring how
        // mzLib's Unimod loader stores N-/C-terminal entries.
        private static Modification MakeTerminalMod(string name, string locationRestriction, double mass,
            string? dbKey = null, string? accession = null)
        {
            ModificationMotif.TryGetMotif("X", out var motif);
            var dr = dbKey == null ? null : new Dictionary<string, IList<string>> { [dbKey] = new List<string> { accession! } };
            return new Modification(_originalId: name, _modificationType: dbKey ?? "testMods", _target: motif,
                _locationRestriction: locationRestriction, _monoisotopicMass: mass, _databaseReference: dr);
        }

        [Test]
        public void Layer2_NullArguments_Throw()
        {
            var allModsKnown = new Dictionary<string, Modification>();
            var term = ProFormaReader.Read("PEPTIDE");

            Assert.That(() => ProFormaConverter.ToModificationDictionary(null!, allModsKnown),
                Throws.TypeOf<ArgumentNullException>());
            Assert.That(() => ProFormaConverter.ToModificationDictionary(term, null!),
                Throws.TypeOf<ArgumentNullException>());
            Assert.That(() => ProFormaConverter.ToProFormaTerm(null!, new Dictionary<int, Modification>()),
                Throws.TypeOf<ArgumentNullException>());
            Assert.That(() => ProFormaConverter.ToProFormaTerm("PEPTIDE", null!),
                Throws.TypeOf<ArgumentNullException>());
        }

        [Test]
        public void Layer2_ResolvesAgainstTheWildcardMotif()
        {
            // Motif "X" means any residue, so it must match a concrete interior residue.
            ModificationMotif.TryGetMotif("X", out var anyResidue);
            var anywhere = new Modification(_originalId: "Label", _modificationType: "testMods", _target: anyResidue,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 8.01);
            var allModsKnown = new Dictionary<string, Modification> { [anywhere.IdWithMotif] = anywhere };

            var term = ProFormaReader.Read("EM[Label]EK");
            var dict = ProFormaConverter.ToModificationDictionary(term, allModsKnown);

            Assert.That(dict[3], Is.SameAs(anywhere));
        }

        [Test]
        public void Layer2_TagIndexOutsideSequence_Throws()
        {
            // A malformed term whose tag points past the end of its own sequence must be rejected,
            // not silently written to a key that means something else in AllModsOneIsNterminus.
            var ox = MakeMod("Oxidation", 'M', 15.99491);
            var allModsKnown = new Dictionary<string, Modification> { [ox.IdWithMotif] = ox };

            var descriptors = new[] { new Tdp.ProFormaDescriptor(Tdp.ProFormaKey.Name, "Oxidation") };
            var term = new Tdp.ProFormaTerm("PEP", new List<Tdp.ProFormaTag> { new(7, descriptors) },
                null, null, null, null, null, null);

            Assert.That(() => ProFormaConverter.ToModificationDictionary(term, allModsKnown),
                Throws.TypeOf<MzLibException>().With.Message.Contains("outside sequence bounds"));
        }

        [Test]
        public void Layer2_TwoTagsOnOneResidue_Throws()
        {
            // ProForma allows several modifications on one residue; mzLib stores one per position.
            // Fail loud rather than let the second silently overwrite the first.
            var ox = MakeMod("Oxidation", 'M', 15.99491);
            var allModsKnown = new Dictionary<string, Modification> { [ox.IdWithMotif] = ox };

            var descriptors = new[] { new Tdp.ProFormaDescriptor(Tdp.ProFormaKey.Name, "Oxidation") };
            var term = new Tdp.ProFormaTerm("EMEK",
                new List<Tdp.ProFormaTag> { new(1, descriptors), new(1, descriptors) },
                null, null, null, null, null, null);

            Assert.That(() => ProFormaConverter.ToModificationDictionary(term, allModsKnown),
                Throws.TypeOf<MzLibException>().With.Message.Contains("Multiple modifications target residue index 1"));
        }

        [Test]
        public void Layer2_ResolvesNameAgainstAContextBearingMotif()
        {
            // The N-glycosylation sequon motif is "Nxs": the modified residue is the upper-case one and
            // the surrounding context is lower-case, so the mod never keys as "Glycan on N". Resolution
            // must match on the motif's modified residue, not on the whole motif string.
            ModificationMotif.TryGetMotif("Nxs", out var sequon);
            var glyco = new Modification(_originalId: "Glycan", _modificationType: "testMods", _target: sequon,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 1234.43);
            var allModsKnown = new Dictionary<string, Modification> { [glyco.IdWithMotif] = glyco };

            var term = ProFormaReader.Read("EN[Glycan]ASK");
            var dict = ProFormaConverter.ToModificationDictionary(term, allModsKnown);

            Assert.That(dict.Keys, Is.EquivalentTo(new[] { 3 }));
            Assert.That(dict[3], Is.SameAs(glyco));
        }

        [Test]
        public void Layer2_ResolvesTerminalNameWhoseMotifIsNotX()
        {
            // Terminal mods usually key as "{name} on X". When a terminal mod is stored against a
            // concrete residue motif instead, fall back to matching the name plus terminus compatibility.
            var acetylOnK = MakeMod("Acetyl", 'K', 42.01057);
            var terminal = new Modification(_originalId: "Acetyl", _modificationType: "testMods",
                _target: acetylOnK.Target, _locationRestriction: "N-terminal.", _monoisotopicMass: 42.01057);
            var allModsKnown = new Dictionary<string, Modification> { [terminal.IdWithMotif] = terminal };

            var term = ProFormaReader.Read("[Acetyl]-KEEP");
            var dict = ProFormaConverter.ToModificationDictionary(term, allModsKnown);

            Assert.That(dict.Keys, Is.EquivalentTo(new[] { 1 }), "N-terminal mods live at key 1");
            Assert.That(dict[1], Is.SameAs(terminal));
        }

        [Test]
        public void Layer2_RoundTrips_PerResidueNameMods()
        {
            var ox = MakeMod("Oxidation", 'M', 15.99491);
            var ph = MakeMod("Phospho", 'S', 79.96633);
            var allModsKnown = new Dictionary<string, Modification> { [ox.IdWithMotif] = ox, [ph.IdWithMotif] = ph };

            // base sequence E M E V E E S P E K -> M is index 1 (key 3), S is index 6 (key 8)
            var term = ProFormaReader.Read("EM[Oxidation]EVEES[Phospho]PEK");
            var dict = ProFormaConverter.ToModificationDictionary(term, allModsKnown);

            Assert.That(dict.Keys, Is.EquivalentTo(new[] { 3, 8 }));
            Assert.That(dict[3], Is.SameAs(ox));
            Assert.That(dict[8], Is.SameAs(ph));

            // inverse: name mods (no accession) write back as names
            var rebuilt = ProFormaConverter.ToProFormaTerm(term.Sequence, dict);
            Assert.That(ProFormaWriter.Write(rebuilt), Is.EqualTo("EM[Oxidation]EVEES[Phospho]PEK"));
        }

        [Test]
        public void Layer2_RoundTrips_UnimodAccessions()
        {
            var ox = MakeMod("Oxidation", 'M', 15.99491, "Unimod", "35");
            var ph = MakeMod("Phospho", 'S', 79.96633, "Unimod", "21");
            var allModsKnown = new Dictionary<string, Modification> { [ox.IdWithMotif] = ox, [ph.IdWithMotif] = ph };

            var term = ProFormaReader.Read("EM[UNIMOD:35]EVEES[UNIMOD:21]PEK");
            var dict = ProFormaConverter.ToModificationDictionary(term, allModsKnown);
            Assert.That(dict[3], Is.SameAs(ox));
            Assert.That(dict[8], Is.SameAs(ph));

            // inverse: accession-bearing mods write back as accessions
            var rebuilt = ProFormaConverter.ToProFormaTerm(term.Sequence, dict);
            Assert.That(ProFormaWriter.Write(rebuilt), Is.EqualTo("EM[UNIMOD:35]EVEES[UNIMOD:21]PEK"));
        }

        [Test]
        public void Layer2_Accession_ResolvesByMotif_WhenAmbiguous()
        {
            // UNIMOD:21 (Phospho) maps to several mods differing only by motif; the residue disambiguates.
            var phS = MakeMod("Phospho", 'S', 79.96633, "Unimod", "21");
            var phT = MakeMod("Phospho", 'T', 79.96633, "Unimod", "21");
            var allModsKnown = new Dictionary<string, Modification> { [phS.IdWithMotif] = phS, [phT.IdWithMotif] = phT };

            // base ASTK -> S index 1 (key 3), T index 2 (key 4)
            var term = ProFormaReader.Read("AS[UNIMOD:21]T[UNIMOD:21]K");
            var dict = ProFormaConverter.ToModificationDictionary(term, allModsKnown);
            Assert.That(dict[3], Is.SameAs(phS));
            Assert.That(dict[4], Is.SameAs(phT));
        }

        [Test]
        public void Layer2_RoundTrips_PsiModAccession()
        {
            var mod = MakeMod("L-methionine sulfoxide", 'M', 15.99491, "PSI-MOD", "00719");
            var allModsKnown = new Dictionary<string, Modification> { [mod.IdWithMotif] = mod };

            var term = ProFormaReader.Read("EM[MOD:00719]EVEESPEK");
            var dict = ProFormaConverter.ToModificationDictionary(term, allModsKnown);
            Assert.That(dict[3], Is.SameAs(mod));

            var rebuilt = ProFormaConverter.ToProFormaTerm(term.Sequence, dict);
            Assert.That(ProFormaWriter.Write(rebuilt), Is.EqualTo("EM[MOD:00719]EVEESPEK"));
        }

        [Test]
        public void Layer2_RoundTrips_TerminalNameMods()
        {
            var nAcetyl = MakeTerminalMod("Acetyl", "N-terminal.", 42.01057);
            var cAmidation = MakeTerminalMod("Amidation", "C-terminal.", -0.98402);
            var ox = MakeMod("Oxidation", 'M', 15.99491);
            var allModsKnown = new Dictionary<string, Modification>
            {
                [nAcetyl.IdWithMotif] = nAcetyl, [cAmidation.IdWithMotif] = cAmidation, [ox.IdWithMotif] = ox
            };

            // base PEMTIDEK -> N=1, M index 2 (key 4), C=10
            var term = ProFormaReader.Read("[Acetyl]-PEM[Oxidation]TIDEK-[Amidation]");
            var dict = ProFormaConverter.ToModificationDictionary(term, allModsKnown);
            Assert.That(dict[1], Is.SameAs(nAcetyl));
            Assert.That(dict[4], Is.SameAs(ox));
            Assert.That(dict[10], Is.SameAs(cAmidation));

            var rebuilt = ProFormaConverter.ToProFormaTerm(term.Sequence, dict);
            Assert.That(ProFormaWriter.Write(rebuilt), Is.EqualTo("[Acetyl]-PEM[Oxidation]TIDEK-[Amidation]"));
        }

        [Test]
        public void Layer2_RoundTrips_TerminalAccession()
        {
            var nMod = MakeTerminalMod("iTRAQ4plex", "N-terminal.", 144.10253, "Unimod", "214");
            var allModsKnown = new Dictionary<string, Modification> { [nMod.IdWithMotif] = nMod };

            var term = ProFormaReader.Read("[UNIMOD:214]-PEPTIDEK");
            var dict = ProFormaConverter.ToModificationDictionary(term, allModsKnown);
            Assert.That(dict[1], Is.SameAs(nMod));

            var rebuilt = ProFormaConverter.ToProFormaTerm(term.Sequence, dict);
            Assert.That(ProFormaWriter.Write(rebuilt), Is.EqualTo("[UNIMOD:214]-PEPTIDEK"));
        }

        [Test]
        public void Layer2_Throws_WhenModNotInDatabase()
        {
            var term = ProFormaReader.Read("EM[Oxidation]EVEESPEK");
            Assert.Throws<MzLibException>(
                () => ProFormaConverter.ToModificationDictionary(term, new Dictionary<string, Modification>()));
        }

        [Test]
        public void Layer2_Throws_OnUnsupportedCrosslink()
        {
            // #XL1 produces a tag group, which is out of scope at Layer 2.
            var term = ProFormaReader.Read("EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]");
            var ex = Assert.Throws<MzLibException>(
                () => ProFormaConverter.ToModificationDictionary(term, new Dictionary<string, Modification>()));
            Assert.That(ex.Message, Does.Contain("tag group"));
        }

        [Test]
        public void Layer2_Throws_OnLabileModifications()
        {
            // {Glycan:Hex} is a labile modification — out of scope at Layer 2 (must fail loud, not drop silently).
            var term = ProFormaReader.Read("{Glycan:Hex}{Glycan:NeuAc}EMEVNESPEK");
            var ex = Assert.Throws<MzLibException>(
                () => ProFormaConverter.ToModificationDictionary(term, new Dictionary<string, Modification>()));
            Assert.That(ex.Message, Does.Contain("labile"));
        }

        [Test]
        public void Layer2_Throws_OnUnlocalizedModifications()
        {
            // [Phospho]? is an unlocalized (unknown-position) mod — out of scope at Layer 2.
            var term = ProFormaReader.Read("[Phospho]?EM[Oxidation]EVTSESPEK");
            var ex = Assert.Throws<MzLibException>(
                () => ProFormaConverter.ToModificationDictionary(term, new Dictionary<string, Modification>()));
            Assert.That(ex.Message, Does.Contain("unlocalized"));
        }

        [Test]
        public void Layer2_Throws_OnGlobalModifications()
        {
            // <13C> is a global isotope label — out of scope at Layer 2.
            var term = ProFormaReader.Read("<13C>ATPEILTVNSIGQLK");
            var ex = Assert.Throws<MzLibException>(
                () => ProFormaConverter.ToModificationDictionary(term, new Dictionary<string, Modification>()));
            Assert.That(ex.Message, Does.Contain("global"));
        }

        [Test]
        public void Layer2_Throws_OnPositionRange()
        {
            // (ESFRMS)[+19.0523] applies one mod across a residue range — out of scope at Layer 2.
            var term = ProFormaReader.Read("PRT(ESFRMS)[+19.0523]ISK");
            var ex = Assert.Throws<MzLibException>(
                () => ProFormaConverter.ToModificationDictionary(term, new Dictionary<string, Modification>()));
            Assert.That(ex.Message, Does.Contain("range"));
        }

        [Test]
        public void Layer2_Throws_OnSequenceAmbiguity()
        {
            // (?N) marks an ambiguous sequence stretch — out of scope at Layer 2.
            var term = ProFormaReader.Read("(?N)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K");
            var ex = Assert.Throws<MzLibException>(
                () => ProFormaConverter.ToModificationDictionary(term, new Dictionary<string, Modification>()));
            Assert.That(ex.Message, Does.Contain("ambiguity"));
        }

        [Test]
        public void Layer2_Throws_WhenTerminalModNotInDatabase()
        {
            // An N- or C-terminal descriptor that resolves to nothing must throw, not drop the mod silently.
            var nTerm = ProFormaReader.Read("[Acetyl]-PEPTIDEK");
            Assert.Throws<MzLibException>(
                () => ProFormaConverter.ToModificationDictionary(nTerm, new Dictionary<string, Modification>()));

            var cTerm = ProFormaReader.Read("PEPTIDEK-[Amidation]");
            Assert.Throws<MzLibException>(
                () => ProFormaConverter.ToModificationDictionary(cTerm, new Dictionary<string, Modification>()));
        }

        [Test]
        public void Layer2_ResolvesTerminalNameMod_WhenMotifIsResidueNotX()
        {
            // Terminal mod stored with a residue motif ("Acetyl on K", N-terminal) rather than the wildcard "X":
            // the "{name} on X" lookup misses and resolution falls back to a same-named, terminus-compatible mod.
            ModificationMotif.TryGetMotif("K", out var motifK);
            var nAcetyl = new Modification(_originalId: "Acetyl", _modificationType: "testMods", _target: motifK,
                _locationRestriction: "N-terminal.", _monoisotopicMass: 42.01057);
            var allModsKnown = new Dictionary<string, Modification> { [nAcetyl.IdWithMotif] = nAcetyl };

            var dict = ProFormaConverter.ToModificationDictionary(ProFormaReader.Read("[Acetyl]-KPEPTIDE"), allModsKnown);
            Assert.That(dict[1], Is.SameAs(nAcetyl));

            // On a peptide that does not start with K it does not fit, and nothing else resolves.
            Assert.That(() => ProFormaConverter.ToModificationDictionary(ProFormaReader.Read("[Acetyl]-PEPTIDEK"), allModsKnown),
                Throws.TypeOf<MzLibException>());
        }

        /// <summary>
        /// Several mods share a name or an accession and differ by residue. At a terminus the one whose
        /// motif fits the terminal residue wins, whatever the order they were loaded in: taking the first
        /// terminus-compatible one put N-terminal myristoylation of G back on C.
        /// </summary>
        [Test]
        public void Layer2_TerminalModIsChosenByTheTerminalResidue()
        {
            ModificationMotif.TryGetMotif("C", out var motifC);
            ModificationMotif.TryGetMotif("G", out var motifG);
            var dr = new Dictionary<string, IList<string>> { ["Unimod"] = new List<string> { "45" } };
            var onC = new Modification(_originalId: "Myristoyl", _modificationType: "Unimod", _target: motifC,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 210.198366, _databaseReference: dr);
            var onG = new Modification(_originalId: "Myristoyl", _modificationType: "Unimod", _target: motifG,
                _locationRestriction: "Peptide N-terminal.", _monoisotopicMass: 210.198366, _databaseReference: dr);
            var allModsKnown = new Dictionary<string, Modification> { [onC.IdWithMotif] = onC, [onG.IdWithMotif] = onG };

            var byAccession = ProFormaConverter.ToModificationDictionary(ProFormaReader.Read("[UNIMOD:45]-GPEPTIDE"), allModsKnown);
            Assert.That(byAccession[1], Is.SameAs(onG));

            var byName = ProFormaConverter.ToModificationDictionary(ProFormaReader.Read("[Myristoyl]-GPEPTIDE"), allModsKnown);
            Assert.That(byName[1], Is.SameAs(onG));

            var onCTerminus = ProFormaConverter.ToModificationDictionary(ProFormaReader.Read("[UNIMOD:45]-CPEPTIDE"), allModsKnown);
            Assert.That(onCTerminus[1], Is.SameAs(onC), "an Anywhere. mod still resolves at a terminus when it fits the residue");

            Assert.That(() => ProFormaConverter.ToModificationDictionary(ProFormaReader.Read("[UNIMOD:45]-APEPTIDE"), allModsKnown),
                Throws.TypeOf<MzLibException>(), "neither fits A");
        }

        /// <summary>
        /// Over the modifications mzLib actually loads: every terminus-restricted mod on a single residue
        /// that is written as an accession must read back as a mod on that same residue. Before the
        /// terminal residue was consulted, 41 of the 121 came back on another residue.
        /// </summary>
        [Test]
        public void Layer2_LoadedTerminalModsReadBackOnTheirOwnResidue()
        {
            var known = Mods.AllModsKnownDictionary;
            int informative = 0;
            var wrongResidue = new List<string>();
            foreach (var mod in known.Values)
            {
                char target = (mod.Target?.ToString() ?? "").FirstOrDefault(char.IsUpper);
                if (target == default || target == 'X')
                    continue;

                bool nTerm = mod.LocationRestriction is "N-terminal." or "Peptide N-terminal.";
                bool cTerm = mod.LocationRestriction is "C-terminal." or "Peptide C-terminal.";
                if (!nTerm && !cTerm)
                    continue;

                string sequence = nTerm ? target + "PEPTIDE" : "PEPTIDE" + target;
                var term = ProFormaConverter.ToProFormaTerm(sequence,
                    new Dictionary<int, Modification> { [nTerm ? 1 : sequence.Length + 2] = mod });
                var descriptors = nTerm ? term.NTerminalDescriptors : term.CTerminalDescriptors;
                if (descriptors[0].Key != Tdp.ProFormaKey.Identifier)
                    continue;

                informative++;
                var back = ProFormaConverter.ToModificationDictionary(term, known).Values.Single();
                char backTarget = (back.Target?.ToString() ?? "").FirstOrDefault(char.IsUpper);
                if (backTarget != target && backTarget != 'X')
                    wrongResidue.Add($"{mod.IdWithMotif} ({mod.LocationRestriction}) read back as {back.IdWithMotif}");
            }

            Assert.That(informative, Is.GreaterThan(100), "the population this test exists for must be present");
            Assert.That(wrongResidue, Is.Empty);
        }

        /// <summary>
        /// A terminal mod restricted to its terminus beats an "Anywhere." one on the same residue, and a
        /// mod on the residue itself beats one on any residue.
        /// </summary>
        [Test]
        public void Layer2_TerminalModPrefersTheTerminusThenTheResidue()
        {
            ModificationMotif.TryGetMotif("K", out var motifK);
            ModificationMotif.TryGetMotif("X", out var motifX);
            var anywhereK = new Modification(_originalId: "Acetyl", _modificationType: "testMods", _target: motifK,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 42.01057);
            var terminalX = new Modification(_originalId: "Acetyl", _modificationType: "testMods", _target: motifX,
                _locationRestriction: "N-terminal.", _monoisotopicMass: 42.01057);
            var allModsKnown = new Dictionary<string, Modification> { [anywhereK.IdWithMotif] = anywhereK, [terminalX.IdWithMotif] = terminalX };
            Assert.That(ProFormaConverter.ToModificationDictionary(ProFormaReader.Read("[Acetyl]-KPEPTIDE"), allModsKnown)[1], Is.SameAs(terminalX));

            var terminalK = new Modification(_originalId: "Acetyl", _modificationType: "testMods", _target: motifK,
                _locationRestriction: "N-terminal.", _monoisotopicMass: 42.01057);
            allModsKnown[terminalK.IdWithMotif] = terminalK;
            Assert.That(ProFormaConverter.ToModificationDictionary(ProFormaReader.Read("[Acetyl]-KPEPTIDE"), allModsKnown)[1], Is.SameAs(terminalK));
        }

        /// <summary>
        /// A term with no residues has no terminal residue to fit, so a terminal mod resolves to the first
        /// terminus-compatible candidate rather than indexing into an empty sequence.
        /// </summary>
        [Test]
        public void Layer2_TerminalModOnAnEmptySequenceTakesTheFirstTerminusCompatibleCandidate()
        {
            ModificationMotif.TryGetMotif("X", out var motifX);
            var terminalX = new Modification(_originalId: "Acetyl", _modificationType: "testMods", _target: motifX,
                _locationRestriction: "N-terminal.", _monoisotopicMass: 42.01057);
            var allModsKnown = new Dictionary<string, Modification> { [terminalX.IdWithMotif] = terminalX };

            var descriptors = new[] { new Tdp.ProFormaDescriptor(Tdp.ProFormaKey.Name, "Acetyl") };
            var term = new Tdp.ProFormaTerm("", null, descriptors, null, null, null, null, null);

            Assert.That(ProFormaConverter.ToModificationDictionary(term, allModsKnown)[1], Is.SameAs(terminalX));
        }

        [Test]
        public void Layer2_ResolvesTerminalMod_WhenRestrictionIsAnywhere()
        {
            // ToProFormaTerm writes any position-1 mod as an N-terminal descriptor regardless of its
            // LocationRestriction, so resolution must accept an "Anywhere." mod at the terminus to keep
            // write and read symmetric — this would previously throw on read.
            var nAcetyl = MakeTerminalMod("Acetyl", "Anywhere.", 42.01057);
            var allModsKnown = new Dictionary<string, Modification> { [nAcetyl.IdWithMotif] = nAcetyl };

            var dict = ProFormaConverter.ToModificationDictionary(ProFormaReader.Read("[Acetyl]-PEPTIDEK"), allModsKnown);
            Assert.That(dict[1], Is.SameAs(nAcetyl));
        }
    }
}
