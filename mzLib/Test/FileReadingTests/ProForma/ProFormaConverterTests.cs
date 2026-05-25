using System.Collections.Generic;
using MzLibUtil;
using NUnit.Framework;
using Omics.Modifications;
using Readers.ProForma;

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
            Assert.Throws<MzLibException>(
                () => ProFormaConverter.ToModificationDictionary(term, new Dictionary<string, Modification>()));
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

            var dict = ProFormaConverter.ToModificationDictionary(ProFormaReader.Read("[Acetyl]-PEPTIDEK"), allModsKnown);
            Assert.That(dict[1], Is.SameAs(nAcetyl));
        }
    }
}
