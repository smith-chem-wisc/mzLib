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
    }
}
