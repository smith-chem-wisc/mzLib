using System.Collections.Generic;
using System.Linq;
using MzLibUtil;
using NUnit.Framework;
using Omics.Modifications;
using Readers.ProForma;

namespace Test.FileReadingTests.ProForma
{
    /// <summary>
    /// Layer-2: ProFormaTerm &lt;-&gt; mzLib (base sequence + AllModsOneIsNterminus). First slice
    /// covers per-residue name modifications and the unsupported-feature guard. Accession-based
    /// resolution and terminal mods are later slices (see known-limitations.md).
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

            // inverse: rebuild a term and confirm it writes back to the same canonical string
            var rebuilt = ProFormaConverter.ToProFormaTerm(term.Sequence, dict);
            Assert.That(ProFormaWriter.Write(rebuilt), Is.EqualTo("EM[Oxidation]EVEES[Phospho]PEK"));
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
