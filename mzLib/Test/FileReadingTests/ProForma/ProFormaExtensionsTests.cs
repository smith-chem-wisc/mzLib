using System.Collections.Generic;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Readers.ProForma;

namespace Test.FileReadingTests.ProForma
{
    /// <summary>
    /// End-to-end bridge: an mzLib peptide (built from its native FullSequence) -&gt; ProForma string.
    /// This is the path MetaMorpheus top-down output will use.
    /// </summary>
    [TestFixture]
    internal class ProFormaExtensionsTests
    {
        [Test]
        public void Peptide_ToProFormaString_NameMod()
        {
            ModificationMotif.TryGetMotif("M", out var motifM);
            var ox = new Modification(_originalId: "Oxidation", _modificationType: "testMods", _target: motifM,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491);
            var allKnown = new Dictionary<string, Modification> { [ox.IdWithMotif] = ox };

            // native mzLib full sequence -> peptide -> ProForma
            var peptide = new PeptideWithSetModifications("PEM[testMods:Oxidation on M]TIDEK", allKnown);
            Assert.That(peptide.ToProFormaString(), Is.EqualTo("PEM[Oxidation]TIDEK"));
        }

        [Test]
        public void Peptide_ToProFormaString_AccessionMod()
        {
            ModificationMotif.TryGetMotif("M", out var motifM);
            var ox = new Modification(_originalId: "Oxidation", _modificationType: "Unimod", _target: motifM,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 15.99491,
                _databaseReference: new Dictionary<string, IList<string>> { ["Unimod"] = new List<string> { "35" } });
            var allKnown = new Dictionary<string, Modification> { [ox.IdWithMotif] = ox };

            var peptide = new PeptideWithSetModifications("PEM[Unimod:Oxidation on M]TIDEK", allKnown);
            // mod carries a Unimod ref, so it emits as an accession descriptor
            Assert.That(peptide.ToProFormaString(), Is.EqualTo("PEM[UNIMOD:35]TIDEK"));
        }
    }
}
