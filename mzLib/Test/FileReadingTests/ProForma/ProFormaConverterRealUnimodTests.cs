using System.Collections.Generic;
using System.IO;
using NUnit.Framework;
using Omics.Modifications;
using Readers.ProForma;
using UsefulProteomicsDatabases;

namespace Test.FileReadingTests.ProForma
{
    /// <summary>
    /// Layer-2 validation against the REAL bundled Unimod database (not hand-built mods), confirming
    /// ProFormaConverter resolves real names and accessions, picks the correct motif amid real
    /// accession ambiguity (e.g. UNIMOD:21 → Phospho on S/T/Y), and round-trips on actual data.
    /// </summary>
    [TestFixture]
    internal class ProFormaConverterRealUnimodTests
    {
        private Dictionary<string, Modification> _allModsKnown = null!;

        [OneTimeSetUp]
        public void LoadRealUnimod()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "unimod_tables2.xml");
            _allModsKnown = new Dictionary<string, Modification>();
            foreach (var mod in Loaders.LoadUnimod(path))
                _allModsKnown.TryAdd(mod.IdWithMotif, mod);
        }

        [Test]
        public void Resolves_RealUnimod_ByName()
        {
            var dict = ProFormaConverter.ToModificationDictionary(
                ProFormaReader.Read("EM[Oxidation]EVEES[Phospho]PEK"), _allModsKnown);

            Assert.That(dict[3].OriginalId, Is.EqualTo("Oxidation"));
            Assert.That(dict[3].Target.ToString(), Is.EqualTo("M"));
            Assert.That(dict[3].DatabaseReference["Unimod"], Does.Contain("35"));
            Assert.That(dict[8].OriginalId, Is.EqualTo("Phospho"));
            Assert.That(dict[8].Target.ToString(), Is.EqualTo("S"));
            Assert.That(dict[8].DatabaseReference["Unimod"], Does.Contain("21"));
        }

        [Test]
        public void Resolves_RealUnimod_ByAccession_MotifAware_AndRoundTrips()
        {
            var term = ProFormaReader.Read("EM[UNIMOD:35]EVEES[UNIMOD:21]PEK");
            var dict = ProFormaConverter.ToModificationDictionary(term, _allModsKnown);

            // UNIMOD:21 is Phospho on S/T/Y; the residue (S) must select the on-S variant.
            Assert.That(dict[3].OriginalId, Is.EqualTo("Oxidation"));
            Assert.That(dict[3].Target.ToString(), Is.EqualTo("M"));
            Assert.That(dict[8].OriginalId, Is.EqualTo("Phospho"));
            Assert.That(dict[8].Target.ToString(), Is.EqualTo("S"));

            // mods carry Unimod refs, so the inverse emits accession descriptors
            var rebuilt = ProFormaConverter.ToProFormaTerm(term.Sequence, dict);
            Assert.That(ProFormaWriter.Write(rebuilt), Is.EqualTo("EM[UNIMOD:35]EVEES[UNIMOD:21]PEK"));
        }
    }
}
