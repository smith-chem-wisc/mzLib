using System.IO;
using System.Text;
using MzIdentML;
using MzLibUtil;
using NUnit.Framework;

namespace Test.MzIdentML
{
    /// <summary>
    /// Correcting the root namespace from ".../mzIdentML/1.1.0" to the schema's ".../mzIdentML/1.1"
    /// fixed the write side, and would have broken the read side for every .mzID mzLib had already
    /// produced: XmlSerializer matches the element namespace exactly, and after the correction no type
    /// in the assembly declared the old value. All three attempts in
    /// <see cref="MzidIdentifications"/> failed, and the innermost has nothing below it to catch, so
    /// the constructor threw on this library's own output.
    ///
    /// These pin both directions: a compliant document reads, and a legacy one reads.
    /// </summary>
    [TestFixture]
    public sealed class MzIdentMlLegacyNamespaceTests
    {
        private const string LegacyNamespace = "http://psidev.info/psi/pi/mzIdentML/1.1.0";
        private const string SchemaNamespace = "http://psidev.info/psi/pi/mzIdentML/1.1";

        /// <summary>
        /// Deliberately carries a nested ParentTolerance three levels below the root. A namespace remap
        /// that only renames the root element deserializes the root and silently drops everything under
        /// it, so a fixture whose assertions stop at root attributes cannot tell the two apart.
        /// </summary>
        private static string Document(string ns) => $@"<?xml version=""1.0"" encoding=""utf-8""?>
<MzIdentML xmlns=""{ns}"" id=""legacy-fixture"" version=""1.1.0"">
  <cvList>
    <cv id=""PSI-MS"" fullName=""PSI-MS"" uri=""https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"" />
  </cvList>
  <AnalysisSoftwareList>
    <AnalysisSoftware id=""AS_metamorpheus"" name=""MetaMorpheus"" />
  </AnalysisSoftwareList>
  <AnalysisProtocolCollection>
    <SpectrumIdentificationProtocol id=""SIP"" analysisSoftware_ref=""AS_metamorpheus"">
      <SearchType>
        <cvParam accession=""MS:1001083"" name=""ms-ms search"" cvRef=""PSI-MS"" />
      </SearchType>
      <ParentTolerance>
        <cvParam accession=""MS:1001412"" name=""search tolerance plus value"" value=""5"" unitName=""parts per million"" cvRef=""PSI-MS"" />
        <cvParam accession=""MS:1001413"" name=""search tolerance minus value"" value=""5"" unitName=""parts per million"" cvRef=""PSI-MS"" />
      </ParentTolerance>
    </SpectrumIdentificationProtocol>
  </AnalysisProtocolCollection>
</MzIdentML>";

        private static string WriteTemp(string contents, string name)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, name);
            File.WriteAllText(path, contents, new UTF8Encoding(false));
            return path;
        }

        [Test]
        [TestCase(SchemaNamespace, "compliant.mzid")]
        [TestCase(LegacyNamespace, "legacy.mzid")]
        public void BothNamespacesAreReadable(string ns, string fileName)
        {
            string path = WriteTemp(Document(ns), fileName);

            try
            {
                // The defect was a throw out of the constructor, so constructing at all is the assertion.
                MzidIdentifications identifications = null;
                Assert.That(() => identifications = new MzidIdentifications(path), Throws.Nothing,
                    $"a document in {ns} could not be read");
                Assert.That(identifications, Is.Not.Null);
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        /// <summary>
        /// The guard against the tempting wrong fix. Handing XmlSerializer an XmlRootAttribute override
        /// for the legacy namespace makes the constructor succeed while leaving every child bound to
        /// ".../1.1", so a file full of identifications reads back empty -- silent loss rather than a
        /// throw. Reaching a nested value is what separates the two.
        /// </summary>
        [Test]
        [TestCase(SchemaNamespace)]
        [TestCase(LegacyNamespace)]
        public void NestedElementsSurviveTheNamespaceRemap(string ns)
        {
            string path = WriteTemp(Document(ns), "nested.mzid");

            try
            {
                var identifications = new MzidIdentifications(path);

                Tolerance parentTolerance = identifications.ParentTolerance;

                Assert.Multiple(() =>
                {
                    Assert.That(parentTolerance, Is.Not.Null,
                        "the root deserialized but AnalysisProtocolCollection did not: the remap is root-only");
                    Assert.That(parentTolerance, Is.TypeOf<PpmTolerance>());
                    Assert.That(parentTolerance.Value, Is.EqualTo(5).Within(1e-9));
                });
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }
    }
}
