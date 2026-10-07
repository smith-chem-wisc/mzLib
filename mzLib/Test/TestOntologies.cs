using System.IO;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Paths to the trimmed ontology fixtures in <c>DataFiles/Ontologies</c>, for tests that need a
    /// realistic modification list but are not themselves testing the download.
    ///
    /// Using these keeps such tests out of the network entirely. <see cref="UsefulProteomicsDatabases.Loaders"/>'s
    /// <c>Load*</c> methods download only when the target file is absent, so a test that points at a
    /// path with no file behind it silently becomes an external-service test — which is how a UniProt
    /// outage once reddened the required CI job through ~150 parser and digestion tests.
    ///
    /// These are NOT a substitute for the live canaries. Tests that genuinely exercise the download
    /// path keep their own <c>*2</c> filenames, which are deliberately absent from the output
    /// directory, and carry <c>[Category("ExternalService")]</c>. See DataFiles/Ontologies/PROVENANCE.md.
    /// </summary>
    public static class TestOntologies
    {
        private static string Path_(string fileName) =>
            Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "Ontologies", fileName);

        /// <summary>UniProt PTM list, complete. Read by <c>Loaders.LoadUniprot</c>.</summary>
        public static string PtmList => Path_("ptmlist.txt");

        /// <summary>
        /// PSI-MOD in XML form, trimmed to the terms carrying a formal charge. Read by
        /// <c>Loaders.LoadPsiMod</c>; yields the same formal-charges dictionary as the full ontology.
        /// </summary>
        public static string PsiModXml => Path_("PSI-MOD.obo.xml");

        /// <summary>
        /// Unimod, trimmed to the Phospho and Oxidation records. Read by <c>Loaders.LoadUnimod</c>.
        /// Not suitable for tests that count modifications.
        /// </summary>
        public static string Unimod => Path_("unimod_tables.xml");
    }
}
