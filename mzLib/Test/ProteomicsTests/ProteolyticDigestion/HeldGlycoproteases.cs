using NUnit.Framework;
using Omics.Digestion;
using Proteomics.ProteolyticDigestion;
using System.IO;

namespace Test.ProteomicsTests.ProteolyticDigestion
{
    /// <summary>
    /// OpeRATOR, IMPa, SmE and their trypsin co-digests: rows written for proteases.tsv and held out of
    /// it until MetaMorpheus can set RespectCleavagePromotingModifications. With the flag off each is a
    /// cut before every Ser and Thr, and MetaMorpheus lists every dictionary entry in its protease
    /// pickers. See TestData/held-glycoproteases.tsv.
    /// </summary>
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal static class HeldGlycoproteases
    {
        private static string HeldPath => Path.Combine(TestContext.CurrentContext.TestDirectory,
            "ProteomicsTests", "ProteolyticDigestion", "TestData", "held-glycoproteases.tsv");

        /// <summary>
        /// Merges the held rows into <see cref="ProteaseDictionary.Dictionary"/> through the same parser
        /// proteases.tsv goes through. Safe to call again, including after another fixture has reset the
        /// dictionary: a name already present is skipped, and its definition is the same.
        /// </summary>
        internal static CustomDigestionAgentLoadResult Register() =>
            ProteaseDictionary.LoadAndMergeCustomProteases(HeldPath);
    }
}
