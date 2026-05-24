using Omics;
using Tdp = TopDownProteomics.ProForma;

namespace Readers.ProForma
{
    /// <summary>
    /// Convenience bridge between mzLib's <see cref="IBioPolymerWithSetMods"/> and ProForma 2.0.
    /// This is the surface MetaMorpheus consumes: a single call turns a scored peptide/proteoform
    /// into its ProForma string. Conversion is Layer-2 (per-residue + terminal mods); see
    /// <see cref="ProFormaConverter"/> for the supported subset.
    /// </summary>
    public static class ProFormaExtensions
    {
        /// <summary>
        /// Builds the ProForma <see cref="Tdp.ProFormaTerm"/> for a biopolymer from its base sequence
        /// and <c>AllModsOneIsNterminus</c>.
        /// </summary>
        public static Tdp.ProFormaTerm ToProFormaTerm(this IBioPolymerWithSetMods bioPolymer)
            => ProFormaConverter.ToProFormaTerm(bioPolymer.BaseSequence, bioPolymer.AllModsOneIsNterminus);

        /// <summary>
        /// Writes the canonical ProForma 2.0 string for a biopolymer (e.g. <c>EM[Oxidation]EVEES[Phospho]PEK</c>).
        /// </summary>
        public static string ToProFormaString(this IBioPolymerWithSetMods bioPolymer)
            => ProFormaWriter.Write(bioPolymer.ToProFormaTerm());
    }
}
