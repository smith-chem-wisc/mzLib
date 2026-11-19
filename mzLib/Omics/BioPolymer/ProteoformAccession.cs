using MzLibUtil;

namespace Omics.BioPolymer
{
    /// <summary>
    /// An accession read back into its parent entry and the variants applied to it. mzLib names a
    /// biopolymer with applied sequence variants "{entry}_{variant}[_{variant}...]", e.g. "P12345_S70N" or
    /// "NP_000537.3_R72P_A80T" (<see cref="VariantApplication.GetAccession"/>); this is the result of
    /// reading such a name back with <see cref="VariantApplication.ParseAccession"/>. Nothing is dropped:
    /// the full accession, the entry and the applied variants are all kept.
    /// </summary>
    /// <param name="Verbatim">Exactly the string given (empty for null).</param>
    /// <param name="Entry">The parent entry; for an accession with no variant suffix, the accession
    /// itself. Unrecognized when neither the UniProt nor the RefSeq grammar matches.</param>
    /// <param name="AppliedVariants">The applied variants as written, without the leading "_" ("S70N" or
    /// "S70N_A80T"), or null when the accession names no applied variant.</param>
    public sealed record ProteoformAccession(string Verbatim, ProteinAccession Entry, string? AppliedVariants)
    {
        /// <summary>True when the accession names a biopolymer with applied variants.</summary>
        public bool HasAppliedVariants => AppliedVariants != null;
    }
}
