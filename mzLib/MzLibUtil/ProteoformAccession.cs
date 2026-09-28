using System.Text.RegularExpressions;

namespace MzLibUtil
{
    /// <summary>
    /// An accession as a search reports it, which may name a proteoform rather than an entry: mzLib
    /// names a protein with applied sequence variants "{entry}_{variant}[_{variant}...]", e.g.
    /// "P12345_S70N" or "NP_000537.3_R72P_A80T" (VariantApplication.GetAccession). This splits that
    /// name into the entry, parsed as a <see cref="ProteinAccession"/>, and the variant suffix, so a
    /// proteoform can be joined to anything keyed by its entry.
    ///
    /// The entry is matched against the full accession grammars, never cut at the first "_": that
    /// would turn "NP_000537" into "NP". A name that only looks like a proteoform stays Unrecognized
    /// and keeps its own key. The load-collision counter "P12345_2" names a DIFFERENT entry whose
    /// accession collided, so mapping it back to P12345 would give it another protein's answers. Decoy,
    /// contaminant and entrapment prefixes stay Unrecognized too.
    /// </summary>
    /// <param name="Verbatim">Exactly the string given (empty for null).</param>
    /// <param name="Entry">The entry the proteoform belongs to; for an accession with no variant
    /// suffix, the accession itself. Unrecognized when neither grammar matches.</param>
    /// <param name="AppliedVariants">The variant suffix without its leading "_" ("S70N" or
    /// "S70N_A80T"), or null when the accession names no variant.</param>
    public sealed record ProteoformAccession(string Verbatim, ProteinAccession Entry, string AppliedVariants)
    {
        /// <summary>
        /// One or more "_{original}{position}{variant}" tokens, as SequenceVariation.SimpleString writes
        /// them: the original residues are never empty, and the variant residues are empty for a deletion.
        /// </summary>
        private static readonly Regex VariantSuffix = new(@"^(.+?)((?:_[A-Z]+\d+[A-Z]*)+)$", RegexOptions.Compiled);

        /// <summary>True when the accession names a proteoform with applied variants.</summary>
        public bool HasAppliedVariants => AppliedVariants != null;

        /// <summary>
        /// Parses, never repairs. An accession that is an entry is returned with no variants; one that is
        /// an entry plus a variant suffix is split; anything else is Unrecognized and kept verbatim.
        /// </summary>
        public static ProteoformAccession Parse(string accession)
        {
            accession ??= "";

            var plain = ProteinAccession.Parse(accession);
            if (plain.Namespace != AccessionNamespace.Unrecognized)
            {
                return new ProteoformAccession(accession, plain, null);
            }

            var m = VariantSuffix.Match(accession);
            if (m.Success)
            {
                var entry = ProteinAccession.Parse(m.Groups[1].Value);
                if (entry.Namespace != AccessionNamespace.Unrecognized)
                {
                    return new ProteoformAccession(accession, entry, m.Groups[2].Value.Substring(1));
                }
            }

            return new ProteoformAccession(accession, plain, null);
        }
    }
}
