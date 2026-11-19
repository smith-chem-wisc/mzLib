using Omics.Modifications;
using Omics.SequenceConversion;

namespace Readers
{
    /// <summary>
    /// Places a modification given at a PSI location into an AllModsOneIsNterminus dictionary. mzIdentML and mzTab
    /// number locations the same way: 0 is the N-terminus, 1 to the peptide length a residue, length + 1 the
    /// C-terminus.
    /// </summary>
    internal static class PsiModificationLocation
    {
        /// <summary>
        /// Null when <paramref name="location"/> is on <paramref name="baseSequence"/>, from the N-terminus (0) to the
        /// C-terminus (length + 1); otherwise why it is not.
        /// </summary>
        internal static string? ReasonNotOnPeptide(string baseSequence, int location) =>
            location < 0 || location > baseSequence.Length + 1
                ? $"modification location {location} is not on the peptide"
                : null;

        /// <summary>
        /// Resolves the Unimod modification <paramref name="unimodAccession"/> at <paramref name="location"/> and adds
        /// it to <paramref name="mods"/>.
        /// </summary>
        /// <param name="mods">Modifications keyed the AllModsOneIsNterminus way. Unchanged when this returns false.</param>
        /// <param name="baseSequence">The unmodified peptide, not empty. A terminal modification's residue is taken
        /// from it.</param>
        /// <param name="location">The PSI location.</param>
        /// <param name="unimodAccession">The accession as written, such as "UNIMOD:35", used in the reason.</param>
        /// <param name="unimodId">The accession's number.</param>
        /// <param name="massDelta">The mass shift the file gives, if any.</param>
        /// <param name="lookup">The Unimod catalog to resolve against.</param>
        /// <param name="reason">Null when the modification was added; otherwise why it was not.</param>
        internal static bool TryAddUnimod(Dictionary<int, Modification> mods, string baseSequence, int location,
            string unimodAccession, int unimodId, double? massDelta, UnimodModificationLookup lookup, out string? reason)
        {
            reason = ReasonNotOnPeptide(baseSequence, location);
            if (reason != null)
            {
                return false;
            }

            // AllModsOneIsNterminus: 1 is the N-terminus, residue r (zero-based) is r + 2, the C-terminus length + 2.
            // Both numberings count the N-terminus as a position, so every key is the location plus one.
            var canonical = location switch
            {
                0 => CanonicalModification.AtNTerminus(unimodAccession, baseSequence[0], massDelta, unimodId: unimodId),
                _ when location == baseSequence.Length + 1 => CanonicalModification.AtCTerminus(unimodAccession, baseSequence[^1], massDelta, unimodId: unimodId),
                _ => CanonicalModification.AtResidue(location - 1, baseSequence[location - 1], unimodAccession, massDelta, unimodId: unimodId),
            };
            int key = location + 1;

            var resolved = lookup.TryResolve(canonical)?.MzLibModification;
            if (resolved == null)
            {
                reason = $"{unimodAccession} at location {location} does not resolve to a Unimod modification";
                return false;
            }

            if (!mods.TryAdd(key, resolved))
            {
                reason = $"more than one modification at location {location}";
                return false;
            }

            return true;
        }
    }
}
