using MzLibUtil;
using Omics.Modifications;
using Tdp = TopDownProteomics.ProForma;

namespace Readers.ProForma
{
    /// <summary>
    /// Layer-2 conversion between a <see cref="Tdp.ProFormaTerm"/> and mzLib's
    /// (base sequence + <c>AllModsOneIsNterminus</c>) representation that
    /// <c>IBioPolymerWithSetMods</c> consumes.
    /// <para>
    /// Lossy by design: only per-residue and N-/C-terminal modifications that resolve to a known
    /// mzLib <see cref="Modification"/> are mapped. Crosslinks, branches, labile mods, ambiguity
    /// groups, position ranges, sequence ambiguity, and global isotope mods are out of scope and
    /// cause <see cref="ToModificationDictionary"/> to throw rather than silently lose data.
    /// </para>
    /// Index convention matches mzLib's <c>AllModsOneIsNterminus</c>: N-terminus = 1,
    /// residue r (0-based) = r + 2, C-terminus = baseSequence.Length + 2 — the same keys used by
    /// <c>IBioPolymerWithSetMods.GetModificationDictionaryFromFullSequence</c>/<c>DetermineFullSequence</c>.
    /// </summary>
    public static class ProFormaConverter
    {
        /// <summary>
        /// Converts a term into the (one-is-N-terminus) modification dictionary mzLib uses.
        /// The base sequence is available from <see cref="Tdp.ProFormaTerm.Sequence"/>.
        /// </summary>
        /// <param name="term">A term containing only Layer-2-supported features.</param>
        /// <param name="allModsKnown">Loaded modifications keyed by <see cref="Modification.IdWithMotif"/>.</param>
        /// <returns>Modifications keyed by one-based position (1 = N-terminus).</returns>
        /// <exception cref="MzLibException">
        /// If the term uses an unsupported feature, or a descriptor cannot be resolved to a known modification.
        /// </exception>
        public static Dictionary<int, Modification> ToModificationDictionary(Tdp.ProFormaTerm term,
            Dictionary<string, Modification> allModsKnown)
        {
            EnsureLayer2Supported(term);

            string sequence = term.Sequence;
            var result = new Dictionary<int, Modification>();

            if (term.NTerminalDescriptors is { Count: > 0 })
                result[1] = Resolve(term.NTerminalDescriptors, residue: null, allModsKnown, term);

            if (term.Tags != null)
            {
                foreach (var tag in term.Tags)
                {
                    int residueIndex = tag.ZeroBasedStartIndex;
                    result[residueIndex + 2] = Resolve(tag.Descriptors, sequence[residueIndex], allModsKnown, term);
                }
            }

            if (term.CTerminalDescriptors is { Count: > 0 })
                result[sequence.Length + 2] = Resolve(term.CTerminalDescriptors, residue: null, allModsKnown, term);

            return result;
        }

        /// <summary>
        /// Builds a term from a base sequence and mzLib modification dictionary (the inverse of
        /// <see cref="ToModificationDictionary"/>). Each modification is emitted as a Name descriptor.
        /// </summary>
        /// <param name="baseSequence">The unmodified residue sequence.</param>
        /// <param name="allModsOneIsNterminus">Modifications keyed by one-based position (1 = N-terminus).</param>
        public static Tdp.ProFormaTerm ToProFormaTerm(string baseSequence,
            IDictionary<int, Modification> allModsOneIsNterminus)
        {
            List<Tdp.ProFormaDescriptor>? nTerm = null;
            List<Tdp.ProFormaDescriptor>? cTerm = null;
            var tags = new List<Tdp.ProFormaTag>();

            foreach (var (position, mod) in allModsOneIsNterminus)
            {
                if (position == 1)
                    (nTerm ??= new()).Add(BuildDescriptor(mod));
                else if (position == baseSequence.Length + 2)
                    (cTerm ??= new()).Add(BuildDescriptor(mod));
                else
                    tags.Add(new Tdp.ProFormaTag(position - 2, new[] { BuildDescriptor(mod) }));
            }

            tags.Sort((a, b) => a.ZeroBasedStartIndex.CompareTo(b.ZeroBasedStartIndex));

            return new Tdp.ProFormaTerm(baseSequence, tags.Count > 0 ? tags : null, nTerm, cTerm,
                labileDescriptors: null, unlocalizedTags: null, tagGroups: null, globalModifications: null);
        }

        /// <summary>Emits a Name descriptor (e.g. <c>[Oxidation]</c>) for a modification.</summary>
        private static Tdp.ProFormaDescriptor BuildDescriptor(Modification mod)
            => new(Tdp.ProFormaKey.Name, mod.OriginalId);

        /// <summary>
        /// Resolves the first descriptor that maps to a known modification. Name descriptors are
        /// matched to <c>"{name} on {residue}"</c> (mzLib's IdWithMotif convention).
        /// </summary>
        private static Modification Resolve(IList<Tdp.ProFormaDescriptor> descriptors, char? residue,
            Dictionary<string, Modification> allModsKnown, Tdp.ProFormaTerm term)
        {
            foreach (var d in descriptors)
            {
                if (d.Key == Tdp.ProFormaKey.Name && residue.HasValue
                    && allModsKnown.TryGetValue($"{d.Value} on {residue}", out var byMotif))
                    return byMotif;

                // bare-name fallback (terminal mods or motif-free ids)
                if (d.Key == Tdp.ProFormaKey.Name && allModsKnown.TryGetValue(d.Value, out var byName))
                    return byName;
            }

            throw new MzLibException(
                $"No known modification resolves descriptor(s) [{string.Join("|", descriptors.Select(d => $"{d.Key}:{d.Value}"))}] " +
                $"on residue '{residue?.ToString() ?? "terminus"}' in ProForma term '{term.Sequence}'.");
        }

        /// <summary>
        /// Throws if the term uses any feature outside the Layer-2 supported subset (per-residue and
        /// terminal mods only). Fails loud so callers never silently drop crosslinks, ranges, etc.
        /// </summary>
        private static void EnsureLayer2Supported(Tdp.ProFormaTerm term)
        {
            if (term.TagGroups is { Count: > 0 })
                throw new MzLibException("ProForma tag groups (crosslinks/localization groups) are not supported at Layer 2.");
            if (term.LabileDescriptors is { Count: > 0 })
                throw new MzLibException("ProForma labile modifications are not supported at Layer 2.");
            if (term.UnlocalizedTags is { Count: > 0 })
                throw new MzLibException("ProForma unlocalized modifications are not supported at Layer 2.");
            if (term.GlobalModifications is { Count: > 0 })
                throw new MzLibException("ProForma global modifications are not supported at Layer 2.");
            if (term.Tags != null)
            {
                foreach (var tag in term.Tags)
                {
                    if (tag.ZeroBasedStartIndex != tag.ZeroBasedEndIndex)
                        throw new MzLibException("ProForma position ranges are not supported at Layer 2.");
                    if (tag.HasAmbiguousSequence)
                        throw new MzLibException("ProForma sequence ambiguity is not supported at Layer 2.");
                }
            }
        }
    }
}
