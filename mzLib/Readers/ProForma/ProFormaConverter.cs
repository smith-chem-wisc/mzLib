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
        private enum Terminus { None, N, C }

        /// <summary>mzLib <see cref="Modification.DatabaseReference"/> key -&gt; ProForma accession prefix.</summary>
        private static readonly Dictionary<string, string> DbKeyToProFormaPrefix = new(StringComparer.OrdinalIgnoreCase)
        {
            ["Unimod"] = "UNIMOD",
            ["PSI-MOD"] = "MOD",
            ["RESID"] = "RESID",
        };

        /// <summary>ProForma accession prefix -&gt; SDK evidence type (for the inverse direction).</summary>
        private static readonly Dictionary<string, Tdp.ProFormaEvidenceType> PrefixToEvidence = new(StringComparer.OrdinalIgnoreCase)
        {
            ["UNIMOD"] = Tdp.ProFormaEvidenceType.Unimod,
            ["MOD"] = Tdp.ProFormaEvidenceType.PsiMod,
            ["RESID"] = Tdp.ProFormaEvidenceType.Resid,
        };

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
            var byAccession = BuildAccessionIndex(allModsKnown.Values);
            var result = new Dictionary<int, Modification>();

            if (term.NTerminalDescriptors is { Count: > 0 })
                result[1] = Resolve(term.NTerminalDescriptors, residue: null, Terminus.N, allModsKnown, byAccession, term);

            if (term.Tags != null)
            {
                foreach (var tag in term.Tags)
                {
                    int residueIndex = tag.ZeroBasedStartIndex;
                    result[residueIndex + 2] = Resolve(tag.Descriptors, sequence[residueIndex], Terminus.None, allModsKnown, byAccession, term);
                }
            }

            if (term.CTerminalDescriptors is { Count: > 0 })
                result[sequence.Length + 2] = Resolve(term.CTerminalDescriptors, residue: null, Terminus.C, allModsKnown, byAccession, term);

            return result;
        }

        /// <summary>
        /// Builds a term from a base sequence and mzLib modification dictionary (the inverse of
        /// <see cref="ToModificationDictionary"/>). Each modification is emitted as an accession
        /// (Identifier) descriptor when it carries a recognized ontology reference, else as a Name.
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

        /// <summary>
        /// Emits an accession (Identifier) descriptor like <c>[UNIMOD:35]</c> when the modification
        /// carries a recognized ontology reference, otherwise a Name descriptor like <c>[Oxidation]</c>.
        /// </summary>
        private static Tdp.ProFormaDescriptor BuildDescriptor(Modification mod)
        {
            if (mod.DatabaseReference != null)
            {
                foreach (var (dbKey, ids) in mod.DatabaseReference)
                {
                    if (DbKeyToProFormaPrefix.TryGetValue(dbKey, out var prefix) && ids is { Count: > 0 })
                        return new Tdp.ProFormaDescriptor(Tdp.ProFormaKey.Identifier, PrefixToEvidence[prefix], $"{prefix}:{ids[0]}");
                }
            }
            return new Tdp.ProFormaDescriptor(Tdp.ProFormaKey.Name, mod.OriginalId);
        }

        /// <summary>
        /// Indexes modifications by their ProForma accession string (e.g. <c>"UNIMOD:35"</c>, upper-cased).
        /// One accession can map to several modifications differing by motif, so values are lists.
        /// </summary>
        private static Dictionary<string, List<Modification>> BuildAccessionIndex(IEnumerable<Modification> mods)
        {
            var index = new Dictionary<string, List<Modification>>(StringComparer.OrdinalIgnoreCase);
            foreach (var mod in mods)
            {
                if (mod.DatabaseReference == null) continue;
                foreach (var (dbKey, ids) in mod.DatabaseReference)
                {
                    if (!DbKeyToProFormaPrefix.TryGetValue(dbKey, out var prefix)) continue;
                    foreach (var id in ids)
                    {
                        string key = $"{prefix}:{id}".ToUpperInvariant();
                        if (!index.TryGetValue(key, out var list))
                            index[key] = list = new List<Modification>();
                        list.Add(mod);
                    }
                }
            }
            return index;
        }

        /// <summary>
        /// Resolves the first descriptor that maps to a known modification. Identifier descriptors are
        /// matched by accession; Name descriptors by mzLib's <c>"{name} on {residue}"</c> IdWithMotif
        /// convention. Interior residues match on target motif; termini match on
        /// <see cref="Modification.LocationRestriction"/>.
        /// </summary>
        private static Modification Resolve(IList<Tdp.ProFormaDescriptor> descriptors, char? residue, Terminus terminus,
            Dictionary<string, Modification> allModsKnown, Dictionary<string, List<Modification>> byAccession,
            Tdp.ProFormaTerm term)
        {
            foreach (var d in descriptors)
            {
                if (d.Key == Tdp.ProFormaKey.Identifier
                    && byAccession.TryGetValue(d.Value, out var candidates)
                    && SelectCandidate(candidates, residue, terminus) is { } byAcc)
                    return byAcc;

                if (d.Key == Tdp.ProFormaKey.Name
                    && ResolveName(d.Value, residue, terminus, allModsKnown) is { } byName)
                    return byName;
            }

            throw new MzLibException(
                $"No known modification resolves descriptor(s) [{string.Join("|", descriptors.Select(d => $"{d.Key}:{d.Value}"))}] " +
                $"at {Where(residue, terminus)} in ProForma term '{term.Sequence}'.");
        }

        /// <summary>Chooses an accession candidate matching the residue motif (interior) or terminus restriction.</summary>
        private static Modification? SelectCandidate(List<Modification> candidates, char? residue, Terminus terminus)
        {
            if (terminus != Terminus.None)
                return candidates.FirstOrDefault(m => IsTerminusCompatible(m, terminus));
            if (residue.HasValue)
                return candidates.FirstOrDefault(m => MotifMatches(m, residue.Value));
            return candidates.Count == 1 ? candidates[0] : null;
        }

        /// <summary>
        /// Resolves a Name descriptor. Interior residues use the <c>"{name} on {residue}"</c> key; termini
        /// try <c>"{name} on X"</c> first (the usual terminal motif) then any same-named terminus-compatible mod.
        /// </summary>
        private static Modification? ResolveName(string name, char? residue, Terminus terminus,
            Dictionary<string, Modification> allModsKnown)
        {
            if (terminus != Terminus.None)
            {
                if (allModsKnown.TryGetValue($"{name} on X", out var onX) && IsTerminusCompatible(onX, terminus))
                    return onX;
                return allModsKnown.Values.FirstOrDefault(m => m.OriginalId == name && IsTerminusCompatible(m, terminus));
            }
            if (residue.HasValue && allModsKnown.TryGetValue($"{name} on {residue}", out var byMotif))
                return byMotif;
            return allModsKnown.TryGetValue(name, out var bare) ? bare : null;
        }

        private static bool MotifMatches(Modification mod, char residue)
            => mod.Target != null && mod.Target.ToString().Equals(residue.ToString(), StringComparison.Ordinal);

        private static bool IsTerminusCompatible(Modification mod, Terminus terminus) => terminus switch
        {
            Terminus.N => mod.LocationRestriction is "N-terminal." or "Peptide N-terminal.",
            Terminus.C => mod.LocationRestriction is "C-terminal." or "Peptide C-terminal.",
            _ => false,
        };

        private static string Where(char? residue, Terminus terminus) => terminus switch
        {
            Terminus.N => "the N-terminus",
            Terminus.C => "the C-terminus",
            _ => $"residue '{residue}'",
        };

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
