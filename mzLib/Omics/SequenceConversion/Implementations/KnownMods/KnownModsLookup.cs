using System.Collections.Generic;

using Omics.Modifications;

namespace Omics.SequenceConversion;

public class KnownModsLookup : ModificationLookupBase
{
    private readonly Dictionary<string, Modification>? _knownMods;

    public KnownModsLookup(Dictionary<string, Modification> knownMods, double massTolerance = 0.001)
        : base(knownMods?.Values, massTolerance)
    {
        _knownMods = knownMods;
    }

    public override string Name => "Known Mods";

    public override CanonicalModification? TryResolve(CanonicalModification mod)
    {
        if (_knownMods != null)
        {
            foreach (var alias in GetLegacyAliasCandidates(mod))
            {
                if (_knownMods.TryGetValue(alias, out var resolved))
                {
                    return mod.WithResolvedModification(resolved, mod.ResidueIndex, mod.PositionType);
                }
            }
        }

        return base.TryResolve(mod);
    }

    private static IEnumerable<string> GetLegacyAliasCandidates(CanonicalModification mod)
    {
        if (!string.IsNullOrWhiteSpace(mod.MzLibId))
        {
            yield return mod.MzLibId;

            var mzLibSuffix = GetSuffixAfterFirstColon(mod.MzLibId);
            if (mzLibSuffix != null)
            {
                yield return mzLibSuffix;
            }
        }

        if (!string.IsNullOrWhiteSpace(mod.OriginalRepresentation))
        {
            yield return mod.OriginalRepresentation;

            var representationSuffix = GetSuffixAfterFirstColon(mod.OriginalRepresentation);
            if (representationSuffix != null)
            {
                yield return representationSuffix;
            }
        }
    }

    private static string? GetSuffixAfterFirstColon(string value)
    {
        var colonIndex = value.IndexOf(':');
        if (colonIndex <= 0 || colonIndex >= value.Length - 1)
        {
            return null;
        }

        return value[(colonIndex + 1)..];
    }
}
