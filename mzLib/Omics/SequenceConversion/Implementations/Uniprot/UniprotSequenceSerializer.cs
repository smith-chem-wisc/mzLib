using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Omics.SequenceConversion;

/// <summary>
/// Serializes canonical sequences to the UniProt modification notation (e.g., [UniProt:Phosphoserine on S]).
/// UniProt has special serialization rules:
/// - "Methionine sulfoxide" REPLACES the M residue (output: [Methionine sulfoxide] without M)
/// - "Phosphoserine/Phosphothreonine" PRESERVES the residue (output: S[Phosphoserine])
/// - "Carbamidomethyl" (fixed modifications) are REMOVED entirely from output
/// </summary>
public class UniProtSequenceSerializer : SequenceSerializerBase
{
    public static UniProtSequenceSerializer Instance { get; } = new();

    /// <summary>
    /// Modification names that REPLACE the residue letter (the mod name represents the entire modified amino acid).
    /// NOTE: Most UniProt modifications DO NOT replace the residue. For example, "Methionine sulfoxide"
    /// is written as M[Methionine sulfoxide], not [Methionine sulfoxide].
    /// This set should only contain modifications where the residue is truly replaced/removed.
    /// Currently empty - we include all residue letters in the output.
    /// </summary>
    private static readonly HashSet<string> ResidueReplacingMods = new(StringComparer.OrdinalIgnoreCase)
    {
        // Currently no known residue-replacing modifications in standard UniProt format
    };

    /// <summary>
    /// Modification types that should be completely removed from output (e.g., fixed modifications).
    /// </summary>
    private static readonly HashSet<string> SuppressedModTypes = new(StringComparer.OrdinalIgnoreCase)
    {
        "Common Fixed"
    };

    /// <summary>
    /// Specific modification names that should be suppressed (regardless of type).
    /// </summary>
    private static readonly HashSet<string> SuppressedModNames = new(StringComparer.OrdinalIgnoreCase)
    {
        "Carbamidomethyl",
        "S-carbamoylmethylcysteine"
    };

    public UniProtSequenceSerializer(SequenceFormatSchema? schema = null, IModificationLookup? lookup = null)
        : base(lookup ?? UniProtModificationLookup.Instance)
    {
        Schema = schema ?? UniProtSequenceSchema.Instance;
    }

    public override string FormatName => Schema.FormatName;

    public override SequenceFormatSchema Schema { get; }

    public override bool CanSerialize(CanonicalSequence sequence) => !string.IsNullOrEmpty(sequence.BaseSequence);

    /// <summary>
    /// Override to set the target residue for N-terminal modifications based on the sequence's first character.
    /// This enables lookup of residue-specific N-terminal mods like "N-acetylproline".
    /// </summary>
    protected override CanonicalSequence EnrichModificationsIfNeeded(CanonicalSequence sequence)
    {
        if (!sequence.HasModifications || string.IsNullOrEmpty(sequence.BaseSequence))
        {
            return base.EnrichModificationsIfNeeded(sequence);
        }

        // Check if we have an N-terminal mod that needs its target residue set
        // TargetResidue may be null OR 'X' (generic N-terminal target) - both need updating
        var nTermMod = sequence.NTerminalModification;
        if (nTermMod.HasValue && 
            (nTermMod.Value.TargetResidue == 'X' || !nTermMod.Value.TargetResidue.HasValue) &&
            sequence.BaseSequence.Length > 0)
        {
            // Set the target residue to the first character of the sequence
            var firstResidue = sequence.BaseSequence[0];
            var updatedNTermMod = nTermMod.Value with { TargetResidue = firstResidue };
            
            // Create a new sequence with the updated N-terminal mod
            var allMods = sequence.Modifications.ToArray();
            for (int i = 0; i < allMods.Length; i++)
            {
                if (allMods[i].PositionType == ModificationPositionType.NTerminus)
                {
                    allMods[i] = updatedNTermMod;
                    break;
                }
            }
            sequence = sequence.WithModifications(allMods);
        }

        // Now call the base implementation to do the actual resolution
        return base.EnrichModificationsIfNeeded(sequence);
    }

    public override bool ShouldResolveMod(CanonicalModification mod)
    {
        var resolved = mod.MzLibModification;
        if (resolved != null)
        {
            // Already resolved - need to re-resolve only if it's NOT already a UniProt mod
            return !string.Equals(resolved.ModificationType, "UniProt", StringComparison.OrdinalIgnoreCase);
        }

        // Not yet resolved - need to resolve via lookup UNLESS it already has inline UniProt representation
        // (i.e., the OriginalRepresentation already contains "UniProt:" which can be used directly)
        return mod.OriginalRepresentation?.IndexOf("UniProt", StringComparison.OrdinalIgnoreCase) < 0;
    }

    protected override string? GetModificationString(CanonicalModification mod, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
    {
        bool writeType = Schema is UniProtSequenceSchema { WriteModType: true };
        bool writeMotif = Schema is UniProtSequenceSchema { WriteMotifs: true };
        var resolved = mod.MzLibModification;
        if (resolved != null &&
            !string.IsNullOrWhiteSpace(resolved.IdWithMotif) &&
            string.Equals(resolved.ModificationType, "UniProt", StringComparison.OrdinalIgnoreCase))
        {
            return writeType switch
            {
                false when !writeMotif => resolved.OriginalId,
                false when writeMotif => resolved.IdWithMotif,
                true when !writeMotif => resolved.ModificationType + ":" + resolved.OriginalId,
                true when writeMotif => resolved.ModificationType + ":" + resolved.IdWithMotif
            };
        }

        var inline = TryGetInlineUniProtRepresentation(mod, writeMotif, writeType);
        if (inline != null)
        {
            return inline;
        }

        warnings.AddIncompatibleItem(mod.ToString());

        if (mode == SequenceConversionHandlingMode.RemoveIncompatibleElements)
        {
            warnings.AddWarning($"Removing incompatible modification without UniProt mapping: {mod}");
            return null;
        }

        if (mode == SequenceConversionHandlingMode.ThrowException)
        {
            throw new SequenceConversionException(
                $"Cannot serialize modification in UniProt format - mapping unavailable: {mod}",
                ConversionFailureReason.IncompatibleModifications,
                new[] { mod.ToString() });
        }

        return null;
    }

    private string? TryGetInlineUniProtRepresentation(CanonicalModification mod, bool writeMotif, bool writeType)
    {
        if (string.IsNullOrWhiteSpace(mod.OriginalRepresentation))
        {
            return null;
        }

        var representation = mod.OriginalRepresentation.Trim().Trim('[', ']');
        if (representation.IndexOf("UniProt", StringComparison.OrdinalIgnoreCase) < 0)
        {
            return null;
        }

        var colonIndex = representation.IndexOf(':');
        var payload = colonIndex >= 0
            ? representation[(colonIndex + 1)..].Trim()
            : representation;

        var motifIndex = payload.IndexOf("on", StringComparison.Ordinal);
        if (motifIndex >= 0)
        {
            var motifPart = payload[(motifIndex + 2)..].Trim();
            payload = payload[..motifIndex].Trim();
            if (writeMotif && !string.IsNullOrEmpty(motifPart))
            {
                payload += $" on {motifPart}";
            }
        }

        if (writeType && colonIndex >= 0)
        {
            var typePart = representation[..colonIndex].Trim();
            if (!string.IsNullOrEmpty(typePart))
            {
                payload = $"{typePart}:{payload}";
            }
        }
        return payload;
    }

    /// <summary>
    /// Override serialization to handle UniProt-specific rules:
    /// - Some mods replace the residue letter
    /// - Some mods (like Carbamidomethyl) are completely suppressed
    /// </summary>
    protected override string? SerializeInternal(
        CanonicalSequence sequence,
        ConversionWarnings warnings,
        SequenceConversionHandlingMode mode)
    {
        var sb = new StringBuilder();

        // Handle N-terminal modification
        var nTermMod = sequence.NTerminalModification;
        if (nTermMod.HasValue)
        {
            if (!ShouldSuppressMod(nTermMod.Value))
            {
                var modString = GetModificationString(nTermMod.Value, warnings, mode);
                if (modString == null && mode == SequenceConversionHandlingMode.ReturnNull)
                    return null;

                if (modString != null)
                {
                    sb.Append(Schema.ModOpenBracket);
                    sb.Append(modString);
                    sb.Append(Schema.ModCloseBracket);

                    // Add N-terminal separator if defined and not empty
                    if (!string.IsNullOrEmpty(Schema.NTermSeparator))
                    {
                        sb.Append(Schema.NTermSeparator);
                    }
                }
            }
        }

        // Handle residue modifications - build a lookup for quick access
        var residueMods = new Dictionary<int, CanonicalModification>();
        foreach (var mod in sequence.ResidueModifications)
        {
            if (mod.ResidueIndex.HasValue)
            {
                residueMods[mod.ResidueIndex.Value] = mod;
            }
        }

        // Write sequence with modifications
        for (int i = 0; i < sequence.BaseSequence.Length; i++)
        {
            var residue = sequence.BaseSequence[i];

            if (residueMods.TryGetValue(i, out var mod))
            {
                // Check if mod should be suppressed entirely (e.g., Carbamidomethyl)
                if (ShouldSuppressMod(mod))
                {
                    // Skip both residue and modification - completely remove from output
                    continue;
                }

                var modString = GetModificationString(mod, warnings, mode);
                if (modString == null && mode == SequenceConversionHandlingMode.ReturnNull)
                    return null;

                if (modString != null)
                {
                    // Check if this mod replaces the residue letter
                    if (!ReplacesResidue(modString))
                    {
                        // Normal case: output residue then mod
                        sb.Append(residue);
                    }
                    // else: Residue-replacing mod - don't output the residue letter

                    sb.Append(Schema.ModOpenBracket);
                    sb.Append(modString);
                    sb.Append(Schema.ModCloseBracket);
                }
                else
                {
                    // No mod string but not suppressed - just output residue
                    sb.Append(residue);
                }
            }
            else
            {
                sb.Append(residue);
            }
        }

        // Handle C-terminal modification
        var cTermMod = sequence.CTerminalModification;
        if (cTermMod.HasValue)
        {
            if (!ShouldSuppressMod(cTermMod.Value))
            {
                var modString = GetModificationString(cTermMod.Value, warnings, mode);
                if (modString == null && mode == SequenceConversionHandlingMode.ReturnNull)
                    return null;

                if (modString != null)
                {
                    // Add C-terminal separator if defined
                    if (!string.IsNullOrEmpty(Schema.CTermSeparator))
                    {
                        sb.Append(Schema.CTermSeparator);
                    }

                    sb.Append(Schema.ModOpenBracket);
                    sb.Append(modString);
                    sb.Append(Schema.ModCloseBracket);
                }
            }
        }

        return sb.ToString();
    }

    /// <summary>
    /// Determines if a modification should be completely suppressed from output.
    /// </summary>
    private static bool ShouldSuppressMod(CanonicalModification mod)
    {
        // Check by modification type (e.g., "Common Fixed")
        var resolved = mod.MzLibModification;
        if (resolved?.ModificationType != null && SuppressedModTypes.Contains(resolved.ModificationType))
        {
            return true;
        }

        // Check original representation for suppressed mod types
        var orig = mod.OriginalRepresentation;
        if (!string.IsNullOrEmpty(orig))
        {
            foreach (var suppressedType in SuppressedModTypes)
            {
                if (orig.Contains(suppressedType, StringComparison.OrdinalIgnoreCase))
                {
                    return true;
                }
            }
        }

        // Check by modification name
        var modName = resolved?.OriginalId ?? resolved?.IdWithMotif ?? ExtractModName(mod.OriginalRepresentation);
        if (!string.IsNullOrEmpty(modName))
        {
            foreach (var suppressed in SuppressedModNames)
            {
                if (modName.Contains(suppressed, StringComparison.OrdinalIgnoreCase))
                {
                    return true;
                }
            }
        }

        return false;
    }

    /// <summary>
    /// Determines if a modification name represents a residue-replacing modification.
    /// </summary>
    private static bool ReplacesResidue(string modString)
    {
        if (string.IsNullOrEmpty(modString))
        {
            return false;
        }

        foreach (var replacingMod in ResidueReplacingMods)
        {
            if (modString.Contains(replacingMod, StringComparison.OrdinalIgnoreCase))
            {
                return true;
            }
        }

        return false;
    }

    /// <summary>
    /// Extracts the modification name from an original representation string.
    /// </summary>
    private static string? ExtractModName(string? representation)
    {
        if (string.IsNullOrEmpty(representation))
        {
            return null;
        }

        // Handle format like "Common Fixed:Carbamidomethyl on C"
        var colonIdx = representation.IndexOf(':');
        if (colonIdx >= 0 && colonIdx < representation.Length - 1)
        {
            var afterColon = representation[(colonIdx + 1)..].Trim();
            var onIdx = afterColon.IndexOf(" on ", StringComparison.OrdinalIgnoreCase);
            if (onIdx > 0)
            {
                return afterColon[..onIdx].Trim();
            }
            return afterColon;
        }

        return representation;
    }
}
