using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text.RegularExpressions;
using Omics.Modifications;
using Omics.SequenceConversion;

namespace PredictionClients.Koina.Util;

/// <summary>
/// Centralizes Koina sequence conversion: parse mzLib, enforce basic constraints, filter modifications, and serialize to UNIMOD.
/// </summary>
internal sealed class KoinaSequenceConverter
{
    private static readonly MzLibSequenceParser Parser = new();
    private static readonly Regex BaseStripper = new(@"\[[^\]]+\]", RegexOptions.Compiled);

    internal sealed record Result(string ApiSequence, WarningException? Warning);

    public Result? Convert(
        string input,
        IReadOnlySet<int> allowedUnimodIds,
        UnimodSequenceFormatSchema schema,
        string allowedAminoAcidPattern,
        int minLength,
        int maxLength,
        SequenceConversionHandlingMode mode)
    {
        var rawBase = BaseStripper.Replace(input, string.Empty);
        if (!Regex.IsMatch(rawBase, allowedAminoAcidPattern))
        {
            HandleFailure(mode, "Invalid base sequence.");
            return null;
        }

        CanonicalSequence canonical;
        try
        {
            var parsed = Parser.Parse(input);
            if (!parsed.HasValue)
            {
                HandleFailure(mode, "Failed to parse sequence.");
                return null;
            }
            canonical = parsed.Value;
        }
        catch (SequenceConversionException ex)
        {
            HandleFailure(mode, $"Failed to parse sequence: {ex.Message}");
            return null;
        }

        if (!IsValidBaseSequence(canonical.BaseSequence, allowedAminoAcidPattern, minLength, maxLength))
        {
            HandleFailure(mode, "Invalid base sequence.");
            return null;
        }

        var lookup = CreateLookup(allowedUnimodIds);
        var filtered = FilterModifications(canonical, lookup, mode, out var policyWarning);
        if (filtered == null)
        {
            return null;
        }

        var warnings = new List<string>();
        if (!string.IsNullOrEmpty(policyWarning))
        {
            warnings.Add(policyWarning!);
        }

        var unimodSerializer = new UnimodSequenceSerializer(schema, lookup);
        var unimodWarnings = new ConversionWarnings();
        var unimodOutput = unimodSerializer.Serialize(filtered.Value, unimodWarnings, mode);
        if (unimodOutput == null)
        {
            if (mode == SequenceConversionHandlingMode.ThrowException)
            {
                throw new ArgumentException("Sequence conversion to UNIMOD format failed.");
            }

            warnings.AddRange(unimodWarnings.Warnings);
            return null;
        }

        warnings.AddRange(unimodWarnings.Warnings);
        var warning = warnings.Count > 0 ? new WarningException(string.Join(" ", warnings)) : null;

        return new Result(unimodOutput, warning);
    }

    public bool IsValidBaseSequence(string baseSequence, string allowedPattern, int minLength, int maxLength)
    {
        return Regex.IsMatch(baseSequence, allowedPattern)
               && baseSequence.Length <= maxLength
               && baseSequence.Length >= minLength;
    }

    private static void HandleFailure(SequenceConversionHandlingMode mode, string message)
    {
        if (mode == SequenceConversionHandlingMode.ThrowException)
        {
            throw new ArgumentException(message);
        }
    }

    private static CanonicalSequence? FilterModifications(
        CanonicalSequence canonical,
        IModificationLookup lookup,
        SequenceConversionHandlingMode mode,
        out string? warning)
    {
        warning = null;

        if (!canonical.HasModifications)
        {
            return canonical;
        }

        if (mode == SequenceConversionHandlingMode.UsePrimarySequence)
        {
            warning = "Sequence modifications were removed for prediction.";
            return canonical.WithModifications(Array.Empty<CanonicalModification>());
        }

        var kept = new List<CanonicalModification>(canonical.ModificationCount);
        var rejected = new List<string>();

        foreach (var mod in canonical.Modifications)
        {
            var resolved = lookup.TryResolve(mod);
            if (resolved.HasValue)
            {
                kept.Add(mod with
                {
                    UnimodId = resolved.Value.UnimodId ?? mod.UnimodId,
                    MonoisotopicMass = mod.MonoisotopicMass ?? resolved.Value.MonoisotopicMass,
                    ChemicalFormula = mod.ChemicalFormula ?? resolved.Value.ChemicalFormula,
                    MzLibId = mod.MzLibId ?? resolved.Value.MzLibId,
                    MzLibModification = resolved.Value.MzLibModification ?? mod.MzLibModification,
                    TargetResidue = mod.TargetResidue ?? resolved.Value.TargetResidue
                });
            }
            else
            {
                rejected.Add(mod.ToString());
            }
        }

        if (rejected.Count > 0)
        {
            var message = $"Sequence contains unsupported modifications: {string.Join(", ", rejected)}";
            switch (mode)
            {
                case SequenceConversionHandlingMode.ThrowException:
                    throw new ArgumentException(message);
                case SequenceConversionHandlingMode.ReturnNull:
                    warning = message;
                    return null;
                case SequenceConversionHandlingMode.RemoveIncompatibleElements:
                    warning = message;
                    break;
            }
        }

        return canonical.WithModifications(kept);
    }

    private static IModificationLookup CreateLookup(IReadOnlySet<int> allowedUnimodIds)
    {
        if (allowedUnimodIds.Count == 0)
        {
            return new UnimodModificationLookup(Enumerable.Empty<Modification>());
        }

        var candidates = Mods.UnimodModifications
            .Where(m => TryGetUnimodId(m, out var id) && allowedUnimodIds.Contains(id))
            .ToList();

        return new UnimodModificationLookup(candidates);
    }

    private static bool TryGetUnimodId(Modification modification, out int id)
    {
        if (modification.Accession?.StartsWith("UNIMOD:", StringComparison.OrdinalIgnoreCase) == true
            && int.TryParse(modification.Accession[7..], out id))
        {
            return true;
        }

        if (modification.DatabaseReference != null)
        {
            foreach (var kvp in modification.DatabaseReference)
            {
                if (!kvp.Key.Equals("UNIMOD", StringComparison.OrdinalIgnoreCase))
                {
                    continue;
                }

                if (kvp.Value.Count > 0)
                {
                    var reference = kvp.Value[0]
                        .Replace("UNIMOD:", string.Empty, StringComparison.OrdinalIgnoreCase)
                        .Replace(":", string.Empty);

                    if (int.TryParse(reference, out id))
                    {
                        return true;
                    }
                }
            }
        }

        id = -1;
        return false;
    }
}
