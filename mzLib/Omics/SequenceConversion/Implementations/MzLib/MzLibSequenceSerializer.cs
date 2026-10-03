using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Serializes a <see cref="CanonicalSequence"/> into mzLib format strings.
/// 
/// Output format examples:
/// - Simple: "PEPTIDE"
/// - With residue modification: "PEP[Oxidation on M]TIDE"
/// - With N-terminal modification: "[Acetyl]PEPTIDE" (no separator)
/// - With C-terminal modification: "PEPTIDE-[Amidated]"
/// 
/// Note: For mass shift notation output, use a separate MassShiftSequenceSerializer instead.
///
/// Readability: a modification resolved from text (a name, UNIMOD id or mass) is written only under a name that
/// mzLib's dictionaries (Mods.AllKnownProteinModsDictionary, Mods.AllKnownRnaModsDictionary) read back as that
/// modification, and where PeptideWithSetModifications or OligoWithSetMods reads it at the same position (a C-terminal
/// one on the last residue reads back at the C-terminus). A written name is held to the same position rule. One it
/// can't write that way fails per the handling mode. A modification that arrives with its own Modification object whose
/// name mzLib's dictionaries don't hold (a search engine's glycan) is written by that name, and reading the output
/// needs a dictionary that has it.
/// </summary>
public class MzLibSequenceSerializer : SequenceSerializerBase
{
    /// <summary>
    /// Singleton instance.
    /// </summary>
    public static MzLibSequenceSerializer Instance { get; } = new();

    /// <summary>
    /// Creates a new MzLibSequenceSerializer.
    /// </summary>
    /// <param name="lookup">Optional modification lookup to resolve modifications.</param>
    public MzLibSequenceSerializer(IModificationLookup? lookup = null)
        : base(lookup ?? GlobalModificationLookup.Instance)
    {
    }

    /// <inheritdoc />
    public override string FormatName => MzLibSequenceFormatSchema.Instance.FormatName;

    /// <inheritdoc />
    public override SequenceFormatSchema Schema => MzLibSequenceFormatSchema.Instance;

    /// <inheritdoc />
    public override bool CanSerialize(CanonicalSequence sequence)
    {
        return !string.IsNullOrEmpty(sequence.BaseSequence);
    }

    /// <inheritdoc />
    public override bool ShouldResolveMod(CanonicalModification mod)
    {
        return !TryGetStrictMzLibToken(mod, out _);
    }

    /// <summary>
    /// Gets the string representation of a modification for serialization: "Type:IdWithMotif" when the name reads
    /// back as the modification (or, for one carrying its own Modification object, when mzLib's dictionaries don't
    /// hold the name at all), otherwise handled per <paramref name="mode"/>.
    /// </summary>
    protected override string? GetModificationString(CanonicalModification mod, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
    {
        return TryGetStrictMzLibToken(mod, out var token) ? token : RejectModification(mod, warnings, mode);
    }

    /// <inheritdoc />
    protected override string? SerializeInternal(CanonicalSequence sequence, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
    {
        var lastResidueIndex = sequence.BaseSequence.Length - 1;
        var misplaced = sequence.Modifications
            .Where(m => TryGetReadBackEntry(m, out var entry) && !IsReadBackWhereWritten(entry, m, lastResidueIndex))
            .ToList();
        if (misplaced.Count == 0)
        {
            return base.SerializeInternal(sequence, warnings, mode);
        }

        foreach (var mod in misplaced)
        {
            RejectModification(mod, warnings, mode);
        }

        return mode == SequenceConversionHandlingMode.ReturnNull
            ? null
            : base.SerializeInternal(sequence.WithModifications(sequence.Modifications.Where(m => !misplaced.Contains(m))), warnings, mode);
    }

    private static string? RejectModification(CanonicalModification mod, ConversionWarnings warnings, SequenceConversionHandlingMode mode)
    {
        // Cannot serialize this modification in mzLib format
        warnings.AddIncompatibleItem(mod.ToString());
        
        if (mode == SequenceConversionHandlingMode.RemoveIncompatibleElements)
        {
            warnings.AddWarning($"Removing incompatible modification: {mod}");
            return null; // Signal to skip this modification
        }

        if (mode == SequenceConversionHandlingMode.ThrowException)
        {
            throw new SequenceConversionException(
                $"Cannot serialize modification in mzLib format: {mod}. Consider using MassShiftSequenceSerializer instead.",
                ConversionFailureReason.IncompatibleModifications,
                new[] { mod.ToString() });
        }

        return null;
    }

    private static bool TryGetStrictMzLibToken(CanonicalModification mod, out string token)
    {
        if (mod.MzLibModification != null)
        {
            var modificationType = mod.MzLibModification.ModificationType;
            var idWithMotif = mod.MzLibModification.IdWithMotif;

            if (!string.IsNullOrWhiteSpace(modificationType) && !string.IsNullOrWhiteSpace(idWithMotif) &&
                ReadsBackAs(idWithMotif, mod.MzLibModification))
            {
                token = $"{modificationType}:{idWithMotif}";
                return true;
            }

            // MzLibId and OriginalRepresentation name this same modification, so they can't stand in for it.
            token = string.Empty;
            return false;
        }

        if (TryNormalizeStrictToken(mod.MzLibId, out token))
        {
            return true;
        }

        if (TryNormalizeStrictToken(mod.OriginalRepresentation, out token))
        {
            return true;
        }

        token = string.Empty;
        return false;
    }

    private static bool TryNormalizeStrictToken(string? value, out string token)
    {
        token = string.Empty;

        if (string.IsNullOrWhiteSpace(value))
        {
            return false;
        }

        var trimmed = value.Trim();
        var separatorIndex = trimmed.IndexOf(':');
        if (separatorIndex <= 0 || separatorIndex >= trimmed.Length - 1)
        {
            return false;
        }

        var modificationType = trimmed.Substring(0, separatorIndex).Trim();
        var idWithMotif = trimmed.Substring(separatorIndex + 1).Trim();
        if (string.IsNullOrEmpty(modificationType) || string.IsNullOrEmpty(idWithMotif) ||
            !TryGetKnownModification(idWithMotif, out _))
        {
            return false;
        }

        token = $"{modificationType}:{idWithMotif}";
        return true;
    }

    // Reading an mzLib sequence keeps only the text after the first colon and looks it up as an IdWithMotif.
    // Mods.AddOrUpdateModification adds to these two dictionaries but not to Mods.AllModsKnownDictionary.
    private static bool TryGetKnownModification(string idWithMotif, out Modification modification) =>
        Mods.AllKnownProteinModsDictionary.TryGetValue(idWithMotif, out modification!) ||
        Mods.AllKnownRnaModsDictionary.TryGetValue(idWithMotif, out modification!);

    // The dictionaries keep one entry per IdWithMotif, and some names cover several modifications ("Methyl on X"
    // is both an N- and a C-terminal entry), so a name must read back with the same mass and terminus. A name with
    // no entry at all belongs to a Modification object that the source (or a custom lookup) brought along; it is
    // written by name, and reading it back needs a dictionary that holds it.
    private static bool ReadsBackAs(string idWithMotif, Modification modification)
    {
        if (!TryGetKnownModification(idWithMotif, out var known) || ReferenceEquals(known, modification))
        {
            return true;
        }

        return known.MonoisotopicMass.HasValue && modification.MonoisotopicMass.HasValue &&
               Math.Abs(known.MonoisotopicMass.Value - modification.MonoisotopicMass.Value) <= 1e-5 &&
               TerminusOf(known) == TerminusOf(modification);
    }

    // The modification a reader gets back for this one: its own object, or the dictionary entry its name reads as.
    private static bool TryGetReadBackEntry(CanonicalModification mod, out Modification entry)
    {
        if (mod.MzLibModification != null)
        {
            entry = mod.MzLibModification;
            return true;
        }

        entry = null!;
        return TryGetStrictMzLibToken(mod, out var token) && TryGetKnownModification(token[(token.IndexOf(':') + 1)..], out entry);
    }

    // Reading a sequence moves a C-terminal modification to the C-terminus wherever it is written. That is right only
    // from the last residue, where digestion puts a protease's C-terminal modification. Any other modification
    // written after the last residue is read as being on that residue.
    private static bool IsReadBackWhereWritten(Modification entry, CanonicalModification mod, int lastResidueIndex) =>
        TerminusOf(entry) == ModificationPositionType.CTerminus
            ? mod.PositionType == ModificationPositionType.CTerminus ||
              (mod.PositionType == ModificationPositionType.Residue && mod.ResidueIndex == lastResidueIndex)
            : mod.PositionType != ModificationPositionType.CTerminus;

    private static ModificationPositionType TerminusOf(Modification modification)
    {
        var restriction = modification.LocationRestriction ?? string.Empty;
        if (restriction.Contains("N-terminal", StringComparison.OrdinalIgnoreCase) ||
            restriction.Contains("5'-terminal", StringComparison.OrdinalIgnoreCase))
        {
            return ModificationPositionType.NTerminus;
        }

        if (restriction.Contains("C-terminal", StringComparison.OrdinalIgnoreCase) ||
            restriction.Contains("3'-terminal", StringComparison.OrdinalIgnoreCase))
        {
            return ModificationPositionType.CTerminus;
        }

        return ModificationPositionType.Residue;
    }
}
