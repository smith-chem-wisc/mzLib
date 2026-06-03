using System.Collections.Immutable;
using Chemistry;
using Omics.Modifications;

namespace Omics.SequenceConversion;

/// <summary>
/// Mutable builder for constructing <see cref="CanonicalSequence"/> instances.
/// This class provides a fluent API for building sequences incrementally,
/// which is useful when parsing sequences from various formats.
/// 
/// Thread Safety: This class is NOT thread-safe. Use separate builders for concurrent parsing.
/// </summary>
public class CanonicalSequenceBuilder
{
    private string _baseSequence = string.Empty;
    private readonly List<CanonicalModification> _modifications = new();
    private string _sourceFormat = "Unknown";

    /// <summary>
    /// Creates a new empty builder.
    /// </summary>
    public CanonicalSequenceBuilder()
    {
    }

    /// <summary>
    /// Creates a builder initialized with the specified base sequence.
    /// </summary>
    public CanonicalSequenceBuilder(string baseSequence)
    {
        _baseSequence = baseSequence ?? string.Empty;
    }

    /// <summary>
    /// Creates a builder initialized from an existing CanonicalSequence.
    /// Useful for creating modified copies of sequences.
    /// </summary>
    public CanonicalSequenceBuilder(CanonicalSequence sequence)
    {
        _baseSequence = sequence.BaseSequence;
        _modifications.AddRange(sequence.Modifications);
        _sourceFormat = sequence.SourceFormat;
    }

    /// <summary>
    /// Gets or sets the base sequence.
    /// </summary>
    public string BaseSequence
    {
        get => _baseSequence;
        set => _baseSequence = value ?? string.Empty;
    }

    /// <summary>
    /// Gets or sets the source format identifier.
    /// </summary>
    public string SourceFormat
    {
        get => _sourceFormat;
        set => _sourceFormat = value ?? "Unknown";
    }

    /// <summary>
    /// Gets the current list of modifications (read-only view).
    /// </summary>
    public IReadOnlyList<CanonicalModification> Modifications => _modifications;

    /// <summary>
    /// Gets the current modification count.
    /// </summary>
    public int ModificationCount => _modifications.Count;

    /// <summary>
    /// Sets the base sequence and returns this builder for chaining.
    /// </summary>
    public CanonicalSequenceBuilder WithBaseSequence(string baseSequence)
    {
        _baseSequence = baseSequence ?? string.Empty;
        return this;
    }

    /// <summary>
    /// Sets the source format and returns this builder for chaining.
    /// </summary>
    public CanonicalSequenceBuilder WithSourceFormat(string sourceFormat)
    {
        _sourceFormat = sourceFormat ?? "Unknown";
        return this;
    }

    /// <summary>
    /// Adds a modification and returns this builder for chaining.
    /// </summary>
    public CanonicalSequenceBuilder AddModification(CanonicalModification modification)
    {
        _modifications.Add(modification);
        return this;
    }

    /// <summary>
    /// Adds multiple modifications and returns this builder for chaining.
    /// </summary>
    public CanonicalSequenceBuilder AddModifications(IEnumerable<CanonicalModification> modifications)
    {
        _modifications.AddRange(modifications);
        return this;
    }

    /// <summary>
    /// Adds a residue modification at the specified zero-based index.
    /// </summary>
    public CanonicalSequenceBuilder AddResidueModification(
        int residueIndex,
        string originalRepresentation,
        double? mass = null,
        ChemicalFormula? formula = null,
        int? unimodId = null,
        string? mzLibId = null,
        Modification? mzLibModification = null)
    {
        if (residueIndex < 0 || residueIndex >= _baseSequence.Length)
        {
            throw new ArgumentOutOfRangeException(nameof(residueIndex),
                $"Residue index {residueIndex} is out of range for sequence of length {_baseSequence.Length}");
        }

        var targetResidue = _baseSequence[residueIndex];
        var mod = CanonicalModification.AtResidue(
            residueIndex, targetResidue, originalRepresentation,
            mass, formula, unimodId, mzLibId, mzLibModification);
        _modifications.Add(mod);
        return this;
    }

    /// <summary>
    /// Adds an N-terminal modification.
    /// </summary>
    public CanonicalSequenceBuilder AddNTerminalModification(
        string originalRepresentation,
        double? mass = null,
        ChemicalFormula? formula = null,
        int? unimodId = null,
        string? mzLibId = null,
        Modification? mzLibModification = null)
    {
        char? targetResidue = _baseSequence.Length > 0 ? _baseSequence[0] : null;
        var mod = CanonicalModification.AtNTerminus(
            originalRepresentation, targetResidue,
            mass, formula, unimodId, mzLibId, mzLibModification);
        _modifications.Add(mod);
        return this;
    }

    /// <summary>
    /// Adds a C-terminal modification.
    /// </summary>
    public CanonicalSequenceBuilder AddCTerminalModification(
        string originalRepresentation,
        double? mass = null,
        ChemicalFormula? formula = null,
        int? unimodId = null,
        string? mzLibId = null,
        Modification? mzLibModification = null)
    {
        char? targetResidue = _baseSequence.Length > 0 ? _baseSequence[^1] : null;
        var mod = CanonicalModification.AtCTerminus(
            originalRepresentation, targetResidue,
            mass, formula, unimodId, mzLibId, mzLibModification);
        _modifications.Add(mod);
        return this;
    }

    /// <summary>
    /// Clears all modifications from the builder.
    /// </summary>
    public CanonicalSequenceBuilder ClearModifications()
    {
        _modifications.Clear();
        return this;
    }

    /// <summary>
    /// Removes a modification at the specified index in the modifications list.
    /// </summary>
    public CanonicalSequenceBuilder RemoveModificationAt(int index)
    {
        if (index >= 0 && index < _modifications.Count)
        {
            _modifications.RemoveAt(index);
        }
        return this;
    }

    /// <summary>
    /// Removes all modifications at a specific residue index.
    /// </summary>
    public CanonicalSequenceBuilder RemoveModificationsAtResidue(int residueIndex)
    {
        _modifications.RemoveAll(m =>
            m.PositionType == ModificationPositionType.Residue &&
            m.ResidueIndex == residueIndex);
        return this;
    }

    /// <summary>
    /// Resets the builder to its initial empty state.
    /// </summary>
    public CanonicalSequenceBuilder Reset()
    {
        _baseSequence = string.Empty;
        _modifications.Clear();
        _sourceFormat = "Unknown";
        return this;
    }

    /// <summary>
    /// Builds the immutable <see cref="CanonicalSequence"/>.
    /// Modifications are sorted by position (N-term, residues in order, C-term).
    /// </summary>
    public CanonicalSequence Build()
    {
        // Sort modifications by position type, then by residue index
        var sortedMods = _modifications
            .OrderBy(m => m.PositionType)
            .ThenBy(m => m.ResidueIndex ?? -1)
            .ToImmutableArray();

        return new CanonicalSequence(_baseSequence, sortedMods, _sourceFormat);
    }

    /// <summary>
    /// Validates the current state and returns any errors found.
    /// </summary>
    public IReadOnlyList<string> Validate()
    {
        var errors = new List<string>();

        if (string.IsNullOrEmpty(_baseSequence))
        {
            errors.Add("Base sequence is empty.");
        }

        foreach (var mod in _modifications)
        {
            if (mod.PositionType == ModificationPositionType.Residue)
            {
                if (!mod.ResidueIndex.HasValue)
                {
                    errors.Add($"Residue modification '{mod.OriginalRepresentation}' has no residue index.");
                }
                else if (mod.ResidueIndex < 0 || mod.ResidueIndex >= _baseSequence.Length)
                {
                    errors.Add($"Residue modification '{mod.OriginalRepresentation}' has invalid index {mod.ResidueIndex} for sequence of length {_baseSequence.Length}.");
                }
                else if (mod.TargetResidue.HasValue && _baseSequence[mod.ResidueIndex.Value] != mod.TargetResidue.Value)
                {
                    errors.Add($"Residue modification '{mod.OriginalRepresentation}' targets '{mod.TargetResidue}' but found '{_baseSequence[mod.ResidueIndex.Value]}' at index {mod.ResidueIndex}.");
                }
            }
        }

        // Check for duplicate positions
        var residueMods = _modifications.Where(m => m.PositionType == ModificationPositionType.Residue && m.ResidueIndex.HasValue);
        var duplicates = residueMods.GroupBy(m => m.ResidueIndex!.Value).Where(g => g.Count() > 1);
        foreach (var dup in duplicates)
        {
            errors.Add($"Multiple modifications at residue index {dup.Key}.");
        }

        var nTermCount = _modifications.Count(m => m.PositionType == ModificationPositionType.NTerminus);
        if (nTermCount > 1)
        {
            errors.Add($"Multiple N-terminal modifications ({nTermCount} found).");
        }

        var cTermCount = _modifications.Count(m => m.PositionType == ModificationPositionType.CTerminus);
        if (cTermCount > 1)
        {
            errors.Add($"Multiple C-terminal modifications ({cTermCount} found).");
        }

        return errors;
    }

    /// <summary>
    /// Returns true if the current state is valid.
    /// </summary>
    public bool IsValid => Validate().Count == 0;
}
