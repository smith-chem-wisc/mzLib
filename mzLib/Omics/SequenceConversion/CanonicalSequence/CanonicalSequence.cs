using Easy.Common.Extensions;
using System.Collections.Immutable;
using System.Linq;
using System.Text;

namespace Omics.SequenceConversion;

/// <summary>
/// Immutable representation of a sequence with modifications in the canonical format.
/// This serves as the universal intermediate representation during sequence conversions.
/// 
/// A CanonicalSequence can be created from any supported input format (mzLib, UNIMOD, 
/// mass-shift, etc.) and serialized to any supported output format. The modifications
/// are stored in a format-agnostic way, with lazy enrichment available to populate
/// additional fields needed for specific output formats.
/// </summary>
/// <param name="BaseSequence">The unmodified amino acid or nucleotide sequence.</param>
/// <param name="Modifications">Immutable list of all modifications on this sequence,
/// including terminal modifications. Ordered by position (N-term first, then residues
/// in order, then C-term).</param>
/// <param name="SourceFormat">Identifier for the format this sequence was parsed from.
/// Used for debugging and to inform serialization decisions.</param>
public readonly record struct CanonicalSequence(
    string BaseSequence,
    ImmutableArray<CanonicalModification> Modifications,
    string SourceFormat)
{
    /// <summary>
    /// Creates an empty canonical sequence (no base sequence, no modifications).
    /// </summary>
    public static readonly CanonicalSequence Empty = new(string.Empty, ImmutableArray<CanonicalModification>.Empty, "Empty");

    /// <summary>
    /// Creates a CanonicalSequence with no modifications.
    /// </summary>
    public static CanonicalSequence Unmodified(string baseSequence, string sourceFormat = "Unknown")
    {
        return new CanonicalSequence(baseSequence, ImmutableArray<CanonicalModification>.Empty, sourceFormat);
    }

    /// <summary>
    /// Gets the N-terminal modification, if present.
    /// </summary>
    public CanonicalModification? NTerminalModification
    {
        get
        {
            var mod = Modifications.FirstOrDefault(m => m.PositionType == ModificationPositionType.NTerminus);

            if (mod.IsDefault())
                return null;
            return mod;
        }
    }

    /// <summary>
    /// Gets the C-terminal modification, if present.
    /// </summary>
    public CanonicalModification? CTerminalModification
    {
        get
        {
            var mod = Modifications.FirstOrDefault(m => m.PositionType == ModificationPositionType.CTerminus);
            if (mod.IsDefault())
                return null;
            return mod;
        }
    }

    /// <summary>
    /// Gets all residue modifications (excludes terminal modifications).
    /// </summary>
    public IEnumerable<CanonicalModification> ResidueModifications =>
        Modifications.Where(m => m.PositionType == ModificationPositionType.Residue);

    /// <summary>
    /// Returns true if this sequence has any modifications.
    /// </summary>
    public bool HasModifications => Modifications.Length > 0;

    /// <summary>
    /// Returns true if this sequence has any terminal modifications.
    /// </summary>
    public bool HasTerminalModifications =>
        NTerminalModification.HasValue || CTerminalModification.HasValue;

    /// <summary>
    /// Gets the number of modifications on this sequence.
    /// </summary>
    public int ModificationCount => Modifications.Length;

    /// <summary>
    /// Gets the length of the base sequence.
    /// </summary>
    public int Length => BaseSequence.Length;

    /// <summary>
    /// Gets the modification at a specific zero-based residue index, if present.
    /// Returns null if no modification exists at that position.
    /// </summary>
    public CanonicalModification? GetModificationAt(int residueIndex)
    {
        foreach (var mod in Modifications)
        {
            if (mod.PositionType == ModificationPositionType.Residue && mod.ResidueIndex == residueIndex)
                return mod;
        }
        return null;
    }

    /// <summary>
    /// Returns true if there is a modification at the specified zero-based residue index.
    /// </summary>
    public bool HasModificationAt(int residueIndex) => GetModificationAt(residueIndex).HasValue;

    /// <summary>
    /// Calculates the total mass shift from all modifications.
    /// Returns null if any modification lacks mass information.
    /// </summary>
    public double? TotalModificationMass
    {
        get
        {
            if (Modifications.Length == 0)
                return 0;

            double total = 0;
            foreach (var mod in Modifications)
            {
                var mass = mod.EffectiveMass;
                if (!mass.HasValue)
                    return null;
                total += mass.Value;
            }
            return total;
        }
    }

    /// <summary>
    /// Returns true if all modifications have mass information available.
    /// </summary>
    public bool AllModificationsHaveMass => Modifications.All(m => m.HasMass);

    /// <summary>
    /// Returns true if all modifications have been resolved to mzLib Modifications.
    /// </summary>
    public bool AllModificationsResolved => Modifications.All(m => m.IsResolved);

    /// <summary>
    /// Creates a new CanonicalSequence with an additional modification.
    /// The modification is inserted in the correct position order.
    /// </summary>
    public CanonicalSequence WithModification(CanonicalModification modification)
    {
        var newMods = Modifications.Add(modification);
        // Sort by position type then by residue index
        newMods = newMods.Sort((a, b) =>
        {
            if (a.PositionType != b.PositionType)
            {
                return a.PositionType.CompareTo(b.PositionType);
            }
            return (a.ResidueIndex ?? -1).CompareTo(b.ResidueIndex ?? -1);
        });
        return this with { Modifications = newMods };
    }

    /// <summary>
    /// Creates a new CanonicalSequence with the specified modifications replacing all existing.
    /// </summary>
    public CanonicalSequence WithModifications(IEnumerable<CanonicalModification> modifications)
    {
        var newMods = modifications.ToImmutableArray();
        newMods = newMods.Sort((a, b) =>
        {
            if (a.PositionType != b.PositionType)
            {
                return a.PositionType.CompareTo(b.PositionType);
            }
            return (a.ResidueIndex ?? -1).CompareTo(b.ResidueIndex ?? -1);
        });
        return this with { Modifications = newMods };
    }

    /// <summary>
    /// Creates a new CanonicalSequence with the specified source format.
    /// </summary>
    public CanonicalSequence WithSourceFormat(string sourceFormat)
    {
        return this with { SourceFormat = sourceFormat };
    }

    /// <summary>
    /// Returns a debug-friendly string representation.
    /// Format: "BaseSequence [mod1, mod2, ...] (SourceFormat)"
    /// </summary>
    public override string ToString()
    {
        if (Modifications.Length == 0)
            return $"{BaseSequence} (unmodified, {SourceFormat})";

        var modStrings = Modifications.Select(m => m.ToString());
        return $"{BaseSequence} [{string.Join(", ", modStrings)}] ({SourceFormat})";
    }

    /// <summary>
    /// Returns a simple representation showing the base sequence with modification
    /// positions marked. Example: "PEP[*]TIDE" where * indicates a modification.
    /// </summary>
    public string ToMarkedSequence()
    {
        if (Modifications.Length == 0)
            return BaseSequence;

        var sb = new StringBuilder();

        // N-terminal mod
        if (NTerminalModification.HasValue)
            sb.Append("[*]-");

        // Residues with modifications marked
        for (int i = 0; i < BaseSequence.Length; i++)
        {
            sb.Append(BaseSequence[i]);
            if (HasModificationAt(i))
                sb.Append("[*]");
        }

        // C-terminal mod
        if (CTerminalModification.HasValue)
            sb.Append("-[*]");

        return sb.ToString();
    }
}
