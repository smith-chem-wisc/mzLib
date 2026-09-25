namespace Omics.BioPolymer;

/// <summary>
/// The value of a results file's "Decoy/Contaminant/Target" column, written and read in one place.
/// </summary>
/// <remarks>
/// <para>Entrapment is a target for the search -- it competes with the real targets and is counted
/// among them for FDR -- so the label marks it rather than reclassifying it: <c>ET</c> is an
/// entrapment target and <c>ED</c> a decoy built from one. An <c>ED</c> is a decoy for every
/// purpose.</para>
/// <para>A PSM whose peptide maps to several parents joins their labels with '|' (<c>T|ET</c>), so
/// every reader tests for a letter rather than for equality: <c>== "D"</c> reads <c>ED</c> as a
/// target.</para>
/// </remarks>
public static class DecoyContaminantTargetLabel
{
    public const string Target = "T";
    public const string Decoy = "D";
    public const string Contaminant = "C";
    public const string EntrapmentTarget = "ET";
    public const string EntrapmentDecoy = "ED";

    /// <summary>The label for one biopolymer, or for a group from its any-member flags.</summary>
    /// <remarks>Precedence is ED, ET, D, C, T. The loaders refuse a biopolymer that is both
    /// entrapment and contaminant, so entrapment outranking contaminant loses nothing.</remarks>
    public static string For(bool isDecoy, bool isContaminant, bool isEntrapment)
    {
        if (isEntrapment)
        {
            return isDecoy ? EntrapmentDecoy : EntrapmentTarget;
        }
        if (isDecoy)
        {
            return Decoy;
        }
        return isContaminant ? Contaminant : Target;
    }

    /// <inheritdoc cref="For(bool, bool, bool)"/>
    public static string For(IBioPolymer bioPolymer) =>
        For(bioPolymer.IsDecoy, bioPolymer.IsContaminant, bioPolymer.IsEntrapment);

    /// <summary>True when any parent named by the label is a decoy, entrapment decoys included.</summary>
    public static bool IsDecoy(string? label) => label?.Contains('D') ?? false;

    /// <summary>True when any parent named by the label is entrapment.</summary>
    public static bool IsEntrapment(string? label) => label?.Contains('E') ?? false;
}
