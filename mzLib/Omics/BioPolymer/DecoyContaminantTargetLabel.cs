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
    /// <remarks>Reads <see cref="Parents"/>, so a value that is not a label is no parent at all.</remarks>
    public static bool IsDecoy(string? label) => AnyParentHas(label, 'D');

    /// <summary>True when any parent named by the label is a contaminant.</summary>
    /// <remarks>Reads <see cref="Parents"/>, so a value that is not a label is no parent at all.</remarks>
    public static bool IsContaminant(string? label) => AnyParentHas(label, 'C');

    /// <summary>True when any parent named by the label is entrapment.</summary>
    /// <remarks>
    /// This answers "does any parent belong to entrapment", not "how much of this PSM is an
    /// entrapment discovery". To count entrapment discoveries, for an FDP estimate, use
    /// <see cref="EntrapmentFraction"/>: a shared <c>T|ET</c> is half a discovery when the two
    /// sequences differ, and none when both proteins carry the same peptide. Reads
    /// <see cref="Parents"/>, so the "E" in "Output too long for Excel" is not entrapment.
    /// </remarks>
    public static bool IsEntrapment(string? label) => AnyParentHas(label, 'E');

    /// <summary>
    /// The label of each parent, in the order the writer joined them; empty when the value is not
    /// a label at all (a blank cell, or MetaMorpheus's "Output too long for Excel").
    /// </summary>
    /// <remarks>A label collapsed to one value (every candidate alike) names a single parent.</remarks>
    public static IReadOnlyList<string> Parents(string? label)
    {
        if (string.IsNullOrWhiteSpace(label))
        {
            return Array.Empty<string>();
        }

        string[] parents = label.Split('|').Select(p => p.Trim()).ToArray();
        return parents.All(IsKnown) ? parents : Array.Empty<string>();
    }

    /// <summary>The share of the PSM that FDR counts as a decoy: decoy parents over all parents.</summary>
    /// <remarks>
    /// MetaMorpheus adds each candidate of a shared PSM as 1/n of a target or a decoy by its own
    /// parent (FdrAnalysisEngine.CalculateQValue), so <c>T|T|D</c> is a third of a decoy. Use
    /// <see cref="IsDecoy"/> for the question MetaMorpheus's SpectralMatch.IsDecoy answers: is any
    /// parent a decoy.
    /// </remarks>
    public static double DecoyFraction(string? label)
    {
        IReadOnlyList<string> parents = Parents(label);
        return parents.Count == 0 ? 0 : parents.Count(p => p.Contains('D')) / (double)parents.Count;
    }

    /// <summary>
    /// The share of the PSM that is an entrapment discovery: entrapment-target candidates whose
    /// full sequence no target or contaminant candidate also carries, over all candidates.
    /// </summary>
    /// <param name="label">The Decoy/Contaminant/Target value.</param>
    /// <param name="fullSequence">The Full Sequence value, '|'-joined in the same order.</param>
    /// <returns>
    /// NaN when the two columns name different numbers of candidates and neither is collapsed, since
    /// they cannot then be lined up.
    /// </returns>
    /// <remarks>
    /// <para>MetaMorpheus searches entrapment as a target, so a peptide present in both a target
    /// and its entrapment partner is written <c>T|ET</c> with a single full sequence. That is the
    /// real peptide, not an ambiguity, and counts as target. Only when the sequences differ is the
    /// PSM truly shared, and then the entrapment candidates count as their share, as decoys do for
    /// FDR (<see cref="DecoyFraction"/>).</para>
    /// <para>Either column collapses to one value when every candidate agrees, and then applies to
    /// all of them. A decoy candidate is never entrapment (<c>ED</c> is a decoy) and claims no
    /// sequence.</para>
    /// </remarks>
    public static double EntrapmentFraction(string? label, string? fullSequence)
    {
        IReadOnlyList<string> parents = Parents(label);
        string[] sequences = (fullSequence ?? string.Empty).Split('|');
        int candidates = Math.Max(parents.Count, sequences.Length);
        if (parents.Count == 0)
        {
            return 0;
        }
        if ((parents.Count != 1 && parents.Count != candidates) || (sequences.Length != 1 && sequences.Length != candidates))
        {
            return double.NaN;
        }

        string ParentOf(int i) => parents.Count == 1 ? parents[0] : parents[i];
        string SequenceOf(int i) => sequences.Length == 1 ? sequences[0] : sequences[i];

        var claimedByRealProteins = new HashSet<string>(Enumerable.Range(0, candidates)
            .Where(i => ParentOf(i) is Target or Contaminant)
            .Select(SequenceOf));
        int entrapment = Enumerable.Range(0, candidates)
            .Count(i => ParentOf(i) == EntrapmentTarget && !claimedByRealProteins.Contains(SequenceOf(i)));
        return entrapment / (double)candidates;
    }

    private static bool AnyParentHas(string? label, char letter) =>
        Parents(label).Any(p => p.Contains(letter));

    private static bool IsKnown(string part) =>
        part is Target or Decoy or Contaminant or EntrapmentTarget or EntrapmentDecoy;
}
