using MassSpectrometry;
using Omics;
using Omics.SpectralMatch;

namespace Omics.BioPolymerGroup;

/// <summary>
/// A group of identified digestion products (peptides or oligonucleotides) sharing one base
/// sequence. The modified forms of that sequence are its members, so modification occupancy is
/// reported against peptide-local positions — position 1 means the first residue of this peptide,
/// not of the protein it came from.
///
/// This is the digestion-product counterpart to <see cref="BioPolymerGroup"/>, not a subclass of
/// it: the two share quantification machinery but almost nothing else. A peptide group has no
/// sequence coverage, no gene or organism, and its identity is a sequence rather than a set of
/// parent accessions.
///
/// Like <see cref="BioPolymerGroup"/>, this class holds data only; rendering it to a results file
/// is the job of <see cref="BioPolymerWithSetModsGroupTsvSchema"/> together with
/// <see cref="TsvWriter"/>.
/// </summary>
public class BioPolymerWithSetModsGroup : IEquatable<BioPolymerWithSetModsGroup>, IHasSampleIntensities
{
    /// <summary>
    /// Creates a group for one base sequence.
    /// </summary>
    /// <param name="baseSequence">The unmodified sequence that identifies this group.</param>
    /// <param name="peptidoforms">The modified forms of that sequence observed in the search.
    /// Every form must carry <paramref name="baseSequence"/>.</param>
    /// <exception cref="ArgumentException">A form's base sequence disagrees with the group's.</exception>
    public BioPolymerWithSetModsGroup(string baseSequence, IEnumerable<IBioPolymerWithSetMods> peptidoforms)
    {
        if (string.IsNullOrEmpty(baseSequence))
            throw new ArgumentException("A digestion product group needs a base sequence.", nameof(baseSequence));

        ArgumentNullException.ThrowIfNull(peptidoforms);

        BaseSequence = baseSequence;

        var supplied = peptidoforms.ToList();

        var mismatched = supplied.FirstOrDefault(p => p.BaseSequence != baseSequence);
        if (mismatched is not null)
        {
            throw new ArgumentException(
                $"All peptidoforms in a group must share its base sequence. Group is '{baseSequence}' " +
                $"but form '{mismatched.FullSequence}' has base sequence '{mismatched.BaseSequence}'.",
                nameof(peptidoforms));
        }

        // Deduplicate explicitly rather than relying on the element type's Equals, which several
        // IBioPolymerWithSetMods implementations define on base sequence alone — that would
        // collapse the very forms this group exists to hold.
        //
        // The key is (form, parent, location). Parent, because a sequence shared across biopolymers
        // has one form per parent with the same FullSequence. Location, because a sequence can occur
        // more than once within a single parent — repeat domains in histones, collagens and mucins
        // are ordinary — and those occurrences are distinct observations that the residue-range
        // columns must report separately.
        Peptidoforms = [.. supplied
            .GroupBy(p => (p.FullSequence, p.Parent?.Accession, p.OneBasedStartResidue, p.OneBasedEndResidue))
            .Select(g => g.First())
            .OrderBy(p => p.FullSequence, StringComparer.Ordinal)
            .ThenBy(p => p.Parent?.Accession, StringComparer.Ordinal)
            .ThenBy(p => p.OneBasedStartResidue)];

        ParentBioPolymers = [.. Peptidoforms.Select(p => p.Parent).Where(p => p is not null)];
        ListOfParentsOrderedByAccession = [.. ParentBioPolymers.OrderBy(p => p.Accession)];

        foreach (var parent in ParentBioPolymers)
        {
            IsDecoy |= parent.IsDecoy;
            IsContaminant |= parent.IsContaminant;
            IsEntrapment |= parent.IsEntrapment;
        }

        AllPsmsBelowOnePercentFDR = [];
    }

    /// <summary>
    /// Buckets PSMs into one group per base sequence, taking each group's forms from the PSMs
    /// assigned to it. PSMs whose sequence was never resolved carry an empty base sequence and are
    /// skipped, since they cannot be attributed to any one sequence.
    /// </summary>
    public static List<BioPolymerWithSetModsGroup> CreateGroups(IEnumerable<ISpectralMatch> psms)
    {
        var groups = new List<BioPolymerWithSetModsGroup>();

        foreach (var bySequence in psms
                     .Where(p => !string.IsNullOrEmpty(p.BaseSequence))
                     .GroupBy(p => p.BaseSequence))
        {
            var psmsForSequence = bySequence.ToList();

            var forms = psmsForSequence
                .SelectMany(p => p.GetIdentifiedBioPolymersWithSetMods())
                .Where(f => f.BaseSequence == bySequence.Key);

            groups.Add(new BioPolymerWithSetModsGroup(bySequence.Key, forms)
            {
                AllPsmsBelowOnePercentFDR = [.. psmsForSequence]
            });
        }

        return groups;
    }

    #region Identity

    /// <summary>
    /// The unmodified sequence shared by every member. Identity key for equality and output.
    /// </summary>
    public string BaseSequence { get; }

    /// <summary>
    /// The modified forms of <see cref="BaseSequence"/> observed for this group, one entry per
    /// distinct (full sequence, parent) pair and ordered for stable output. A sequence found in
    /// several biopolymers therefore appears once per parent; the report deduplicates to distinct
    /// full sequences where that is what a column means.
    /// </summary>
    public IReadOnlyList<IBioPolymerWithSetMods> Peptidoforms { get; }

    /// <summary>
    /// Biopolymers this sequence was found in. Reported for provenance but not part of identity,
    /// so one group covers a sequence shared across several parents.
    /// </summary>
    public HashSet<IBioPolymer> ParentBioPolymers { get; }

    /// <summary>
    /// <see cref="ParentBioPolymers"/> in accession order, for deterministic output.
    /// </summary>
    public List<IBioPolymer> ListOfParentsOrderedByAccession { get; }

    /// <summary>True if any parent is a decoy, used for FDR estimation.</summary>
    public bool IsDecoy { get; }

    /// <summary>True if any parent is marked as a contaminant.</summary>
    public bool IsContaminant { get; }

    /// <summary>True if any parent is marked as an entrapment sequence.</summary>
    public bool IsEntrapment { get; }

    #endregion

    #region Identification statistics

    /// <summary>
    /// PSMs for this sequence that pass the 1% FDR threshold. Setting this invalidates
    /// <see cref="SampleGroupResults"/>.
    /// </summary>
    /// <exception cref="ArgumentException">A PSM's base sequence disagrees with the group's.</exception>
    private HashSet<ISpectralMatch> _allPsmsBelowOnePercentFDR = null!;
    public HashSet<ISpectralMatch> AllPsmsBelowOnePercentFDR
    {
        get => _allPsmsBelowOnePercentFDR;
        set
        {
            // Assigning null would defer the failure to whichever later read dereferences it, so
            // reject it here where the caller's mistake is still visible.
            ArgumentNullException.ThrowIfNull(value);

            // Occupancy is computed in this sequence's coordinate space, so a PSM for a different
            // sequence would silently land on the wrong positions. Reject it at the boundary
            // instead. A PSM whose sequence was never resolved carries an empty base sequence
            // rather than null; those are allowed here and skipped when occupancy is computed.
            var mismatched = value.FirstOrDefault(
                p => !string.IsNullOrEmpty(p.BaseSequence) && p.BaseSequence != BaseSequence);
            if (mismatched is not null)
            {
                throw new ArgumentException(
                    $"All PSMs in a group must share its base sequence. Group is '{BaseSequence}' " +
                    $"but a PSM from '{mismatched.FullFilePath}' has base sequence '{mismatched.BaseSequence}'.",
                    nameof(value));
            }

            // Copy rather than alias. Validating the caller's set and then holding a reference to it
            // leaves the guard bypassable: the caller keeps its own handle and can add a foreign PSM
            // afterwards, which then gets filed under this group's base sequence and reports a
            // modification at a residue that sequence does not have.
            _allPsmsBelowOnePercentFDR = [.. value];
            SampleGroupResults = null;
        }
    }

    /// <summary>The q-value for this group. Lower is more confident (0.01 = 1% FDR).</summary>
    public double QValue { get; set; }

    /// <summary>Best (highest) score among the PSMs in this group.</summary>
    public double BestPsmScore { get; set; }

    /// <summary>Best (lowest) q-value among the PSMs in this group.</summary>
    public double BestPsmQValue { get; set; }

    /// <summary>Cumulative target count at this group's rank, for FDR calculation.</summary>
    public int CumulativeTarget { get; set; }

    /// <summary>Cumulative decoy count at this group's rank, for FDR calculation.</summary>
    public int CumulativeDecoy { get; set; }

    /// <summary>
    /// Sets <see cref="BestPsmScore"/> to the highest score among the group's PSMs, or 0 when there
    /// are none.
    /// </summary>
    public void Score()
    {
        BestPsmScore = AllPsmsBelowOnePercentFDR.Count == 0
            ? 0
            : AllPsmsBelowOnePercentFDR.Max(p => p.Score);
    }

    #endregion

    #region Quantification

    /// <summary>
    /// Samples contributing quantification. Setting this invalidates <see cref="SampleGroupResults"/>.
    /// </summary>
    private List<ISampleInfo>? _samplesForQuantification;
    public List<ISampleInfo>? SamplesForQuantification
    {
        get => _samplesForQuantification;
        set
        {
            _samplesForQuantification = value;
            SampleGroupResults = null;
        }
    }

    /// <summary>
    /// Measured intensities keyed by sample. Setting this invalidates <see cref="SampleGroupResults"/>.
    /// </summary>
    private Dictionary<ISampleInfo, double>? _intensitiesBySample;
    public Dictionary<ISampleInfo, double>? IntensitiesBySample
    {
        get => _intensitiesBySample;
        set
        {
            _intensitiesBySample = value;
            SampleGroupResults = null;
        }
    }

    /// <summary>
    /// Per-sample-group quantification and occupancy, built by <see cref="PopulateSampleGroupResults"/>.
    /// </summary>
    public List<SampleGroupResult>? SampleGroupResults { get; set; }

    /// <summary>
    /// Buckets this group's PSMs and intensities by experimental design and computes modification
    /// occupancy for each bucket in peptide-local coordinates.
    /// </summary>
    public void PopulateSampleGroupResults()
    {
        SampleGroupResults = SampleGroupBuilder.Build(
            SamplesForQuantification,
            IntensitiesBySample,
            AllPsmsBelowOnePercentFDR,
            PopulateOccupancy);
    }

    /// <summary>
    /// Attaches occupancy for one sample group, keyed by this group's base sequence. Positions use
    /// the AllModsOneIsNterminus convention relative to the peptide, not the parent biopolymer.
    /// </summary>
    private void PopulateOccupancy(SampleGroupResult result, List<ISpectralMatch> psms)
    {
        if (psms.Count == 0)
            return;

        var occupancy = ModificationOccupancyCalculator.CalculateDigestionProductLevelOccupancy(psms);

        if (occupancy.Count > 0)
            result.DigestionProductOccupancy[BaseSequence] = occupancy;
    }

    #endregion

    #region Equality

    /// <summary>Two groups are equal when they describe the same base sequence.</summary>
    public bool Equals(BioPolymerWithSetModsGroup? other)
    {
        if (other is null) return false;
        if (ReferenceEquals(this, other)) return true;
        return BaseSequence == other.BaseSequence;
    }

    public override bool Equals(object? obj) => obj is BioPolymerWithSetModsGroup other && Equals(other);

    public override int GetHashCode() => BaseSequence.GetHashCode();

    #endregion
}
