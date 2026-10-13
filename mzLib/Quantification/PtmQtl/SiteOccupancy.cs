using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using MzLibUtil;
using Omics.Modifications;

namespace Quantification.PtmQtl;

/// <summary>
/// One peptidoform identified by MS/MS in one run, mapped to one protein. A peptide shared by several
/// proteins is one row per protein. Plain values only, so any caller (a result file, a stored catalog, a
/// language binding) can supply it.
/// </summary>
/// <param name="Run">Run (raw file) identifier.</param>
/// <param name="FullSequence">
/// The peptidoform in MetaMorpheus full-sequence notation, modifications as <c>[Category:IdWithMotif]</c>
/// after the residue they sit on, a protein N-terminal modification before the first residue.
/// </param>
/// <param name="ProteinAccession">Protein the peptide maps to.</param>
/// <param name="StartResidue">First residue of the peptide in the protein, 1-based.</param>
/// <param name="EndResidue">Last residue of the peptide in the protein, 1-based.</param>
/// <param name="Intensity">
/// The peptidoform's MS/MS-quantified intensity in this run (e.g. FlashLFQ apex, <c>DEF-PEP-INT</c>), or NaN
/// when it was identified by MS/MS in the run but not quantified. Match-between-runs transfers must NOT be
/// supplied: intensity-based occupancy is MS/MS-only (<c>DEF-OCC-INT</c>).
/// </param>
/// <param name="PreviousResidue">
/// The protein residue before the peptide (MetaMorpheus's "Previous Residue"), or null when unknown. A peptide
/// starting at residue 2 is the protein N-terminus only after removal of an initiator Met, so when this is
/// given and is not 'M', an N-terminal modification on such a peptide is not a protein N-terminal site. When
/// null, the caller guarantees that residue 1 is Met for every peptide starting at residue 2.
/// </param>
public sealed record PeptidoformObservation(string Run, string FullSequence, string ProteinAccession,
    int StartResidue, int EndResidue, double Intensity, char? PreviousResidue = null);

/// <summary>What an occupancy value is, for one site in one run.</summary>
public enum OccupancyState
{
    /// <summary>The modified form was quantified: a measurement in (0, 1].</summary>
    Quantified,
    /// <summary>
    /// The modified form was identified but not quantified while forms covering the site were: the
    /// occupancy is low but NOT zero, a floor (<c>DEF-OCC-INT-ZERO</c>). Censored, never averaged in as 0.
    /// </summary>
    Floor,
    /// <summary>The modified form was identified, but nothing covering the site was quantified (<c>DEF-OCC-COUNTONLY</c>).</summary>
    CountOnly,
    /// <summary>
    /// The site was covered by quantified peptides, but no modified form was identified in this run. The
    /// modification may be absent or below detection; occupancy is unknown, not 0 (<c>DEF-OCC-ABSENT</c>).
    /// </summary>
    NotDetected,
}

/// <summary>A modification site on a protein, keyed as the stored catalog keys it.</summary>
/// <param name="ProteinAccession">Protein accession.</param>
/// <param name="Position">1-based residue position; a protein N-terminal modification is at the residue it sits on.</param>
/// <param name="Residue">One-letter residue.</param>
/// <param name="Modification">The modification as the engine names it, <c>Category:IdWithMotif</c>.</param>
public sealed record ModificationSite(string ProteinAccession, int Position, char Residue, string Modification)
{
    /// <summary>Stable text key: <c>accession:residue position:modification</c>.</summary>
    public string Key => $"{ProteinAccession}:{Residue}{Position}:{Modification}";

    /// <summary>The modification's category, the text before the first ':' (e.g. <c>Common Biological</c>).</summary>
    public string Category => Modification.Contains(':') ? Modification[..Modification.IndexOf(':')] : "";
}

/// <summary>Occupancy of one site in one run.</summary>
/// <remarks>
/// A site seen only in modified form in a run has occupancy 1: every quantified form covering it carries the
/// modification. <see cref="UnmodifiedQuantified"/> is false there, so a caller can treat that 1 as a ceiling
/// (the unmodified form may be present below detection), as a floor is treated at the other end. Despite its
/// name, <see cref="SiteOccupancyCalculator.Calculate(IEnumerable{PeptidoformObservation}, Func{string, bool}?)"/> sets <see cref="UnmodifiedQuantified"/> true whenever some quantified form covering the position does not
/// carry this modification there: the unmodified form, or a form with another modification at the same position
/// (including one excluded from the sites, e.g. <c>Common Variable</c>).
/// <para>
/// Occupancy read from a stored table (e.g. MetaMorpheus's <c>IntensityOccupancy_</c> cells, <c>DEF-OCC-CELL</c>)
/// sets <see cref="ReportedFraction"/>: there the fraction is written exactly while the intensities are rounded,
/// so their ratio is not the value that was measured.
/// </para>
/// </remarks>
public sealed record SiteRunOccupancy(ModificationSite Site, string Run, OccupancyState State,
    double ModifiedIntensity, double CoveringIntensity, bool UnmodifiedQuantified)
{
    /// <summary>The fraction as its producer reported it, when that is more exact than the intensity ratio; otherwise null.</summary>
    public double? ReportedFraction { get; init; }

    /// <summary>
    /// When <see cref="State"/> is Quantified: <see cref="ReportedFraction"/> if set, else ModifiedIntensity / CoveringIntensity.
    /// Otherwise NaN.
    /// </summary>
    public double Fraction => State == OccupancyState.Quantified ? ReportedFraction ?? ModifiedIntensity / CoveringIntensity : double.NaN;
}

/// <summary>
/// Intensity-based site occupancy per run from plain peptidoform tables: QuantProject's <c>DEF-OCC-INT</c>
/// (the estimator MetaMorpheus writes as <c>IntensityOccupancy_</c> when no experimental design is used),
/// computed from stored per-run intensities instead of from live search objects.
/// </summary>
/// <remarks>
/// <para>
/// For a site (protein, position, modification) in a run: the numerator is the summed intensity of the
/// quantified peptidoforms that carry the modification at that position, and the denominator is the summed
/// intensity of every quantified peptidoform covering the position, modified or not. MetaMorpheus splits a
/// peptidoform's intensity evenly over its PSMs in the file and then sums over PSMs, which gives each
/// peptidoform its intensity once; this sums peptidoforms directly, with the same result.
/// </para>
/// <para>
/// As in <c>DEF-OCC-PSMS</c>: modifications of category <c>Common Variable</c> and <c>Common Fixed</c> are
/// excluded, as are peptide-terminal modifications. A protein N-terminal modification counts only on a
/// peptide starting at residue 1, or at residue 2 after removal of the initiator methionine (see
/// <see cref="PeptidoformObservation.PreviousResidue"/>). Whether a modification written before the first
/// residue is peptide- or protein-terminal is read from its location restriction in the mzLib modification
/// registry (<see cref="Mods"/>), matched on <c>Category:IdWithMotif</c>, else on <c>IdWithMotif</c>; a
/// modification the registry does not know is treated as protein N-terminal. Two modifications at one
/// position are two sites sharing one denominator.
/// </para>
/// <para>
/// Known differences from MetaMorpheus's output: PSMs whose full sequence was never resolved add to
/// MetaMorpheus's denominator but cannot be supplied here; and MetaMorpheus computes per protein group,
/// giving a shared peptide a proportional share, where this gives it to every protein it maps to.
/// <see cref="ModificationOccupancyCalculator"/> in Omics is the same estimator over live search objects.
/// </para>
/// </remarks>
public static class SiteOccupancyCalculator
{
    private static readonly string[] ExcludedCategories = ["Common Variable", "Common Fixed"];
    private static readonly ConcurrentDictionary<string, bool> PeptideNTerminalCache = new(StringComparer.Ordinal);

    /// <summary>
    /// Computes every site's occupancy state in every run where it was identified or covered.
    /// </summary>
    /// <param name="observations">
    /// MS/MS identifications, one row per (run, peptidoform, protein); a repeated row (e.g. one per PSM) is refused,
    /// since it would count the form's intensity several times.
    /// </param>
    /// <param name="includeModification">
    /// Optional filter on the modification name (<c>Category:IdWithMotif</c>), e.g. biological modifications
    /// only. Applied after the <c>DEF-OCC-PSMS</c> exclusions.
    /// </param>
    public static IReadOnlyList<SiteRunOccupancy> Calculate(IEnumerable<PeptidoformObservation> observations,
        Func<string, bool>? includeModification = null)
    {
        var parsed = Parse(observations, includeModification);

        var result = new List<SiteRunOccupancy>();
        // Sites of one protein are only ever covered by that protein's rows.
        foreach (var protein in parsed.GroupBy(p => p.obs.ProteinAccession, StringComparer.Ordinal)
                                      .OrderBy(g => g.Key, StringComparer.Ordinal))
        {
            var sites = protein.SelectMany(p => p.sites).Distinct()
                .OrderBy(s => s.Position).ThenBy(s => s.Modification, StringComparer.Ordinal).ToList();
            if (sites.Count == 0) continue;
            foreach (var run in protein.GroupBy(p => p.obs.Run, StringComparer.Ordinal).OrderBy(g => g.Key, StringComparer.Ordinal))
            {
                foreach (var site in sites)
                {
                    double modified = 0, covering = 0, unmodified = 0;
                    bool identified = false, covered = false;
                    foreach (var (obs, _, obsSites) in run)
                    {
                        if (site.Position < obs.StartResidue || site.Position > obs.EndResidue) continue;
                        bool carries = obsSites.Contains(site);
                        if (carries) identified = true;
                        if (double.IsNaN(obs.Intensity)) continue;
                        covered = true;
                        covering += obs.Intensity;
                        if (carries) modified += obs.Intensity; else unmodified += obs.Intensity;
                    }
                    OccupancyState? state =
                        identified && covered && covering > 0 ? (modified > 0 ? OccupancyState.Quantified : OccupancyState.Floor)
                        : identified ? OccupancyState.CountOnly
                        : covered && covering > 0 ? OccupancyState.NotDetected
                        : null;
                    if (state is OccupancyState s)
                        result.Add(new SiteRunOccupancy(site, run.Key, s, modified, covering, unmodified > 0));
                }
            }
        }
        return result;
    }

    /// <summary>
    /// Computes occupancy per sample rather than per run, for a sample measured as several runs (fractions).
    /// Each peptidoform's intensity is summed over the sample's runs, then <see cref="Calculate(IEnumerable{PeptidoformObservation}, Func{string, bool}?)"/>
    /// runs on the sums, so a site's numerator and denominator are both summed over the fractions.
    /// </summary>
    /// <remarks>
    /// Fractions are summed, as FlashLFQ and <c>CollapseFractions</c> sum them. Map only fractions together:
    /// technical replicates (repeat injections) are separate measurements of one sample and should keep their
    /// own key. A peptidoform identified but not quantified in every run of a sample stays unquantified (NaN);
    /// quantified in any run, it carries the sum of the quantified runs. So a modified form identified only in
    /// a fraction where it was not quantified, with the site covered in another fraction, is a Floor for the sample.
    /// </remarks>
    /// <param name="observations">As for the per-run overload, validated per run before anything is summed.</param>
    /// <param name="runToSample">Every run's sample; the sample becomes <see cref="SiteRunOccupancy.Run"/> of the result.</param>
    /// <param name="includeModification">As for the per-run overload.</param>
    /// <exception cref="ArgumentException">A run is missing from <paramref name="runToSample"/>, or maps to an empty sample.</exception>
    public static IReadOnlyList<SiteRunOccupancy> Calculate(IEnumerable<PeptidoformObservation> observations,
        IReadOnlyDictionary<string, string> runToSample, Func<string, bool>? includeModification = null)
        => Calculate(CombineObservations(observations, runToSample), includeModification);

    /// <summary>
    /// Sums each peptidoform's intensity over a sample's runs (fractions), giving one observation per (sample,
    /// peptidoform, protein) with <see cref="PeptidoformObservation.Run"/> set to the sample. Feed the result to
    /// <see cref="Calculate(IEnumerable{PeptidoformObservation}, Func{string, bool}?)"/> or
    /// <see cref="PtmPairEngine.Physical"/>, so that both see samples rather than fractions.
    /// </summary>
    /// <remarks>
    /// A peptidoform unquantified (NaN) in every run of a sample stays NaN; quantified in any, it carries the sum of
    /// the quantified runs. Map only fractions together; technical replicates keep their own key.
    /// </remarks>
    /// <param name="observations">Validated per run first, as for <see cref="Calculate(IEnumerable{PeptidoformObservation}, Func{string, bool}?)"/>.</param>
    /// <param name="runToSample">Every run's sample.</param>
    /// <exception cref="ArgumentException">A run is missing from <paramref name="runToSample"/>, or maps to an empty sample.</exception>
    public static IReadOnlyList<PeptidoformObservation> CombineObservations(IEnumerable<PeptidoformObservation> observations,
        IReadOnlyDictionary<string, string> runToSample)
    {
        ArgumentNullException.ThrowIfNull(runToSample);
        return Parse(observations, null)
            .GroupBy(p => (sample: SampleOf(p.obs.Run, runToSample), p.obs.FullSequence, p.obs.ProteinAccession, p.obs.StartResidue))
            .Select(g =>
            {
                var first = g.First().obs;
                var quantified = g.Select(p => p.obs.Intensity).Where(i => !double.IsNaN(i)).ToList();
                return first with { Run = g.Key.sample, Intensity = quantified.Count > 0 ? quantified.Sum() : double.NaN };
            })
            .ToList();
    }

    /// <summary>
    /// Combines stored per-run occupancy cells into one cell per (site, sample), for a sample measured as several
    /// runs (fractions). The catalog-side counterpart of the per-sample <see cref="Calculate(IEnumerable{PeptidoformObservation}, IReadOnlyDictionary{string, string}, Func{string, bool}?)"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Per (site, sample): the numerator is the sum over Quantified cells of <see cref="SiteRunOccupancy.Fraction"/>
    /// × <see cref="SiteRunOccupancy.CoveringIntensity"/>, so a cell's <see cref="SiteRunOccupancy.ReportedFraction"/>
    /// is respected; the denominator is the summed covering intensity of every cell that has one (Quantified,
    /// Floor, NotDetected). The state follows <see cref="Calculate(IEnumerable{PeptidoformObservation}, Func{string, bool}?)"/>:
    /// Quantified when any run is; otherwise Floor when the modified form was identified in some run and the site
    /// was covered in some run; otherwise CountOnly when identified; otherwise NotDetected.
    /// <see cref="SiteRunOccupancy.UnmodifiedQuantified"/> is true when it is true in any run.
    /// </para>
    /// <para>
    /// As with the observations, map only fractions together, never technical replicates.
    /// </para>
    /// </remarks>
    /// <param name="occupancy">One cell per (site, run).</param>
    /// <param name="runToSample">Every run's sample; the sample becomes <see cref="SiteRunOccupancy.Run"/> of the result.</param>
    /// <exception cref="ArgumentException">
    /// A run is missing from <paramref name="runToSample"/> or maps to an empty sample; a (site, run) repeats; or a
    /// Quantified cell has no finite, positive covering intensity to weight its fraction by.
    /// </exception>
    public static IReadOnlyList<SiteRunOccupancy> CombineRuns(IEnumerable<SiteRunOccupancy> occupancy,
        IReadOnlyDictionary<string, string> runToSample)
    {
        ArgumentNullException.ThrowIfNull(occupancy);
        ArgumentNullException.ThrowIfNull(runToSample);
        var cells = new List<SiteRunOccupancy>();
        var seen = new HashSet<(ModificationSite, string)>();
        foreach (var c in occupancy)
        {
            if (c is null) throw new ArgumentException("An occupancy cell is null.", nameof(occupancy));
            if (!seen.Add((c.Site, c.Run)))
                throw new ArgumentException($"{c.Site.Key} appears twice in {c.Run}; supply one cell per (site, run).", nameof(occupancy));
            if (c.State == OccupancyState.Quantified && !(double.IsFinite(c.CoveringIntensity) && c.CoveringIntensity > 0))
                throw new ArgumentException(
                    $"{c.Site.Key} in {c.Run} is Quantified but its covering intensity is {c.CoveringIntensity}; its fraction cannot be weighted.",
                    nameof(occupancy));
            cells.Add(c);
        }

        var result = new List<SiteRunOccupancy>();
        foreach (var g in cells.GroupBy(c => (c.Site, sample: SampleOf(c.Run, runToSample)))
                     .OrderBy(g => g.Key.Site.Key, StringComparer.Ordinal).ThenBy(g => g.Key.sample, StringComparer.Ordinal))
        {
            double modified = g.Where(c => c.State == OccupancyState.Quantified).Sum(c => c.Fraction * c.CoveringIntensity);
            double covering = g.Where(c => c.State != OccupancyState.CountOnly && double.IsFinite(c.CoveringIntensity))
                .Sum(c => c.CoveringIntensity);
            bool identified = g.Any(c => c.State is OccupancyState.Quantified or OccupancyState.Floor or OccupancyState.CountOnly);
            var state =
                g.Any(c => c.State == OccupancyState.Quantified) ? OccupancyState.Quantified
                : identified && covering > 0 ? OccupancyState.Floor
                : identified ? OccupancyState.CountOnly
                : OccupancyState.NotDetected;
            result.Add(new SiteRunOccupancy(g.Key.Site, g.Key.sample, state, modified, covering, g.Any(c => c.UnmodifiedQuantified)));
        }
        return result;
    }

    private static string SampleOf(string run, IReadOnlyDictionary<string, string> runToSample) =>
        runToSample.TryGetValue(run, out var sample)
            ? string.IsNullOrEmpty(sample)
                ? throw new ArgumentException($"Run {run} maps to an empty sample.", nameof(runToSample))
                : sample
            : throw new ArgumentException($"Run {run} has no sample in the run-to-sample map.", nameof(runToSample));

    /// <summary>
    /// Validates the observations (no null rows, coordinates matching the sequence, intensity finite and
    /// non-negative or NaN, one row per (run, peptidoform, protein, start)) and parses each one's sites.
    /// </summary>
    internal static List<(PeptidoformObservation obs, string baseSeq, List<ModificationSite> sites)> Parse(
        IEnumerable<PeptidoformObservation> observations, Func<string, bool>? includeModification)
    {
        ArgumentNullException.ThrowIfNull(observations);
        var parsed = new List<(PeptidoformObservation obs, string baseSeq, List<ModificationSite> sites)>();
        var seen = new HashSet<(string, string, string, int)>();
        foreach (var o in observations)
        {
            if (o is null) throw new ArgumentException("An observation is null.", nameof(observations));
            if (string.IsNullOrEmpty(o.Run) || string.IsNullOrEmpty(o.FullSequence) || string.IsNullOrEmpty(o.ProteinAccession))
                throw new ArgumentException("Every observation needs a run, a full sequence and a protein accession.", nameof(observations));
            string baseSeq = o.FullSequence.GetBaseSequenceFromFullSequence();
            if (o.StartResidue < 1 || o.EndResidue - o.StartResidue + 1 != baseSeq.Length)
                throw new ArgumentException(
                    $"'{o.FullSequence}' has {baseSeq.Length} residues but spans {o.StartResidue}-{o.EndResidue} of {o.ProteinAccession}.",
                    nameof(observations));
            if (double.IsInfinity(o.Intensity) || o.Intensity < 0)
                throw new ArgumentException($"Intensity of '{o.FullSequence}' in {o.Run} must be finite and non-negative, or NaN.", nameof(observations));
            if (!seen.Add((o.Run, o.FullSequence, o.ProteinAccession, o.StartResidue)))
                throw new ArgumentException(
                    $"'{o.FullSequence}' at {o.ProteinAccession} {o.StartResidue} appears twice in {o.Run}; supply one row per (run, peptidoform, protein), not one per PSM.",
                    nameof(observations));
            parsed.Add((o, baseSeq, SitesOf(o, baseSeq, includeModification)));
        }
        return parsed;
    }

    /// <summary>The sites a peptidoform carries, in protein coordinates, after the exclusions.</summary>
    internal static List<ModificationSite> SitesOf(PeptidoformObservation o, string baseSeq, Func<string, bool>? include)
    {
        var sites = new List<ModificationSite>();
        foreach (var (key, mod) in o.FullSequence.ParseModifications())
        {
            if (ExcludedCategories.Any(c => mod.StartsWith(c + ":", StringComparison.Ordinal))) continue;
            if (include != null && !include(mod)) continue;
            int position;
            char residue;
            if (key == 0)
            {
                // N-terminal: a protein N-terminal modification only when the peptide is the protein's N-terminus.
                if (o.StartResidue > 2) continue;
                if (o.StartResidue == 2 && o.PreviousResidue is char previous && previous != 'M') continue;
                if (IsPeptideNTerminal(mod)) continue;
                position = o.StartResidue;
                residue = baseSeq[0];
            }
            else if (key > baseSeq.Length)
            {
                continue; // C-terminal: protein length is not supplied, so protein and peptide C-termini cannot be told apart.
            }
            else
            {
                position = o.StartResidue + key - 1;
                residue = baseSeq[key - 1];
            }
            sites.Add(new ModificationSite(o.ProteinAccession, position, residue, mod));
        }
        return sites;
    }

    /// <summary>
    /// True when the registry knows <paramref name="mod"/> (<c>Category:IdWithMotif</c>) only as a peptide
    /// N-terminal modification: matched on the full name first, else on <c>IdWithMotif</c> alone.
    /// </summary>
    internal static bool IsPeptideNTerminal(string mod) => PeptideNTerminalCache.GetOrAdd(mod, name =>
    {
        var candidates = Mods.AllProteinModsList.Where(m => $"{m.ModificationType}:{m.IdWithMotif}" == name).ToList();
        if (candidates.Count == 0)
        {
            string idWithMotif = name.Contains(':') ? name[(name.IndexOf(':') + 1)..] : name;
            candidates = Mods.AllProteinModsList.Where(m => m.IdWithMotif == idWithMotif).ToList();
        }
        return candidates.Count > 0 && candidates.All(m => m.LocationRestriction == "Peptide N-terminal.");
    });
}
