using System;
using System.Collections.Generic;
using System.Linq;
using MzLibUtil;

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
public sealed record PeptidoformObservation(string Run, string FullSequence, string ProteinAccession,
    int StartResidue, int EndResidue, double Intensity);

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
/// (the unmodified form may be present below detection), as a floor is treated at the other end.
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
/// peptide starting at residue 1 or 2 (after removal of the initiator methionine). Two modifications at one
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

    /// <summary>
    /// Computes every site's occupancy state in every run where it was identified or covered.
    /// </summary>
    /// <param name="observations">MS/MS identifications, one row per (run, peptidoform, protein).</param>
    /// <param name="includeModification">
    /// Optional filter on the modification name (<c>Category:IdWithMotif</c>), e.g. biological modifications
    /// only. Applied after the <c>DEF-OCC-PSMS</c> exclusions.
    /// </param>
    public static IReadOnlyList<SiteRunOccupancy> Calculate(IEnumerable<PeptidoformObservation> observations,
        Func<string, bool>? includeModification = null)
    {
        ArgumentNullException.ThrowIfNull(observations);
        var parsed = new List<(PeptidoformObservation obs, string baseSeq, List<ModificationSite> sites)>();
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
            parsed.Add((o, baseSeq, SitesOf(o, baseSeq, includeModification)));
        }

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
}
