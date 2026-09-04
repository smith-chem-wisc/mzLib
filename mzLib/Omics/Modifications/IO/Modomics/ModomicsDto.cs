namespace Omics.Modifications.IO.Modomics;

/// <summary>
/// A single entry of the MODOMICS RNA modifications database, as exported by its REST API
/// (https://genesilico.pl/modomics/api/) and embedded as Resources/modomicsmods.json.
/// Field semantics follow the MODOMICS help page (https://genesilico.pl/modomics/help); see the
/// 2025 database update (doi:10.1093/nar/gkaf1284) for the citable reference.
/// Every published field is loaded, including those the loader does not currently consume,
/// so future features can draw on them.
/// </summary>
public class ModomicsDto
{
    /// <summary>MODOMICS numeric database id (the JSON key of this entry).</summary>
    public required int Id { get; init; }

    /// <summary>Single-character Unicode abbreviation used in MODOMICS sequence notation.</summary>
    public required string Abbrev { get; init; } = string.Empty;

    /// <summary>
    /// Sum formula of the whole moiety, verbatim from the source; what it describes depends on
    /// <see cref="MoietyType"/>. Cationic entries are suffixed with '+' and nucleotide entries are
    /// published as phosphate anions — consumers must handle both.
    /// </summary>
    public required string Formula { get; init; } = string.Empty;

    /// <summary>LC elution order/characteristics relative to canonical nucleosides — free text.</summary>
    public string? LcElutionComment { get; init; }

    /// <summary>
    /// Normalized LC retention time (reversed-phase C18, normalized to guanosine), often with a
    /// literature citation — free text, not a number.
    /// </summary>
    public string? LcElutionTime { get; init; }

    /// <summary>Average mass as reported by MODOMICS (may be null, or inconsistent with <see cref="Formula"/>).</summary>
    public double? MassAvg { get; init; }

    /// <summary>
    /// Neutral monoisotopic mass as reported by MODOMICS. Unreliable: the source is internally
    /// inconsistent for cationic entries (e.g. m3C) and some caps; the loader always computes mass
    /// from <see cref="Formula"/> instead.
    /// </summary>
    public double? MassMonoiso { get; init; }

    /// <summary>[M+H]+ protonated mass as reported by MODOMICS.</summary>
    public double? MassProt { get; init; }

    /// <summary>Full chemical name (also used as the modification's OriginalId by the loader).</summary>
    public required string Name { get; init; } = string.Empty;

    /// <summary>
    /// Product ions: protonated fragment ions generated from the precursor [M+H]+, showing the typical
    /// neutral losses for the nucleoside. Published as nominal integer m/z, '/'-separated; predominantly
    /// the protonated (modified) base ion plus secondary fragments.
    /// </summary>
    public string? ProductIons { get; init; }

    /// <summary>
    /// The original nucleobase(s) this modification derives from: A, C, G, or U. Generic/novel bases use
    /// "X"; one entry (peroxywybutosine) uses "QtRNA" for a queuosine-containing tRNA context; caps and
    /// terminal structures list all four bases.
    /// </summary>
    public required List<string> ReferenceMoieties { get; init; } = [];

    /// <summary>Common short name / acronym (e.g. "m6A"); unique per database entry.</summary>
    public required string ShortName { get; init; } = string.Empty;

    /// <summary>SMILES structure of the whole moiety.</summary>
    public string? Smile { get; init; }

    /// <summary>
    /// Moiety type, joined from the CSV export rather than the JSON: "nucleoside" (neutral, whole
    /// nucleoside), "nucleotide" (whole monophosphate, published as an anion), or "base" (base only,
    /// requires a new residue definition). Empty when the CSV has no row for the short name.
    /// </summary>
    public required string MoietyType { get; init; } = string.Empty;
}
