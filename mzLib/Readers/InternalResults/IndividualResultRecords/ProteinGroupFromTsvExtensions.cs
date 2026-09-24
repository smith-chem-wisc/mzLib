using UsefulProteomicsDatabases.GeneOntology;

namespace Readers;

/// <summary>
/// Turns stored MetaMorpheus protein-group rows into the inputs of other mzLib engines, so they run on
/// results already on disk rather than inside a search.
/// </summary>
public static class ProteinGroupFromTsvExtensions
{
    /// <summary>
    /// The group as <see cref="GoGroupAnnotator"/> takes it: its name, every member accession (no member
    /// is privileged), the decoy and contaminant marks, and the group's q-value.
    /// </summary>
    /// <remarks>The flags are the reader's: <see cref="ProteinGroupFromTsv.IsDecoy"/> is true for D and
    /// for an entrapment decoy (ED); <see cref="ProteinGroupFromTsv.IsContaminant"/> only for C. An
    /// entrapment target (ET) is neither, and the descriptor carries no entrapment mark.</remarks>
    /// <exception cref="ArgumentNullException">row is null.</exception>
    public static GoAnnotationGroup ToGoAnnotationGroup(this ProteinGroupFromTsv row)
    {
        ArgumentNullException.ThrowIfNull(row);
        return new GoAnnotationGroup(row.ProteinGroupName, row.Accessions, row.IsDecoy, row.IsContaminant, row.QValue);
    }

    /// <summary>
    /// Every non-decoy group, in file order, contaminants included -- the population the GO annotation
    /// file is defined over. Decoys are skipped here because the annotator refuses them. Rows are not
    /// filtered by q-value: every group is kept and the consumer filters.
    /// </summary>
    /// <exception cref="ArgumentNullException">rows is null.</exception>
    public static IEnumerable<GoAnnotationGroup> ToGoAnnotationGroups(this IEnumerable<ProteinGroupFromTsv> rows)
    {
        ArgumentNullException.ThrowIfNull(rows);
        return rows.Where(r => !r.IsDecoy).Select(r => r.ToGoAnnotationGroup());
    }
}
