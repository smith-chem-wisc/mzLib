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
    /// for an entrapment decoy (ED); <see cref="ProteinGroupFromTsv.IsContaminant"/> only for C. The
    /// descriptor carries no entrapment mark because the label cannot be trusted for one: MetaMorpheus
    /// writes T for an entrapment group. The annotator finds entrapment members from their accessions
    /// and names them on every row.</remarks>
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
    ///
    /// MetaMorpheus sometimes writes the same group twice, identically. Each group is returned once, at
    /// its first position, so it is annotated and counted once. Two rows with the same group name that
    /// disagree on members, marks or q-value cannot be told apart and are refused.
    /// </summary>
    /// <exception cref="ArgumentNullException">rows is null.</exception>
    /// <exception cref="InvalidDataException">Two rows share a group name but differ, raised when enumerated.</exception>
    public static IEnumerable<GoAnnotationGroup> ToGoAnnotationGroups(this IEnumerable<ProteinGroupFromTsv> rows)
    {
        ArgumentNullException.ThrowIfNull(rows);
        return Distinct(rows.Where(r => !r.IsDecoy).Select(r => r.ToGoAnnotationGroup()));
    }

    private static IEnumerable<GoAnnotationGroup> Distinct(IEnumerable<GoAnnotationGroup> groups)
    {
        var seen = new Dictionary<string, GoAnnotationGroup>(StringComparer.Ordinal);
        foreach (var group in groups)
        {
            if (seen.TryGetValue(group.Name, out var first))
            {
                if (!first.MemberAccessions.SequenceEqual(group.MemberAccessions, StringComparer.Ordinal)
                    || first.IsDecoy != group.IsDecoy || first.IsContaminant != group.IsContaminant
                    || !first.QValue.Equals(group.QValue))
                {
                    throw new InvalidDataException(
                        $"Protein group '{group.Name}' appears twice with different members, marks or q-value.");
                }
                continue;
            }
            seen.Add(group.Name, group);
            yield return group;
        }
    }
}
