namespace Omics.BioPolymerGroup;

/// <summary>
/// Builds the per-sample-group quantification columns shared by every group-level report:
/// spectral count and count-based occupancy for each sample group, plus intensity and
/// intensity-based occupancy for any sample group carrying intensity data.
///
/// The sample-group labels are the union across the whole dataset, so a record that was not
/// quantified in some condition contributes empty fields rather than a narrower row. Each report
/// supplies only the part that differs — how to render its own occupancy.
/// </summary>
public static class SampleGroupColumnBuilder
{
    /// <summary>
    /// </summary>
    /// <param name="records">Every record destined for the file.</param>
    /// <param name="sampleGroups">Returns a record's sample group results, populating them if needed.</param>
    /// <param name="countOccupancy">Formats count-based occupancy for one record and sample group.</param>
    /// <param name="intensityOccupancy">Formats intensity-based occupancy for one record and sample group.</param>
    public static IEnumerable<TsvColumn<T>> Build<T>(
        IReadOnlyCollection<T> records,
        Func<T, IReadOnlyList<SampleGroupResult>> sampleGroups,
        Func<T, SampleGroupResult, string> countOccupancy,
        Func<T, SampleGroupResult, string> intensityOccupancy)
        where T : IHasSampleIntensities
    {
        // Union of sample groups in first-seen order, plus which of them carry intensity anywhere in
        // the dataset. Taking the union is what keeps every row the header's width.
        //
        // Sample groups are matched across records by Identity — the files they cover — never by
        // Label and never by position. Labels are not unique (SampleGroupBuilder names a sample
        // group after its file whenever there is no design to name it by, so same-named files in
        // different directories collide), and position is not stable across records because a
        // record only has sample groups for the files it was observed in — so the same index means
        // different files for different records, which files one record's counts under another's
        // column.
        var identities = new List<string>();
        var seen = new HashSet<string>();

        // Whether the intensity columns exist is a property of the run, not of a sample group: an
        // engine either populated these records or it did not. Deliberately not
        // SampleGroupResult.HasIntensityData, which answers "did this sample group get a value" --
        // using that to choose columns made a record whose samples were all unobserved describe a
        // narrower table than its neighbour. The T constraint is what makes the predicate reachable,
        // and is why a record type must opt in via IHasSampleIntensities to be rendered here.
        //
        // Distinct from the TotalIntensity filter in FormatOccupancy, which stops a site printing
        // "0.0000(0/0)" inside a column that already exists.
        bool quantified = records.Any(r => r.SamplesForQuantification is { Count: > 0 }
                                        && r.IntensitiesBySample is not null);
        var labelByIdentity = new Dictionary<string, (string Label, string? LabelSourcePath)>();

        foreach (var record in records)
        {
            foreach (var result in sampleGroups(record))
            {
                if (seen.Add(result.Identity))
                {
                    identities.Add(result.Identity);
                    labelByIdentity[result.Identity] = (result.Label, result.LabelSourcePath);
                }

            }
        }

        var displayLabels = SampleGroupLabels.Disambiguate(identities, labelByIdentity);

        var index = new SampleGroupIndex<T>(sampleGroups);
        var columns = new List<TsvColumn<T>>();

        foreach (var identity in identities)
        {
            string thisIdentity = identity;
            string display = displayLabels[thisIdentity];

            columns.Add(new TsvColumn<T>($"SpectralCount_{display}",
                r => index.ResultFor(r, thisIdentity)?.SpectralCount.ToString() ?? string.Empty));

            if (quantified)
                columns.Add(new TsvColumn<T>($"Intensity_{display}",
                    r => index.ResultFor(r, thisIdentity) is { HasIntensityData: true } result
                        ? result.Intensity.ToString()
                        : string.Empty));

            columns.Add(new TsvColumn<T>($"CountOccupancy_{display}",
                r => index.ResultFor(r, thisIdentity) is { } result
                    ? countOccupancy(r, result)
                    : string.Empty));

            if (quantified)
                columns.Add(new TsvColumn<T>($"IntensityOccupancy_{display}",
                    r => index.ResultFor(r, thisIdentity) is { } result
                        ? intensityOccupancy(r, result)
                        : string.Empty));
        }

        return columns;
    }

    /// <summary>
    /// Finds a record's sample group by identity, in place of scanning its list once per cell.
    ///
    /// Decoupling a column from its position in the record's list is what fixes the ragged-row
    /// defect, but it also means a column can no longer read its value off the list in order. Left
    /// as a scan, each of the ~4 columns per sample group re-walks the whole list, so writing costs
    /// O(samples² × records) where rendering a row used to cost O(samples).
    ///
    /// A row is rendered one record at a time across every column, so remembering the most recent
    /// record's lookup restores O(samples × records).
    ///
    /// Keyed on the sample-group list instance rather than on the record: every invalidation path
    /// nulls SampleGroupResults, and repopulating assigns a fresh list, so a record whose
    /// quantification changed mid-write rebuilds its index instead of serving a stale one.
    /// Not thread-safe, in common with the lazy population it wraps.
    /// </summary>
    private sealed class SampleGroupIndex<T>(Func<T, IReadOnlyList<SampleGroupResult>> sampleGroups)
    {
        private IReadOnlyList<SampleGroupResult>? _source;
        private Dictionary<string, SampleGroupResult> _byIdentity = [];

        public SampleGroupResult? ResultFor(T record, string identity)
        {
            var results = sampleGroups(record);

            if (!ReferenceEquals(results, _source))
            {
                _source = results;
                _byIdentity = new Dictionary<string, SampleGroupResult>(results.Count);

                // First wins, matching the FirstOrDefault this replaces — identities are unique
                // within a record by construction, but TryAdd keeps a hand-built list from throwing.
                foreach (var result in results)
                    _byIdentity.TryAdd(result.Identity, result);
            }

            return _byIdentity.GetValueOrDefault(identity);
        }
    }
}
