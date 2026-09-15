namespace Omics.BioPolymerGroup;

/// <summary>
/// Turns sample group labels into column names that are unique within one output file.
///
/// A label is only a display name and is not unique — <see cref="SampleGroupBuilder"/> names a
/// sample group after its file whenever there is no experimental design to name it by, so
/// <c>Rep1\sample.raw</c> and <c>Rep2\sample.raw</c> both label <c>sample</c>. Emitting two columns
/// with the same name leaves a reader unable to tell which file a value came from, so colliding
/// labels are widened with as much of their parent path as it takes to separate them.
/// </summary>
public static class SampleGroupLabels
{
    /// <summary>
    /// Maps each sample group identity to a column name unique across <paramref name="identities"/>.
    /// Labels that do not collide are returned unchanged, so output is unaffected for the datasets
    /// that were already unambiguous.
    /// </summary>
    /// <param name="identities">Sample group identities, in the order their columns appear.</param>
    /// <param name="labels">Label and originating file path for each identity.</param>
    public static Dictionary<string, string> Disambiguate(
        IReadOnlyList<string> identities,
        IReadOnlyDictionary<string, (string Label, string? LabelSourcePath)> labels)
    {
        // Widen only what still clashes, and re-check against the whole set each round rather than
        // within one label's members: widening is itself a source of collisions, because the name a
        // widened column takes ("run1_sample") can be the plain label of some other file that was
        // never in that collision group ("run1_sample.raw").
        var depth = identities.ToDictionary(id => id, _ => 0);

        // Bounded by the deepest path; the no-progress branch below is the real exit.
        for (int round = 0; round <= MaxWideningRounds; round++)
        {
            var display = identities.ToDictionary(
                id => id,
                id => WidenPath(labels[id].LabelSourcePath, labels[id].Label, depth[id]));

            var clashing = display
                .GroupBy(entry => entry.Value)
                .Where(sameName => sameName.Count() > 1)
                .SelectMany(sameName => sameName.Select(entry => entry.Key))
                .ToList();

            if (clashing.Count == 0)
                return display;

            bool widened = false;
            foreach (var id in clashing)
            {
                if (CanWiden(labels[id].LabelSourcePath, depth[id] + 1))
                {
                    depth[id]++;
                    widened = true;
                }
            }

            // Nothing left to widen with — paths are exhausted, or the labels came from the
            // experimental design and have no path at all. Separate them by ordinal instead.
            if (!widened)
                return AppendOrdinals(identities, display);
        }

        return AppendOrdinals(
            identities,
            identities.ToDictionary(id => id, id => WidenPath(labels[id].LabelSourcePath, labels[id].Label, depth[id])));
    }

    private const int MaxWideningRounds = 64;

    /// <summary>
    /// Last resort when no more path is available: keeps the first occurrence and suffixes the rest,
    /// so the names are at least distinct even though they no longer say which file they came from.
    /// </summary>
    private static Dictionary<string, string> AppendOrdinals(
        IReadOnlyList<string> identities, Dictionary<string, string> display)
    {
        var used = new HashSet<string>();
        var resolved = new Dictionary<string, string>();

        foreach (var id in identities)
        {
            string name = display[id];

            if (used.Add(name))
            {
                resolved[id] = name;
                continue;
            }

            int ordinal = 2;
            while (!used.Add($"{name}_{ordinal}"))
                ordinal++;

            resolved[id] = $"{name}_{ordinal}";
        }

        return resolved;
    }

    /// <summary>
    /// Prefixes <paramref name="label"/> with up to <paramref name="depth"/> parent directory names
    /// from <paramref name="path"/>, nearest first — <c>Rep1_sample</c>, then <c>ExpA_Rep1_sample</c>.
    /// </summary>
    private static string WidenPath(string? path, string label, int depth)
    {
        if (string.IsNullOrEmpty(path))
            return label;

        var parents = ParentDirectories(path).Take(depth).Reverse();
        var parts = parents.Append(label);

        return string.Join("_", parts);
    }

    /// <summary>True when <paramref name="path"/> has at least <paramref name="depth"/> ancestors to widen with.</summary>
    private static bool CanWiden(string? path, int depth)
        => !string.IsNullOrEmpty(path) && ParentDirectories(path).Skip(depth - 1).Any();

    /// <summary>
    /// Directory names containing <paramref name="path"/>, nearest ancestor first.
    /// </summary>
    private static IEnumerable<string> ParentDirectories(string path)
    {
        var directory = Path.GetDirectoryName(path);

        while (!string.IsNullOrEmpty(directory))
        {
            string name = Path.GetFileName(directory);

            // A root such as "C:\" has no name of its own; stop rather than emit an empty segment.
            if (string.IsNullOrEmpty(name))
                yield break;

            yield return name;
            directory = Path.GetDirectoryName(directory);
        }
    }
}
