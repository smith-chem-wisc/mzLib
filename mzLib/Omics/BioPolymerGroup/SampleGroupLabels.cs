using MassSpectrometry;

namespace Omics.BioPolymerGroup;

/// <summary>
/// Builds the label a sample's quantification columns are named by, and turns those labels into
/// column names that are unique within one output file.
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
    /// The label for one sample's columns. The grouped protein table and the quantification matrices
    /// both start from it, so a channel whose label is unique within its output is named the same way
    /// in both.
    ///
    /// Labels that collide are still made unique by each writer separately, and not alike: the grouped
    /// table widens them with <see cref="Disambiguate"/> (directory, then ordinals), while
    /// <c>QuantificationWriter.UniqueColumnLabels</c> only appends ordinals. A collided channel can
    /// therefore carry a different header in each file.
    ///
    /// An isobaric channel is labelled <c>{sample}_{file}_{channel}</c> when the design names the
    /// sample in it, and <c>{file}_{channel}</c> otherwise. The file stays in both forms: the fractions
    /// of one plex share their sample and differ only by file, so a label without it would name six
    /// fractions alike and leave <see cref="Disambiguate"/> nothing to separate them by but ordinals.
    ///
    /// A label-free sample is labelled by its file name when <paramref name="labelFreeByFileName"/>
    /// is set and <c>{condition}_{biorep + 1}</c> otherwise.
    /// </summary>
    /// <remarks>
    /// An unnamed channel keeps exactly the label it had before sample names existed, rather than
    /// falling back to its condition and replicate. That is what makes adding the name additive: a
    /// caller that does not yet supply <see cref="IsobaricQuantSampleInfo.SampleName"/> sees no
    /// column change at all, so it can adopt this release first and opt into the new labels when it
    /// starts passing names.
    ///
    /// The parts are joined with <c>_</c>, which a sample name or file stem may itself contain, so
    /// different channels can share a label: <c>S_1</c> in <c>run.raw</c> and <c>S</c> in
    /// <c>1_run.raw</c> both give <c>S_1_run_126</c>, and a channel named <c>plex1</c> in <c>x.raw</c>
    /// takes the label of the same, unnamed channel of <c>plex1_x.raw</c>. Such labels collide and are
    /// made unique as above. Values stay on their own channels, because columns are keyed on the
    /// sample, never on its label.
    /// </remarks>
    /// <param name="sample">The sample to label.</param>
    /// <param name="labelFreeByFileName">For a label-free sample, whether to label it by file name
    /// rather than by condition and replicate. Ignored for an isobaric sample.</param>
    public static string ForSample(ISampleInfo sample, bool labelFreeByFileName = false)
    {
        if (sample is IsobaricQuantSampleInfo isobaric)
        {
            string fileAndChannel = $"{Path.GetFileNameWithoutExtension(isobaric.FullFilePathWithExtension)}_{isobaric.ChannelLabel}";

            return string.IsNullOrWhiteSpace(isobaric.SampleName)
                ? fileAndChannel
                : $"{isobaric.SampleName}_{fileAndChannel}";
        }

        return labelFreeByFileName
            ? sample.FilenameWithoutExtension
            : $"{sample.Condition}_{sample.BiologicalReplicate + 1}";
    }

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
