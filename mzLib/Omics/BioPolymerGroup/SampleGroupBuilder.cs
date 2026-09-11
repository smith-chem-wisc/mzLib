using MassSpectrometry;
using Omics.SpectralMatch;

namespace Omics.BioPolymerGroup;

/// <summary>
/// Buckets PSMs and intensities into <see cref="SampleGroupResult"/>s according to the
/// experimental design. This is the part of quantification that does not care what is being
/// grouped, so parent-level and digestion-product-level groups share one copy; each supplies its
/// own occupancy calculation through a callback.
/// </summary>
public static class SampleGroupBuilder
{
    /// <summary>
    /// Groups by (Condition × BiologicalReplicate) for label-free samples, by (File × Channel)
    /// for isobaric samples, and by PSM source file when no experimental design is supplied.
    ///
    /// Isobaric results additionally declare a <see cref="SampleGroupResult.CountIdentity"/> of the
    /// file they came from, because a spectral count is a property of the acquired file while an
    /// intensity is a property of the channel. The other two designs leave it at its default, where
    /// the two spaces coincide.
    /// </summary>
    /// <param name="samples">Samples contributing quantification, or null when no design exists.</param>
    /// <param name="intensitiesBySample">Measured intensities keyed by sample, or null when unavailable.
    /// When null, results carry no intensity data and only spectral counts are reported.</param>
    /// <param name="psms">The PSMs to distribute across the resulting sample groups.</param>
    /// <param name="populateOccupancy">Invoked once per result with the PSMs belonging to it, so the
    /// caller can attach modification occupancy in whichever coordinate space it reports.</param>
    public static List<SampleGroupResult> Build(
        IReadOnlyList<ISampleInfo>? samples,
        IReadOnlyDictionary<ISampleInfo, double>? intensitiesBySample,
        IReadOnlyCollection<ISpectralMatch> psms,
        Action<SampleGroupResult, List<ISpectralMatch>> populateOccupancy)
    {
        var results = new List<SampleGroupResult>();

        var spectraFiles = samples?.OfType<SpectraFileInfo>().ToList() ?? [];
        var isobaricSamples = samples?.OfType<IsobaricQuantSampleInfo>().ToList() ?? [];

        if (spectraFiles.Count > 0)
        {
            bool unfractionated = spectraFiles.Select(p => p.Fraction).Distinct().Count() == 1;
            bool conditionsUndefined = spectraFiles.All(p => string.IsNullOrEmpty(p.Condition));
            bool silacExperimentalDesign = spectraFiles.Any(p => !File.Exists(p.FullFilePathWithExtension));

            foreach (var conditionGroup in spectraFiles.GroupBy(p => p.Condition))
            {
                foreach (var bioRepGroup in conditionGroup.GroupBy(p => p.BiologicalReplicate).OrderBy(p => p.Key))
                {
                    var filesInGroup = bioRepGroup.ToList();
                    bool labelFromFileName = (conditionsUndefined && unfractionated) || silacExperimentalDesign;
                    string label = labelFromFileName
                        ? filesInGroup.First().FilenameWithoutExtension
                        : $"{conditionGroup.Key}_{bioRepGroup.Key + 1}";

                    var filePaths = new HashSet<string>(filesInGroup.Select(f => f.FullFilePathWithExtension));
                    var psmsInGroup = psms.Where(p => filePaths.Contains(p.FullFilePath)).ToList();

                    // With intensities available, carry the per-file values; otherwise leave them
                    // empty so HasIntensityData stays false and only spectral counts are reported.
                    var groupIntensities = new Dictionary<string, double>();
                    if (intensitiesBySample != null)
                    {
                        foreach (var file in filesInGroup)
                        {
                            // Keyed by full path, not file name: a sample group spans fractions and
                            // technical replicates, which are routinely stored one per directory
                            // under the same name. Keying by name silently summed the wrong pair
                            // here and threw outright in FilesInGroup below.
                            if (intensitiesBySample.TryGetValue(file, out var fileIntensity))
                                groupIntensities[file.FullFilePathWithExtension] = fileIntensity;
                        }
                    }

                    var result = new SampleGroupResult(conditionGroup.Key, bioRepGroup.Key)
                    {
                        Label = label,
                        // (Condition, BiologicalReplicate) identifies the sample group, not the label:
                        // two files with the same name in different directories share a label but are
                        // distinct samples. These two values are unique here by construction — the
                        // enclosing GroupBy pair produces exactly one bucket per (condition, replicate) —
                        // and the replicate is a trailing integer, so the join cannot be ambiguous.
                        Identity = $"{conditionGroup.Key}|{bioRepGroup.Key}",
                        LabelSourcePath = labelFromFileName ? filesInGroup.First().FullFilePathWithExtension : null,
                        SpectralCount = psmsInGroup.Count,
                        FilesInGroup = filesInGroup.ToDictionary(kvp => kvp.FullFilePathWithExtension, kvp => (ISampleInfo)kvp),
                        IntensitiesBySample = groupIntensities
                    };

                    populateOccupancy(result, psmsInGroup);
                    results.Add(result);
                }
            }
        }
        else if (isobaricSamples.Count > 0)
        {
            foreach (var fileGroup in isobaricSamples.GroupBy(p => p.FullFilePathWithExtension).OrderBy(g => g.Key))
            {
                var psmsInFile = psms.Where(p => p.FullFilePath.Equals(fileGroup.Key)).ToList();

                // Ascending reporter-ion m/z, which is the order mzLib already declares canonical:
                // SampleExperimentalDesign.FromSamples says to "pass isobaric channels already
                // sorted the way the search writes them (ascending reporter m/z, for MetaMorpheus)",
                // and QuantificationWriter's positional Reporter_n columns follow that same design
                // array. Ordering by ChannelLabel instead sorts "127C" before "127N" in every N/C
                // pair -- C precedes N ordinally, but the C reporter is the heavier ion -- so the
                // grouped table disagreed with the raw one about which 127 came first.
                //
                // Sorting on an explicit key rather than OrderBy(p => p): IsobaricQuantSampleInfo's
                // CompareTo does not look at ChannelLabel or ReporterIonMz at all, and returns 1 in
                // both directions for two channels of one file, which is not a total order.
                //
                // Label is the tiebreak so the order stays deterministic when a design leaves
                // ReporterIonMz at its default -- nothing validates it -- rather than permuting
                // arbitrarily.
                foreach (var sample in fileGroup.OrderBy(p => p.ReporterIonMz)
                                                .ThenBy(p => p.ChannelLabel, StringComparer.Ordinal))
                {
                    string label = $"{Path.GetFileNameWithoutExtension(sample.FullFilePathWithExtension)}_{sample.ChannelLabel}";

                    var channelIntensities = new Dictionary<string, double>();
                    if (intensitiesBySample != null && intensitiesBySample.TryGetValue(sample, out var channelIntensity))
                        channelIntensities[label] = channelIntensity;

                    var result = new SampleGroupResult(sample.Condition, sample.BiologicalReplicate)
                    {
                        Label = label,
                        Identity = $"{sample.FullFilePathWithExtension}|{sample.ChannelLabel}",
                        LabelSourcePath = sample.FullFilePathWithExtension,

                        // The count belongs to the file, not to the channel: psmsInFile is one
                        // file's PSMs and every channel of that file is handed the same list. Each
                        // channel result still carries the value, because a result is what a column
                        // reads from, but they all carry the SAME value and the schema renders it
                        // once per file instead of once per channel.
                        CountIdentity = fileGroup.Key,
                        CountLabel = Path.GetFileNameWithoutExtension(fileGroup.Key),

                        SpectralCount = psmsInFile.Count,
                        FilesInGroup = new Dictionary<string, ISampleInfo> { { label, sample } },
                        IntensitiesBySample = channelIntensities
                    };

                    populateOccupancy(result, psmsInFile);
                    results.Add(result);
                }
            }
        }
        else
        {
            // No experimental design — group PSMs by source file for count-only results
            foreach (var fileGroup in psms.GroupBy(p => p.FullFilePath).OrderBy(g => g.Key))
            {
                var psmsInFile = fileGroup.ToList();

                var result = new SampleGroupResult(string.Empty, 0)
                {
                    Label = Path.GetFileNameWithoutExtension(fileGroup.Key),
                    Identity = fileGroup.Key,
                    LabelSourcePath = fileGroup.Key,
                    SpectralCount = psmsInFile.Count
                    // FilesInGroup and IntensitiesBySample left empty → HasIntensityData = false
                };

                populateOccupancy(result, psmsInFile);
                results.Add(result);
            }
        }

        return results;
    }
}
