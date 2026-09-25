using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace MassSpectrometry
{
    /// <summary>
    /// The general-purpose <see cref="IExperimentalDesign"/>: which samples were measured in which file.
    /// Label-free designs hold one <see cref="ISampleInfo"/> per file; isobaric designs hold one per
    /// channel, and the per-file order is the contract that aligns them with
    /// <c>ISpectralMatch.Intensities</c>.
    ///
    /// Named for the samples it holds rather than simply "ExperimentalDesign", because MetaMorpheus
    /// already has an <c>EngineLayer.ExperimentalDesign</c> -- the reader for its ExperimentalDesign.tsv
    /// -- and files there import both namespaces. An unqualified "ExperimentalDesign" in
    /// <c>MassSpectrometry</c> makes every one of those references ambiguous.
    ///
    /// Before this type, the only implementations of the interface lived in the Development project and
    /// in test fixtures, so every consumer outside MetaMorpheus had to write its own.
    ///
    /// Files are keyed by name-with-extension, because that is what
    /// <c>QuantificationEngine</c> looks up: it reduces each spectral match's <c>FullFilePath</c> with
    /// <see cref="Path.GetFileName(string)"/>. <see cref="Add"/> therefore accepts a path or a bare name
    /// and stores the name either way, so a design built from full paths still resolves.
    ///
    /// Lookups are case-insensitive. A design that names <c>Sample1.raw</c> should match data at
    /// <c>sample1.raw</c>, which it would not under the default ordinal comparer — and a case-only
    /// collision is rejected at <see cref="Add"/> rather than silently keeping one of the two.
    /// </summary>
    public class SampleExperimentalDesign : IExperimentalDesign
    {
        /// <inheritdoc />
        public Dictionary<string, ISampleInfo[]> FileNameSampleInfoDictionary { get; }

        /// <summary>
        /// Creates an empty design. Populate it with <see cref="Add"/>.
        /// </summary>
        public SampleExperimentalDesign()
        {
            FileNameSampleInfoDictionary =
                new Dictionary<string, ISampleInfo[]>(StringComparer.OrdinalIgnoreCase);
        }

        /// <summary>
        /// Adds the samples measured in one file.
        /// </summary>
        /// <param name="fileNameOrPath">
        /// The file the samples were measured in, as either a full path or a bare name. Stored as the
        /// name with its extension.
        /// </param>
        /// <param name="samples">
        /// The samples, in the order their intensities appear in each spectral match. For an isobaric
        /// design this is the channel order; getting it wrong shifts every later channel, so the caller
        /// owns it. Must contain at least one sample and no nulls.
        /// </param>
        /// <exception cref="ArgumentException">
        /// The file name is empty, the sample array is empty or contains a null, the array lists one
        /// sample more than once (see <see cref="DescribeRepeatedSample"/>), or the file has already
        /// been added (including under different casing).
        /// </exception>
        public void Add(string fileNameOrPath, params ISampleInfo[] samples)
        {
            if (string.IsNullOrWhiteSpace(fileNameOrPath))
            {
                throw new ArgumentException("A file name is required.", nameof(fileNameOrPath));
            }

            if (samples == null || samples.Length == 0)
            {
                throw new ArgumentException(
                    $"File '{fileNameOrPath}' was added with no samples. A file with no samples cannot be quantified.",
                    nameof(samples));
            }

            if (samples.Any(s => s == null))
            {
                throw new ArgumentException(
                    $"File '{fileNameOrPath}' was added with a null sample. Intensities map to samples by position, " +
                    "so a missing channel has to be described rather than omitted.",
                    nameof(samples));
            }

            string fileName = Path.GetFileName(fileNameOrPath);
            if (string.IsNullOrWhiteSpace(fileName))
            {
                throw new ArgumentException(
                    $"'{fileNameOrPath}' has no file name component.", nameof(fileNameOrPath));
            }

            if (FileNameSampleInfoDictionary.ContainsKey(fileName))
            {
                throw new ArgumentException(
                    $"File '{fileName}' is already in this design. Add all of a file's samples in one call.",
                    nameof(fileNameOrPath));
            }

            string? repeated = DescribeRepeatedSample(samples);
            if (repeated != null)
            {
                throw new ArgumentException(
                    $"File '{fileNameOrPath}': {repeated}. Quantification indexes columns by sample, so " +
                    "the repeats would merge into one column under one of their names. List each sample " +
                    "once.",
                    nameof(samples));
            }

            FileNameSampleInfoDictionary[fileName] = samples.ToArray();
        }

        /// <summary>
        /// Describes the first sample listed more than once among <paramref name="samples"/>, or returns
        /// null when every sample is distinct. Null entries are skipped.
        /// </summary>
        /// <remarks>
        /// "The same sample" means equal, which for an isobaric channel is the same file and channel
        /// label — deliberately not its <see cref="IsobaricQuantSampleInfo.SampleName"/> — and for a
        /// label-free <see cref="SpectraFileInfo"/> the same file, condition and replicates. Two such
        /// entries are one column to quantification, whose matrices index columns by sample: the later
        /// silently takes the earlier one's place, and the merged column carries only one of their
        /// names. Public SDRF has exactly this shape — <c>PXD040455</c>, a TMT × SILAC design, lists a
        /// light and a heavy sample under one reporter channel of one file 551 times.
        ///
        /// Shared by <see cref="Add"/>, which refuses such a design one file at a time as it is built, and
        /// by the quantification engine, which asks it of every file's samples at once, from any
        /// <see cref="IExperimentalDesign"/> implementation. Not every design is built through this
        /// class, and the engine merges all files' columns into one matrix, so a sample listed under two
        /// file keys collides there as surely as one listed twice under one.
        /// </remarks>
        public static string? DescribeRepeatedSample(IEnumerable<ISampleInfo> samples)
        {
            var repeated = samples
                .Where(sample => sample != null)
                .GroupBy(sample => sample)
                .FirstOrDefault(sameSample => sameSample.Count() > 1);

            if (repeated == null)
            {
                return null;
            }

            int count = repeated.Count();

            if (repeated.Key is IsobaricQuantSampleInfo channel)
            {
                var names = repeated
                    .OfType<IsobaricQuantSampleInfo>()
                    .Select(c => c.SampleName)
                    .Where(name => !string.IsNullOrWhiteSpace(name))
                    .Distinct(StringComparer.Ordinal)
                    // Sorted, so the message does not depend on the order a caller's collection enumerates
                    // in; the engine hands over a Dictionary's values, whose order is not guaranteed.
                    .OrderBy(name => name, StringComparer.Ordinal)
                    .Select(name => $"'{name}'")
                    .ToList();

                string asNamed = names.Count == 0 ? string.Empty : $", as {string.Join(" and ", names)}";
                return $"channel {channel.ChannelLabel} of '{Path.GetFileName(channel.FullFilePathWithExtension)}' " +
                       $"is listed {count} times{asNamed}";
            }

            return $"sample '{repeated.Key}' is listed {count} times";
        }

        /// <summary>
        /// Builds a design from a flat sequence of samples, grouping them by
        /// <see cref="ISampleInfo.FullFilePathWithExtension"/>. Works for both modalities: label-free
        /// samples each name a different file, isobaric channels share one.
        ///
        /// Order within a file is the order of the input sequence, and that is the order the engine will
        /// align intensities to — so pass isobaric channels already sorted the way the search writes
        /// them (ascending reporter m/z, for MetaMorpheus).
        ///
        /// Samples that share a file name but came from different directories are rejected rather than
        /// grouped together. Isobaric channels sharing one file is the intended case; two different
        /// files that merely happen to be called the same thing is not, and merging them would put one
        /// run's channels in another run's row. The design cannot tell them apart afterwards either --
        /// <see cref="IExperimentalDesign"/> is keyed by file name because that is what the engine looks
        /// up -- so this is caught where the information still exists.
        /// </summary>
        /// <exception cref="ArgumentException">
        /// A sample is null, names no file, or shares a file name with a sample from a different path.
        /// </exception>
        public static SampleExperimentalDesign FromSamples(IEnumerable<ISampleInfo> samples)
        {
            if (samples == null)
            {
                throw new ArgumentNullException(nameof(samples));
            }

            var design = new SampleExperimentalDesign();

            var byFile = samples
                .Select((sample, index) => (sample, index))
                .GroupBy(t =>
                {
                    if (t.sample == null)
                    {
                        throw new ArgumentException(
                            $"The sample at index {t.index} is null.", nameof(samples));
                    }

                    string path = t.sample.FullFilePathWithExtension;
                    if (string.IsNullOrWhiteSpace(path))
                    {
                        throw new ArgumentException(
                            $"The sample at index {t.index} names no file, so it cannot be grouped into a design.",
                            nameof(samples));
                    }

                    return Path.GetFileName(path);
                }, StringComparer.OrdinalIgnoreCase);

            foreach (var fileGroup in byFile)
            {
                // Channels of one run share a path; two different paths reaching the same key are two
                // different files whose samples would otherwise be silently interleaved into one row.
                var distinctPaths = fileGroup
                    .Select(t => t.sample.FullFilePathWithExtension)
                    .Distinct(StringComparer.OrdinalIgnoreCase)
                    .ToList();

                if (distinctPaths.Count > 1)
                {
                    throw new ArgumentException(
                        $"'{fileGroup.Key}' names {distinctPaths.Count} different files: " +
                        string.Join(", ", distinctPaths) +
                        ". An experimental design is keyed by file name, so these cannot be told apart " +
                        "once it is built.",
                        nameof(samples));
                }

                design.Add(fileGroup.Key, fileGroup.Select(t => t.sample).ToArray());
            }

            return design;
        }

        /// <summary>
        /// Builds a label-free design: one sample per file.
        /// </summary>
        /// <exception cref="ArgumentException">
        /// A file is null, names no path, or appears twice. Label-free measures a file once; a repeat is
        /// a caller mistake rather than a second channel.
        /// </exception>
        public static SampleExperimentalDesign LabelFree(IEnumerable<SpectraFileInfo> files)
        {
            if (files == null)
            {
                throw new ArgumentNullException(nameof(files));
            }

            var design = new SampleExperimentalDesign();

            foreach (var file in files)
            {
                if (file == null)
                {
                    throw new ArgumentException("A null file cannot be added to a design.", nameof(files));
                }

                if (string.IsNullOrWhiteSpace(file.FullFilePathWithExtension))
                {
                    throw new ArgumentException(
                        "A file with no path cannot be added to a design.", nameof(files));
                }

                design.Add(file.FullFilePathWithExtension, file);
            }

            return design;
        }
    }
}
