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
    public class ExperimentalDesign : IExperimentalDesign
    {
        /// <inheritdoc />
        public Dictionary<string, ISampleInfo[]> FileNameSampleInfoDictionary { get; }

        /// <summary>
        /// Creates an empty design. Populate it with <see cref="Add"/>.
        /// </summary>
        public ExperimentalDesign()
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
        /// The file name is empty, the sample array is empty or contains a null, or the file has already
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

            FileNameSampleInfoDictionary[fileName] = samples.ToArray();
        }

        /// <summary>
        /// Builds a design from a flat sequence of samples, grouping them by
        /// <see cref="ISampleInfo.FullFilePathWithExtension"/>. Works for both modalities: label-free
        /// samples each name a different file, isobaric channels share one.
        ///
        /// Order within a file is the order of the input sequence, and that is the order the engine will
        /// align intensities to — so pass isobaric channels already sorted the way the search writes
        /// them (ascending reporter m/z, for MetaMorpheus).
        /// </summary>
        /// <exception cref="ArgumentException">A sample is null or names no file.</exception>
        public static ExperimentalDesign FromSamples(IEnumerable<ISampleInfo> samples)
        {
            if (samples == null)
            {
                throw new ArgumentNullException(nameof(samples));
            }

            var design = new ExperimentalDesign();

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
        public static ExperimentalDesign LabelFree(IEnumerable<SpectraFileInfo> files)
        {
            if (files == null)
            {
                throw new ArgumentNullException(nameof(files));
            }

            var design = new ExperimentalDesign();

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
