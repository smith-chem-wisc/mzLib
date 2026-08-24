using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using MassSpectrometry;
using Omics;
using Omics.BioPolymerGroup;
using Omics.SpectralMatch;

namespace Quantification
{
    /// <summary>
    /// Writes quantification output as tab-separated text.
    ///
    /// Three files are produced, each optional and controlled by the corresponding
    /// <see cref="QuantificationParameters"/> flag:
    ///
    /// <list type="bullet">
    /// <item><c>RawQuantification.tsv</c> — the PSM-level snapshot, taken before any normalization or
    /// roll-up. This is the reloadable record: re-running with different strategies does not require
    /// re-searching.</item>
    /// <item><c>PeptideQuantification.tsv</c> — the final peptide matrix.</item>
    /// <item><c>ProteinGroupQuantification.tsv</c> — the final protein-group matrix.</item>
    /// </list>
    ///
    /// Each method returns the full path it wrote, or null when there was nothing to write.
    /// </summary>
    public static class QuantificationWriter
    {
        public const string RawFileName = "RawQuantification.tsv";
        public const string PeptideFileName = "PeptideQuantification.tsv";
        public const string ProteinGroupFileName = "ProteinGroupQuantification.tsv";

        /// <summary>
        /// Writes the PSM-level intensities exactly as they arrived, before normalization or roll-up.
        /// One row per spectral match; one column per entry in <see cref="ISpectralMatch.Intensities"/>.
        /// Matches with a null <c>Intensities</c> array are skipped, as they carry no quantification.
        /// </summary>
        /// <returns>The path written, or null if <paramref name="matches"/> held nothing quantified.</returns>
        public static string WriteRawData(List<ISpectralMatch> matches, string outputDirectory)
        {
            var quantified = matches?
                .Where(m => m?.Intensities != null)
                .OrderBy(m => m.FullFilePath, StringComparer.Ordinal)
                .ThenBy(m => m.OneBasedScanNumber)
                .ToList();

            if (quantified is not { Count: > 0 })
                return null;

            int channelCount = quantified.Max(m => m.Intensities.Length);
            string path = ResolvePath(outputDirectory, RawFileName);

            using var output = new StreamWriter(path);

            var header = new StringBuilder("File\tScan\tBaseSequence\tFullSequence\tAccession\tIsDecoy\tScore");
            for (int i = 0; i < channelCount; i++)
                header.Append("\tIntensity_").Append(i + 1);
            output.WriteLine(header.ToString());

            var line = new StringBuilder();
            foreach (var match in quantified)
            {
                line.Clear();
                line.Append(match.FullFilePath).Append('\t')
                    .Append(match.OneBasedScanNumber.ToString(CultureInfo.InvariantCulture)).Append('\t')
                    .Append(match.BaseSequence).Append('\t')
                    .Append(match.FullSequence).Append('\t')
                    .Append(match.Accession).Append('\t')
                    .Append(match.IsDecoy ? "D" : "T").Append('\t')
                    .Append(Format(match.Score));

                for (int i = 0; i < channelCount; i++)
                {
                    line.Append('\t');
                    line.Append(i < match.Intensities.Length ? Format(match.Intensities[i]) : string.Empty);
                }
                output.WriteLine(line.ToString());
            }

            return path;
        }

        /// <summary>
        /// Writes the peptide matrix: one row per peptide, one column per sample.
        /// </summary>
        /// <returns>The path written, or null if the matrix was empty.</returns>
        public static string WritePeptideMatrix(QuantMatrix<IBioPolymerWithSetMods> peptideMatrix, string outputDirectory) =>
            WriteMatrix(peptideMatrix, outputDirectory, PeptideFileName, "Peptide", p => p.FullSequence);

        /// <summary>
        /// Writes the protein-group matrix: one row per group, one column per sample.
        /// </summary>
        /// <returns>The path written, or null if the matrix was empty.</returns>
        public static string WriteProteinGroupMatrix(QuantMatrix<IBioPolymerGroup> proteinMatrix, string outputDirectory) =>
            WriteMatrix(proteinMatrix, outputDirectory, ProteinGroupFileName, "Protein Group", g => g.BioPolymerGroupName);

        /// <summary>
        /// Shared row-per-entity, column-per-sample writer. Column headers come from
        /// <see cref="SampleColumnLabel"/> and are made unique by suffixing any repeats, so a
        /// malformed experimental design produces an awkward header rather than silently
        /// overwritten columns.
        /// </summary>
        private static string WriteMatrix<T>(
            QuantMatrix<T> matrix,
            string outputDirectory,
            string fileName,
            string rowHeader,
            Func<T, string> rowLabel) where T : IEquatable<T>
        {
            if (matrix == null || matrix.RowCount == 0)
                return null;

            string path = ResolvePath(outputDirectory, fileName);
            using var output = new StreamWriter(path);

            output.WriteLine(rowHeader + "\t" + string.Join("\t", UniqueColumnLabels(matrix.ColumnKeys)));

            var line = new StringBuilder();
            for (int row = 0; row < matrix.RowCount; row++)
            {
                line.Clear();
                line.Append(rowLabel(matrix.RowKeys[row]));
                for (int col = 0; col < matrix.ColumnCount; col++)
                    line.Append('\t').Append(Format(matrix.Matrix[row, col]));
                output.WriteLine(line.ToString());
            }

            return path;
        }

        /// <summary>
        /// The column header for one sample. Isobaric samples are labelled
        /// <c>{file}_{channel}</c>; label-free samples by file name. Falls back to
        /// <see cref="object.ToString"/> for other <see cref="ISampleInfo"/> implementations.
        /// </summary>
        internal static string SampleColumnLabel(ISampleInfo sample)
        {
            if (sample is IsobaricQuantSampleInfo isobaric)
                return isobaric.ToString();

            string name = sample.FilenameWithoutExtension;
            return string.IsNullOrEmpty(name) ? sample.ToString() : name;
        }

        /// <summary>
        /// Produces one label per column, appending "_2", "_3", … to any repeated label so that every
        /// column header is distinct.
        /// </summary>
        internal static List<string> UniqueColumnLabels(IEnumerable<ISampleInfo> samples)
        {
            var seen = new Dictionary<string, int>(StringComparer.Ordinal);
            var labels = new List<string>();

            foreach (var sample in samples)
            {
                string label = SampleColumnLabel(sample);
                if (seen.TryGetValue(label, out int count))
                {
                    seen[label] = count + 1;
                    label = label + "_" + (count + 1).ToString(CultureInfo.InvariantCulture);
                }
                else
                {
                    seen[label] = 1;
                }
                labels.Add(label);
            }

            return labels;
        }

        private static string ResolvePath(string outputDirectory, string fileName)
        {
            if (string.IsNullOrEmpty(outputDirectory))
                return fileName;

            Directory.CreateDirectory(outputDirectory);
            return Path.Combine(outputDirectory, fileName);
        }

        private static string Format(double value) => value.ToString("G17", CultureInfo.InvariantCulture);
    }
}
