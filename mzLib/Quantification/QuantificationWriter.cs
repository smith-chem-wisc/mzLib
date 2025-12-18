using Omics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Omics.BioPolymerGroup;

namespace Quantification
{
    /// <summary>
    /// This will need to be implemented or changed or relocated later. Right now it's just a placeholder.
    /// </summary>
    public static class QuantificationWriter
    {
        private const char Delimiter = '\t';

        public static void WriteRawData(List<ISpectralMatch> matches, string outputDirectory)
        {
            if (matches == null || matches.Count == 0)
                return;

            string outputPath = Path.Combine(outputDirectory, "RawQuantification.tsv");

            using var writer = new StreamWriter(outputPath);

            // Determine max reporter channels from data (for TMT support)
            int maxChannels = matches
                .Where(m => m.QuantValues != null)
                .Select(m => m.QuantValues!.Length)
                .DefaultIfEmpty(1)
                .Max();

            // Write header
            var header = new StringBuilder();
            header.Append("FileName");
            header.Append(Delimiter);
            header.Append("Ms2ScanNumber");
            header.Append(Delimiter);
            header.Append("FullSequence");
            header.Append(Delimiter);
            header.Append("BaseSequence");

            if (maxChannels == 1)
            {
                header.Append(Delimiter);
                header.Append("Intensity");
            }
            else
            {
                // TMT reporters
                for (int i = 0; i < maxChannels; i++)
                {
                    header.Append(Delimiter);
                    header.Append($"Reporter_{i + 1}"); // Placeholder naming - will be updated with actual TMT channel names
                }
            }

            writer.WriteLine(header);

            // Write data rows
            foreach (var match in matches.OrderBy(m => m.FullFilePath).ThenBy(m => m.OneBasedScanNumber))
            {
                var row = new StringBuilder();
                row.Append(Path.GetFileNameWithoutExtension(match.FullFilePath));
                row.Append(Delimiter);
                row.Append(match.OneBasedScanNumber);
                row.Append(Delimiter);
                row.Append(match.FullSequence);
                row.Append(Delimiter);
                row.Append(match.BaseSequence);

                if (match.QuantValues == null || match.QuantValues.Length == 0)
                {
                    // No quant data - write empty columns
                    for (int i = 0; i < maxChannels; i++)
                    {
                        row.Append(Delimiter);
                    }
                }
                else
                {
                    foreach (var intensity in match.QuantValues)
                    {
                        row.Append(Delimiter);
                        row.Append(intensity);
                    }
                    // Pad remaining columns if this PSM has fewer channels
                    for (int i = match.QuantValues.Length; i < maxChannels; i++)
                    {
                        row.Append(Delimiter);
                    }
                }

                writer.WriteLine(row);
            }
        }

        public static void WritePeptideMatrix(PeptideMatrix peptideMatrix, string outputDirectory)
        {
            throw new NotImplementedException();
        }

        public static void WriteProteinGroupMatrix(ProteinMatrix proteinMatrix, string outputDirectory)
        {
            throw new NotImplementedException();
        }
    }
}