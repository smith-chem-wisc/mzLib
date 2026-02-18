using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;

namespace Readers.InternalIons
{
    public static class InternalFragmentTsvWriter
    {
        private const char Delimiter = '\t';

        public static void WriteToTsv(List<InternalFragmentIon> ions, string outputPath)
        {
            if (ions == null)
                throw new ArgumentNullException(nameof(ions));
            if (string.IsNullOrWhiteSpace(outputPath))
                throw new ArgumentException("Output path cannot be null or empty.", nameof(outputPath));

            using var writer = new StreamWriter(outputPath);
            writer.WriteLine(string.Join(Delimiter, InternalFragmentIon.GetHeaderNames()));

            foreach (var ion in ions)
            {
                writer.WriteLine(string.Join(Delimiter, ion.GetValues()));
            }
        }

        public static List<InternalFragmentIon> ReadFromTsv(string path)
        {
            if (!File.Exists(path))
                throw new FileNotFoundException($"File not found: {path}");

            var results = new List<InternalFragmentIon>();
            var lines = File.ReadAllLines(path);

            if (lines.Length == 0)
                return results;

            var headerColumns = lines[0].Split(Delimiter);
            var columnMap = new Dictionary<string, int>();
            for (int i = 0; i < headerColumns.Length; i++)
            {
                columnMap[headerColumns[i]] = i;
            }

            for (int lineIdx = 1; lineIdx < lines.Length; lineIdx++)
            {
                var values = lines[lineIdx].Split(Delimiter);
                if (values.Length == 0)
                    continue;

                try
                {
                    var ion = ParseInternalFragmentIon(values, columnMap);
                    results.Add(ion);
                }
                catch
                {
                    continue;
                }
            }

            return results;
        }

        private static InternalFragmentIon ParseInternalFragmentIon(
            string[] values,
            Dictionary<string, int> columnMap)
        {
            return new InternalFragmentIon
            {
                PeptideSequence = GetValue(values, columnMap, nameof(InternalFragmentIon.PeptideSequence)),
                InternalSequence = GetValue(values, columnMap, nameof(InternalFragmentIon.InternalSequence)),
                StartResidue = ParseInt(GetValue(values, columnMap, nameof(InternalFragmentIon.StartResidue))),
                EndResidue = ParseInt(GetValue(values, columnMap, nameof(InternalFragmentIon.EndResidue))),
                TheoreticalMass = ParseDouble(GetValue(values, columnMap, nameof(InternalFragmentIon.TheoreticalMass))),
                ObservedMass = ParseDouble(GetValue(values, columnMap, nameof(InternalFragmentIon.ObservedMass))),
                NormalizedIntensity = ParseDouble(GetValue(values, columnMap, nameof(InternalFragmentIon.NormalizedIntensity))),
                LocalIntensityRank = ParseDouble(GetValue(values, columnMap, nameof(InternalFragmentIon.LocalIntensityRank))),
                PrecursorCharge = ParseInt(GetValue(values, columnMap, nameof(InternalFragmentIon.PrecursorCharge))),
                CollisionEnergy = ParseDouble(GetValue(values, columnMap, nameof(InternalFragmentIon.CollisionEnergy))),
                PeptidePEP = ParseDouble(GetValue(values, columnMap, nameof(InternalFragmentIon.PeptidePEP))),
                PeptideScore = ParseDouble(GetValue(values, columnMap, nameof(InternalFragmentIon.PeptideScore))),
                NTerminalFlankingResidue = ParseChar(GetValue(values, columnMap, nameof(InternalFragmentIon.NTerminalFlankingResidue))),
                CTerminalFlankingResidue = ParseChar(GetValue(values, columnMap, nameof(InternalFragmentIon.CTerminalFlankingResidue))),
                IsDecoy = ParseBool(GetValue(values, columnMap, nameof(InternalFragmentIon.IsDecoy))),
                SourceFile = GetValue(values, columnMap, nameof(InternalFragmentIon.SourceFile)),
                ScanNumber = GetValue(values, columnMap, nameof(InternalFragmentIon.ScanNumber)),
                IsIsobaricAmbiguous = ParseBool(GetValue(values, columnMap, nameof(InternalFragmentIon.IsIsobaricAmbiguous))),
                FullModifiedSequence = GetValue(values, columnMap, nameof(InternalFragmentIon.FullModifiedSequence)),
                ModificationsInInternalFragment = GetValue(values, columnMap, nameof(InternalFragmentIon.ModificationsInInternalFragment)),
                BIonIntensityAtNTerm = ParseDouble(GetValue(values, columnMap, nameof(InternalFragmentIon.BIonIntensityAtNTerm))),
                YIonIntensityAtCTerm = ParseDouble(GetValue(values, columnMap, nameof(InternalFragmentIon.YIonIntensityAtCTerm))),
                HasMatchedBIonAtNTerm = ParseBool(GetValue(values, columnMap, nameof(InternalFragmentIon.HasMatchedBIonAtNTerm))),
                HasMatchedYIonAtCTerm = ParseBool(GetValue(values, columnMap, nameof(InternalFragmentIon.HasMatchedYIonAtCTerm))),
                BYProductScore = ParseDouble(GetValue(values, columnMap, nameof(InternalFragmentIon.BYProductScore)))
            };
        }

        private static string GetValue(string[] values, Dictionary<string, int> columnMap, string columnName)
        {
            if (!columnMap.TryGetValue(columnName, out int idx) || idx >= values.Length)
                return string.Empty;
            return values[idx];
        }

        private static int ParseInt(string value) =>
            int.TryParse(value, out int result) ? result : 0;

        private static double ParseDouble(string value)
        {
            if (string.Equals(value, "NaN", StringComparison.OrdinalIgnoreCase))
                return double.NaN;
            return double.TryParse(value, NumberStyles.Any, CultureInfo.InvariantCulture, out double result)
                ? result : double.NaN;
        }

        private static char ParseChar(string value) =>
            string.IsNullOrEmpty(value) ? '-' : value[0];

        private static bool ParseBool(string value) =>
            bool.TryParse(value, out bool result) && result;
    }
}