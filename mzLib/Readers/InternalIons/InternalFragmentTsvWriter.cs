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
            if (ions == null) throw new ArgumentNullException(nameof(ions));
            if (string.IsNullOrWhiteSpace(outputPath))
                throw new ArgumentException("Output path cannot be null or empty.", nameof(outputPath));

            using var writer = new StreamWriter(outputPath);
            writer.WriteLine(string.Join(Delimiter, InternalFragmentIon.GetHeaderNames()));
            foreach (var ion in ions)
                writer.WriteLine(string.Join(Delimiter, ion.GetValues()));
        }

        public static List<InternalFragmentIon> ReadFromTsv(string path)
        {
            if (!File.Exists(path)) throw new FileNotFoundException($"File not found: {path}");
            var results = new List<InternalFragmentIon>();
            var lines = File.ReadAllLines(path);
            if (lines.Length == 0) return results;

            var columnMap = new Dictionary<string, int>();
            var headers = lines[0].Split(Delimiter);
            for (int i = 0; i < headers.Length; i++) columnMap[headers[i]] = i;

            for (int idx = 1; idx < lines.Length; idx++)
            {
                var vals = lines[idx].Split(Delimiter);
                if (vals.Length == 0) continue;
                try { results.Add(Parse(vals, columnMap)); } catch { }
            }
            return results;
        }

        private static InternalFragmentIon Parse(string[] v, Dictionary<string, int> m)
        {
            return new InternalFragmentIon
            {
                PeptideSequence = Get(v, m, nameof(InternalFragmentIon.PeptideSequence)),
                InternalSequence = Get(v, m, nameof(InternalFragmentIon.InternalSequence)),
                StartResidue = Int(Get(v, m, nameof(InternalFragmentIon.StartResidue))),
                EndResidue = Int(Get(v, m, nameof(InternalFragmentIon.EndResidue))),
                TheoreticalMass = Dbl(Get(v, m, nameof(InternalFragmentIon.TheoreticalMass))),
                ObservedMass = Dbl(Get(v, m, nameof(InternalFragmentIon.ObservedMass))),
                NormalizedIntensity = Dbl(Get(v, m, nameof(InternalFragmentIon.NormalizedIntensity))),
                TicNormalizedIntensity = Dbl(Get(v, m, nameof(InternalFragmentIon.TicNormalizedIntensity))),
                TotalIonCurrent = Dbl(Get(v, m, nameof(InternalFragmentIon.TotalIonCurrent))),
                LocalIntensityRank = Dbl(Get(v, m, nameof(InternalFragmentIon.LocalIntensityRank))),
                PrecursorCharge = Int(Get(v, m, nameof(InternalFragmentIon.PrecursorCharge))),
                CollisionEnergy = Dbl(Get(v, m, nameof(InternalFragmentIon.CollisionEnergy))),
                PeptidePEP = Dbl(Get(v, m, nameof(InternalFragmentIon.PeptidePEP))),
                PeptideScore = Dbl(Get(v, m, nameof(InternalFragmentIon.PeptideScore))),
                NTerminalFlankingResidue = Chr(Get(v, m, nameof(InternalFragmentIon.NTerminalFlankingResidue))),
                CTerminalFlankingResidue = Chr(Get(v, m, nameof(InternalFragmentIon.CTerminalFlankingResidue))),
                IsDecoy = Bool(Get(v, m, nameof(InternalFragmentIon.IsDecoy))),
                SourceFile = Get(v, m, nameof(InternalFragmentIon.SourceFile)),
                ScanNumber = Get(v, m, nameof(InternalFragmentIon.ScanNumber)),
                IsIsobaricAmbiguous = Bool(Get(v, m, nameof(InternalFragmentIon.IsIsobaricAmbiguous))),
                FullModifiedSequence = Get(v, m, nameof(InternalFragmentIon.FullModifiedSequence)),
                ModificationsInInternalFragment = Get(v, m, nameof(InternalFragmentIon.ModificationsInInternalFragment)),
                BIonIntensityAtNTerm = Dbl(Get(v, m, nameof(InternalFragmentIon.BIonIntensityAtNTerm))),
                YIonIntensityAtCTerm = Dbl(Get(v, m, nameof(InternalFragmentIon.YIonIntensityAtCTerm))),
                HasMatchedBIonAtNTerm = Bool(Get(v, m, nameof(InternalFragmentIon.HasMatchedBIonAtNTerm))),
                HasMatchedYIonAtCTerm = Bool(Get(v, m, nameof(InternalFragmentIon.HasMatchedYIonAtCTerm))),
                BYProductScore = Dbl(Get(v, m, nameof(InternalFragmentIon.BYProductScore))),
                DistanceFromCTerm = Int(Get(v, m, nameof(InternalFragmentIon.DistanceFromCTerm))),
                MaxTerminalIonIntensity = Dbl(Get(v, m, nameof(InternalFragmentIon.MaxTerminalIonIntensity))),
                HasBothTerminalIons = Bool(Get(v, m, nameof(InternalFragmentIon.HasBothTerminalIons))),
                BasicResiduesInBIonSpan = Int(Get(v, m, nameof(InternalFragmentIon.BasicResiduesInBIonSpan))),
                BasicResiduesInYIonSpan = Int(Get(v, m, nameof(InternalFragmentIon.BasicResiduesInYIonSpan)))
            };
        }

        private static string Get(string[] v, Dictionary<string, int> m, string col) =>
            m.TryGetValue(col, out int i) && i < v.Length ? v[i] : "";
        private static int Int(string s) => int.TryParse(s, out int r) ? r : 0;
        private static double Dbl(string s) => s.Equals("NaN", StringComparison.OrdinalIgnoreCase) ? double.NaN :
            double.TryParse(s, NumberStyles.Any, CultureInfo.InvariantCulture, out double r) ? r : double.NaN;
        private static char Chr(string s) => string.IsNullOrEmpty(s) ? '-' : s[0];
        private static bool Bool(string s) => bool.TryParse(s, out bool r) && r;
    }
}