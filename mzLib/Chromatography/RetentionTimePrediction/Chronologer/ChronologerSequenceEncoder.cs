using System.Text.RegularExpressions;
using TorchSharp;
using static TorchSharp.torch;

namespace Chromatography.RetentionTimePrediction.Chronologer
{
    /// <summary>
    /// Encodes peptide sequences with mass shifts into Chronologer's internal representation.
    /// This is a stateless, thread-safe encoder that performs the transformation from
    /// mass-shift annotated sequences to the character-based format required by Chronologer.
    /// </summary>
    public static class ChronologerSequenceEncoder
    {
        private const int MaxPeptideLength = 50;
        private const int EncodedLength = MaxPeptideLength + 2; // +2 for N/C termini tokens

        // Chronologer alphabet
        // 20 canonical (1..20) + 17 modified (21..37) + 7 N/C states (38..44) + 10 user slots (45..54)
        private static readonly char[] Residues = (
            "ACDEFGHIKLMNPQRSTVWY" +     // 1-20: canonical amino acids
            "cmdestyabunopqrxz" +         // 21-37: modified residues
            "-^()&*_" +                   // 38-44: N/C terminus states
            "0123456789"                  // 45-54: user-defined slots
        ).ToCharArray();

        private static readonly Dictionary<char, int> CodeToInt =
            Residues.Select((c, i) => (c, i + 1)).ToDictionary(t => t.c, t => t.Item2);

        // Compiled regex patterns for performance
        private static readonly (Regex pattern, string replacement)[] ModificationPatterns = new[]
        {
            (new Regex(@"M\[\+15\.99\d*\]", RegexOptions.Compiled), "m"),  // Oxidation on M
            (new Regex(@"C\[\+57\.02\d*\]", RegexOptions.Compiled), "c"),  // Carbamidomethyl on C
            (new Regex(@"C\[\+39\.99\d*\]", RegexOptions.Compiled), "d"),  // Alternative C mod
            (new Regex(@"\[\-18\.01\d*\]E", RegexOptions.Compiled), "e"), // PyroGlu from E (prefix)
            (new Regex(@"E\[\-18\.01\d*\]", RegexOptions.Compiled), "e"),  // PyroGlu from E (suffix)
            (new Regex(@"\[\-17\.02\d*\]Q", RegexOptions.Compiled), "e"), // PyroGlu from Q (prefix)
            (new Regex(@"Q\[\-17\.02\d*\]", RegexOptions.Compiled), "e"),  // PyroGlu from Q (suffix)
            (new Regex(@"S\[\+79\.96\d*\]", RegexOptions.Compiled), "s"),  // Phosphorylation on S
            (new Regex(@"T\[\+79\.96\d*\]", RegexOptions.Compiled), "t"),  // Phosphorylation on T
            (new Regex(@"Y\[\+79\.96\d*\]", RegexOptions.Compiled), "y"),  // Phosphorylation on Y
            (new Regex(@"K\[\+42\.01\d*\]", RegexOptions.Compiled), "a"),  // Acetylation on K
            (new Regex(@"K\[\+100\.0\d*\]", RegexOptions.Compiled), "b"),  // Succinylation on K
            (new Regex(@"K\[\+114\.0\d*\]", RegexOptions.Compiled), "u"),  // Ubiquitination on K
            (new Regex(@"K\[\+14\.01\d*\]", RegexOptions.Compiled), "n"),  // Methylation on K
            (new Regex(@"K\[\+28\.03\d*\]", RegexOptions.Compiled), "o"),  // Dimethylation on K
            (new Regex(@"K\[\+42\.04\d*\]", RegexOptions.Compiled), "p"),  // Trimethylation on K
            (new Regex(@"R\[\+14\.01\d*\]", RegexOptions.Compiled), "q"),  // Methylation on R
            (new Regex(@"R\[\+28\.03\d*\]", RegexOptions.Compiled), "r"),  // Dimethylation on R
            (new Regex(@"K\[\+224\.1\d*\]", RegexOptions.Compiled), "z"),  // GlyGly on K
            (new Regex(@"K\[\+229\.1\d*\]", RegexOptions.Compiled), "x"),  // Heavy GlyGly on K
        };

        // N-terminus modification codes
        private static readonly Dictionary<string, char> NTerminusCodes = new()
        {
            { "+42.01", '^' },  // N-term acetylation
            { "+224.1", '&' },  // N-term GlyGly
            { "+229.1", '*' }   // N-term heavy GlyGly
        };

        /// <summary>
        /// Encodes a peptide sequence with mass shifts into a tensor for Chronologer prediction.
        /// Returns null if the sequence cannot be encoded.
        /// </summary>
        /// <param name="sequenceWithMassShifts">Sequence with mass shift annotations, e.g., "PEPTIDE[+15.995]K"</param>
        /// <returns>Encoded tensor of shape (1, 52), or null if sequence cannot be encoded</returns>
        public static Tensor? EncodeTensor(string sequenceWithMassShifts)
        {
            if (string.IsNullOrEmpty(sequenceWithMassShifts))
                return null;

            string? encoded = EncodeToChronologerFormat(sequenceWithMassShifts);
            if (encoded == null)
                return null;

            return ConvertToTensor(encoded);
        }

        /// <summary>
        /// Converts a sequence with mass shifts to Chronologer's character-based format.
        /// Returns null if conversion fails.
        /// </summary>
        public static string? EncodeToChronologerFormat(string sequenceWithMassShifts)
        {
            if (string.IsNullOrEmpty(sequenceWithMassShifts))
                return null;

            // Step 1: Apply modification pattern replacements
            string seq = sequenceWithMassShifts;
            foreach (var (pattern, replacement) in ModificationPatterns)
            {
                seq = pattern.Replace(seq, replacement);
            }

            if (seq.Length == 0)
                return null;

            // Step 2: Add N-terminus token
            seq = AddNTerminusToken(seq);
            if (seq == null)
                return null;

            // Step 3: Add C-terminus token
            seq += "_";

            // Step 4: Validate - no unhandled modifications should remain
            if (seq.Contains('['))
                return null; // Unhandled modification present

            // Step 5: Length check
            if (seq.Length > EncodedLength)
                return null;

            return seq;
        }

        /// <summary>
        /// Adds appropriate N-terminus token based on modifications or default state.
        /// </summary>
        private static string? AddNTerminusToken(string seq)
        {
            // Check for PyroGlu at N-terminus
            if (seq[0] == 'd') // pyroGlu at first position
                return ")" + seq;

            if (seq[0] == 'e') // cyclized CAM-Cys at first
                return "(" + seq;

            // Check for N-terminal mass modification
            if (seq[0] == '[')
            {
                // grab [+xx.xx]
                int close = seq.IndexOf(']');
                if (close < 0 || close < 7)
                    return null; // Invalid format

                string key = seq.Substring(1, 6);
                if (!NTerminusCodes.TryGetValue(key, out char nterm))
                    return null; // Unsupported N-terminal modification

                return nterm + seq.Substring(close + 1);
            }

            // Free N-terminus
            return "-" + seq;
        }

        /// <summary>
        /// Converts Chronologer-formatted sequence string to tensor.
        /// </summary>
        private static Tensor? ConvertToTensor(string chronologerSequence)
        {
            var ids = new long[EncodedLength]; // Zero-padded

            for (int i = 0; i < chronologerSequence.Length; i++)
            {
                if (!CodeToInt.TryGetValue(chronologerSequence[i], out int v))
                    return null; // Invalid character

                ids[i] = v;
            }

            // Output shape: [1, MaxPepLen+2], dtype int64
            return tensor(ids, dtype: ScalarType.Int64).reshape(1, EncodedLength);
        }

        /// <summary>
        /// Validates if a sequence can be encoded by Chronologer.
        /// Provides detailed failure reason for diagnostics.
        /// </summary>
        public static bool CanEncode(string baseSequence, string? sequenceWithMassShifts, out string? failureReason)
        {
            // Basic checks
            if (string.IsNullOrEmpty(baseSequence))
            {
                failureReason = "Empty sequence";
                return false;
            }

            if (baseSequence.Contains('U'))
            {
                failureReason = "Chronologer does not support selenocysteine (U)";
                return false;
            }

            if (baseSequence.Length > MaxPeptideLength)
            {
                failureReason = $"Sequence length ({baseSequence.Length}) exceeds maximum ({MaxPeptideLength})";
                return false;
            }

            // Try encoding to check for unsupported modifications
            if (!string.IsNullOrEmpty(sequenceWithMassShifts))
            {
                string? encoded = EncodeToChronologerFormat(sequenceWithMassShifts);
                if (encoded == null)
                {
                    failureReason = "Sequence contains unsupported modifications or invalid format";
                    return false;
                }
            }

            failureReason = null;
            return true;
        }
    }
}