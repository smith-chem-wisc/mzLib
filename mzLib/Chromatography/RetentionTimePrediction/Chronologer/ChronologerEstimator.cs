using System.Text.RegularExpressions;
using TorchSharp;
using static TorchSharp.torch;

namespace Chromatography.RetentionTimePrediction.Chronologer
{
    public static class ChronologerEstimator
    {
        private static Chronologer ChronologerModel = new Chronologer();

        private const int MaxPepLen = 50; 

        private static readonly char[] Residues = (
            "ACDEFGHIKLMNPQRSTVWY" +     // 20 canonical (1..20)
            "cmdestyabunopqrxz" +      // 17 modified (21..37)
            "-^()&*_" +                   // 7 N/C states (38..44)
            "0123456789"                  // 10 user slots (45..54)
        ).ToCharArray();

        private static readonly Dictionary<char, int> CodeToInt =
            Residues.Select((c, i) => (c, i + 1)).ToDictionary(t => t.c, t => t.Item2);

        // 2) Regex map exactly like Python
        private static readonly (string pattern, string repl)[] ModRegex = new (string, string)[] {
            (@"M\[\+15\.99.{0,6}\]", "m"),
            (@"C\[\+57\.02.{0,6}\]", "c"),
            (@"C\[\+39\.99.{0,6}\]", "d"),
            (@"\[\-18\.01.{0,6}\]E", "e"), (@"E\[\-18\.01.{0,6}\]", "e"),
            (@"\[\-17\.02.{0,6}\]Q", "e"), (@"Q\[\-17\.02.{0,6}\]", "e"),
            (@"S\[\+79\.96.{0,6}\]", "s"), (@"T\[\+79\.96.{0,6}\]", "t"), (@"Y\[\+79\.96.{0,6}\]", "y"),
            (@"K\[\+42\.01.{0,6}\]", "a"), (@"K\[\+100\.0.{0,6}\]", "b"), (@"K\[\+114\.0.{0,6}\]", "u"),
            (@"K\[\+14\.01.{0,6}\]", "n"), (@"K\[\+28\.03.{0,6}\]", "o"), (@"K\[\+42\.04.{0,6}\]", "p"),
            (@"R\[\+14\.01.{0,6}\]", "q"), (@"R\[\+28\.03.{0,6}\]", "r"),
            (@"K\[\+224\.1.{0,6}\]", "z"), (@"K\[\+229\.1.{0,6}\]", "x"),
        };

        private static readonly string[] ModsToSanitize =
        [
            "[+43.005814]-", // Carbamyl- N-Term
            "[+43.005814]", // Carbamyl
            "[+0.984016]", // Deamidation
        ];

        private static readonly Dictionary<string, char> NTermKeys = new()
        {
            { "+42.01", '^' }, { "+224.1", '&' }, { "+229.1", '*' }
        };

        /// <summary>
        /// Uses the Chronologer model to predict C18 retention times (reported in % ACN).
        /// Only modifications present in the Chronologer dictionary are supported.
        /// Returns null if the sequence is not valid.
        /// <code>
        /// "Carbamidomethyl on C"
        /// "Oxidation on M"
        /// "Glu to PyroGlu"
        /// "Phosphorylation on S"
        /// "Phosphorylation on T"
        /// "Phosphorylation on Y"
        /// "Accetylation on K"
        /// "Succinylation on K"
        /// "Ubiquitination on K"
        /// "Methylation on K"
        /// "Dimethylation on K"
        /// "Trimethylation on K"
        /// "Methylation on R"
        /// "Dimethylation on R"
        /// </code>
        /// </summary>
        /// <param name="baseSequence"></param>
        /// <param name="fullSequence"></param>
        /// <returns></returns>
        public static double? PredictRetentionTime(string baseSequence, string fullSequence)
        {
            if (fullSequence.Contains("[Metal"))
                return null;
            if (baseSequence.Contains('U'))
                return null;

            var tensor = Tensorize(fullSequence);
            if (tensor is null)
                return null;

            var prediction = ChronologerModel.Predict(tensor);
            return prediction[0].ToDouble();
        }

        public static torch.Tensor Tensorize(string fullSequence)
        {

            // Step A: substitute mods → single-letter codes
            string seq = fullSequence;
            foreach (var (pattern, repl) in ModRegex)
                seq = Regex.Replace(seq, pattern, repl);

            // Step A2: Sanitize mods that are unsupported THIS IS NOT GOOD CODE
            if (ModsToSanitize.Any(p => seq.Contains(p)))
                foreach (var mod in ModsToSanitize)
                    if (seq.Contains(mod))
                        seq = seq.Replace(mod, "");

            // Step B: N-term state prefix
            if (seq.Length == 0) 
                return null!;
            if (seq[0] == 'd') // pyroGlu at first position
                seq = ")" + seq;         
            else if (seq[0] == 'e') // cyclized CAM-Cys at first
                seq = "(" + seq;         
            else if (seq[0] == '[') // raw N-term mass present
            {
                // grab [+xx.xx]
                int close = seq.IndexOf(']');
                if (close < 0 || close < 7) 
                    return null!;
                string key = seq.Substring(1, 6);
                if (!NTermKeys.TryGetValue(key, out char nterm)) 
                    return null!;
                seq = nterm + seq.Substring(close + 1);
            }
            else
            {
                seq = "-" + seq; // free N-term
            }

            // Step C: append free C-term
            seq += "_";

            // Step D: reject if any brackets remain (means unhandled mod)
            if (seq.IndexOf('[') >= 0) 
                return null!;

            // Step E: length guard (without padding, excluding padding zeros)
            // Python allows length up to MaxPepLen+2 (N & C tokens included)
            if (seq.Length > MaxPepLen + 2) 
                return null!;

            // Step F: map chars → ints
            var ids = new long[MaxPepLen + 2]; // zero-padded (0) to the right
            for (int i = 0; i < seq.Length; i++)
            {
                if (!CodeToInt.TryGetValue(seq[i], out int v)) 
                    return null!;
                ids[i] = v;
            }

            // Output shape: [1, MaxPepLen+2], dtype int64
            var t = torch.tensor(ids, dtype: ScalarType.Int64).reshape(1, MaxPepLen + 2);
            return t;
        }

    }
}
