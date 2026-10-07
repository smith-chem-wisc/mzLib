using System.Collections.Generic;
using System.Linq;

namespace Quantification
{
    /// <summary>
    /// The confidence value <see cref="QuantifiedPsmRule"/> filters a search's matches on. One tier is chosen for the
    /// whole search, the first that the search can support, in this order.
    /// </summary>
    public enum QuantifiedPsmTier
    {
        /// <summary>The PEP q-value, when PEP was trained for the search.</summary>
        PepQValue,

        /// <summary>The notch q-value, when PEP was not trained but the search reports notch q-values.</summary>
        QValueNotch,

        /// <summary>The q-value, otherwise.</summary>
        QValue,
    }

    /// <summary>
    /// Which spectral matches are quantified. One rule, so that every caller (MetaMorpheus, FlashLFQ) follows the same
    /// order of filters, and falls back the same way when a value is missing:
    /// <list type="bullet">
    /// <item>one confidence value decides, chosen once per search as a fallback (<see cref="ChooseTier"/>): the PEP
    /// q-value when PEP was trained, otherwise the notch q-value when the search reports one, otherwise the q-value;</item>
    /// <item>strictly below the threshold: a match at exactly the threshold is not quantified;</item>
    /// <item>never an ambiguous match, one that names more than one sequence or modified form.</item>
    /// </list>
    /// Which tier applies belongs to the whole search, not to one match, so a caller decides it once with
    /// <see cref="ChooseTier"/> and passes the answer to <see cref="PassesConfidence"/> for every match.
    /// </summary>
    public static class QuantifiedPsmRule
    {
        public const double DefaultThreshold = 0.01;

        /// <summary>
        /// The tier a search's matches are filtered on: <see cref="QuantifiedPsmTier.PepQValue"/> when
        /// <see cref="PepIsUsable"/>, else <see cref="QuantifiedPsmTier.QValueNotch"/> when <see cref="NotchIsUsable"/>,
        /// else <see cref="QuantifiedPsmTier.QValue"/>. Pass the values of the search's target matches.
        /// </summary>
        public static QuantifiedPsmTier ChooseTier(IEnumerable<double> pepQValues, IEnumerable<double> notchQValues,
            IEnumerable<double>? peps = null)
        {
            if (PepIsUsable(pepQValues, peps))
            {
                return QuantifiedPsmTier.PepQValue;
            }

            return NotchIsUsable(notchQValues) ? QuantifiedPsmTier.QValueNotch : QuantifiedPsmTier.QValue;
        }

        /// <summary>
        /// True when a search's PEP was trained, so its PEP q-values can be filtered on:
        /// <list type="bullet">
        /// <item>at least one PEP q-value is a real value in [0, 1]. MetaMorpheus writes 2 for every match when PEP was
        /// not trained, and a reader fills NaN when the column is absent;</item>
        /// <item>and, when <paramref name="peps"/> are given, they are not all the same value. When PEP training fails,
        /// MetaMorpheus leaves PEP at 0 for every match and still writes PEP q-values, ordered by score, which look
        /// usable but are not. A search with a single PEP value cannot show that it was trained, so it is not usable
        /// either. Without <paramref name="peps"/> (a reader that does not expose PEP), this failure cannot be
        /// seen.</item>
        /// </list>
        /// </summary>
        public static bool PepIsUsable(IEnumerable<double> pepQValues, IEnumerable<double>? peps = null)
        {
            if (!pepQValues.Any(IsProbability))
            {
                return false;
            }

            if (peps is null)
            {
                return true;
            }

            List<double> realPeps = peps.Where(double.IsFinite).ToList();
            return realPeps.Count == 0 || realPeps.Distinct().Skip(1).Any();
        }

        /// <summary>True when a search reports notch q-values: at least one is a real value in [0, 1].</summary>
        public static bool NotchIsUsable(IEnumerable<double> notchQValues) => notchQValues.Any(IsProbability);

        /// <summary>
        /// Whether a match is confident enough to quantify: the value <paramref name="tier"/> names is strictly below
        /// the threshold. A match with no value for that tier is not quantified. The other values are ignored.
        /// </summary>
        public static bool PassesConfidence(QuantifiedPsmTier tier, double qValue, double? notchQValue, double? pepQValue,
            double threshold = DefaultThreshold)
        {
            double? value = tier switch
            {
                QuantifiedPsmTier.PepQValue => pepQValue,
                QuantifiedPsmTier.QValueNotch => notchQValue,
                _ => qValue,
            };

            return value is double v && v < threshold;
        }

        /// <summary>
        /// True when a match does not name one sequence and one modified form: either is missing, or either lists
        /// alternatives joined with '|', as MetaMorpheus writes them.
        /// </summary>
        public static bool IsAmbiguous(string baseSequence, string fullSequence) =>
            string.IsNullOrEmpty(baseSequence) || string.IsNullOrEmpty(fullSequence)
            || baseSequence.Contains('|') || fullSequence.Contains('|');

        private static bool IsProbability(double value) => value >= 0 && value <= 1;
    }
}
