using System.Collections.Generic;
using System.Linq;

namespace Quantification
{
    /// <summary>
    /// Which spectral matches are quantified. One rule, so that MetaMorpheus and FlashLFQ quantify the same matches
    /// from the same search:
    /// <list type="bullet">
    /// <item>the PEP q-value is below the threshold, when PEP was trained for the search;</item>
    /// <item>otherwise the q-value AND the notch q-value are below it (the notch only where one is reported);</item>
    /// <item>strictly below: a match at exactly the threshold is not quantified;</item>
    /// <item>never an ambiguous match, one that names more than one sequence or modified form.</item>
    /// </list>
    /// Whether PEP was trained belongs to the whole search, not to one match, so a caller decides it once with
    /// <see cref="PepQValueIsUsable"/> and passes the answer to <see cref="PassesConfidence"/> for every match.
    /// </summary>
    public static class QuantifiedPsmRule
    {
        public const double DefaultThreshold = 0.01;

        /// <summary>
        /// True when a search's PEP q-values can be filtered on: at least one is a real value in [0, 1]. MetaMorpheus
        /// writes 2 for every match when PEP was not trained, and a reader fills NaN when the column is absent; both
        /// mean there is no PEP tier, and filtering on them would drop every match.
        /// </summary>
        public static bool PepQValueIsUsable(IEnumerable<double> pepQValues) => pepQValues.Any(v => v >= 0 && v <= 1);

        /// <summary>
        /// Whether a match is confident enough to quantify. With <paramref name="usePepQValue"/> only the PEP q-value
        /// counts; without it the q-value and, when reported, the notch q-value must both be below the threshold.
        /// </summary>
        public static bool PassesConfidence(double qValue, double? notchQValue, double? pepQValue, bool usePepQValue,
            double threshold = DefaultThreshold)
        {
            if (usePepQValue)
            {
                return pepQValue is double pep && pep < threshold;
            }

            return qValue < threshold && (notchQValue is not double notch || notch < threshold);
        }

        /// <summary>
        /// True when a match does not name one sequence and one modified form: either is missing, or either lists
        /// alternatives joined with '|', as MetaMorpheus writes them.
        /// </summary>
        public static bool IsAmbiguous(string baseSequence, string fullSequence) =>
            string.IsNullOrEmpty(baseSequence) || string.IsNullOrEmpty(fullSequence)
            || baseSequence.Contains('|') || fullSequence.Contains('|');
    }
}
