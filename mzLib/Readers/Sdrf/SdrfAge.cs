using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Text.RegularExpressions;

namespace Readers
{
    /// <summary>
    /// How much of an age an <see cref="SdrfAge"/> cell actually pins down.
    /// </summary>
    public enum SdrfAgePrecision
    {
        /// <summary>One age: <c>58Y</c>, <c>30Y6M</c>, <c>1 week</c>.</summary>
        Exact,

        /// <summary>
        /// Both ends known: <c>40Y-85Y</c>, <c>6-8 weeks</c>. <see cref="SdrfAge.Years"/> is the
        /// midpoint, which for a wide cohort range is a poor stand-in for any one sample.
        /// </summary>
        Range,

        /// <summary>
        /// At least this old: <c>&gt;=90Y</c>, the de-identification cap most human cohorts use.
        /// <see cref="SdrfAge.Years"/> is the bound itself, and <see cref="SdrfAge.MaxYears"/> is
        /// infinite.
        /// </summary>
        LowerBound,

        /// <summary>At most this old: <c>&lt;1Y</c>. <see cref="SdrfAge.MinYears"/> is zero.</summary>
        UpperBound
    }

    /// <summary>
    /// A <c>characteristics[age]</c> cell read into years.
    ///
    /// The confidence is two facts rather than one score, because they fail independently.
    /// <see cref="Precision"/> says how much of the age the cell pins down. <see cref="FollowsSpecification"/>
    /// says whether the cell was written in the spec's own grammar (<c>nYnMnD</c>, <c>nW</c>, ranges
    /// and bounds of those) or in words whose unit is still unambiguous (<c>3 year</c>,
    /// <c>6-8 weeks</c>, <c>4 hour</c>). A caller that only trusts the spec, or only point ages,
    /// filters on the one it cares about.
    ///
    /// What is refused matters more than what is read. A bare number -- <c>63</c>, 11% of the real age
    /// cells in the curated corpus -- has no unit, and a unit cannot be recovered from the cell: 63
    /// years and 63 days are both plausible in one study. So <see cref="TryParse"/> returns false for
    /// it rather than assume years, and the same for reserved words, free text, and anything else it
    /// cannot read without guessing.
    /// </summary>
    /// <param name="Years">
    /// The single figure to place the sample on an age axis: the age, a range's midpoint, or a bound.
    /// </param>
    /// <param name="MinYears">The youngest the cell allows; zero for an upper bound.</param>
    /// <param name="MaxYears">The oldest the cell allows; positive infinity for a lower bound.</param>
    public sealed record SdrfAge(
        double Years,
        double MinYears,
        double MaxYears,
        SdrfAgePrecision Precision,
        bool FollowsSpecification)
    {
        private const double DaysPerYear = 365.25;

        // The spec's grammar: one or more number+unit parts, largest unit first, each at most once
        // (58Y, 30Y6M, 1Y2M3D, 16W). Case-insensitive, because "54y" is a slip of the shift key, not
        // an ambiguity.
        private static readonly Regex SpecPart =
            new(@"(?<n>\d+(?:\.\d+)?)(?<u>[YMWD])", RegexOptions.IgnoreCase | RegexOptions.Compiled);

        // Words whose unit is unambiguous even though the spec would not write them: "3 year",
        // "51 years", "1 week", "4 hour", "6-8 weeks". Hours appear only in developmental time
        // courses (0-2 hour ... 83-85 hour in one Drosophila study).
        private static readonly Regex WordValue = new(
            @"^(?<n>\d+(?:\.\d+)?)\s*(?<u>years?|yrs?|months?|weeks?|wks?|days?|hours?|hrs?)(?:\s+old)?$",
            RegexOptions.IgnoreCase | RegexOptions.Compiled);

        private static readonly Regex Bound =
            new(@"^(?<op>>=|<=|>|<|≥|≤)\s*(?<v>.+)$", RegexOptions.Compiled);

        /// <summary>
        /// Reads one cell. Never throws: false means the cell does not state an age this can read
        /// without guessing, and the caller should treat the sample as having no age.
        /// </summary>
        public static bool TryParse(string? cell, [MaybeNullWhen(false)] out SdrfAge age)
        {
            age = null;
            if (string.IsNullOrWhiteSpace(cell)) return false;
            string text = cell.Trim();

            if (SdrfValidator.ReservedWords.Any(w => string.Equals(text, w, StringComparison.OrdinalIgnoreCase)))
                return false;

            var bound = Bound.Match(text);
            if (bound.Success)
            {
                if (!TryParseValue(bound.Groups["v"].Value, out double years, out bool spec)) return false;
                bool lower = bound.Groups["op"].Value is ">=" or ">" or "≥";
                age = lower
                    ? new SdrfAge(years, years, double.PositiveInfinity, SdrfAgePrecision.LowerBound, spec)
                    : new SdrfAge(years, 0d, years, SdrfAgePrecision.UpperBound, spec);
                return true;
            }

            if (TryParseValue(text, out double exact, out bool exactSpec))
            {
                age = new SdrfAge(exact, exact, exact, SdrfAgePrecision.Exact, exactSpec);
                return true;
            }

            return TryParseRange(text, out age);
        }

        /// <summary>
        /// "40Y-85Y" (spec) or "6-8 weeks" (words, the unit written once and shared). A range whose
        /// left end has no unit borrows the right end's only when the right end is in words; "40-85Y"
        /// is not spec grammar, so it is read as words-with-a-shared-unit, not as spec.
        /// </summary>
        private static bool TryParseRange(string text, [MaybeNullWhen(false)] out SdrfAge age)
        {
            age = null;
            int dash = text.IndexOf('-', 1);
            if (dash < 0) return false;

            string left = text[..dash].Trim();
            string right = text[(dash + 1)..].Trim();
            if (!TryParseValue(right, out double max, out bool rightSpec)) return false;

            double min;
            bool spec;
            if (TryParseValue(left, out min, out bool leftSpec))
            {
                spec = leftSpec && rightSpec;
            }
            else if (double.TryParse(left, NumberStyles.AllowDecimalPoint, CultureInfo.InvariantCulture, out double bare))
            {
                // The left end is a bare number: it takes the right end's unit, scaled the same way.
                string unit = Regex.Match(right, @"[A-Za-z]+(?:\s+old)?$").Value;
                if (!TryParseValue(bare.ToString(CultureInfo.InvariantCulture) + " " + unit, out min, out _)
                    && !TryParseValue(bare.ToString(CultureInfo.InvariantCulture) + unit, out min, out _))
                    return false;
                spec = false;
            }
            else
            {
                return false;
            }

            if (min > max) return false;
            age = new SdrfAge((min + max) / 2d, min, max, SdrfAgePrecision.Range, spec);
            return true;
        }

        /// <summary>One age, in the spec's grammar or in unambiguous words.</summary>
        private static bool TryParseValue(string text, out double years, out bool followsSpecification)
        {
            years = 0d;
            followsSpecification = false;
            text = text.Trim();

            var word = WordValue.Match(text);
            if (word.Success)
            {
                years = ToYears(double.Parse(word.Groups["n"].Value, CultureInfo.InvariantCulture),
                    char.ToUpperInvariant(word.Groups["u"].Value[0]));
                return true;
            }

            // Spec grammar: the parts must cover the whole cell, in Y > M > W > D order, each once.
            // "not availableY" matches nothing; "6M30Y" is out of order; "63" has no unit.
            var parts = SpecPart.Matches(text);
            if (parts.Count == 0 || parts.Sum(p => p.Length) != text.Length) return false;

            const string order = "YMWD";
            int last = -1;
            foreach (Match part in parts)
            {
                int rank = order.IndexOf(char.ToUpperInvariant(part.Groups["u"].Value[0]));
                if (rank <= last) return false;
                last = rank;
                years += ToYears(double.Parse(part.Groups["n"].Value, CultureInfo.InvariantCulture), order[rank]);
            }

            followsSpecification = true;
            return true;
        }

        // A month is a twelfth of a year and a year is 365.25 days: the same conventions the spec's
        // own nYnMnD form implies, and precise far beyond what any age cell claims.
        private static double ToYears(double value, char unit) => unit switch
        {
            'Y' => value,
            'M' => value / 12d,
            'W' => value * 7d / DaysPerYear,
            'D' => value / DaysPerYear,
            'H' => value / (DaysPerYear * 24d),
            _ => throw new ArgumentOutOfRangeException(nameof(unit), unit, "not an age unit")
        };
    }
}
