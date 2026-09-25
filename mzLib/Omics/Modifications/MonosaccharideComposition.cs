using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace Omics.Modifications
{
    /// <summary>
    /// How many of each monosaccharide a glycan is built from, with no claim about how they are joined.
    /// </summary>
    /// <remarks>
    /// <para><b>Composition, deliberately, and not structure.</b> It is the part of a glycan's identity
    /// that is portable, unambiguous and cheap, and it separates the cases enzymology most often turns
    /// on: OpeRATOR refuses the Tn antigen (HexNAc alone), works on core 1 (Hex + HexNAc), and is blocked
    /// again by the second HexNAc that makes core 2.</para>
    ///
    /// <para><b>What it cannot do, stated so nobody expects it to.</b> Two glycans with the same
    /// composition can behave completely differently. Alpha-2,6-sialyl core 1 blocks OpeRATOR while
    /// alpha-2,3 merely slows it, and those are the same counts differing only in linkage; beta-O-GlcNAc
    /// and alpha-O-GalNAc are both one HexNAc, yet StcE accepts one and refuses the other. Anything
    /// turning on linkage, anomeric configuration, branching topology or which hexosamine it is needs a
    /// real structure model. This is not one, and a rule written against it must not pretend otherwise.
    /// </para>
    ///
    /// <para>Not a <see cref="Chemistry.ChemicalFormula"/>: that is keyed on <c>Element</c> and its
    /// vocabulary is the periodic table, so it is structurally closed to a Hex/HexNAc key. The two are
    /// complementary -- a composition says which sugars, a formula says which atoms.</para>
    ///
    /// <para>The monosaccharide names and one-letter codes match the slots the glycan databases already
    /// use, so a composition can be read from one without translation.</para>
    /// </remarks>
    public sealed class MonosaccharideComposition
    {
        // Canonical name and its one-letter database code. A private table rather than a public
        // vocabulary: Omics expresses fixed vocabularies as private static lookups or enums, and there is
        // no static readonly IReadOnlyList<T> vocabulary anywhere in the repo to follow.
        private static readonly (string Name, char Code)[] Slots =
        {
            ("Hex", 'H'), ("HexNAc", 'N'), ("NeuAc", 'A'), ("NeuGc", 'G'), ("Fuc", 'F'),
            ("Phospho", 'P'), ("Sulfo", 'S'), ("Na", 'Y'), ("Ac", 'C'), ("Xylose", 'X'), ("Kdn", 'K'),
        };

        // Both spellings resolve to one canonical name, so "Hex1HexNAc1" and "H1N1" are one composition.
        private static readonly Dictionary<string, string> CanonicalNames = BuildNameLookup();

        private static Dictionary<string, string> BuildNameLookup()
        {
            var lookup = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
            foreach ((string name, char code) in Slots)
            {
                lookup[name] = name;
                lookup[code.ToString()] = name;
            }

            return lookup;
        }

        private static readonly Regex Term = new(@"([A-Za-z]+)\(?(\d+)\)?", RegexOptions.Compiled);

        private readonly SortedDictionary<string, int> _counts;

        private MonosaccharideComposition(SortedDictionary<string, int> counts) => _counts = counts;

        /// <summary>How many of this monosaccharide the glycan carries; zero when it carries none.</summary>
        public int this[string monosaccharide] =>
            monosaccharide is not null
            && CanonicalNames.TryGetValue(monosaccharide, out string canonical)
            && _counts.TryGetValue(canonical, out int count)
                ? count
                : 0;

        /// <summary>The monosaccharides actually present, canonically cased.</summary>
        public IEnumerable<string> Monosaccharides => _counts.Keys;

        /// <summary>The total number of monosaccharide units.</summary>
        public int TotalUnits => _counts.Values.Sum();

        /// <summary>
        /// True when this composition holds at least as many of every monosaccharide as
        /// <paramref name="other"/> does -- "is at least this large, in every component".
        /// </summary>
        /// <remarks>
        /// Named for <see cref="Chemistry.ChemicalFormula.IsSupersetOf"/>, which is the same question
        /// asked of elements. A component-wise floor rather than a total, because that is the shape
        /// enzymology takes: "core 1 or bigger" is Hex(1)HexNAc(1) and up, and a second HexNAc is what
        /// makes core 2 regardless of what else is attached. A null argument is no constraint at all.
        /// </remarks>
        public bool IsSupersetOf(MonosaccharideComposition other)
        {
            if (other is null)
            {
                return true;
            }

            foreach (KeyValuePair<string, int> required in other._counts)
            {
                if (this[required.Key] < required.Value)
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Reads a composition written as monosaccharide-and-count pairs, in either the long spelling or
        /// the one-letter database codes, with or without parentheses: <c>HexNAc1Hex1</c>,
        /// <c>HexNAc(1)Hex(1)</c> and <c>N1H1</c> are all the same composition.
        /// </summary>
        /// <exception cref="MzLibException">The text is not a composition this can read.</exception>
        public static MonosaccharideComposition Parse(string text)
        {
            if (!TryParse(text, out MonosaccharideComposition composition, out string error))
            {
                throw new MzLibException(error);
            }

            return composition;
        }

        public static bool TryParse(string text, [MaybeNullWhen(false)] out MonosaccharideComposition composition) =>
            TryParse(text, out composition, out _);

        private static bool TryParse(string text, [MaybeNullWhen(false)] out MonosaccharideComposition composition,
            out string error)
        {
            composition = null;
            error = null;

            if (string.IsNullOrWhiteSpace(text))
            {
                error = "A glycan composition cannot be empty. Write it as monosaccharide-and-count pairs, "
                        + "for example HexNAc1Hex1.";
                return false;
            }

            string trimmed = text.Replace(" ", string.Empty);
            var counts = new SortedDictionary<string, int>(StringComparer.Ordinal);
            int consumed = 0;

            foreach (Match match in Term.Matches(trimmed))
            {
                // Anchored to the end of the previous term, so trailing or interleaved junk cannot be
                // silently skipped -- a composition that half-parsed would quietly relax an enzyme's rule.
                if (match.Index != consumed)
                {
                    break;
                }

                consumed = match.Index + match.Length;

                string name = match.Groups[1].Value;
                if (!CanonicalNames.TryGetValue(name, out string canonical))
                {
                    error = "Unrecognized monosaccharide '" + name + "' in glycan composition '" + text
                            + "'. Known monosaccharides are "
                            + string.Join(", ", Slots.Select(s => s.Name + " (" + s.Code + ")")) + ".";
                    return false;
                }

                int count = int.Parse(match.Groups[2].Value);
                counts[canonical] = counts.TryGetValue(canonical, out int already) ? already + count : count;
            }

            if (consumed != trimmed.Length || counts.Count == 0)
            {
                error = "Unrecognized glycan composition '" + text + "'. Expected monosaccharide-and-count "
                        + "pairs such as HexNAc1Hex1, HexNAc(1)Hex(1) or N1H1.";
                return false;
            }

            composition = new MonosaccharideComposition(counts);
            return true;
        }

        /// <summary>
        /// The canonical long spelling, which <see cref="Parse"/> reads back unchanged.
        /// </summary>
        public override string ToString()
        {
            var text = new StringBuilder();
            foreach (KeyValuePair<string, int> entry in _counts)
            {
                text.Append(entry.Key).Append('(').Append(entry.Value).Append(')');
            }

            return text.ToString();
        }

        // Value equality, because this is a value. Both halves are implemented together and over the same
        // field, unlike Modification, whose Equals and GetHashCode read different fields. No operator==:
        // no type in Omics or Chemistry defines one.
        //
        // Equal means each covers the other: IsSupersetOf is a component-wise floor, so one direction
        // alone makes Hex2 "equal" Hex1. A zero count is the same as an absent monosaccharide -- the
        // indexer answers 0 for both -- so Hex0HexNAc1 equals HexNAc1, and the hash skips zeros to agree.
        public override bool Equals(object obj) =>
            obj is MonosaccharideComposition other
            && IsSupersetOf(other)
            && other.IsSupersetOf(this);

        public override int GetHashCode()
        {
            int hash = 17;
            foreach (KeyValuePair<string, int> entry in _counts)
            {
                if (entry.Value == 0)
                {
                    continue;
                }

                hash = hash * 31 + entry.Key.GetHashCode();
                hash = hash * 31 + entry.Value;
            }

            return hash;
        }
    }
}
