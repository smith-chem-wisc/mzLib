using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using Omics.Modifications;

namespace Readers
{
    /// <summary>What a text says about isobaric labelling: whether it names TMT or iTRAQ at all (<see cref="Labelled"/>), and
    /// the kit when it names exactly one with its size (<see cref="Kit"/>), in the text's own words (<see cref="Quote"/>).</summary>
    internal readonly record struct SdrfKitReading(IsobaricMassTagType? Kit, bool Labelled, string Quote);

    /// <summary>
    /// Reads an isobaric kit from free text: a repository record's title, description and protocols, or the value of an
    /// evidence claim (sdrf PXR-F15 P1). Measured first on the 41 curated isobaric deposits with a cached PRIDE record
    /// (sdrf results/isobaric-p0-2026-11-16):
    /// <list type="bullet">
    /// <item>a record names the kit with its size in 17 and names another kit in 1;</item>
    /// <item>PRIDE's quantification method names only the family (<c>TMT</c>, <c>iTRAQ</c>), never the size;</item>
    /// <item><c>TMT6plex</c> is also Unimod's name for the 6-, 10- and 11-plex reagent, and a record whose only kit word was
    /// TMT6plex was a 10-plex in 3 of 12. So no spelling of TMT 6-plex states a kit here: it only says the runs are
    /// labelled. TMTpro without a number is 16 or 18, and is family only too.</item>
    /// </list>
    /// A text that names two sized kits names neither: it may describe two experiments. So does a 6-plex written apart
    /// (<c>TMT 10-plex/ 6-plex</c>, PXD009243) beside another kit, and so does a text that also names another way of
    /// quantifying (label-free, SILAC, 18O, dimethyl): measured on the 1,148 cached deposits, those deposits hold runs of
    /// both kinds (PXD014565 has label-free, SILAC and TMT runs), and a kit for the whole deposit would label all of them.
    /// </summary>
    internal static class SdrfIsobaricKit
    {
        private const RegexOptions Options = RegexOptions.IgnoreCase | RegexOptions.CultureInvariant | RegexOptions.Compiled;

        // A size written beside its reagent, either way round: TMT10plex, TMT 10-plex, 10-plex TMT, iTRAQ Reagent-8plex.
        private static Regex Sized(string reagent, string size) => new(
            $@"\b{reagent}(?:\s+reagents?)?\W{{0,3}}{size}(?:\W{{0,2}}plex)?(?!\d)|(?<!\d){size}\W{{0,2}}plex\W{{0,3}}{reagent}\b", Options);

        private static readonly (IsobaricMassTagType Kit, Regex Pattern)[] Kits =
        {
            (IsobaricMassTagType.TMT18, Sized(@"tmt\s*pro", "18")),
            (IsobaricMassTagType.TMT18, Sized("tmt", "18")),
            (IsobaricMassTagType.TMT16, Sized(@"tmt\s*pro", "16")),
            (IsobaricMassTagType.TMT16, Sized("tmt", "16")),
            (IsobaricMassTagType.TMT11, Sized("tmt", "11")),
            (IsobaricMassTagType.TMT10, Sized("tmt", "10")),
            (IsobaricMassTagType.iTRAQ8, Sized("itraq", "8")),
            (IsobaricMassTagType.iTRAQ4, Sized("itraq", "4")),
        };

        // A 6-plex written apart from its "plex" is a kit; TMT6plex written as one word may be the modification's name.
        private static readonly Regex SixPlexApart = new(@"(?<!\d)6\W{1,2}plex\b|\btmt\s*sixplex\b|\bsix-?plex\b", Options);

        // Another way of quantifying named in the same text: the deposit may hold runs that are not labelled with the kit.
        // Whole words only: dimethylsulfoxide is a solvent (PXD017201), and a label-free buffer is a reagent (PXD013753).
        private static readonly Regex OtherQuantification = new(@"\blabel[\s-]*free\b(?!\s+buffer)|\blfq\b|\bsilac\b|\b18o\b|\bdimethyl\b", Options);

        // The whole word, so a reason can quote what the text said (TMT6plex, TMTpro), not only "TMT".
        private static readonly Regex Family = new(@"\b(?:tmt|itraq)\w*|tandem\s+mass\s+tags?", Options);

        /// <summary>Reads <paramref name="text"/>. Never throws; null or empty text is not labelled.</summary>
        internal static SdrfKitReading Read(string? text)
        {
            if (string.IsNullOrWhiteSpace(text)) return default;
            var family = Family.Match(text);
            if (!family.Success) return default;
            var found = new List<(IsobaricMassTagType Kit, string Quote)>();
            foreach (var (kit, pattern) in Kits)
            {
                var m = pattern.Match(text);
                if (m.Success && found.All(f => f.Kit != kit)) found.Add((kit, m.Value));
            }
            if (found.Count == 1 && SixPlexApart.Match(text) is { Success: true } six)
                return new SdrfKitReading(null, true, $"{found[0].Quote}', '{six.Value}");
            if (found.Count == 1 && OtherQuantification.Match(text) is { Success: true } other)
                return new SdrfKitReading(null, true, $"{found[0].Quote}', '{other.Value}");
            return found.Count == 1
                ? new SdrfKitReading(found[0].Kit, true, found[0].Quote)
                : new SdrfKitReading(null, true, found.Count > 1 ? string.Join("', '", found.Select(f => f.Quote)) : family.Value);
        }
    }
}
