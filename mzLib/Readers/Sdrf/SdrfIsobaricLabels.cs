using System;
using System.Linq;
using System.Text.RegularExpressions;
using Omics.Modifications;

namespace Readers
{
    /// <summary>
    /// The one translation between an isobaric kit's channel as mzLib spells it
    /// (<see cref="IsobaricMassTag.ChannelLabels"/>, which is also the spelling of <c>TmtDesign.txt</c>) and the SDRF
    /// <c>comment[label]</c> that names it. A writer of SDRFs and <see cref="SdrfIsobaricDesign"/> both go through it, so
    /// a design written as an SDRF reads back to the same channels (sdrf QP-S26).
    /// <para>
    /// <b>The labels are PRIDE's, as sdrf-pipelines reads them:</b> <c>TMT</c> or <c>ITRAQ</c> and the channel, except
    /// that TMT10's last channel, 131N to mzLib, is <c>TMT131</c> (PRIDE:0000529; the <c>tmt10plex</c> table of
    /// sdrf-pipelines 0.1.6; of the TMT10 SDRFs in bigbio/sdrf-annotated-datasets that label it, 23 write <c>TMT131</c>
    /// and 2 <c>TMT131N</c>). TMT11 and the TMTpro kits
    /// keep <c>TMT131N</c>. DiLeu has no PRIDE term; its label is <c>DiLeu</c> and the channel, the form this class reads.
    /// </para>
    /// <para>
    /// <b>Reading accepts more than writing produces:</b> the label with or without its family, in any case, bare or
    /// accessioned (<c>NT=TMT131;AC=PRIDE:0000529</c>), and the kit's own channel name, so TMT10 also reads
    /// <c>TMT131N</c>, as PXD008841 writes it. A family must name the kit: <c>TMTpro</c> fits TMT16 and TMT18 only, and
    /// <c>ITRAQ115</c> is never DiLeu4's 115 (MAP-06).
    /// </para>
    /// </summary>
    public static class SdrfIsobaricLabels
    {
        // MAP-06 / N4: an optional family, then the channel as the kit or the SDRF spells it.
        private static readonly Regex LabelName =
            new(@"^(?<family>tmtpro|tmt|itraq|dileu)?[-_ ]?(?<channel>\d{3}[a-z]?)$",
                RegexOptions.IgnoreCase | RegexOptions.CultureInvariant);

        /// <summary>
        /// The SDRF <c>comment[label]</c> name of <paramref name="channel"/> in the kit of <paramref name="tagType"/>,
        /// e.g. <c>TMT131</c> for TMT10's 131N and <c>TMT127N</c> for 127N. Matching the channel ignores case.
        /// </summary>
        /// <exception cref="ArgumentException">The kit has no channel <paramref name="channel"/>.</exception>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="tagType"/> is not a defined kit.</exception>
        public static string GetSdrfLabel(IsobaricMassTagType tagType, string channel)
        {
            string? kitChannel = IsobaricMassTag.GetReporterIonLabels(tagType)
                .FirstOrDefault(c => string.Equals(c, channel?.Trim(), StringComparison.OrdinalIgnoreCase));
            if (kitChannel == null)
            {
                throw new ArgumentException($"'{channel}' is not a channel of {tagType}. Its channels are " +
                    $"{string.Join(", ", IsobaricMassTag.GetReporterIonLabels(tagType))}.", nameof(channel));
            }

            return Family(tagType) + SdrfChannel(tagType, kitChannel);
        }

        /// <summary>
        /// The channel of the kit of <paramref name="tagType"/> that a <c>comment[label]</c> cell names, in the kit's own
        /// spelling: TMT10's <c>TMT131</c> reads as <c>131N</c>. The cell may be bare or accessioned. False, with
        /// <paramref name="channel"/> null, for an empty cell or one that names no channel of the kit.
        /// </summary>
        /// <exception cref="ArgumentOutOfRangeException"><paramref name="tagType"/> is not a defined kit.</exception>
        public static bool TryGetChannel(IsobaricMassTagType tagType, string? label, out string? channel)
        {
            channel = null;
            var channels = IsobaricMassTag.GetReporterIonLabels(tagType);
            if (string.IsNullOrWhiteSpace(label))
            {
                return false;
            }

            var match = LabelName.Match(SdrfQuantAuditor.ReadLabel(label).Name);
            if (!match.Success || (match.Groups["family"].Success && !FamilyFits(match.Groups["family"].Value, tagType)))
            {
                return false;
            }

            string written = match.Groups["channel"].Value;
            channel = channels.FirstOrDefault(c => string.Equals(c, written, StringComparison.OrdinalIgnoreCase)
                                                  || string.Equals(SdrfChannel(tagType, c), written, StringComparison.OrdinalIgnoreCase));
            return channel != null;
        }

        // The one channel whose SDRF spelling differs from mzLib's: TMT10's 131N is PRIDE's TMT131 (PRIDE:0000529).
        private static string SdrfChannel(IsobaricMassTagType tagType, string channel) =>
            tagType == IsobaricMassTagType.TMT10 && channel == "131N" ? "131" : channel;

        private static string Family(IsobaricMassTagType tagType) => tagType switch
        {
            IsobaricMassTagType.TMT6 or IsobaricMassTagType.TMT10 or IsobaricMassTagType.TMT11
                or IsobaricMassTagType.TMT16 or IsobaricMassTagType.TMT18 => "TMT",
            IsobaricMassTagType.iTRAQ4 or IsobaricMassTagType.iTRAQ8 => "ITRAQ",
            IsobaricMassTagType.diLeu4 or IsobaricMassTagType.diLeu12 => "DiLeu",
            _ => throw new ArgumentOutOfRangeException(nameof(tagType), tagType, null)
        };

        private static bool FamilyFits(string family, IsobaricMassTagType tagType)
        {
            string kit = tagType.ToString();
            if (family.Equals("tmtpro", StringComparison.OrdinalIgnoreCase))
                return tagType is IsobaricMassTagType.TMT16 or IsobaricMassTagType.TMT18;
            return kit.StartsWith(family, StringComparison.OrdinalIgnoreCase);
        }
    }
}
