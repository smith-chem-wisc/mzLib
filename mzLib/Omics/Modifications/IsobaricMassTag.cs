using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;

namespace Omics.Modifications
{
    /// <summary>
    /// The isobaric labelling kits whose reporter ions this library can name. TMT16 is TMTpro 16-plex,
    /// the lowest sixteen channels of the TMTpro (TMT18) set; it has no modification entry of its own.
    /// </summary>
    public enum IsobaricMassTagType
    {
        TMT6,
        TMT10,
        TMT11,
        TMT16,
        TMT18,
        iTRAQ4,
        iTRAQ8,
        diLeu4,
        diLeu12
    }

    /// <summary>
    /// An isobaric mass tag (TMT, iTRAQ, DiLeu): its channel labels and the theoretical m/z of each
    /// channel's reporter ion, in ascending m/z, plus the lookup that reads reporter intensities out of
    /// a spectrum. It holds no intensities itself.
    ///
    /// <para>The m/z are derived, never typed in: each is a <c>DI HCD:</c> line of the tag's
    /// <c>Multiplex Label</c> modification (<c>Omics/Resources/TMT.txt</c> by default) plus one proton.
    /// The labels are the only data held here. <see cref="ChannelLabels"/>[i] names the channel whose
    /// reporter ion is <see cref="ReporterIonMzs"/>[i].</para>
    /// </summary>
    public class IsobaricMassTag
    {
        /// <summary>
        /// Tolerance in daltons for matching an observed peak to a reporter ion.
        /// Taken from https://doi.org/10.1021/acs.jproteome.1c00168
        /// </summary>
        public const double AbsoluteToleranceValue = 0.003;

        private const string MultiplexLabelType = "Multiplex Label";

        public IsobaricMassTagType TagType { get; }

        /// <summary>
        /// One label per channel, in the same order as <see cref="ReporterIonMzs"/>.
        /// </summary>
        public IReadOnlyList<string> ChannelLabels { get; }

        /// <summary>
        /// Theoretical m/z (charge 1) of each channel's reporter ion, ascending.
        /// </summary>
        public double[] ReporterIonMzs { get; }

        public DoubleRange[] ReporterIonMzRanges { get; }

        private IsobaricMassTag(IsobaricMassTagType tagType, IReadOnlyList<string> channelLabels, double[] reporterIonMzs)
        {
            TagType = tagType;
            ChannelLabels = channelLabels;
            ReporterIonMzs = reporterIonMzs;
            var tolerance = new AbsoluteTolerance(AbsoluteToleranceValue);
            ReporterIonMzRanges = reporterIonMzs
                .Select(mz => new DoubleRange(tolerance.GetMinimumValue(mz), tolerance.GetMaximumValue(mz)))
                .ToArray();
        }

        /// <summary>
        /// Resolves a tag from a modification ID, as MetaMorpheus stores it (<c>"TMT10"</c>,
        /// <c>"TMT6-plex"</c>, <c>"iTRAQ-4plex on K"</c>), reading reporter ions from
        /// <see cref="Mods.IsobaricLabelModifications"/>.
        /// </summary>
        public static bool TryGetIsobaricMassTag(string? modificationId, out IsobaricMassTag? tag) =>
            TryGetIsobaricMassTag(modificationId, Mods.IsobaricLabelModifications, out tag);

        /// <summary>
        /// Resolves a tag from a modification ID, reading reporter ions from <paramref name="knownModifications"/>,
        /// so a user-defined <c>Multiplex Label</c> modification is found as well as the embedded ones.
        /// </summary>
        public static bool TryGetIsobaricMassTag(string? modificationId, IEnumerable<Modification> knownModifications,
            out IsobaricMassTag? tag)
        {
            tag = null;
            return TryGetTagType(modificationId, out var tagType)
                && TryGetIsobaricMassTag(tagType, knownModifications, out tag);
        }

        /// <summary>
        /// Resolves a tag of the given type, reading reporter ions from <see cref="Mods.IsobaricLabelModifications"/>.
        /// </summary>
        public static bool TryGetIsobaricMassTag(IsobaricMassTagType tagType, out IsobaricMassTag? tag) =>
            TryGetIsobaricMassTag(tagType, Mods.IsobaricLabelModifications, out tag);

        /// <summary>
        /// Resolves a tag of the given type, reading reporter ions from <paramref name="knownModifications"/>.
        /// False when no <c>Multiplex Label</c> modification for the type is known, when it has no HCD
        /// diagnostic ions, or when its reporter-ion count is not the kit's channel count, since the labels
        /// and ions are paired by position and a mismatch would mislabel every channel after it.
        /// </summary>
        public static bool TryGetIsobaricMassTag(IsobaricMassTagType tagType, IEnumerable<Modification> knownModifications,
            out IsobaricMassTag? tag)
        {
            tag = null;
            string sourceId = SourceModificationId(tagType);
            Modification? modification = knownModifications?.FirstOrDefault(m =>
                m.ModificationType == MultiplexLabelType
                && string.Equals(m.OriginalId, sourceId, StringComparison.OrdinalIgnoreCase));

            if (modification?.DiagnosticIons == null
                || !modification.DiagnosticIons.TryGetValue(DissociationType.HCD, out var reporterIonMasses)
                || reporterIonMasses == null)
            {
                return false;
            }

            // Diagnostic ions are neutral masses in file order, which is not channel order: sort, then
            // pair by position. TMT16 is the lowest sixteen of the TMTpro set, so its source carries 18.
            IReadOnlyList<string> labels = GetReporterIonLabels(tagType);
            int expectedIonCount = tagType == IsobaricMassTagType.TMT16
                ? GetReporterIonLabels(IsobaricMassTagType.TMT18).Count
                : labels.Count;
            if (reporterIonMasses.Count != expectedIonCount)
            {
                return false;
            }

            double[] reporterIonMzs = reporterIonMasses
                .Select(mass => mass.ToMz(1))
                .OrderBy(mz => mz)
                .Take(labels.Count)
                .ToArray();

            tag = new IsobaricMassTag(tagType, labels, reporterIonMzs);
            return true;
        }

        /// <summary>
        /// The channel labels of a kit, in ascending reporter-ion m/z.
        /// </summary>
        public static IReadOnlyList<string> GetReporterIonLabels(IsobaricMassTagType tagType)
        {
            return tagType switch
            {
                IsobaricMassTagType.TMT6 => new[] { "126", "127", "128", "129", "130", "131" },
                IsobaricMassTagType.TMT10 => new[] { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N" },
                IsobaricMassTagType.TMT11 => new[] { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C" },
                IsobaricMassTagType.TMT16 => new[] { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N" },
                IsobaricMassTagType.TMT18 => new[] { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N" },
                IsobaricMassTagType.iTRAQ4 => new[] { "114", "115", "116", "117" },
                // 8-plex has no 120 channel (the phenylalanine immonium ion sits at 120.081); the eighth
                // reagent is 121, whatever the "# 120" comment on its DI line in TMT.txt says.
                IsobaricMassTagType.iTRAQ8 => new[] { "113", "114", "115", "116", "117", "118", "119", "121" },
                IsobaricMassTagType.diLeu4 => new[] { "115", "116", "117", "118" },
                IsobaricMassTagType.diLeu12 => new[] { "115a", "115b", "116a", "116b", "116c", "117a", "117b", "117c", "118a", "118b", "118c", "118d" },
                _ => throw new ArgumentOutOfRangeException(nameof(tagType), tagType, null)
            };
        }

        /// <summary>
        /// Maps a modification ID to a kit. Matches the whole name, case-insensitively, optionally
        /// followed by <c>" on &lt;motif&gt;"</c>; never a substring, so an unknown kit is not mistaken
        /// for a known one whose name it contains.
        /// </summary>
        public static bool TryGetTagType(string? modificationId, out IsobaricMassTagType tagType)
        {
            tagType = default;
            if (string.IsNullOrWhiteSpace(modificationId))
            {
                return false;
            }

            string name = modificationId.Trim();
            int on = name.IndexOf(" on ", StringComparison.OrdinalIgnoreCase);
            if (on >= 0)
            {
                name = name.Substring(0, on).TrimEnd();
            }

            return NamesByTagType.TryGetValue(name, out tagType);
        }

        /// <summary>
        /// Every name <see cref="TryGetTagType"/> accepts. The first spelling of each kit is its
        /// modification's ID in TMT.txt; the rest are the short forms MetaMorpheus has accepted.
        /// </summary>
        private static readonly Dictionary<string, IsobaricMassTagType> NamesByTagType =
            new(StringComparer.OrdinalIgnoreCase)
            {
                ["TMT6-plex"] = IsobaricMassTagType.TMT6,
                ["TMT6"] = IsobaricMassTagType.TMT6,
                ["TMT10"] = IsobaricMassTagType.TMT10,
                ["TMT11"] = IsobaricMassTagType.TMT11,
                ["TMT16"] = IsobaricMassTagType.TMT16,
                ["TMT18"] = IsobaricMassTagType.TMT18,
                ["iTRAQ-4plex"] = IsobaricMassTagType.iTRAQ4,
                ["iTRAQ4"] = IsobaricMassTagType.iTRAQ4,
                ["iTRAQ-8plex"] = IsobaricMassTagType.iTRAQ8,
                ["iTRAQ8"] = IsobaricMassTagType.iTRAQ8,
                ["DiLeu-4plex"] = IsobaricMassTagType.diLeu4,
                ["DiLeu4"] = IsobaricMassTagType.diLeu4,
                ["DiLeu-12plex"] = IsobaricMassTagType.diLeu12,
                ["DiLeu12"] = IsobaricMassTagType.diLeu12,
            };

        /// <summary>
        /// The <c>OriginalId</c> of the modification that carries a kit's reporter ions.
        /// </summary>
        private static string SourceModificationId(IsobaricMassTagType tagType)
        {
            return tagType switch
            {
                IsobaricMassTagType.TMT6 => "TMT6-plex",
                IsobaricMassTagType.TMT10 => "TMT10",
                IsobaricMassTagType.TMT11 => "TMT11",
                IsobaricMassTagType.TMT16 => "TMT18",
                IsobaricMassTagType.TMT18 => "TMT18",
                IsobaricMassTagType.iTRAQ4 => "iTRAQ-4plex",
                IsobaricMassTagType.iTRAQ8 => "iTRAQ-8plex",
                IsobaricMassTagType.diLeu4 => "DiLeu-4plex",
                IsobaricMassTagType.diLeu12 => "DiLeu-12plex",
                _ => throw new ArgumentOutOfRangeException(nameof(tagType), tagType, null)
            };
        }

        /// <summary>
        /// Finds each reporter ion in <paramref name="spectrum"/> and returns its intensity, in the order of
        /// <see cref="ReporterIonMzs"/>. A reporter ion with no peak within tolerance gets zero; where several
        /// peaks fall within tolerance, the most intense wins. Null for a null or empty spectrum.
        /// </summary>
        public double[]? GetReporterIonIntensities(MzSpectrum? spectrum)
        {
            if (spectrum == null || spectrum.Size < 1) return null;
            double[] reporterIonIntensities = new double[ReporterIonMzs.Length];
            if (spectrum.XArray[0] > ReporterIonMzRanges[^1].Maximum)
            {
                return reporterIonIntensities;
            }

            // Two pointers: the for loop walks the theoretical ions, the while loops walk the peaks.
            int spectrumIndex = 0;
            for (int theoreticalIonIndex = 0; theoreticalIonIndex < reporterIonIntensities.Length; theoreticalIonIndex++)
            {
                double minMz = ReporterIonMzRanges[theoreticalIonIndex].Minimum;
                double maxMz = ReporterIonMzRanges[theoreticalIonIndex].Maximum;
                double maxIntensity = 0;
                while (spectrumIndex < spectrum.Size
                    && spectrum.XArray[spectrumIndex] < minMz)
                {
                    spectrumIndex++;
                }
                if (spectrumIndex >= spectrum.Size)
                {
                    break;
                }
                while (spectrumIndex < spectrum.Size
                    && minMz <= spectrum.XArray[spectrumIndex]
                    && spectrum.XArray[spectrumIndex] <= maxMz)
                {
                    if (spectrum.YArray[spectrumIndex] > maxIntensity)
                    {
                        maxIntensity = spectrum.YArray[spectrumIndex];
                    }
                    spectrumIndex++;
                }
                reporterIonIntensities[theoreticalIonIndex] = maxIntensity;
            }

            return reporterIonIntensities;
        }
    }
}
