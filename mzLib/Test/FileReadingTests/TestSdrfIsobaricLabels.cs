using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Omics.Modifications;
using Readers;
using UsefulProteomicsDatabases;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Covers <see cref="SdrfIsobaricLabels"/>, the one translation between a kit's channel as mzLib spells it and the
    /// SDRF <c>comment[label]</c> that names it (sdrf QP-S26). The labels it writes are checked against sdrf-pipelines'
    /// own table and against the PRIDE vocabulary, not against mzLib.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfIsobaricLabels
    {
        /// <summary>
        /// sdrf-pipelines 0.1.6, <c>converters/channel_map.yaml</c> (schema_version 2): the labels its converters read for
        /// each kit, in channel order. Independent of mzLib. DiLeu has no table there and no PRIDE term.
        /// </summary>
        private static readonly Dictionary<IsobaricMassTagType, string[]> SdrfPipelinesLabels = new()
        {
            [IsobaricMassTagType.TMT6] = new[] { "TMT126", "TMT127", "TMT128", "TMT129", "TMT130", "TMT131" },
            [IsobaricMassTagType.TMT10] = new[]
            {
                "TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131"
            },
            [IsobaricMassTagType.TMT11] = new[]
            {
                "TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131N",
                "TMT131C"
            },
            [IsobaricMassTagType.TMT16] = new[]
            {
                "TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131N",
                "TMT131C", "TMT132N", "TMT132C", "TMT133N", "TMT133C", "TMT134N"
            },
            [IsobaricMassTagType.TMT18] = new[]
            {
                "TMT126", "TMT127N", "TMT127C", "TMT128N", "TMT128C", "TMT129N", "TMT129C", "TMT130N", "TMT130C", "TMT131N",
                "TMT131C", "TMT132N", "TMT132C", "TMT133N", "TMT133C", "TMT134N", "TMT134C", "TMT135N"
            },
            [IsobaricMassTagType.iTRAQ4] = new[] { "ITRAQ114", "ITRAQ115", "ITRAQ116", "ITRAQ117" },
            [IsobaricMassTagType.iTRAQ8] = new[]
            {
                "ITRAQ113", "ITRAQ114", "ITRAQ115", "ITRAQ116", "ITRAQ117", "ITRAQ118", "ITRAQ119", "ITRAQ121"
            },
        };

        private static IEnumerable<IsobaricMassTagType> AllTypes() =>
            Enum.GetValues(typeof(IsobaricMassTagType)).Cast<IsobaricMassTagType>();

        private static string[] Written(IsobaricMassTagType type) =>
            IsobaricMassTag.GetReporterIonLabels(type).Select(c => SdrfIsobaricLabels.GetSdrfLabel(type, c)).ToArray();

        [Test]
        public void EveryLabelIsTheOneSdrfPipelinesReads()
        {
            var disagreements = SdrfPipelinesLabels
                .Where(kit => !Written(kit.Key).SequenceEqual(kit.Value))
                .Select(kit => $"{kit.Key}: {string.Join(", ", Written(kit.Key))}")
                .ToList();

            Assert.That(disagreements, Is.Empty);
            Assert.That(AllTypes().Except(SdrfPipelinesLabels.Keys),
                Is.EquivalentTo(new[] { IsobaricMassTagType.diLeu4, IsobaricMassTagType.diLeu12 }),
                "a new kit needs its labels pinned here");
        }

        [Test]
        public void Tmt10sLastChannelIsPridesTmt131AndTmt11KeepsTmt131N()
        {
            Assert.That(SdrfIsobaricLabels.GetSdrfLabel(IsobaricMassTagType.TMT10, "131N"), Is.EqualTo("TMT131"));
            Assert.That(ControlledVocabulary.Pride.TryGetByName("TMT131", out var term), Is.True);
            Assert.That(term.Accession, Is.EqualTo("PRIDE:0000529"));
            Assert.That(SdrfIsobaricLabels.GetSdrfLabel(IsobaricMassTagType.TMT11, "131N"), Is.EqualTo("TMT131N"));
        }

        [Test]
        public void EveryTmtAndItraqLabelIsAPrideTerm()
        {
            var missing = SdrfPipelinesLabels.Keys
                .SelectMany(type => Written(type).Select(label => (type, label)))
                .Where(x => !ControlledVocabulary.Pride.TryGetByName(x.label, out _))
                .Select(x => $"{x.type} {x.label}")
                .ToList();

            Assert.That(missing, Is.Empty);
        }

        [Test]
        public void EveryChannelOfEveryKitSurvivesTheRoundTrip()
        {
            var disagreements = new List<string>();
            foreach (var type in AllTypes())
            {
                foreach (string channel in IsobaricMassTag.GetReporterIonLabels(type))
                {
                    string label = SdrfIsobaricLabels.GetSdrfLabel(type, channel);
                    var forms = new List<string> { label, label.ToLowerInvariant() };
                    if (ControlledVocabulary.Pride.TryGetByName(label, out var term))
                        forms.Add($"NT={label};AC={term.Accession}");

                    foreach (string form in forms)
                    {
                        if (!SdrfIsobaricLabels.TryGetChannel(type, form, out string? read) || read != channel)
                            disagreements.Add($"{type} {channel}: '{form}' read as '{read}'");
                    }
                }
            }

            Assert.That(disagreements, Is.Empty);
        }

        [TestCase(IsobaricMassTagType.TMT10, "131C")]
        [TestCase(IsobaricMassTagType.TMT6, "131N")]
        [TestCase(IsobaricMassTagType.iTRAQ4, "113")]
        public void GetSdrfLabelRefusesAChannelTheKitDoesNotHave(IsobaricMassTagType type, string channel)
        {
            Assert.That(() => SdrfIsobaricLabels.GetSdrfLabel(type, channel),
                Throws.ArgumentException.With.Message.Contains($"not a channel of {type}"));
        }

        // TMT11's 131N is TMT131N, as sdrf-pipelines' tmt11plex table has it; TMT131C is no TMT10 channel; 115 is iTRAQ's
        // and DiLeu's, so only the family tells them apart.
        [TestCase(IsobaricMassTagType.TMT11, "TMT131")]
        [TestCase(IsobaricMassTagType.TMT10, "TMT131C")]
        [TestCase(IsobaricMassTagType.diLeu4, "ITRAQ115")]
        [TestCase(IsobaricMassTagType.TMT10, "TMTpro126")]
        [TestCase(IsobaricMassTagType.TMT10, "TMT10plex")]
        [TestCase(IsobaricMassTagType.TMT10, " ")]
        [TestCase(IsobaricMassTagType.TMT10, null)]
        public void TryGetChannelRefusesWhatIsNotALabelOfTheKit(IsobaricMassTagType type, string? label)
        {
            Assert.That(SdrfIsobaricLabels.TryGetChannel(type, label, out string? channel), Is.False);
            Assert.That(channel, Is.Null);
        }

        [Test]
        public void AnUndefinedKitThrows()
        {
            Assert.That(() => SdrfIsobaricLabels.TryGetChannel((IsobaricMassTagType)99, "TMT126", out _),
                Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => SdrfIsobaricLabels.GetSdrfLabel((IsobaricMassTagType)99, "126"),
                Throws.TypeOf<ArgumentOutOfRangeException>());
        }
    }
}
