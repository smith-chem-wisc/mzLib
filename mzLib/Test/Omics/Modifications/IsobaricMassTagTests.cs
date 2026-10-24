using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Modifications;

namespace Test.Omics.Modifications
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class IsobaricMassTagTests
    {
        /// <summary>
        /// Published reporter-ion m/z for the TMT/TMTpro isotopologues (Thermo TMTpro user guide;
        /// Unimod 737 and 2016). Independent of TMT.txt, which the tags are derived from.
        /// </summary>
        private static readonly Dictionary<string, double> PublishedTmtReporterMz = new()
        {
            ["126"] = 126.127726, ["127N"] = 127.124761, ["127C"] = 127.131081,
            ["128N"] = 128.128116, ["128C"] = 128.134436, ["129N"] = 129.131471,
            ["129C"] = 129.137790, ["130N"] = 130.134825, ["130C"] = 130.141145,
            ["131N"] = 131.138180, ["131C"] = 131.144500, ["132N"] = 132.141535,
            ["132C"] = 132.147855, ["133N"] = 133.144890, ["133C"] = 133.151210,
            ["134N"] = 134.148245, ["134C"] = 134.154565, ["135N"] = 135.151600,
        };

        /// <summary>
        /// TMT6 channels carry bare nominal labels; these are the isotopologues the kit uses.
        /// </summary>
        private static readonly Dictionary<string, string> Tmt6Isotopologues = new()
        {
            ["126"] = "126", ["127"] = "127N", ["128"] = "128C",
            ["129"] = "129N", ["130"] = "130C", ["131"] = "131N",
        };

        private static IsobaricMassTag Resolve(IsobaricMassTagType type)
        {
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(type, out var tag), Is.True, type.ToString());
            return tag!;
        }

        private static IEnumerable<IsobaricMassTagType> AllTypes() =>
            Enum.GetValues(typeof(IsobaricMassTagType)).Cast<IsobaricMassTagType>();

        [Test]
        public void EveryKitResolvesWithItsChannelCount()
        {
            var expected = new Dictionary<IsobaricMassTagType, int>
            {
                [IsobaricMassTagType.TMT6] = 6,
                [IsobaricMassTagType.TMT10] = 10,
                [IsobaricMassTagType.TMT11] = 11,
                [IsobaricMassTagType.TMT16] = 16,
                [IsobaricMassTagType.TMT18] = 18,
                [IsobaricMassTagType.iTRAQ4] = 4,
                [IsobaricMassTagType.iTRAQ8] = 8,
                [IsobaricMassTagType.diLeu4] = 4,
                [IsobaricMassTagType.diLeu12] = 12,
            };
            Assert.That(AllTypes(), Is.EquivalentTo(expected.Keys), "a new kit needs its count pinned here");

            foreach (var type in AllTypes())
            {
                var tag = Resolve(type);
                Assert.That(tag.TagType, Is.EqualTo(type));
                Assert.That(tag.ChannelLabels, Has.Count.EqualTo(expected[type]), type.ToString());
                Assert.That(tag.ReporterIonMzs, Has.Length.EqualTo(expected[type]), type.ToString());
                Assert.That(tag.ReporterIonMzRanges, Has.Length.EqualTo(expected[type]), type.ToString());
                Assert.That(tag.ChannelLabels, Is.EqualTo(IsobaricMassTag.GetReporterIonLabels(type)));
            }
        }

        [Test]
        public void ReporterIonsAscendStrictly()
        {
            foreach (var type in AllTypes())
            {
                var mzs = Resolve(type).ReporterIonMzs;
                for (int i = 1; i < mzs.Length; i++)
                {
                    Assert.That(mzs[i], Is.GreaterThan(mzs[i - 1]), $"{type} at index {i}");
                }
            }
        }

        /// <summary>
        /// The label at index i must name the ion at index i. Checked against published m/z, not
        /// against TMT.txt: the N and C forms of one nominal mass differ by 0.006, which nominal
        /// mass cannot see.
        /// </summary>
        [Test]
        [TestCase(IsobaricMassTagType.TMT6)]
        [TestCase(IsobaricMassTagType.TMT10)]
        [TestCase(IsobaricMassTagType.TMT11)]
        [TestCase(IsobaricMassTagType.TMT16)]
        [TestCase(IsobaricMassTagType.TMT18)]
        public void TmtChannelsMatchThePublishedReporterMz(IsobaricMassTagType type)
        {
            var tag = Resolve(type);
            var disagreements = new List<string>();
            for (int i = 0; i < tag.ChannelLabels.Count; i++)
            {
                string label = tag.ChannelLabels[i];
                string isotopologue = type == IsobaricMassTagType.TMT6 ? Tmt6Isotopologues[label] : label;
                double published = PublishedTmtReporterMz[isotopologue];
                if (Math.Abs(tag.ReporterIonMzs[i] - published) > 1e-5)
                {
                    disagreements.Add($"{label} at {i}: {tag.ReporterIonMzs[i]:F6}, published {published:F6}");
                }
            }
            Assert.That(disagreements, Is.Empty);
        }

        /// <summary>
        /// Every label begins with its channel's nominal mass. This is the only check iTRAQ and DiLeu
        /// get, and it is the one that caught iTRAQ8's eighth channel named 120 rather than 121.
        /// </summary>
        [Test]
        public void EveryLabelNamesTheNominalMassAtItsOwnIndex()
        {
            var disagreements = new List<string>();
            foreach (var type in AllTypes())
            {
                var tag = Resolve(type);
                for (int i = 0; i < tag.ChannelLabels.Count; i++)
                {
                    string digits = new string(tag.ChannelLabels[i].TakeWhile(char.IsDigit).ToArray());
                    int observed = (int)Math.Round(tag.ReporterIonMzs[i]);
                    if (digits != observed.ToString())
                    {
                        disagreements.Add($"{type} '{tag.ChannelLabels[i]}' at {i} is {tag.ReporterIonMzs[i]:F4}");
                    }
                }
            }
            Assert.That(disagreements, Is.Empty);
        }

        [Test]
        public void Tmt16IsTheLowestSixteenOfTmt18()
        {
            var tmt16 = Resolve(IsobaricMassTagType.TMT16);
            var tmt18 = Resolve(IsobaricMassTagType.TMT18);
            Assert.That(tmt16.ChannelLabels, Is.EqualTo(tmt18.ChannelLabels.Take(16)));
            Assert.That(tmt16.ReporterIonMzs, Is.EqualTo(tmt18.ReporterIonMzs.Take(16)));
        }

        /// <summary>
        /// A kit added to TMT.txt without a label list here would be invisible; this makes it fail.
        /// </summary>
        [Test]
        public void EveryEmbeddedMultiplexLabelResolvesToAKit()
        {
            var multiplexMods = Mods.IsobaricLabelModifications
                .Where(m => m.ModificationType == "Multiplex Label")
                .ToList();
            Assert.That(multiplexMods.Select(m => m.OriginalId).Distinct().Count(), Is.EqualTo(8));

            foreach (var mod in multiplexMods)
            {
                Assert.That(IsobaricMassTag.TryGetTagType(mod.OriginalId, out _), Is.True, mod.OriginalId);
                Assert.That(IsobaricMassTag.TryGetTagType(mod.IdWithMotif, out _), Is.True, mod.IdWithMotif);
                Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(mod.IdWithMotif, out _), Is.True, mod.IdWithMotif);
            }
        }

        [Test]
        [TestCase("TMT6-plex", IsobaricMassTagType.TMT6)]
        [TestCase("TMT6 on K", IsobaricMassTagType.TMT6)]
        [TestCase("TMT10", IsobaricMassTagType.TMT10)]
        [TestCase("TMT10 on K", IsobaricMassTagType.TMT10)]
        [TestCase("TMT11 on X", IsobaricMassTagType.TMT11)]
        [TestCase("TMT16", IsobaricMassTagType.TMT16)]
        [TestCase("TMT18 on X", IsobaricMassTagType.TMT18)]
        [TestCase("iTRAQ-4plex on K", IsobaricMassTagType.iTRAQ4)]
        [TestCase("iTRAQ4 on K", IsobaricMassTagType.iTRAQ4)]
        [TestCase("iTRAQ-8plex", IsobaricMassTagType.iTRAQ8)]
        [TestCase("iTRAQ8 on K", IsobaricMassTagType.iTRAQ8)]
        [TestCase("DiLeu-4plex on K", IsobaricMassTagType.diLeu4)]
        [TestCase("DiLeu4 on K", IsobaricMassTagType.diLeu4)]
        [TestCase("DiLeu-12plex on X", IsobaricMassTagType.diLeu12)]
        [TestCase("DiLeu12 on X", IsobaricMassTagType.diLeu12)]
        [TestCase("tmt10 on k", IsobaricMassTagType.TMT10)]
        [TestCase("ITRAQ-4PLEX ON K", IsobaricMassTagType.iTRAQ4)]
        [TestCase("  TMT11  ", IsobaricMassTagType.TMT11)]
        public void ModificationIdsResolveToTheirKit(string modificationId, IsobaricMassTagType expected)
        {
            Assert.That(IsobaricMassTag.TryGetTagType(modificationId, out var type), Is.True);
            Assert.That(type, Is.EqualTo(expected));
        }

        /// <summary>
        /// Whole-name matching: a name that merely contains a known kit's name is not that kit.
        /// </summary>
        [Test]
        [TestCase(null)]
        [TestCase("")]
        [TestCase("   ")]
        [TestCase("InvalidModification")]
        [TestCase("TMT1")]
        [TestCase("TMT100")]
        [TestCase("TMT64 on K")]
        [TestCase("TMT10plex")]
        [TestCase("xTMT10")]
        [TestCase("TMTpro")]
        public void UnknownModificationIdsDoNotResolve(string? modificationId)
        {
            Assert.That(IsobaricMassTag.TryGetTagType(modificationId, out _), Is.False);
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(modificationId, out var tag), Is.False);
            Assert.That(tag, Is.Null);
        }

        [Test]
        public void KnownModificationsAreAParameter()
        {
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(IsobaricMassTagType.TMT10, new List<Modification>(), out var none), Is.False);
            Assert.That(none, Is.Null);

            // A user-defined TMT10 modification is read in place of the embedded one.
            var embedded = Resolve(IsobaricMassTagType.TMT10).ReporterIonMzs;
            var custom = MultiplexMod("TMT10", embedded.Select(mz => mz + 0.5).ToList());
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag("TMT10 on K", new[] { custom }, out var tag), Is.True);
            Assert.That(tag!.ReporterIonMzs, Is.EqualTo(embedded.Select(mz => mz + 0.5)).Within(1e-9));
        }

        [Test]
        public void AModificationWithTheWrongIonCountIsRefused()
        {
            var nine = Resolve(IsobaricMassTagType.TMT10).ReporterIonMzs.Take(9).ToList();
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(IsobaricMassTagType.TMT10, new[] { MultiplexMod("TMT10", nine) }, out _), Is.False);

            var eleven = Resolve(IsobaricMassTagType.TMT11).ReporterIonMzs.ToList();
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(IsobaricMassTagType.TMT10, new[] { MultiplexMod("TMT10", eleven) }, out _), Is.False);
        }

        [Test]
        public void AModificationWithoutHcdDiagnosticIonsIsRefused()
        {
            var noIons = new Modification("TMT10", _modificationType: "Multiplex Label", _monoisotopicMass: 229.162932);
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(IsobaricMassTagType.TMT10, new[] { noIons }, out _), Is.False);

            var cidOnly = new Modification("TMT10", _modificationType: "Multiplex Label", _monoisotopicMass: 229.162932,
                _diagnosticIons: new Dictionary<DissociationType, List<double>> { [DissociationType.CID] = new() { 125.12 } });
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(IsobaricMassTagType.TMT10, new[] { cidOnly }, out _), Is.False);

            var notMultiplex = MultiplexMod("TMT10", Resolve(IsobaricMassTagType.TMT10).ReporterIonMzs.ToList(), "Common Fixed");
            Assert.That(IsobaricMassTag.TryGetIsobaricMassTag(IsobaricMassTagType.TMT10, new[] { notMultiplex }, out _), Is.False);
        }

        [Test]
        public void ToleranceIsThreeMillidaltons()
        {
            var tag = Resolve(IsobaricMassTagType.TMT10);
            Assert.That(IsobaricMassTag.AbsoluteToleranceValue, Is.EqualTo(0.003));
            Assert.That(tag.ReporterIonMzRanges[0].Minimum, Is.EqualTo(tag.ReporterIonMzs[0] - 0.003).Within(1e-12));
            Assert.That(tag.ReporterIonMzRanges[0].Maximum, Is.EqualTo(tag.ReporterIonMzs[0] + 0.003).Within(1e-12));
        }

        [Test]
        public void IntensitiesAreReadAtEveryReporterIon()
        {
            var tag = Resolve(IsobaricMassTagType.TMT10);
            double[] intensities = { 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 };
            var spectrum = new MzSpectrum(tag.ReporterIonMzs.ToArray(), intensities, false);
            Assert.That(tag.GetReporterIonIntensities(spectrum), Is.EqualTo(intensities));
        }

        [Test]
        public void AMissingReporterIonReadsZero()
        {
            var tag = Resolve(IsobaricMassTagType.TMT10);
            var spectrum = new MzSpectrum(
                new[] { tag.ReporterIonMzs[0], tag.ReporterIonMzs[2], tag.ReporterIonMzs[5] },
                new double[] { 100, 300, 600 }, false);
            Assert.That(tag.GetReporterIonIntensities(spectrum),
                Is.EqualTo(new double[] { 100, 0, 300, 0, 0, 600, 0, 0, 0, 0 }));
        }

        [Test]
        public void AnEmptyOrNullSpectrumReadsNull()
        {
            var tag = Resolve(IsobaricMassTagType.TMT10);
            Assert.That(tag.GetReporterIonIntensities(new MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), false)), Is.Null);
            Assert.That(tag.GetReporterIonIntensities(null), Is.Null);
        }

        [Test]
        public void ASpectrumStartingAfterTheReporterIonsReadsZeros()
        {
            var tag = Resolve(IsobaricMassTagType.TMT10);
            var spectrum = new MzSpectrum(new double[] { 200, 300, 400, 500 }, new double[] { 100, 200, 300, 400 }, false);
            Assert.That(tag.GetReporterIonIntensities(spectrum), Is.EqualTo(new double[10]));
        }

        [Test]
        public void ASpectrumEndingPartWayReadsZerosAfterItsLastPeak()
        {
            var tag = Resolve(IsobaricMassTagType.TMT10);
            var spectrum = new MzSpectrum(new[] { tag.ReporterIonMzs[0], tag.ReporterIonMzs[1] }, new double[] { 100, 200 }, false);
            Assert.That(tag.GetReporterIonIntensities(spectrum),
                Is.EqualTo(new double[] { 100, 200, 0, 0, 0, 0, 0, 0, 0, 0 }));
        }

        [Test]
        public void PeaksWithinToleranceCountAndPeaksBeyondItDoNot()
        {
            var tag = Resolve(IsobaricMassTagType.TMT10);
            double[] intensities = { 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000 };

            var inside = new MzSpectrum(tag.ReporterIonMzs.Select(mz => mz + 0.002).ToArray(), intensities, false);
            Assert.That(tag.GetReporterIonIntensities(inside), Is.EqualTo(intensities));

            var outside = new MzSpectrum(tag.ReporterIonMzs.Select(mz => mz + 0.00301).ToArray(), intensities, false);
            Assert.That(tag.GetReporterIonIntensities(outside), Is.EqualTo(new double[10]));
        }

        [Test]
        public void TheMostIntensePeakWithinToleranceWins()
        {
            var tag = Resolve(IsobaricMassTagType.TMT10);
            var mzs = new List<double> { tag.ReporterIonMzs[0] - 0.002 };
            var intensities = new List<double> { 50 };
            for (int i = 0; i < tag.ReporterIonMzs.Length; i++)
            {
                mzs.Add(tag.ReporterIonMzs[i]);
                intensities.Add((i + 1) * 100);
                mzs.Add(tag.ReporterIonMzs[i] + 0.001);
                intensities.Add(25);
            }
            mzs.Add(tag.ReporterIonMzs[^1] + 0.002);
            intensities.Add(75);

            var result = tag.GetReporterIonIntensities(new MzSpectrum(mzs.ToArray(), intensities.ToArray(), false));
            Assert.That(result, Is.EqualTo(Enumerable.Range(1, 10).Select(i => i * 100.0).ToArray()));
        }

        [Test]
        public void GetReporterIonLabelsRefusesAnUndefinedKit()
        {
            Assert.That(() => IsobaricMassTag.GetReporterIonLabels((IsobaricMassTagType)99), Throws.TypeOf<ArgumentOutOfRangeException>());
            Assert.That(() => IsobaricMassTag.TryGetIsobaricMassTag((IsobaricMassTagType)99, out _), Throws.TypeOf<ArgumentOutOfRangeException>());
        }

        /// <summary>
        /// A Multiplex Label modification whose HCD diagnostic ions put reporter ions at <paramref name="reporterMzs"/>.
        /// </summary>
        private static Modification MultiplexMod(string id, List<double> reporterMzs, string modificationType = "Multiplex Label")
        {
            var neutralMasses = reporterMzs.Select(mz => mz - 1.00727646688).ToList();
            return new Modification(id, _modificationType: modificationType, _monoisotopicMass: 229.162932,
                _diagnosticIons: new Dictionary<DissociationType, List<double>> { [DissociationType.HCD] = neutralMasses });
        }
    }
}
