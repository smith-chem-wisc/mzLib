using FlashLFQ;
using MassSpectrometry;
using NUnit.Framework;
using Quantification;
using Readers;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using Test.FlashLFQ;
using System.IO;
using System.Linq;

namespace Test.Quantification
{
    /// <summary>
    /// <see cref="QuantifiedPsmRule"/>, and <see cref="MzLibExtensions.MakeQuantifiedIdentifications(IQuantifiableResultFile, List{SpectraFileInfo}, out QuantifiedPsmTier, double)"/>
    /// applying it to result files. The file tests rewrite columns of a small MetaMorpheus search
    /// (BottomUpExample.psmtsv, eight target PSMs that all pass, PEP trained) so that each case picks one tier and
    /// fails it in exactly one way.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class QuantifiedPsmRuleTests
    {
        #region Rule

        private static readonly double[] Trained = { 0.0003, 0.006 };
        private static readonly double[] TrainedPeps = { 2.4e-7, 0.065 };
        private static readonly double[] Notch = { 0.0001 };
        private static readonly double[] None = Array.Empty<double>();

        [Test]
        public static void TrainedPepChoosesThePepTier()
        {
            Assert.That(QuantifiedPsmRule.ChooseTier(Trained, Notch, TrainedPeps), Is.EqualTo(QuantifiedPsmTier.PepQValue));
            Assert.That(QuantifiedPsmRule.ChooseTier(Trained, None, TrainedPeps), Is.EqualTo(QuantifiedPsmTier.PepQValue));
        }

        [Test]
        public static void UntrainedPepFallsBackToTheNotch()
        {
            Assert.That(QuantifiedPsmRule.ChooseTier(new[] { 2.0, 2.0 }, Notch, TrainedPeps), Is.EqualTo(QuantifiedPsmTier.QValueNotch),
                "MetaMorpheus writes 2 when PEP is untrained");
            Assert.That(QuantifiedPsmRule.ChooseTier(new[] { double.NaN, double.NaN }, Notch), Is.EqualTo(QuantifiedPsmTier.QValueNotch),
                "a reader fills NaN when the column is absent");
            Assert.That(QuantifiedPsmRule.ChooseTier(None, Notch), Is.EqualTo(QuantifiedPsmTier.QValueNotch));
        }

        /// <summary>
        /// When PEP training fails, MetaMorpheus leaves PEP 0 for every match but still writes PEP q-values in [0, 1].
        /// Those q-values look usable; the identical PEPs show they are not.
        /// </summary>
        [Test]
        public static void FailedPepTrainingFallsBackToTheNotch()
        {
            Assert.That(QuantifiedPsmRule.ChooseTier(Trained, Notch, new[] { 0.0, 0.0, 0.0 }), Is.EqualTo(QuantifiedPsmTier.QValueNotch));
            Assert.That(QuantifiedPsmRule.ChooseTier(Trained, Notch, new[] { 0.3 }), Is.EqualTo(QuantifiedPsmTier.QValueNotch),
                "a single PEP value cannot show training");
            Assert.That(QuantifiedPsmRule.ChooseTier(Trained, Notch, new[] { double.NaN, double.NaN }), Is.EqualTo(QuantifiedPsmTier.PepQValue),
                "no real PEP to check: the PEP q-values decide");
            Assert.That(QuantifiedPsmRule.ChooseTier(Trained, Notch), Is.EqualTo(QuantifiedPsmTier.PepQValue),
                "a reader that exposes no PEP cannot show the failure");
        }

        [Test]
        public static void NoUsableNotchFallsBackToTheQValue()
        {
            Assert.That(QuantifiedPsmRule.ChooseTier(new[] { 2.0 }, None), Is.EqualTo(QuantifiedPsmTier.QValue));
            Assert.That(QuantifiedPsmRule.ChooseTier(new[] { 2.0 }, new[] { double.NaN, 2.0 }), Is.EqualTo(QuantifiedPsmTier.QValue));
            Assert.That(QuantifiedPsmRule.ChooseTier(None, None), Is.EqualTo(QuantifiedPsmTier.QValue));
        }

        [Test]
        public static void UsabilityNeedsOneValueInTheUnitInterval()
        {
            Assert.That(QuantifiedPsmRule.PepIsUsable(new[] { 2.0, 0.3 }), Is.True);
            Assert.That(QuantifiedPsmRule.PepIsUsable(new[] { 2.0, 0.0 }, new[] { 0.1, 0.2 }), Is.True, "0 and 1 are in range");
            Assert.That(QuantifiedPsmRule.PepIsUsable(new[] { -0.1, 1.1 }), Is.False);
            Assert.That(QuantifiedPsmRule.NotchIsUsable(new[] { 2.0, 1.0 }), Is.True);
            Assert.That(QuantifiedPsmRule.NotchIsUsable(new[] { 2.0, double.NaN }), Is.False);
        }

        /// <summary>Each tier reads its own value only, strictly below the threshold.</summary>
        [TestCase(QuantifiedPsmTier.PepQValue, 0.5, 0.5, 0.009, true)]
        [TestCase(QuantifiedPsmTier.PepQValue, 0.0, 0.0, 0.01, false)]
        [TestCase(QuantifiedPsmTier.QValueNotch, 0.5, 0.009, 0.5, true)]
        [TestCase(QuantifiedPsmTier.QValueNotch, 0.0, 0.01, 0.0, false)]
        [TestCase(QuantifiedPsmTier.QValue, 0.009, 0.5, 0.5, true)]
        [TestCase(QuantifiedPsmTier.QValue, 0.01, 0.0, 0.0, false)]
        public static void EachTierReadsOnlyItsValue(QuantifiedPsmTier tier, double qValue, double notchQValue, double pepQValue, bool passes)
        {
            Assert.That(QuantifiedPsmRule.PassesConfidence(tier, qValue, notchQValue, pepQValue), Is.EqualTo(passes));
        }

        [Test]
        public static void AMatchWithNoValueForItsTierIsDropped()
        {
            Assert.That(QuantifiedPsmRule.PassesConfidence(QuantifiedPsmTier.PepQValue, 0.0, 0.0, pepQValue: null), Is.False);
            Assert.That(QuantifiedPsmRule.PassesConfidence(QuantifiedPsmTier.PepQValue, 0.0, 0.0, pepQValue: double.NaN), Is.False);
            Assert.That(QuantifiedPsmRule.PassesConfidence(QuantifiedPsmTier.QValueNotch, 0.0, notchQValue: null, 0.0), Is.False);
            Assert.That(QuantifiedPsmRule.PassesConfidence(QuantifiedPsmTier.QValue, double.NaN, 0.0, 0.0), Is.False);
        }

        [TestCase("PEPTIDE", "PEPTIDE", false)]
        [TestCase("PEPTIDE", "PEPT[Common Biological:Phosphorylation on T]IDE|PEPTIDE", true)]
        [TestCase("PEPTIDE|PEPTLDE", "PEPTIDE|PEPTLDE", true)]
        [TestCase(null, "PEPTIDE", true)]
        [TestCase("PEPTIDE", "", true)]
        public static void AmbiguousMeansMoreThanOneSequenceOrForm(string baseSequence, string fullSequence, bool ambiguous)
        {
            Assert.That(QuantifiedPsmRule.IsAmbiguous(baseSequence, fullSequence), Is.EqualTo(ambiguous));
        }

        #endregion

        #region Result files

        private static string Source => Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\SearchResults", "BottomUpExample.psmtsv");

        private static List<SpectraFileInfo> BothFiles => new()
        {
            new SpectraFileInfo("04-30-13_CAST_Frac5_4uL.raw", "A", 0, 0, 0),
            new SpectraFileInfo("04-30-13_CAST_Frac4_6uL.raw", "A", 1, 0, 0),
        };

        private delegate void Edit(Func<string[], string, string> get, Action<string[], string, string> set, List<string[]> rows);

        private static readonly Func<string, IQuantifiableResultFile> Full = FileReader.ReadQuantifiableResultFile;
        private static readonly Func<string, IQuantifiableResultFile> Lightweight = path => new LightWeightSpectralMatchFile(path);

        /// <summary>
        /// Writes BottomUpExample.psmtsv with <paramref name="edit"/> applied to its rows (column by header name) and
        /// the <paramref name="drop"/> columns removed, opens it with <paramref name="open"/>, and returns the
        /// quantified identifications made from it and the tier chosen.
        /// </summary>
        private static List<Identification> Quantify(string name, Edit edit, out QuantifiedPsmTier tier,
            Func<string, IQuantifiableResultFile> open = null, List<SpectraFileInfo> spectraFiles = null, params string[] drop)
        {
            string path = Write(name, edit, drop);
            try
            {
                return (open ?? Full)(path).MakeQuantifiedIdentifications(spectraFiles ?? BothFiles, out tier);
            }
            finally
            {
                File.Delete(path);
            }
        }

        private static string Write(string name, Edit edit, params string[] drop)
        {
            string[] lines = File.ReadAllLines(Source);
            string[] header = lines[0].Split('\t');
            var rows = lines.Skip(1).Select(line => line.Split('\t')).ToList();
            edit((row, column) => row[Array.IndexOf(header, column)], (row, column, value) => row[Array.IndexOf(header, column)] = value, rows);

            var kept = Enumerable.Range(0, header.Length).Where(i => !drop.Contains(header[i])).ToList();
            string Join(string[] cells) => string.Join('\t', kept.Select(i => cells[i]));

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "QuantifiedPsmRule_" + name + ".psmtsv");
            File.WriteAllLines(path, new[] { Join(header) }.Concat(rows.Select(Join)));
            return path;
        }

        private static readonly Edit NoEdit = (get, set, rows) => { };

        /// <summary>PEP untrained: MetaMorpheus writes 2 for every PEP q-value.</summary>
        private static readonly Edit Untrained = (get, set, rows) => rows.ForEach(row => set(row, "PEP_QValue", "2"));

        [TestCase("Full")]
        [TestCase("Lightweight")]
        public static void ASearchWhoseMatchesAllPassKeepsThemAllOnThePepTier(string reader)
        {
            var ids = Quantify("allPass" + reader, NoEdit, out var tier, Open(reader));

            Assert.That(tier, Is.EqualTo(QuantifiedPsmTier.PepQValue));
            Assert.That(ids.Count, Is.EqualTo(8));
        }

        /// <summary>With PEP trained, the PEP q-value alone decides: a poor q-value or notch does not drop a match.</summary>
        [TestCase("Full")]
        [TestCase("Lightweight")]
        public static void WithPepTrainedThePepQValueDecides(string reader)
        {
            var ids = Quantify("trained" + reader, (get, set, rows) =>
            {
                set(rows[0], "PEP_QValue", "0.01"); // at the threshold: dropped
                set(rows[1], "QValue", "0.5");      // ignored in the PEP tier: kept
                set(rows[1], "QValue Notch", "0.5");
            }, out var tier, Open(reader));

            Assert.That(tier, Is.EqualTo(QuantifiedPsmTier.PepQValue));
            Assert.That(ids.Count, Is.EqualTo(7));
        }

        /// <summary>
        /// With PEP untrained, the notch q-value alone decides, strictly below 0.01, and an ambiguous match is dropped.
        /// The same file gives the same set whichever reader opens it. Filtering on the untrained PEP would have
        /// dropped all eight.
        /// </summary>
        [TestCase("Full")]
        [TestCase("Lightweight")]
        public static void WithPepUntrainedTheNotchDecides(string reader)
        {
            var ids = Quantify("untrained" + reader, (get, set, rows) =>
            {
                Untrained(get, set, rows);
                set(rows[0], "QValue", "0.5");          // ignored in the notch tier: kept
                set(rows[1], "QValue Notch", "0.01");   // at the threshold: dropped
                set(rows[2], "Full Sequence", get(rows[2], "Full Sequence") + "|" + get(rows[2], "Full Sequence")); // ambiguous: dropped
            }, out var tier, Open(reader));

            Assert.That(tier, Is.EqualTo(QuantifiedPsmTier.QValueNotch));
            Assert.That(ids.Count, Is.EqualTo(6));
        }

        /// <summary>
        /// PEP training failed: every PEP is 0, yet the PEP q-values are in range. Both readers see it from PEP and
        /// fall back to the notch.
        /// </summary>
        [TestCase("Full")]
        [TestCase("Lightweight")]
        public static void WithPepTrainingFailedTheNotchDecides(string reader)
        {
            var ids = Quantify("failed" + reader, (get, set, rows) =>
            {
                rows.ForEach(row => set(row, "PEP", "0"));
                set(rows[0], "PEP_QValue", "0.5");      // ignored in the notch tier: kept
                set(rows[1], "QValue Notch", "0.05");   // notch fails: dropped
            }, out var tier, Open(reader));

            Assert.That(tier, Is.EqualTo(QuantifiedPsmTier.QValueNotch));
            Assert.That(ids.Count, Is.EqualTo(7));
        }

        /// <summary>With PEP untrained and no notch column, the q-value decides, strictly below 0.01.</summary>
        [TestCase("Full")]
        [TestCase("Lightweight")]
        public static void WithNoNotchTheQValueDecides(string reader)
        {
            var ids = Quantify("noNotch" + reader, (get, set, rows) =>
            {
                Untrained(get, set, rows);
                set(rows[0], "QValue", "0.01");    // at the threshold: dropped
                set(rows[1], "QValue", "0.0099");  // just below: kept
            }, out var tier, Open(reader), drop: "QValue Notch");

            Assert.That(tier, Is.EqualTo(QuantifiedPsmTier.QValue));
            Assert.That(ids.Count, Is.EqualTo(7));
        }

        /// <summary>
        /// The tier is decided per result file: the same match (a notch of 0.05) is kept in a file whose PEP was trained
        /// and dropped in one whose PEP was not.
        /// </summary>
        [Test]
        public static void EachResultFileChoosesItsOwnTier()
        {
            Edit notchFails = (get, set, rows) => set(rows[1], "QValue Notch", "0.05");
            var trained = Quantify("mixedTrained", notchFails, out var trainedTier);
            var untrained = Quantify("mixedUntrained", (get, set, rows) => { Untrained(get, set, rows); notchFails(get, set, rows); }, out var untrainedTier);

            Assert.That((trainedTier, trained.Count), Is.EqualTo((QuantifiedPsmTier.PepQValue, 8)));
            Assert.That((untrainedTier, untrained.Count), Is.EqualTo((QuantifiedPsmTier.QValueNotch, 7)));
        }

        /// <summary>
        /// The tier is chosen from target matches only: a decoy's PEP q-value in range does not make an untrained
        /// PEP usable. The decoy is still kept, for match-between-runs.
        /// </summary>
        [Test]
        public static void TheTierIsChosenFromTargetsOnly()
        {
            var ids = Quantify("decoyPep", (get, set, rows) =>
            {
                Untrained(get, set, rows);
                set(rows[7], "Decoy/Contaminant/Target", "D");
                set(rows[7], "PEP_QValue", "0.001");
            }, out var tier);

            Assert.That(tier, Is.EqualTo(QuantifiedPsmTier.QValueNotch));
            Assert.That(ids.Count(id => id.IsDecoy), Is.EqualTo(1));
        }

        /// <summary>Each identification carries the value its tier filtered on.</summary>
        [Test]
        public static void AnIdentificationCarriesItsTiersValue()
        {
            var pep = Quantify("valuePep", NoEdit, out _);
            var notch = Quantify("valueNotch", Untrained, out _);
            var q = Quantify("valueQ", Untrained, out _, drop: "QValue Notch");

            Assert.That(pep[0].QValue, Is.EqualTo(0.000314465));
            Assert.That(notch[0].QValue, Is.EqualTo(0.000124));
            Assert.That(q[0].QValue, Is.EqualTo(0.000121));
        }

        [Test]
        public static void TheLightweightReaderCarriesQValueAndScore()
        {
            var ids = Quantify("trainedLight", NoEdit, out _, Lightweight);
            var heavy = Quantify("trainedHeavy", NoEdit, out _);

            Assert.That(ids.Select(id => (id.QValue, id.PsmScore)), Is.EqualTo(heavy.Select(id => (id.QValue, id.PsmScore))));
            Assert.That(ids.All(id => id.PsmScore > 0), Is.True);
        }

        [Test]
        public static void TheLightweightReaderReadsPep()
        {
            string path = Write("lightPep", NoEdit);
            try
            {
                var light = new LightWeightSpectralMatchFile(path).Cast<LightWeightSpectralMatch>().ToList();
                Assert.That(light.Select(psm => psm.Pep).Distinct().Count(), Is.GreaterThan(1));
                Assert.That(light[1].Pep, Is.EqualTo(0.06503135));
            }
            finally
            {
                File.Delete(path);
            }

            path = Write("lightNoPep", NoEdit, "PEP");
            try
            {
                Assert.That(new LightWeightSpectralMatchFile(path).Cast<LightWeightSpectralMatch>().All(psm => double.IsNaN(psm.Pep)), Is.True);
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>A record type that reports no q-value is kept: its tool filtered it.</summary>
        [Test]
        public static void ARecordWithNoQValueIsKept()
        {
            var record = new MockQuantifiableRecord
            {
                BaseSequence = "PEPTIDE",
                FullSequence = "PEPTIDE",
                RetentionTime = 5.0,
                MonoisotopicMass = 800.0,
                ChargeState = 2,
                FileName = "file1.mzML",
                ProteinGroupInfos = new List<(string, string, string)> { ("P1", "Gene1", "Organism1") },
            };
            var file = new MockQuantifiableResultFile(new List<IQuantifiableRecord> { record });

            var ids = file.MakeQuantifiedIdentifications(new List<SpectraFileInfo> { new SpectraFileInfo("file1.mzML", "A", 0, 0, 0) });

            Assert.That(ids.Count, Is.EqualTo(1));
        }

        /// <summary>
        /// A MetaMorpheus file with no q-value column, an untrained PEP and no notch would have every match dropped in
        /// silence; it is refused instead. With a notch it is not: the notch decides and the q-value is not read.
        /// (Only the lightweight reader accepts a .psmtsv with no q-value column; the full reader requires one.)
        /// </summary>
        [Test]
        public static void AFileWithNoQValueColumnIsRefusedOnlyOnTheQValueTier()
        {
            Assert.That(() => Quantify("noQ", Untrained, out _, Lightweight, drop: new[] { "QValue", "QValue Notch" }),
                Throws.TypeOf<MzLibUtil.MzLibException>());

            var ids = Quantify("noQNotch", Untrained, out var tier, Lightweight, drop: "QValue");
            Assert.That((tier, ids.Count), Is.EqualTo((QuantifiedPsmTier.QValueNotch, 8)));
        }

        [Test]
        public static void MatchesFromASpectraFileNotSuppliedAreSkipped()
        {
            var ids = Quantify("oneFile", NoEdit, out _,
                spectraFiles: new List<SpectraFileInfo> { new SpectraFileInfo("04-30-13_CAST_Frac4_6uL.raw", "A", 0, 0, 0) });

            Assert.That(ids.Count, Is.EqualTo(2));
        }

        /// <summary>
        /// DIA-NN reports no PEP q-value and no notch, so Global.Q.Value is the only tier; its raw PEP is never read
        /// as a q-value.
        /// </summary>
        [Test]
        public static void ADiaNnReportIsFilteredOnGlobalQValue()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, @"FileReadingTests\ExternalFileTypes\DiaNn_LongFormat_report.tsv");
            var file = new DiaNnReportFile(path);
            var spectraFiles = file.Select(p => p.SpectraFilePath).Distinct()
                .Select((filePath, i) => new SpectraFileInfo(filePath, "PBMC", i, 0, 0)).ToList();

            var quantified = file.MakeQuantifiedIdentifications(spectraFiles, out var tier);
            int expected = file.Count(p => p.GlobalQValue < QuantifiedPsmRule.DefaultThreshold);

            Assert.That(expected, Is.GreaterThan(0).And.LessThan(file.Count()), "premise: the global q-value keeps some rows and drops some");
            Assert.That(tier, Is.EqualTo(QuantifiedPsmTier.QValue));
            Assert.That(quantified.Count, Is.EqualTo(expected));
        }

        /// <summary>The existing reader is unchanged: it filters nothing.</summary>
        [Test]
        public static void MakeIdentificationsStillKeepsEveryMatch()
        {
            IQuantifiableResultFile file = FileReader.ReadQuantifiableResultFile(Source);
            Assert.That(file.MakeIdentifications(BothFiles).Count, Is.EqualTo(8));
        }

        private static Func<string, IQuantifiableResultFile> Open(string reader) => reader == "Lightweight" ? Lightweight : Full;

        #endregion
    }
}
