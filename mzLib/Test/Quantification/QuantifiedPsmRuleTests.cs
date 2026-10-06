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
    /// <see cref="QuantifiedPsmRule"/>, and <see cref="MzLibExtensions.MakeQuantifiedIdentifications"/> applying it to
    /// result files. The file tests rewrite one column of a small MetaMorpheus search (BottomUpExample.psmtsv, eight
    /// PSMs that all pass) so that each case fails the rule in exactly one way.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class QuantifiedPsmRuleTests
    {
        #region Rule

        [TestCase(0.009, true)]
        [TestCase(0.01, false)] // strictly below
        [TestCase(0.02, false)]
        public static void ThePepTierReadsOnlyThePepQValue(double pepQValue, bool passes)
        {
            Assert.That(QuantifiedPsmRule.PassesConfidence(qValue: 0.5, notchQValue: 0.5, pepQValue, usePepQValue: true), Is.EqualTo(passes));
        }

        [TestCase(0.009, 0.009, true)]
        [TestCase(0.009, 0.01, false)] // the notch must pass too
        [TestCase(0.01, 0.009, false)] // strictly below
        [TestCase(0.009, null, true)]  // no notch reported: the q-value alone
        public static void WithoutPepTheQValueAndNotchMustBothPass(double qValue, double? notchQValue, bool passes)
        {
            Assert.That(QuantifiedPsmRule.PassesConfidence(qValue, notchQValue, pepQValue: 0.0, usePepQValue: false), Is.EqualTo(passes));
        }

        [Test]
        public static void APepTierWithNoPepQValueDropsTheMatch()
        {
            Assert.That(QuantifiedPsmRule.PassesConfidence(0.0, 0.0, pepQValue: null, usePepQValue: true), Is.False);
            Assert.That(QuantifiedPsmRule.PassesConfidence(0.0, 0.0, pepQValue: double.NaN, usePepQValue: true), Is.False);
        }

        [Test]
        public static void PepIsUsableOnlyWhenItWasTrained()
        {
            Assert.That(QuantifiedPsmRule.PepQValueIsUsable(new[] { 2.0, 2.0 }), Is.False, "MetaMorpheus writes 2 when PEP is untrained");
            Assert.That(QuantifiedPsmRule.PepQValueIsUsable(new[] { double.NaN, double.NaN }), Is.False, "a reader fills NaN when the column is absent");
            Assert.That(QuantifiedPsmRule.PepQValueIsUsable(Array.Empty<double>()), Is.False);
            Assert.That(QuantifiedPsmRule.PepQValueIsUsable(new[] { 2.0, 0.3 }), Is.True);
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

        /// <summary>
        /// Writes BottomUpExample.psmtsv with <paramref name="edit"/> applied to its rows, column by header name, and
        /// returns the quantified identifications made from it.
        /// </summary>
        private static List<Identification> Quantify(string name, Action<Func<string[], string, string>, Action<string[], string, string>, List<string[]>> edit,
            List<SpectraFileInfo> spectraFiles = null, Func<string, IQuantifiableResultFile> open = null)
        {
            string[] lines = File.ReadAllLines(Source);
            string[] header = lines[0].Split('\t');
            var rows = lines.Skip(1).Select(line => line.Split('\t')).ToList();
            edit((row, column) => row[Array.IndexOf(header, column)], (row, column, value) => row[Array.IndexOf(header, column)] = value, rows);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "QuantifiedPsmRule_" + name + ".psmtsv");
            File.WriteAllLines(path, new[] { lines[0] }.Concat(rows.Select(row => string.Join('\t', row))));
            try
            {
                return (open ?? FileReader.ReadQuantifiableResultFile)(path).MakeQuantifiedIdentifications(spectraFiles ?? BothFiles);
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public static void ASearchWhoseMatchesAllPassKeepsThemAll()
        {
            Assert.That(Quantify("allPass", (get, set, rows) => { }).Count, Is.EqualTo(8));
        }

        /// <summary>
        /// With PEP untrained (every PEP_QValue 2), the q-value and notch decide, strictly below 0.01, and an ambiguous
        /// match is dropped. Filtering on the untrained PEP would have dropped all eight.
        /// </summary>
        [Test]
        public static void WithPepUntrainedTheQValueAndNotchDecide()
        {
            var ids = Quantify("untrained", UntrainedEdit);

            Assert.That(ids.Count, Is.EqualTo(5));
        }

        private static readonly Action<Func<string[], string, string>, Action<string[], string, string>, List<string[]>> UntrainedEdit =
            (get, set, rows) =>
            {
                rows.ForEach(row => set(row, "PEP_QValue", "2"));
                set(rows[0], "QValue", "0.01");                         // at the threshold: dropped
                set(rows[1], "QValue Notch", "0.05");                   // notch fails: dropped
                set(rows[2], "Full Sequence", get(rows[2], "Full Sequence") + "|" + get(rows[2], "Full Sequence")); // ambiguous: dropped
            };

        /// <summary>
        /// The same file gives the same quantified set whichever reader opens it: the lightweight reader reads the
        /// notch q-value too, so the row whose notch fails is dropped there as well.
        /// </summary>
        [Test]
        public static void TheLightweightReaderQuantifiesTheSameMatches()
        {
            var ids = Quantify("untrainedLight", UntrainedEdit, open: path => new LightWeightSpectralMatchFile(path));

            Assert.That(ids.Count, Is.EqualTo(5));
        }

        [Test]
        public static void TheLightweightReaderCarriesQValueAndScore()
        {
            var ids = Quantify("trainedLight", (get, set, rows) => { }, open: path => new LightWeightSpectralMatchFile(path));
            var heavy = Quantify("trainedHeavy", (get, set, rows) => { });

            Assert.That(ids.Select(id => (id.QValue, id.PsmScore)), Is.EqualTo(heavy.Select(id => (id.QValue, id.PsmScore))));
            Assert.That(ids.All(id => id.PsmScore > 0), Is.True);
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
        /// A MetaMorpheus file with no q-value column and an untrained PEP would have every match dropped in
        /// silence; it is refused instead.
        /// </summary>
        [Test]
        public static void AFileWithNoQValueColumnIsRefused()
        {
            string[] lines = File.ReadAllLines(Source);
            int column = Array.IndexOf(lines[0].Split('	'), "QValue");
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "QuantifiedPsmRule_noQValue.psmtsv");
            File.WriteAllLines(path, lines.Select(line =>
            {
                var cells = line.Split('	').ToList();
                cells.RemoveAt(column);
                return string.Join('	', cells);
            }).Select((line, i) => i == 0 ? line.Replace("PEP_QValue", "PEP_QValue_untrained") : line));
            try
            {
                var file = new LightWeightSpectralMatchFile(path);
                Assert.That(() => file.MakeQuantifiedIdentifications(BothFiles), Throws.TypeOf<MzLibUtil.MzLibException>());
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>With PEP trained, the PEP q-value alone decides: a poor q-value does not drop a match.</summary>
        [Test]
        public static void WithPepTrainedThePepQValueDecides()
        {
            var ids = Quantify("trained", (get, set, rows) =>
            {
                set(rows[0], "PEP_QValue", "0.01"); // at the threshold: dropped
                set(rows[1], "QValue", "0.5");      // ignored in the PEP tier: kept
                set(rows[1], "QValue Notch", "0.5");
            });

            Assert.That(ids.Count, Is.EqualTo(7));
        }

        [Test]
        public static void MatchesFromASpectraFileNotSuppliedAreSkipped()
        {
            var ids = Quantify("oneFile", (get, set, rows) => { },
                new List<SpectraFileInfo> { new SpectraFileInfo("04-30-13_CAST_Frac4_6uL.raw", "A", 0, 0, 0) });

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

            var quantified = file.MakeQuantifiedIdentifications(spectraFiles);
            int expected = file.Count(p => p.GlobalQValue < QuantifiedPsmRule.DefaultThreshold);

            Assert.That(expected, Is.GreaterThan(0).And.LessThan(file.Count()), "premise: the global q-value keeps some rows and drops some");
            Assert.That(quantified.Count, Is.EqualTo(expected));
        }

        /// <summary>The existing reader is unchanged: it filters nothing.</summary>
        [Test]
        public static void MakeIdentificationsStillKeepsEveryMatch()
        {
            IQuantifiableResultFile file = FileReader.ReadQuantifiableResultFile(Source);
            Assert.That(file.MakeIdentifications(BothFiles).Count, Is.EqualTo(8));
        }

        #endregion
    }
}
