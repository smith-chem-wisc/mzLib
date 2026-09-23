using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Tests for finding a deposit's conditions as file-name tokens its own PRIDE record also uses. The
    /// cases come from the blind benchmark's fresh set (sdrf project, 2026-09-23), where the draft lost on
    /// factors because the conditions spanned several name shapes, and over-read animal IDs as factors.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestSdrfRecordAnchors
    {
        private static string[] Record(params string[] text) => text;

        [Test]
        public void ConditionsTheRecordNamesAreFoundAcrossDifferentlyShapedNames()
        {
            var names = new List<string>();
            foreach (var arm in new[] { "CTRL", "IL1BETA", "LPS", "TNF" })
                names.AddRange(Enumerable.Range(1, 3).Select(i => $"Astro_{arm}_R{i}.raw"));
            names.AddRange(new[] { "PHA_CTRL_R1.raw", "PHA_CTRL_R2.raw", "PHA_LPS_R1.raw", "PHA_LPS_R2.raw" });
            var record = Record("Astrocytes were exposed to TNF, IL-1β, and LPS or left as control (CTRL).",
                "Primary human astrocytes (PHA) were treated the same way.");

            var a = SdrfRecordAnchors.Read(names, record);

            var treatment = a.Factors.Single(f => f.Levels.Contains("LPS", StringComparer.OrdinalIgnoreCase));
            Assert.That(treatment.Levels.Select(l => l.ToUpperInvariant()), Is.EquivalentTo(new[] { "CTRL", "IL1BETA", "LPS", "TNF" }));
            int k = a.Factors.ToList().IndexOf(treatment);
            Assert.That(a.LevelsByFile["PHA_LPS_R2.raw"][k], Is.EqualTo("LPS").IgnoreCase, "one factor across both name shapes");
            Assert.That(a.LevelsByFile["Astro_IL1BETA_R3.raw"][k], Is.EqualTo("IL1BETA").IgnoreCase, "IL-1β in the text anchors IL1BETA");
            Assert.That(treatment.Evidence, Does.Contain("record"));
        }

        [Test]
        public void AnAnimalIdTheRecordNeverMentionsIsNotACondition()
        {
            var names = new[] { "sham_0316_1", "sham_0321_1", "TAC_0317_1", "TAC_0322_1", "sham_0316_2", "TAC_0317_2" };
            var record = Record("Proteomes of healthy and diseased murine hearts after TAC or sham surgery.");

            var a = SdrfRecordAnchors.Read(names, record);

            Assert.That(a.Factors.Count, Is.EqualTo(1));
            Assert.That(a.Factors[0].Levels, Is.EquivalentTo(new[] { "sham", "TAC" }));
        }

        [Test]
        public void CommonWordsAndInstrumentWordsAreNeverAnchors()
        {
            var names = new[] { "as_01", "as_02", "c_01", "c_02", "DIA_x_1", "DIA_x_2", "DDA_x_1", "DDA_x_2" };
            var record = Record("Samples were prepared as described and acquired in DIA and DDA mode on an Orbitrap.");

            var a = SdrfRecordAnchors.Read(names, record);

            Assert.That(a.Factors, Is.Empty);
        }

        [Test]
        public void ALevelOnOneFileOnlyOrOnEveryFileIsNotACondition()
        {
            var names = new[] { "HeLa_LPS_1", "HeLa_LPS_2", "HeLa_TNF_1", "HeLa_TNF_2", "HeLa_IL6_1" };
            var record = Record("HeLa cells were treated with LPS, TNF or IL6.");

            var a = SdrfRecordAnchors.Read(names, record);

            var f = a.Factors.Single();
            Assert.That(f.Levels, Is.EquivalentTo(new[] { "LPS", "TNF" }), "IL6 names one file; HeLa names every file");
            Assert.That(a.LevelsByFile["HeLa_IL6_1"][0], Is.Empty, "a file with no level of a factor is left blank, not guessed");
        }

        [Test]
        public void NothingAnchoredMeansNoFactors()
        {
            var a = SdrfRecordAnchors.Read(new[] { "run_01", "run_02" }, Record("A proteome."));

            Assert.That(a.Factors, Is.Empty);
            Assert.That(a.LevelsByFile["run_01"], Is.Empty);
        }

        [Test]
        public void MalformedArgumentsThrow()
        {
            Assert.Throws<ArgumentNullException>(() => SdrfRecordAnchors.Read(null!, Record("x")));
            Assert.Throws<ArgumentNullException>(() => SdrfRecordAnchors.Read(new[] { "a" }, null!));
        }
    }
}
