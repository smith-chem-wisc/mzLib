using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;

namespace Test.FileReadingTests
{
    /// <summary>
    /// Critical round-trip tests for <see cref="MsAlignWriter"/>.
    /// Each test writes synthetic data, reads it back through the existing
    /// <see cref="Ms1Align"/> or <see cref="Ms2Align"/> readers, and asserts
    /// that values are preserved.  This verifies that the writer output is
    /// compatible with the readers without modifying either reader.
    /// </summary>
    [TestFixture]
    public sealed class TestMsAlignWriter
    {
        private string _tempDir;

        [OneTimeSetUp]
        public void SetUp()
        {
            _tempDir = Path.Combine(
                TestContext.CurrentContext.TestDirectory,
                "MsAlignWriterTemp");
            Directory.CreateDirectory(_tempDir);
        }

        [OneTimeTearDown]
        public void TearDown()
        {
            if (Directory.Exists(_tempDir))
                Directory.Delete(_tempDir, recursive: true);
        }

        // ── helpers ───────────────────────────────────────────────────────────

        /// <summary>
        /// Builds a minimal MS1 <see cref="MsDataScan"/> with a trivial spectrum.
        /// </summary>
        private static MsDataScan MakeMs1Scan(int scanNumber, double retentionTimeMinutes)
        {
            var spectrum = new MzSpectrum(
                new double[] { 500.0 },
                new double[] { 1000.0 },
                shouldCopy: false);

            return new MsDataScan(
                spectrum,
                scanNumber, 1, true, Polarity.Positive,
                retentionTimeMinutes,
                new MzRange(400, 1600), null,
                MZAnalyzerType.Orbitrap, spectrum.SumOfAllY,
                null, null, $"scan={scanNumber}");
        }

        /// <summary>
        /// Builds a minimal MS2 <see cref="MsDataScan"/> with precursor metadata.
        /// </summary>
        private static MsDataScan MakeMs2Scan(
            int scanNumber,
            double retentionTimeMinutes,
            int precursorScanNumber,
            double precursorMz,
            int precursorCharge,
            double precursorIntensity)
        {
            var spectrum = new MzSpectrum(
                new double[] { 200.0, 300.0 },
                new double[] { 500.0, 250.0 },
                shouldCopy: false);

            return new MsDataScan(
                spectrum,
                scanNumber, 2, true, Polarity.Positive,
                retentionTimeMinutes,
                new MzRange(100, 2000), null,
                MZAnalyzerType.Orbitrap, spectrum.SumOfAllY,
                null, null, $"scan={scanNumber}",
                selectedIonMz: precursorMz,
                selectedIonChargeStateGuess: precursorCharge,
                selectedIonIntensity: precursorIntensity,
                dissociationType: DissociationType.HCD,
                oneBasedPrecursorScanNumber: precursorScanNumber,
                selectedIonMonoisotopicGuessMz: precursorMz);
        }

        /// <summary>
        /// Creates an <see cref="IsotopicEnvelope"/> using the pre-scored
        /// constructor so that the test controls all field values precisely.
        /// </summary>
        private static IsotopicEnvelope MakeEnvelope(
            double monoMass,
            int charge,
            double intensity)
        {
            double mz = monoMass.ToMz(Math.Abs(charge));
            var peaks = new List<(double mz, double intensity)> { (mz, intensity) };
            return new IsotopicEnvelope(0, peaks, monoMass, charge, intensity, 1.0);
        }

        // ══════════════════════════════════════════════════════════════════════
        // 1. MS1 — file is created and readable
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void WriteMs1Align_FileIsCreated()
        {
            var path = Path.Combine(_tempDir, "created_ms1.msalign");
            var scan = MakeMs1Scan(1, 1.0);
            var envelope = MakeEnvelope(5000.0, 5, 1000.0);

            MsAlignWriter.WriteMs1Align(
                new[] { (scan, (IEnumerable<IsotopicEnvelope>)new[] { envelope }) },
                path);

            Assert.That(File.Exists(path), Is.True);
        }

        [Test]
        public void WriteMs1Align_ExtensionAppendedWhenMissing()
        {
            var pathWithoutExt = Path.Combine(_tempDir, "noext_ms1");
            var scan = MakeMs1Scan(1, 1.0);
            var envelope = MakeEnvelope(5000.0, 5, 1000.0);

            MsAlignWriter.WriteMs1Align(
                new[] { (scan, (IEnumerable<IsotopicEnvelope>)new[] { envelope }) },
                pathWithoutExt);

            Assert.That(File.Exists(pathWithoutExt + ".msalign"), Is.True);
        }

        // ══════════════════════════════════════════════════════════════════════
        // 2. MS1 — scan count round-trip
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void WriteMs1Align_ThreeScans_ReadBackCountIsThree()
        {
            var path = Path.Combine(_tempDir, "count_ms1.msalign");

            var pairs = Enumerable.Range(1, 3).Select(i =>
            {
                var scan = MakeMs1Scan(i, i * 1.0);
                IEnumerable<IsotopicEnvelope> envs = new[] { MakeEnvelope(5000.0 * i, 5, 1000.0) };
                return (scan, envs);
            });

            MsAlignWriter.WriteMs1Align(pairs, path);

            var read = new Ms1Align(path);
            read.LoadAllStaticData();

            Assert.That(read.Scans.Length, Is.EqualTo(3));
        }

        // ══════════════════════════════════════════════════════════════════════
        // 3. MS1 — peak values survive the round-trip
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void WriteMs1Align_PeakValues_RoundTripCorrectly()
        {
            var path = Path.Combine(_tempDir, "peaks_ms1.msalign");

            double[] masses = { 5000.0, 10000.0, 15000.0 };
            int[] charges = { 5, 10, 15 };
            double[] intensities = { 1000.0, 2000.0, 3000.0 };

            var scan = MakeMs1Scan(1, 1.0);
            var envs = masses
                .Select((m, i) => MakeEnvelope(m, charges[i], intensities[i]))
                .ToArray();

            MsAlignWriter.WriteMs1Align(
                new[] { (scan, (IEnumerable<IsotopicEnvelope>)envs) },
                path);

            var read = new Ms1Align(path);
            read.LoadAllStaticData();

            var readScan = read.GetOneBasedScan(1);
            var neutral = (NeutralMassSpectrum)readScan.MassSpectrum;

            // Peaks are re-sorted by mass on read; sort expected values the same way
            var sorted = masses
                .Select((m, i) => (Mass: m, Charge: charges[i], Intensity: intensities[i]))
                .OrderBy(x => x.Mass)
                .ToArray();

            Assert.That(neutral.Size, Is.EqualTo(3));

            for (int i = 0; i < 3; i++)
            {
                Assert.That(neutral.XArray[i], Is.EqualTo(sorted[i].Mass).Within(0.001));
                Assert.That(neutral.Charges[i], Is.EqualTo(sorted[i].Charge));
                Assert.That(neutral.YArray[i], Is.EqualTo(sorted[i].Intensity).Within(0.1));
            }
        }

        // ══════════════════════════════════════════════════════════════════════
        // 4. MS1 — retention time converts minutes → seconds → minutes
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void WriteMs1Align_RetentionTime_RoundTripsCorrectly()
        {
            var path = Path.Combine(_tempDir, "rt_ms1.msalign");
            const double rtMinutes = 1.5;

            var scan = MakeMs1Scan(1, rtMinutes);
            MsAlignWriter.WriteMs1Align(
                new[] { (scan, (IEnumerable<IsotopicEnvelope>)new[] { MakeEnvelope(5000.0, 5, 1000.0) }) },
                path);

            var read = new Ms1Align(path);
            read.LoadAllStaticData();

            // Writer: minutes * 60 → seconds.  Reader: seconds / 60 → minutes.
            Assert.That(
                read.GetOneBasedScan(1).RetentionTime,
                Is.EqualTo(rtMinutes).Within(0.001));
        }

        // ══════════════════════════════════════════════════════════════════════
        // 5. MS1 — scans with no envelopes are silently skipped
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void WriteMs1Align_EmptyEnvelopes_ScanIsSkipped()
        {
            var path = Path.Combine(_tempDir, "skip_ms1.msalign");

            var pairs = new (MsDataScan, IEnumerable<IsotopicEnvelope>)[]
            {
                (MakeMs1Scan(1, 1.0), new[] { MakeEnvelope(5000.0, 5, 1000.0) }),
                (MakeMs1Scan(2, 2.0), Array.Empty<IsotopicEnvelope>()),          // skipped
                (MakeMs1Scan(3, 3.0), new[] { MakeEnvelope(6000.0, 6, 2000.0) }),
            };

            MsAlignWriter.WriteMs1Align(pairs, path);

            // Raw text check: exactly two BEGIN IONS blocks
            var content = File.ReadAllText(path);
            var beginCount = CountOccurrences(content, "BEGIN IONS");
            Assert.That(beginCount, Is.EqualTo(2));

            // Reader check: two scans loaded
            var read = new Ms1Align(path);
            read.LoadAllStaticData();
            Assert.That(read.Scans.Length, Is.EqualTo(2));
        }

        // ══════════════════════════════════════════════════════════════════════
        // 6. MS1 — negative-mode charges are written as positive integers
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void WriteMs1Align_NegativeCharge_WrittenAsPositiveInteger()
        {
            var path = Path.Combine(_tempDir, "negcharge_ms1.msalign");

            // IsotopicEnvelope pre-scored constructor: charge = -5 (negative mode)
            double monoMass = 5000.0;
            int negCharge = -5;
            double mz = monoMass.ToMz(Math.Abs(negCharge));
            var peaks = new List<(double mz, double intensity)> { (mz, 1000.0) };
            var negEnv = new IsotopicEnvelope(0, peaks, monoMass, negCharge, 1000.0, 1.0);

            var scan = MakeMs1Scan(1, 1.0);
            MsAlignWriter.WriteMs1Align(
                new[] { (scan, (IEnumerable<IsotopicEnvelope>)new[] { negEnv }) },
                path);

            // The file must contain "5000" and "5", never "-5"
            var content = File.ReadAllText(path);
            Assert.That(content, Does.Not.Contain("\t-5"));

            // Reader stores charges as positive ints in NeutralMassSpectrum
            var read = new Ms1Align(path);
            read.LoadAllStaticData();
            var neutral = (NeutralMassSpectrum)read.GetOneBasedScan(1).MassSpectrum;
            Assert.That(neutral.Charges[0], Is.EqualTo(5));
        }

        // ══════════════════════════════════════════════════════════════════════
        // 7. MS2 — precursor metadata round-trip
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void WriteMs2Align_PrecursorMetadata_RoundTripsCorrectly()
        {
            var path = Path.Combine(_tempDir, "precursor_ms2.msalign");

            const double precursorMz = 500.25;
            const int precursorCharge = 3;
            const double precursorIntensity = 12345.67;
            const int precursorScanNum = 10;

            var scan = MakeMs2Scan(11, 2.5, precursorScanNum,
                precursorMz, precursorCharge, precursorIntensity);

            MsAlignWriter.WriteMs2Align(
                new[] { (scan, (IEnumerable<IsotopicEnvelope>)new[]
                    { MakeEnvelope(1497.74, 3, precursorIntensity) }) },
                path);

            var read = new Ms2Align(path);
            read.LoadAllStaticData();

            var readScan = read.GetOneBasedScan(11);

            Assert.That(readScan.SelectedIonMZ,
                Is.EqualTo(precursorMz).Within(0.001));
            Assert.That(readScan.SelectedIonChargeStateGuess,
                Is.EqualTo(precursorCharge));
            Assert.That(readScan.OneBasedPrecursorScanNumber,
                Is.EqualTo(precursorScanNum));
            Assert.That(readScan.DissociationType,
                Is.EqualTo(DissociationType.HCD));
        }

        // ══════════════════════════════════════════════════════════════════════
        // 8. Integration — NeutralMassSpectrum short-circuit after write/read
        // ══════════════════════════════════════════════════════════════════════

        /// <summary>
        /// After writing and re-reading an msAlign file, passing the resulting
        /// <see cref="MsDataScan"/> through <see cref="Deconvoluter"/> with
        /// <see cref="FLASHDeconvolutionParameters"/> must hit the
        /// NeutralMassSpectrum short-circuit and return one envelope per peak —
        /// no re-deconvolution occurs.
        /// </summary>
        [Test]
        public void WriteReadMs1Align_ThenDeconvolute_UsesNeutralMassShortCircuit()
        {
            var path = Path.Combine(_tempDir, "shortcircuit_ms1.msalign");

            double[] masses = { 5000.0, 10000.0 };
            int[] charges = { 5, 10 };
            double[] intensities = { 1000.0, 2000.0 };

            var scan = MakeMs1Scan(1, 1.0);
            var envs = masses
                .Select((m, i) => MakeEnvelope(m, charges[i], intensities[i]))
                .ToArray();

            MsAlignWriter.WriteMs1Align(
                new[] { (scan, (IEnumerable<IsotopicEnvelope>)envs) },
                path);

            // Re-read
            var read = new Ms1Align(path);
            read.LoadAllStaticData();
            var readScan = read.GetOneBasedScan(1);

            // Deconvolute — NeutralMassSpectrum short-circuits, ignoring algorithm
            var result = Deconvoluter
                .Deconvolute(readScan, new FLASHDeconvolutionParameters())
                .ToList();

            // Expect exactly as many envelopes as peaks written
            Assert.That(result.Count, Is.EqualTo(masses.Length));

            var sortedMasses = masses.OrderBy(m => m).ToArray();
            var resultMasses = result.Select(e => e.MonoisotopicMass).OrderBy(m => m).ToArray();
            for (int i = 0; i < sortedMasses.Length; i++)
                Assert.That(resultMasses[i], Is.EqualTo(sortedMasses[i]).Within(0.001));
        }

        // ── private utility ───────────────────────────────────────────────────

        private static int CountOccurrences(string source, string pattern)
        {
            int count = 0;
            int index = 0;
            while ((index = source.IndexOf(pattern, index, StringComparison.Ordinal)) >= 0)
            {
                count++;
                index += pattern.Length;
            }
            return count;
        }
    }
}