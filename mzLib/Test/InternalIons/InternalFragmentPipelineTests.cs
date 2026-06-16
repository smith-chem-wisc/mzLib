using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using Readers.InternalIons;

namespace Test.InternalIons
{
    /// <summary>
    /// Coverage for the observed-internal-fragment analysis pipeline
    /// (InternalFragmentIon, InternalFragmentTsvWriter, InternalFragmentFeatureExtractor,
    /// InternalFragmentAnalysisRunner). Uses the committed internalIons.psmtsv and the
    /// null-tolerant MsDataFile path so no raw spectra are required.
    /// </summary>
    [TestFixture]
    public static class InternalFragmentPipelineTests
    {
        private static string PsmTsvPath => Path.Combine(
            TestContext.CurrentContext.TestDirectory,
            "FileReadingTests", "SearchResults", "internalIons.psmtsv");

        private static InternalFragmentIon SampleIon() => new()
        {
            PeptideSequence = "PEPDKTIDE",
            InternalSequence = "PDK",
            StartResidue = 3,
            EndResidue = 5,
            TheoreticalMass = 400.0,
            ObservedMass = 400.001,
            NormalizedIntensity = 0.5,
            TicNormalizedIntensity = 0.02,
            TotalIonCurrent = 1234.5,
            PrecursorCharge = 2,
            CollisionEnergy = 42.0,
            NTerminalFlankingResidue = 'E',
            CTerminalFlankingResidue = 'T',
            SourceFile = "file1",
            ScanNumber = "100",
            ModificationsInInternalFragment = "Common Fixed:Carbamidomethyl on C at position 4",
            MetalModCount = 0,
            CommonBiologicalModCount = 1,
        };

        // ---- InternalFragmentIon ----------------------------------------------------------

        [Test]
        public static void InternalFragmentIon_HeaderAndValueArity_Match()
        {
            var headers = InternalFragmentIon.GetHeaderNames();
            var values = SampleIon().GetValues();
            Assert.That(values.Length, Is.EqualTo(headers.Length),
                "GetValues must emit exactly one value per GetHeaderNames column");
            Assert.That(headers, Is.Unique);
        }

        [Test]
        public static void InternalFragmentIon_ComputedProperties_DeriveFromFields()
        {
            var ion = SampleIon();
            Assert.That(ion.FragmentLength, Is.EqualTo(3));                 // 5 - 3 + 1
            Assert.That(ion.MassError, Is.EqualTo(0.001).Within(1e-9));
            Assert.That(ion.MassErrorPpm, Is.EqualTo(2.5).Within(0.01));    // 0.001 / 400 * 1e6
            Assert.That(ion.InternalNTerminalAA, Is.EqualTo('P'));
            Assert.That(ion.InternalCTerminalAA, Is.EqualTo('K'));
            Assert.That(ion.HasProlineAtEitherTerminus, Is.True);
            Assert.That(ion.HasAspartateAtEitherTerminus, Is.False);       // D is interior, not terminal
            Assert.That(ion.NumberOfBasicResidues, Is.EqualTo(1));         // K
            Assert.That(ion.HasModifiedResidue, Is.True);
            Assert.That(ion.PassesMassAccuracyFilter, Is.True);            // |2.5| ppm < 5
        }

        [Test]
        public static void InternalFragmentIon_PassesMassAccuracyFilter_FalseForLargeError()
        {
            var ion = new InternalFragmentIon
            {
                InternalSequence = "AAA",
                TheoreticalMass = 400.0,
                ObservedMass = 400.02,   // 50 ppm
                NormalizedIntensity = 0.0,
            };
            Assert.That(ion.PassesMassAccuracyFilter, Is.False);
        }

        // ---- InternalFragmentTsvWriter ----------------------------------------------------

        [Test]
        public static void InternalFragmentTsvWriter_RoundTrip_PreservesFields()
        {
            var path = Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"internalFragmentRoundTrip_{Guid.NewGuid():N}.tsv");
            try
            {
                var original = SampleIon();
                InternalFragmentTsvWriter.WriteToTsv(new List<InternalFragmentIon> { original }, path);

                var roundTripped = InternalFragmentTsvWriter.ReadFromTsv(path);
                Assert.That(roundTripped.Count, Is.EqualTo(1));
                var r = roundTripped[0];
                Assert.That(r.PeptideSequence, Is.EqualTo(original.PeptideSequence));
                Assert.That(r.InternalSequence, Is.EqualTo(original.InternalSequence));
                Assert.That(r.StartResidue, Is.EqualTo(original.StartResidue));
                Assert.That(r.EndResidue, Is.EqualTo(original.EndResidue));
                Assert.That(r.PrecursorCharge, Is.EqualTo(original.PrecursorCharge));
                Assert.That(r.TheoreticalMass, Is.EqualTo(original.TheoreticalMass).Within(1e-6));
                Assert.That(r.NTerminalFlankingResidue, Is.EqualTo(original.NTerminalFlankingResidue));
                Assert.That(r.ModificationsInInternalFragment, Is.EqualTo(original.ModificationsInInternalFragment));
                Assert.That(r.CommonBiologicalModCount, Is.EqualTo(original.CommonBiologicalModCount));
            }
            finally
            {
                if (File.Exists(path)) File.Delete(path);
            }
        }

        [Test]
        public static void InternalFragmentTsvWriter_ReadFromTsv_MissingFile_Throws()
        {
            var missing = Path.Combine(TestContext.CurrentContext.TestDirectory,
                $"does_not_exist_{Guid.NewGuid():N}.tsv");
            Assert.That(() => InternalFragmentTsvWriter.ReadFromTsv(missing),
                Throws.TypeOf<FileNotFoundException>());
        }

        // ---- InternalFragmentFeatureExtractor ---------------------------------------------

        [Test]
        public static void InternalFragmentFeatureExtractor_ExtractFromPsms_EmptyInput_ReturnsEmpty()
        {
            var result = InternalFragmentFeatureExtractor.ExtractFromPsms(
                new List<PsmFromTsv>(), msDataFile: null, defaultCollisionEnergy: 42.0);
            Assert.That(result, Is.Empty);
        }

        [Test]
        public static void InternalFragmentFeatureExtractor_ExtractFromPsms_PopulatesInternalFragments()
        {
            Assert.That(File.Exists(PsmTsvPath), Is.True, $"Missing test data: {PsmTsvPath}");

            var psms = SpectrumMatchTsvReader.ReadPsmTsv(PsmTsvPath, out _);
            Assert.That(psms, Is.Not.Empty, "internalIons.psmtsv should contain PSMs");

            // Null MsDataFile is tolerated (BuildScanLookup short-circuits); extraction runs on the
            // PSMs' matched-ion data with TIC/CE fallbacks.
            var ions = InternalFragmentFeatureExtractor.ExtractFromPsms(psms, null, defaultCollisionEnergy: 42.0);

            Assert.That(ions, Is.Not.Empty, "expected internal fragment ions from the PSM matched ions");
            var ion = ions[0];
            Assert.That(ion.StartResidue, Is.GreaterThan(0));
            Assert.That(ion.EndResidue, Is.GreaterThanOrEqualTo(ion.StartResidue));
            Assert.That(ion.InternalSequence, Is.Not.Empty);
            Assert.That(ion.CollisionEnergy, Is.EqualTo(42.0));        // fallback used (no scan)
            Assert.That(ions.All(i => i.PeptideSequence.Length > 0), Is.True);
        }

        // ---- InternalFragmentAnalysisRunner -----------------------------------------------

        [Test]
        public static void InternalFragmentAnalysisRunner_Run_MissingPsmTsv_Throws()
        {
            var missing = Path.Combine(TestContext.CurrentContext.TestDirectory, "no_such.psmtsv");
            var outDir = Path.Combine(TestContext.CurrentContext.TestDirectory, $"ifrun_out_{Guid.NewGuid():N}");
            Assert.That(() => InternalFragmentAnalysisRunner.Run(missing, TestContext.CurrentContext.TestDirectory, outDir),
                Throws.TypeOf<FileNotFoundException>());
        }

        [Test]
        public static void InternalFragmentAnalysisRunner_Run_WritesTsvAndExtractsIons()
        {
            Assert.That(File.Exists(PsmTsvPath), Is.True, $"Missing test data: {PsmTsvPath}");

            var sourceMzml = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", "sliced_ethcd.mzML");
            Assert.That(File.Exists(sourceMzml), Is.True, $"Missing test data: {sourceMzml}");

            var rawDir = Path.Combine(TestContext.CurrentContext.TestDirectory, $"ifrun_raw_{Guid.NewGuid():N}");
            var outDir = Path.Combine(TestContext.CurrentContext.TestDirectory, $"ifrun_out_{Guid.NewGuid():N}");
            Directory.CreateDirectory(rawDir);
            try
            {
                // The runner only loads files whose name matches a PSM's source file; copy the small
                // mzML under that name. Its scans won't match the PSM scan numbers, which the extractor
                // tolerates (falls back to matched-ion TIC), so ions still come from the PSM data.
                var psms = SpectrumMatchTsvReader.ReadPsmTsv(PsmTsvPath, out _);
                var sourceName = psms.First().FileNameWithoutExtension;
                File.Copy(sourceMzml, Path.Combine(rawDir, sourceName + ".mzML"));

                var ions = InternalFragmentAnalysisRunner.Run(PsmTsvPath, rawDir, outDir, defaultCollisionEnergy: 42.0);

                Assert.That(ions, Is.Not.Empty);
                Assert.That(File.Exists(Path.Combine(outDir, "InternalFragmentIons.tsv")), Is.True);
            }
            finally
            {
                if (Directory.Exists(rawDir)) Directory.Delete(rawDir, true);
                if (Directory.Exists(outDir)) Directory.Delete(outDir, true);
            }
        }

        [Test]
        public static void InternalFragmentAnalysisRunner_Run_LoadsMgfSpectrumFile()
        {
            // Exercises the .mgf arm of LoadRawFiles' format switch.
            RunWithRenamedSpectrumFile("tester.mgf", expectTsv: true);
        }

        [Test]
        public static void InternalFragmentAnalysisRunner_Run_MissingRawFolder_Throws()
        {
            var outDir = Path.Combine(TestContext.CurrentContext.TestDirectory, $"ifrun_out_{Guid.NewGuid():N}");
            var missingRaw = Path.Combine(TestContext.CurrentContext.TestDirectory, $"ifrun_noraw_{Guid.NewGuid():N}");
            Assert.That(() => InternalFragmentAnalysisRunner.Run(PsmTsvPath, missingRaw, outDir, 42.0),
                Throws.TypeOf<DirectoryNotFoundException>());
        }

        [Test]
        public static void InternalFragmentAnalysisRunner_Run_UnreadableFile_WarnsAndContinues()
        {
            // An unparseable spectrum file makes the loader throw; LoadRawFiles catches it, warns,
            // and continues, so Run still completes and writes the (empty) output TSV.
            Assert.That(File.Exists(PsmTsvPath), Is.True, $"Missing test data: {PsmTsvPath}");

            var rawDir = Path.Combine(TestContext.CurrentContext.TestDirectory, $"ifrun_raw_{Guid.NewGuid():N}");
            var outDir = Path.Combine(TestContext.CurrentContext.TestDirectory, $"ifrun_out_{Guid.NewGuid():N}");
            Directory.CreateDirectory(rawDir);
            try
            {
                var psms = SpectrumMatchTsvReader.ReadPsmTsv(PsmTsvPath, out _);
                var sourceName = psms.First().FileNameWithoutExtension;
                // Garbage content under a supported extension -> the mzML loader throws -> caught.
                File.WriteAllText(Path.Combine(rawDir, sourceName + ".mzML"), "this is not valid mzML xml");

                Assert.That(() => InternalFragmentAnalysisRunner.Run(PsmTsvPath, rawDir, outDir, 42.0),
                    Throws.Nothing);
                Assert.That(File.Exists(Path.Combine(outDir, "InternalFragmentIons.tsv")), Is.True);
            }
            finally
            {
                if (Directory.Exists(rawDir)) Directory.Delete(rawDir, true);
                if (Directory.Exists(outDir)) Directory.Delete(outDir, true);
            }
        }

        // Copies <DataFiles/spectrumFileName> into a temp folder under the PSM's source-file name so
        // the runner picks it up, runs the full pipeline, and asserts the output TSV was written.
        private static void RunWithRenamedSpectrumFile(string spectrumFileName, bool expectTsv)
        {
            Assert.That(File.Exists(PsmTsvPath), Is.True, $"Missing test data: {PsmTsvPath}");
            var source = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", spectrumFileName);
            Assert.That(File.Exists(source), Is.True, $"Missing test data: {source}");

            var rawDir = Path.Combine(TestContext.CurrentContext.TestDirectory, $"ifrun_raw_{Guid.NewGuid():N}");
            var outDir = Path.Combine(TestContext.CurrentContext.TestDirectory, $"ifrun_out_{Guid.NewGuid():N}");
            Directory.CreateDirectory(rawDir);
            try
            {
                var psms = SpectrumMatchTsvReader.ReadPsmTsv(PsmTsvPath, out _);
                var sourceName = psms.First().FileNameWithoutExtension;
                File.Copy(source, Path.Combine(rawDir, sourceName + Path.GetExtension(spectrumFileName)));

                Assert.That(() => InternalFragmentAnalysisRunner.Run(PsmTsvPath, rawDir, outDir, 42.0),
                    Throws.Nothing);
                if (expectTsv)
                    Assert.That(File.Exists(Path.Combine(outDir, "InternalFragmentIons.tsv")), Is.True);
            }
            finally
            {
                if (Directory.Exists(rawDir)) Directory.Delete(rawDir, true);
                if (Directory.Exists(outDir)) Directory.Delete(outDir, true);
            }
        }
    }
}
