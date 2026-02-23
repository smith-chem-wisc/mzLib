using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.LocalModels;
using PredictionClients.MixedModels;
using PredictionClients.MixedModels.Components;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using CategoryAttribute = NUnit.Framework.CategoryAttribute;

namespace Test.MixedModelTests
{
    /// <summary>
    /// Tests for the MixedModels infrastructure.
    ///
    /// TEST ORGANISATION
    /// -----------------
    /// 1. LibrarySpectrumMerger unit tests  — pure logic, no models, no files
    /// 2. CombinedLibraryModel construction — validates the factory and component wiring
    /// 3. Integration tests                 — require Koina + ONNX model, tagged accordingly
    ///
    /// [Category("RequiresKoina")]      — needs a live Koina endpoint
    /// [Category("RequiresOnnxModel")]  — needs the ONNX file on disk
    /// [Category("Integration")]        — needs both
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class CombinedLibraryModelTests
    {
        private static readonly string OnnxModelPath =
            Environment.GetEnvironmentVariable("INTERNAL_FRAGMENT_ONNX_PATH")
            ?? Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"LocalModels\TestData\InternalFragmentScorer_v3_AllProteases.onnx");

        // ════════════════════════════════════════════════════════════════════════
        // 1. LibrarySpectrumMerger — pure unit tests (no models, no files)
        // ════════════════════════════════════════════════════════════════════════

        /// <summary>
        /// Helper to build a minimal LibrarySpectrum with one fragment ion.
        /// Avoids needing a real model for merger unit tests.
        /// </summary>
        private static LibrarySpectrum MakeSpectrum(
            string sequence, int charge, double rt,
            ProductType ionType, int fragNum, double mz, double intensity,
            ProductType? secondaryType = null, int secondaryFragNum = 0)
        {
            var product = new Product(
                ionType,
                secondaryType == null ? FragmentationTerminus.N : FragmentationTerminus.None,
                mz,
                fragNum,
                fragNum,
                neutralLoss: 0,
                secondaryType,
                secondaryFragNum);

            var ion = new MatchedFragmentIon(product, mz, intensity, charge: 1);
            return new LibrarySpectrum(sequence, mz * charge, charge,
                new List<MatchedFragmentIon> { ion }, rt);
        }

        [Test]
        public static void Merger_TwoComplementaryComponents_UnionsFragmentIons()
        {
            // Primary: one b-ion for PEPTIDEK/2
            var primarySpectrum = MakeSpectrum("PEPTIDEK", 2, 45.0, ProductType.b, 3, 300.15, 0.9);
            // Internal: one internal ion for PEPTIDEK/2
            var internalSpectrum = MakeSpectrum("PEPTIDEK", 2, 0.0,
                ProductType.b, 2, 270.13, 0.5, secondaryType: ProductType.b, secondaryFragNum: 4);

            var results = new List<MixedModelResult>
            {
                MixedModelResult.FromSpectra("Prosit", ContributionType.PrimaryFragmentIntensities,
                    new[] { primarySpectrum }),
                MixedModelResult.FromSpectra("Internal", ContributionType.InternalFragmentIntensities,
                    new[] { internalSpectrum }),
            };

            var merged = LibrarySpectrumMerger.Merge(results, out var warnings);

            Assert.That(merged.ContainsKey("PEPTIDEK/2"), Is.True);
            Assert.That(merged["PEPTIDEK/2"].MatchedFragmentIons.Count, Is.EqualTo(2),
                "Merged spectrum should contain both the primary b-ion and the internal ion");
            Assert.That(warnings, Is.Null);
        }

        [Test]
        public static void Merger_PrimaryRtPreferred_WhenNoRtComponent()
        {
            var primarySpectrum = MakeSpectrum("PEPTIDEK", 2, 45.0, ProductType.b, 3, 300.15, 0.9);
            var internalSpectrum = MakeSpectrum("PEPTIDEK", 2, 0.0,
                ProductType.b, 2, 270.13, 0.5, ProductType.b, 4);

            var results = new List<MixedModelResult>
            {
                MixedModelResult.FromSpectra("Prosit",   ContributionType.PrimaryFragmentIntensities, new[] { primarySpectrum }),
                MixedModelResult.FromSpectra("Internal", ContributionType.InternalFragmentIntensities, new[] { internalSpectrum }),
            };

            var merged = LibrarySpectrumMerger.Merge(results, out _);

            Assert.That(merged["PEPTIDEK/2"].RetentionTime, Is.EqualTo(45.0).Within(1e-9),
                "RT should come from the primary spectrum when no dedicated RT component is present");
        }

        [Test]
        public static void Merger_DedicatedRtComponent_OverridesPrimaryRt()
        {
            var primarySpectrum = MakeSpectrum("PEPTIDEK", 2, 45.0, ProductType.b, 3, 300.15, 0.9);

            var rtData = new Dictionary<string, double> { ["PEPTIDEK/2"] = 99.7 };

            var results = new List<MixedModelResult>
            {
                MixedModelResult.FromSpectra("Prosit", ContributionType.PrimaryFragmentIntensities,
                    new[] { primarySpectrum }),
                MixedModelResult.FromScalars("AlphaPeptRT", ContributionType.RetentionTime,
                    rtData),
            };

            var merged = LibrarySpectrumMerger.Merge(results, out _);

            Assert.That(merged["PEPTIDEK/2"].RetentionTime, Is.EqualTo(99.7).Within(1e-9),
                "Dedicated RT component should override the primary spectrum's RT");
        }

        [Test]
        public static void Merger_InternalOnlyComponent_ProducesValidSpectrum()
        {
            // No primary component — internal-only should still produce a merged spectrum
            var internalSpectrum = MakeSpectrum("PEPTIDEK", 2, 0.0,
                ProductType.b, 2, 270.13, 0.5, ProductType.b, 4);

            var results = new List<MixedModelResult>
            {
                MixedModelResult.FromSpectra("Internal", ContributionType.InternalFragmentIntensities,
                    new[] { internalSpectrum }),
            };

            var merged = LibrarySpectrumMerger.Merge(results, out var warnings);

            Assert.That(merged.ContainsKey("PEPTIDEK/2"), Is.True);
            Assert.That(merged["PEPTIDEK/2"].MatchedFragmentIons.Count, Is.EqualTo(1));
            Assert.That(warnings, Is.Null);
        }

        [Test]
        public static void Merger_FailedComponent_SkippedWithWarning()
        {
            var primarySpectrum = MakeSpectrum("PEPTIDEK", 2, 45.0, ProductType.b, 3, 300.15, 0.9);

            var results = new List<MixedModelResult>
            {
                MixedModelResult.FromSpectra("Prosit", ContributionType.PrimaryFragmentIntensities,
                    new[] { primarySpectrum }),
                MixedModelResult.FromError("BrokenModel", ContributionType.InternalFragmentIntensities,
                    new Exception("ONNX file not found")),
            };

            var merged = LibrarySpectrumMerger.Merge(results, out var warnings);

            // Merge still succeeds with just the primary ions
            Assert.That(merged.ContainsKey("PEPTIDEK/2"), Is.True);
            Assert.That(merged["PEPTIDEK/2"].MatchedFragmentIons.Count, Is.EqualTo(1),
                "Only primary ions should be present when internal component failed");
            Assert.That(warnings, Is.Not.Null,
                "Failed component should generate a warning");
            Assert.That(warnings!.Message, Does.Contain("BrokenModel"),
                "Warning should identify the failed component by name");
        }

        [Test]
        public static void Merger_DuplicateFragmentIons_HigherIntensityKept()
        {
            // Both components claim an ion at m/z 300.15 with the same annotation — simulates a collision
            var ion1 = new MatchedFragmentIon(
                new Product(ProductType.b, FragmentationTerminus.N, 300.15, 3, 3, 0),
                300.15, 0.4, 1);
            var ion2 = new MatchedFragmentIon(
                new Product(ProductType.b, FragmentationTerminus.N, 300.15, 3, 3, 0),
                300.15, 0.9, 1);

            var s1 = new LibrarySpectrum("PEPTIDEK", 500.0, 2, new List<MatchedFragmentIon> { ion1 }, 45.0);
            var s2 = new LibrarySpectrum("PEPTIDEK", 500.0, 2, new List<MatchedFragmentIon> { ion2 }, 0.0);

            var results = new List<MixedModelResult>
            {
                MixedModelResult.FromSpectra("ModelA", ContributionType.PrimaryFragmentIntensities,  new[] { s1 }),
                MixedModelResult.FromSpectra("ModelB", ContributionType.InternalFragmentIntensities, new[] { s2 }),
            };

            var merged = LibrarySpectrumMerger.Merge(results, out var warnings);

            Assert.That(merged["PEPTIDEK/2"].MatchedFragmentIons.Count, Is.EqualTo(1),
                "Duplicate ions should be collapsed to one");
            Assert.That(merged["PEPTIDEK/2"].MatchedFragmentIons[0].Intensity,
                Is.EqualTo(0.9).Within(1e-9),
                "Higher intensity ion should be kept");
            Assert.That(warnings, Is.Not.Null,
                "Duplicate resolution should be noted in warnings");
        }

        [Test]
        public static void Merger_MultiplePeptides_AllPresent()
        {
            var s1 = MakeSpectrum("PEPTIDEK", 2, 45.0, ProductType.b, 3, 300.15, 0.9);
            var s2 = MakeSpectrum("ELVISLIVESK", 2, 60.0, ProductType.y, 5, 600.32, 0.8);
            var i1 = MakeSpectrum("PEPTIDEK", 2, 0.0, ProductType.b, 2, 270.13, 0.5, ProductType.b, 4);

            var results = new List<MixedModelResult>
            {
                MixedModelResult.FromSpectra("Prosit",   ContributionType.PrimaryFragmentIntensities,
                    new[] { s1, s2 }),
                MixedModelResult.FromSpectra("Internal", ContributionType.InternalFragmentIntensities,
                    new[] { i1 }),  // only PEPTIDEK has internal ion predictions here
            };

            var merged = LibrarySpectrumMerger.Merge(results, out _);

            Assert.That(merged.Count, Is.EqualTo(2));
            Assert.That(merged["PEPTIDEK/2"].MatchedFragmentIons.Count, Is.EqualTo(2)); // b-ion + internal
            Assert.That(merged["ELVISLIVESK/2"].MatchedFragmentIons.Count, Is.EqualTo(1)); // y-ion only
        }

        [Test]
        public static void Merger_ComponentWarning_SurfacedInOutput()
        {
            var spectrum = MakeSpectrum("PEPTIDEK", 2, 45.0, ProductType.b, 3, 300.15, 0.9);
            var componentWarning = new WarningException("2 duplicate spectra removed");

            var results = new List<MixedModelResult>
            {
                MixedModelResult.FromSpectra("Prosit", ContributionType.PrimaryFragmentIntensities,
                    new[] { spectrum }, componentWarning),
            };

            LibrarySpectrumMerger.Merge(results, out var warnings);

            Assert.That(warnings, Is.Not.Null);
            Assert.That(warnings!.Message, Does.Contain("2 duplicate spectra removed"),
                "Per-component warnings should be forwarded to the caller");
        }

        // ════════════════════════════════════════════════════════════════════════
        // 2. CombinedLibraryModel construction
        // ════════════════════════════════════════════════════════════════════════

        [Test]
        public static void Constructor_NoComponents_ThrowsArgumentException()
        {
            Assert.Throws<ArgumentException>(() =>
                new CombinedLibraryModel(new List<IMixedModelComponent>()));
        }

        [Test]
        public static void Constructor_NullComponents_ThrowsArgumentException()
        {
            Assert.Throws<ArgumentException>(() =>
                new CombinedLibraryModel(null!));
        }

        [Test]
        public static void WithPrimaryAndInternalFragments_Factory_BuildsCorrectly()
        {
            // Just verify the factory constructs without throwing (no inference run)
            var peptides = new List<string> { "PEPTIDEK" };
            var charges = new List<int> { 2 };
            var energies = new List<int> { 35 };
            var rts = new List<double?> { null };

            var primary = new Prosit2020IntensityHCD(peptides, charges, energies, rts, out _);
            var internal_ = new InternalFragmentIntensityModel(
                peptides, charges, rts, out _, onnxModelPath: OnnxModelPath);

            Assert.DoesNotThrow(() =>
                CombinedLibraryModel.WithPrimaryAndInternalFragments(primary, internal_));
        }

        // ════════════════════════════════════════════════════════════════════════
        // 3. Integration tests (require both Koina and ONNX model)
        // ════════════════════════════════════════════════════════════════════════

        [Test, NUnit.Framework.Category("Integration"), NUnit.Framework.Category("RequiresKoina"), NUnit.Framework.Category("RequiresOnnxModel")]
        public static async Task RunAsync_CombinedSpectra_ContainBothPrimaryAndInternalIons()
        {
            var peptides = new List<string> { "PEPTIDEK", "ELVISLIVESK" };
            var charges = new List<int> { 2, 2 };
            var energies = new List<int> { 35, 35 };
            var rts = new List<double?> { 100.0, 200.0 };

            var primary = new Prosit2020IntensityHCD(peptides, charges, energies, rts, out _);
            var internal_ = new InternalFragmentIntensityModel(
                peptides, charges, rts, out _, onnxModelPath: OnnxModelPath);

            var combined = CombinedLibraryModel.WithPrimaryAndInternalFragments(primary, internal_);
            var warning = await combined.RunAsync();

            Assert.That(combined.PredictedSpectra.Count, Is.EqualTo(2));

            foreach (var spectrum in combined.PredictedSpectra)
            {
                var primaryIons = spectrum.MatchedFragmentIons.Where(f => !f.IsInternalFragment).ToList();
                var internalIons = spectrum.MatchedFragmentIons.Where(f => f.IsInternalFragment).ToList();

                Assert.That(primaryIons.Count, Is.GreaterThan(0),
                    $"{spectrum.Name}: should have primary (b/y) ions from Prosit");
                Assert.That(internalIons.Count, Is.GreaterThan(0),
                    $"{spectrum.Name}: should have internal fragment ions from local model");
            }
        }

        [Test, Category("Integration"), Category("RequiresKoina"), Category("RequiresOnnxModel")]
        public static async Task RunAsync_AllIonMzValues_AreChemicallyReasonable()
        {
            var peptides = new List<string> { "PEPTIDEK" };
            var charges = new List<int> { 2 };
            var energies = new List<int> { 35 };
            var rts = new List<double?> { null };

            var combined = CombinedLibraryModel.WithPrimaryAndInternalFragments(
                new Prosit2020IntensityHCD(peptides, charges, energies, rts, out _),
                new InternalFragmentIntensityModel(peptides, charges, rts, out _,
                    onnxModelPath: OnnxModelPath));

            await combined.RunAsync();

            foreach (var ion in combined.PredictedSpectra.SelectMany(s => s.MatchedFragmentIons))
            {
                Assert.That(ion.Mz, Is.GreaterThan(0).And.LessThan(5000));
                Assert.That(ion.Intensity, Is.GreaterThanOrEqualTo(0));
                Assert.That(ion.Charge, Is.GreaterThan(0));
            }
        }

        [Test, Category("Integration"), Category("RequiresKoina"), Category("RequiresOnnxModel")]
        public static async Task RunAsync_RetentionTimeTakenFromPrimary_WhenNoRtComponent()
        {
            var peptides = new List<string> { "PEPTIDEK" };
            var charges = new List<int> { 2 };
            var energies = new List<int> { 35 };
            var rts = new List<double?> { 42.5 };

            var combined = CombinedLibraryModel.WithPrimaryAndInternalFragments(
                new Prosit2020IntensityHCD(peptides, charges, energies, rts, out _),
                new InternalFragmentIntensityModel(peptides, charges, rts, out _,
                    onnxModelPath: OnnxModelPath));

            await combined.RunAsync();

            Assume.That(combined.PredictedSpectra.Count, Is.EqualTo(1));
            Assert.That(combined.PredictedSpectra[0].RetentionTime, Is.EqualTo(42.5).Within(1e-6),
                "RT from primary model input should flow through to merged spectrum");
        }

        [Test, Category("Integration"), Category("RequiresKoina"), Category("RequiresOnnxModel")]
        public static async Task RunAsync_MspRoundTrip_ParsedSpectraHaveBothIonTypes()
        {
            var outPath = Path.Combine(
                TestContext.CurrentContext.TestDirectory,
                "combinedLibraryRoundTripTest.msp");

            SpectralLibrary? savedLib = null;
            try
            {
                var peptides = new List<string> { "PEPTIDEK", "SAMPLER" };
                var charges = new List<int> { 2, 2 };
                var energies = new List<int> { 35, 35 };
                var rts = new List<double?> { 100.0, 200.0 };

                var combined = CombinedLibraryModel.WithPrimaryAndInternalFragments(
                    new Prosit2020IntensityHCD(peptides, charges, energies, rts, out _),
                    new InternalFragmentIntensityModel(peptides, charges, rts, out _,
                        onnxModelPath: OnnxModelPath),
                    spectralLibrarySavePath: outPath);

                await combined.RunAsync();

                Assert.That(File.Exists(outPath), Is.True);

                savedLib = new SpectralLibrary(new List<string> { outPath });
                var savedSpectra = savedLib.GetAllLibrarySpectra().ToList();

                Assert.That(savedSpectra.Count, Is.EqualTo(combined.PredictedSpectra.Count));

                foreach (var spectrum in savedSpectra)
                {
                    // After round-trip through MSP, IsInternalFragment should still be correct
                    var internal_ = spectrum.MatchedFragmentIons.Where(f => f.IsInternalFragment).ToList();
                    var primary_ = spectrum.MatchedFragmentIons.Where(f => !f.IsInternalFragment).ToList();

                    Assert.That(primary_.Count, Is.GreaterThan(0),
                        $"{spectrum.Name}: primary ions should survive MSP round-trip");
                    Assert.That(internal_.Count, Is.GreaterThan(0),
                        $"{spectrum.Name}: internal ions should survive MSP round-trip");
                }
            }
            finally
            {
                savedLib?.CloseConnections();
                if (File.Exists(outPath)) File.Delete(outPath);
            }
        }

        [Test, Category("Integration"), Category("RequiresKoina"), Category("RequiresOnnxModel")]
        public static async Task RunAsync_KoinaFails_InternalOnlyLibraryStillProduced()
        {
            // Simulate Koina being unavailable by using a broken primary component
            var brokenComponent = new BrokenComponentStub(
                ContributionType.PrimaryFragmentIntensities, "Prosit (unreachable)");

            var peptides = new List<string> { "PEPTIDEK" };
            var charges = new List<int> { 2 };
            var rts = new List<double?> { null };

            var internal_ = new InternalFragmentIntensityModel(
                peptides, charges, rts, out _, onnxModelPath: OnnxModelPath);

            var combined = new CombinedLibraryModel(new List<IMixedModelComponent>
            {
                brokenComponent,
                new InternalIntensityComponent(internal_),
            });

            var warning = await combined.RunAsync();

            Assert.That(warning, Is.Not.Null,
                "Broken component should produce a warning");
            Assert.That(warning!.Message, Does.Contain("Prosit (unreachable)"),
                "Warning should name the failed component");

            // Internal-only spectra should still be present
            Assert.That(combined.PredictedSpectra.Count, Is.GreaterThan(0),
                "Internal-only library should be produced even when Koina fails");

            var allInternal = combined.PredictedSpectra
                .All(s => s.MatchedFragmentIons.All(f => f.IsInternalFragment));
            Assert.That(allInternal, Is.True,
                "All ions should be internal when primary component failed");
        }
    }

    // ── Test helper ─────────────────────────────────────────────────────────────

    /// <summary>
    /// A stub component that always fails — used to test graceful degradation.
    /// </summary>
    internal class BrokenComponentStub : IMixedModelComponent
    {
        public string ComponentName { get; }
        public ContributionType ContributionType { get; }

        public BrokenComponentStub(ContributionType type, string name)
        {
            ContributionType = type;
            ComponentName = name;
        }

        public Task<MixedModelResult> RunAsync()
            => Task.FromResult(MixedModelResult.FromError(
                ComponentName, ContributionType,
                new Exception("Simulated component failure")));
    }
}
