using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using PredictionClients.Koina.AbstractClasses;
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
using System.Net.Http;
using System.Net.Sockets;
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
    /// 3. Component failure handling       — offline, stub components and the local ONNX model
    ///
    /// Nothing here calls Koina. The live end-to-end tests are in CombinedLibraryModelLiveTests,
    /// tagged ExternalService + Koina so the required CI job never runs them.
    ///
    /// [Category("RequiresOnnxModel")]  — needs the ONNX file on disk (it ships with PredictionClients)
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class CombinedLibraryModelTests
    {
        private static readonly string OnnxModelPath =
            Environment.GetEnvironmentVariable("INTERNAL_FRAGMENT_ONNX_PATH")
            ?? InternalFragmentIntensityModel.DefaultOnnxModelPath;

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
            var rts = new List<double?> { null };

            var primary = new Prosit2020IntensityHCD();
            var internal_ = new InternalFragmentIntensityModel(
                peptides, charges, rts, out _, onnxModelPath: OnnxModelPath);

            Assert.DoesNotThrow(() =>
                CombinedLibraryModel.WithPrimaryAndInternalFragments(primary, internal_, collisionEnergy: 35));
        }

        // ════════════════════════════════════════════════════════════════════════
        // 3. Component failure handling — offline (stubs + local ONNX only)
        //    The live Koina tests are in CombinedLibraryModelLiveTests.
        // ════════════════════════════════════════════════════════════════════════

        [Test, Category("RequiresOnnxModel")]
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

        [Test]
        public static void PrimaryIntensityComponent_ModelThrows_ReturnsFailedResultHoldingTheException()
        {
            var fault = new HttpRequestException("Koina unreachable (simulated)");
            var inputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 35, null, null),
            };
            var component = new PrimaryIntensityComponent(
                new ThrowingProsit2020IntensityHCD(fault), inputs, new double?[] { null });

            MixedModelResult result = null!;
            Assert.DoesNotThrowAsync(async () => result = await component.RunAsync(),
                "A component failure is recorded in the result, not thrown");

            Assert.That(result.Succeeded, Is.False);
            Assert.That(result.Error, Is.SameAs(fault));
            Assert.That(result.ContributionType, Is.EqualTo(ContributionType.PrimaryFragmentIntensities));
            Assert.That(result.Spectra, Is.Empty);
        }

        [Test, Category("RequiresOnnxModel")]
        public static async Task InternalIntensityComponent_ModelThrows_ReturnsFailedResultHoldingTheException()
        {
            var model = new InternalFragmentIntensityModel(
                new List<string> { "PEPTIDEK" }, new List<int> { 2 }, new List<double?> { null },
                out _, onnxModelPath: OnnxModelPath);
            model.Dispose(); // RunInferenceAsync on a disposed model throws ObjectDisposedException

            var result = await new InternalIntensityComponent(model).RunAsync();

            Assert.That(result.Succeeded, Is.False);
            Assert.That(result.Error, Is.TypeOf<ObjectDisposedException>());
            Assert.That(result.ContributionType, Is.EqualTo(ContributionType.InternalFragmentIntensities));
        }

        [Test]
        public static async Task RunAsync_ComponentResults_KeepTheCapturedException()
        {
            var fault = new InvalidOperationException("Simulated component failure");
            var combined = new CombinedLibraryModel(new List<IMixedModelComponent>
            {
                new BrokenComponentStub(ContributionType.PrimaryFragmentIntensities, "Broken", fault),
                new SucceedingComponentStub(ContributionType.InternalFragmentIntensities, "Internal", InternalOnlySpectrum()),
            });

            Assert.That(combined.ComponentResults, Is.Empty, "Empty until RunAsync is called");

            await combined.RunAsync();

            Assert.That(combined.ComponentResults.Select(r => r.ComponentName),
                Is.EqualTo(new[] { "Broken", "Internal" }), "One result per component, in component order");
            Assert.That(combined.ComponentResults[0].Succeeded, Is.False);
            Assert.That(combined.ComponentResults[0].Error, Is.SameAs(fault));
            Assert.That(combined.ComponentResults[1].Succeeded, Is.True);
        }

        // ── ThrowIfAnyComponentFailed: what the live tests do with a captured error ──
        // A transport fault must SKIP under ExternalServiceTestHelper.RunAsync; anything else must FAIL.

        private static IEnumerable<TestCaseData> TransportFaults()
        {
            yield return new TestCaseData(new HttpRequestException("Request failed with status 503 Service Unavailable"))
                .SetArgDisplayNames("HttpRequestException");
            yield return new TestCaseData(new TaskCanceledException("The request was canceled due to the configured HttpClient Timeout"))
                .SetArgDisplayNames("TaskCanceledException");
            yield return new TestCaseData(new SocketException((int)SocketError.ConnectionRefused))
                .SetArgDisplayNames("SocketException");
        }

        private static async Task<CombinedLibraryModel> RunWithPrimaryFailing(Exception fault)
        {
            var combined = new CombinedLibraryModel(new List<IMixedModelComponent>
            {
                new BrokenComponentStub(ContributionType.PrimaryFragmentIntensities, "Prosit", fault),
                new SucceedingComponentStub(ContributionType.InternalFragmentIntensities, "Internal", InternalOnlySpectrum()),
            });
            await combined.RunAsync();
            return combined;
        }

        private static Task RunUnderExternalServiceHelper(CombinedLibraryModel combined)
            => ExternalServiceTestHelper.RunAsync("Koina", () =>
            {
                CombinedLibraryModelLiveTests.ThrowIfAnyComponentFailed(combined);
                return Task.CompletedTask;
            });

        [Test]
        public static async Task ThrowIfAnyComponentFailed_NoFailure_DoesNotThrow()
        {
            var combined = new CombinedLibraryModel(new List<IMixedModelComponent>
            {
                new SucceedingComponentStub(ContributionType.InternalFragmentIntensities, "Internal", InternalOnlySpectrum()),
            });
            await combined.RunAsync();

            Assert.DoesNotThrow(() => CombinedLibraryModelLiveTests.ThrowIfAnyComponentFailed(combined));
        }

        [TestCaseSource(nameof(TransportFaults))]
        public static async Task ThrowIfAnyComponentFailed_TransportFault_IsSkippedByExternalServiceHelper(Exception fault)
        {
            var combined = await RunWithPrimaryFailing(fault);

            Assert.ThrowsAsync<IgnoreException>(() => RunUnderExternalServiceHelper(combined));
        }

        [Test]
        public static async Task ThrowIfAnyComponentFailed_ContractBreak_FailsUnderExternalServiceHelper()
        {
            // What FragmentIntensityModel.ResponseToPredictions throws when Koina answers with something it cannot parse
            var fault = new Exception("Something went wrong during deserialization of responses.");
            var combined = await RunWithPrimaryFailing(fault);

            var thrown = Assert.ThrowsAsync<Exception>(() => RunUnderExternalServiceHelper(combined));
            Assert.That(thrown, Is.SameAs(fault));
        }

        [Test]
        public static async Task ThrowIfAnyComponentFailed_TwoFailures_FailsEvenWhenOneIsATransportFault()
        {
            var combined = new CombinedLibraryModel(new List<IMixedModelComponent>
            {
                new BrokenComponentStub(ContributionType.PrimaryFragmentIntensities, "Prosit",
                    new HttpRequestException("Koina unreachable (simulated)")),
                new BrokenComponentStub(ContributionType.InternalFragmentIntensities, "Internal",
                    new InvalidOperationException("ONNX session broke")),
            });
            await combined.RunAsync();

            var thrown = Assert.ThrowsAsync<AggregateException>(() => RunUnderExternalServiceHelper(combined));
            Assert.That(thrown!.InnerExceptions, Has.Count.EqualTo(2));
        }

        [Test]
        public static async Task ThrowIfAnyComponentFailed_FailureWithNoException_Fails()
        {
            var combined = new CombinedLibraryModel(new List<IMixedModelComponent>
            {
                new SucceedingComponentStub(ContributionType.PrimaryFragmentIntensities, "Prosit",
                    new MixedModelResult { ComponentName = "Prosit", Succeeded = false }),
            });
            await combined.RunAsync();

            var thrown = Assert.ThrowsAsync<InvalidOperationException>(() => RunUnderExternalServiceHelper(combined));
            Assert.That(thrown!.Message, Does.Contain("Prosit"));
        }

        private static MixedModelResult InternalOnlySpectrum() =>
            MixedModelResult.FromSpectra("Internal", ContributionType.InternalFragmentIntensities, new[]
            {
                MakeSpectrum("PEPTIDEK", 2, 10.0, ProductType.b, 3, 300.0, 1.0,
                    secondaryType: ProductType.b, secondaryFragNum: 5),
            });
    }

    // ── Test helper ─────────────────────────────────────────────────────────────

    /// <summary>
    /// A stub component that always fails — used to test graceful degradation.
    /// </summary>
    internal class BrokenComponentStub : IMixedModelComponent
    {
        public string ComponentName { get; }
        public ContributionType ContributionType { get; }

        private readonly Exception _error;

        public BrokenComponentStub(ContributionType type, string name, Exception? error = null)
        {
            ContributionType = type;
            ComponentName = name;
            _error = error ?? new Exception("Simulated component failure");
        }

        public Task<MixedModelResult> RunAsync()
            => Task.FromResult(MixedModelResult.FromError(
                ComponentName, ContributionType, _error));
    }

    /// <summary>
    /// A stub component that returns a fixed result, so failure handling can be tested without a model.
    /// </summary>
    internal class SucceedingComponentStub : IMixedModelComponent
    {
        private readonly MixedModelResult _result;
        public string ComponentName { get; }
        public ContributionType ContributionType { get; }

        public SucceedingComponentStub(ContributionType type, string name, MixedModelResult result)
        {
            ContributionType = type;
            ComponentName = name;
            _result = result;
        }

        public Task<MixedModelResult> RunAsync() => Task.FromResult(_result);
    }

    /// <summary>
    /// Prosit 2020 HCD with the Koina call replaced by a fault, so PrimaryIntensityComponent's
    /// catch can be exercised offline.
    /// </summary>
    internal class ThrowingProsit2020IntensityHCD : Prosit2020IntensityHCD
    {
        private readonly Exception _fault;

        public ThrowingProsit2020IntensityHCD(Exception fault) => _fault = fault;

        protected override Task<List<PeptideFragmentIntensityPrediction>> AsyncThrottledPredictor(
            List<FragmentIntensityPredictionInput> modelInputs)
            => Task.FromException<List<PeptideFragmentIntensityPrediction>>(_fault);
    }
}
