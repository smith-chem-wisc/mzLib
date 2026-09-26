using NUnit.Framework;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.LocalModels;
using PredictionClients.MixedModels;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Runtime.ExceptionServices;
using System.Threading.Tasks;

namespace Test.MixedModelTests
{
    /// <summary>
    /// CombinedLibraryModel end to end: Prosit 2020 HCD on the live Koina server merged with the
    /// local internal-fragment ONNX model.
    ///
    /// Tagged ExternalService + Koina so the required (coverage-collecting) CI job never runs it.
    ///
    /// The Koina call happens inside PrimaryIntensityComponent, which catches every exception into
    /// MixedModelResult.FromError. Left alone, an outage would surface as an empty or internal-only
    /// library and a failed assertion, and a contract break could pass with nothing to assert on.
    /// So every test rethrows the captured exception (<see cref="ThrowIfAnyComponentFailed"/>)
    /// inside ExternalServiceTestHelper.RunAsync, which skips a transport fault and fails anything else.
    ///
    /// TODO once master is merged into this branch: derive from KoinaLiveTestFixture and drop the
    /// [OneTimeSetUp] below. It has the same reachability probe plus the narrower fault classifier.
    /// </summary>
    [TestFixture]
    [Category("ExternalService")]
    [Category("Koina")]
    [ExcludeFromCodeCoverage]
    public class CombinedLibraryModelLiveTests
    {
        private const string KoinaReadyUrl = "https://koina.wilhelmlab.org:443/v2/health/ready";

        private static readonly string OnnxModelPath =
            Environment.GetEnvironmentVariable("INTERNAL_FRAGMENT_ONNX_PATH")
            ?? InternalFragmentIntensityModel.DefaultOnnxModelPath;

        [OneTimeSetUp]
        public void EnsureKoinaReachable()
        {
            ExternalServiceTestHelper.EnsureReachable("Koina", KoinaReadyUrl);
        }

        /// <summary>
        /// Rethrows the exception a component captured during RunAsync, with its original stack trace,
        /// so ExternalServiceTestHelper.RunAsync can classify it. A single failure is rethrown as is.
        /// More than one is thrown as an AggregateException, which always fails: only the primary
        /// component touches the network, so two failures at once are not an outage.
        /// </summary>
        internal static void ThrowIfAnyComponentFailed(CombinedLibraryModel combined)
        {
            var errors = combined.ComponentResults
                .Where(r => !r.Succeeded)
                .Select(r => r.Error ?? new InvalidOperationException(
                    $"Component '{r.ComponentName}' failed without recording an exception"))
                .ToList();

            if (errors.Count == 1)
                ExceptionDispatchInfo.Capture(errors[0]).Throw();
            if (errors.Count > 1)
                throw new AggregateException("More than one mixed-model component failed", errors);
        }

        [Test]
        public static async Task RunAsync_CombinedSpectra_ContainBothPrimaryAndInternalIons()
        {
            Assume.That(File.Exists(OnnxModelPath), Is.True,
                $"ONNX model not found at: {OnnxModelPath}");

            await ExternalServiceTestHelper.RunAsync("Koina", async () =>
            {
                var peptides = new List<string> { "PEPTIDEK", "ELVISLIVESK" };
                var charges = new List<int> { 2, 2 };
                var rts = new List<double?> { 100.0, 200.0 };

                var primary = new Prosit2020IntensityHCD();
                var internal_ = new InternalFragmentIntensityModel(
                    peptides, charges, rts, out _, onnxModelPath: OnnxModelPath);

                var combined = CombinedLibraryModel.WithPrimaryAndInternalFragments(primary, internal_, collisionEnergy: 35);
                await combined.RunAsync();
                ThrowIfAnyComponentFailed(combined);

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
            });
        }

        [Test]
        public static async Task RunAsync_AllIonMzValues_AreChemicallyReasonable()
        {
            await ExternalServiceTestHelper.RunAsync("Koina", async () =>
            {
                var peptides = new List<string> { "PEPTIDEK" };
                var charges = new List<int> { 2 };
                var rts = new List<double?> { null };

                var combined = CombinedLibraryModel.WithPrimaryAndInternalFragments(
                    new Prosit2020IntensityHCD(),
                    new InternalFragmentIntensityModel(peptides, charges, rts, out _,
                        onnxModelPath: OnnxModelPath),
                    collisionEnergy: 35);

                await combined.RunAsync();
                ThrowIfAnyComponentFailed(combined);

                var ions = combined.PredictedSpectra.SelectMany(s => s.MatchedFragmentIons).ToList();
                Assert.That(ions, Is.Not.Empty, "Nothing to check: the combined library has no fragment ions");

                foreach (var ion in ions)
                {
                    Assert.That(ion.Mz, Is.GreaterThan(0).And.LessThan(5000));
                    Assert.That(ion.Intensity, Is.GreaterThanOrEqualTo(0));
                    Assert.That(ion.Charge, Is.GreaterThan(0));
                }
            });
        }

        [Test]
        public static async Task RunAsync_RetentionTimeTakenFromPrimary_WhenNoRtComponent()
        {
            await ExternalServiceTestHelper.RunAsync("Koina", async () =>
            {
                var peptides = new List<string> { "PEPTIDEK" };
                var charges = new List<int> { 2 };
                var rts = new List<double?> { 42.5 };

                var combined = CombinedLibraryModel.WithPrimaryAndInternalFragments(
                    new Prosit2020IntensityHCD(),
                    new InternalFragmentIntensityModel(peptides, charges, rts, out _,
                        onnxModelPath: OnnxModelPath),
                    collisionEnergy: 35);

                await combined.RunAsync();
                ThrowIfAnyComponentFailed(combined);

                Assert.That(combined.PredictedSpectra.Count, Is.EqualTo(1));
                Assert.That(combined.PredictedSpectra[0].RetentionTime, Is.EqualTo(42.5).Within(1e-6),
                    "RT from primary model input should flow through to merged spectrum");
            });
        }

        [Test]
        public static async Task RunAsync_MspRoundTrip_ParsedSpectraHaveBothIonTypes()
        {
            var outPath = Path.Combine(
                TestContext.CurrentContext.TestDirectory,
                "combinedLibraryRoundTripTest.msp");

            SpectralLibrary? savedLib = null;
            try
            {
                await ExternalServiceTestHelper.RunAsync("Koina", async () =>
                {
                    var peptides = new List<string> { "PEPTIDEK", "SAMPLER" };
                    var charges = new List<int> { 2, 2 };
                    var rts = new List<double?> { 100.0, 200.0 };

                    var combined = CombinedLibraryModel.WithPrimaryAndInternalFragments(
                        new Prosit2020IntensityHCD(),
                        new InternalFragmentIntensityModel(peptides, charges, rts, out _,
                            onnxModelPath: OnnxModelPath),
                        collisionEnergy: 35,
                        spectralLibrarySavePath: outPath);

                    await combined.RunAsync();
                    ThrowIfAnyComponentFailed(combined);

                    Assert.That(File.Exists(outPath), Is.True);

                    savedLib = new SpectralLibrary(new List<string> { outPath });
                    var savedSpectra = savedLib.GetAllLibrarySpectra().ToList();

                    Assert.That(savedSpectra, Is.Not.Empty);
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
                });
            }
            finally
            {
                savedLib?.CloseConnections();
                if (File.Exists(outPath)) File.Delete(outPath);
            }
        }
    }
}
