using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;

namespace Test.KoinaTests
{
    /// <summary>
    /// Drives each concrete model's request-building code (ToBatchedRequests and the
    /// model-specific TryCleanSequence/Validate overrides) without network access, so the Koina
    /// request-construction paths count toward coverage. ToBatchedRequests does no validation —
    /// it only reads the Validated* fields — so one fully-populated input per family is enough.
    /// </summary>
    [TestFixture]
    public class KoinaRequestBuildingTests
    {
        private static readonly Assembly KoinaAssembly = typeof(FragmentIntensityModel).Assembly;

        private static IEnumerable<Type> Concrete<TBase>() =>
            KoinaAssembly.GetTypes().Where(t => !t.IsAbstract && typeof(TBase).IsAssignableFrom(t));

        public static IEnumerable<Type> FragmentModels() => Concrete<FragmentIntensityModel>();
        public static IEnumerable<Type> RtModels() => Concrete<RetentionTimeModel>();
        public static IEnumerable<Type> CcsModels() => Concrete<CollisionalCrossSectionModel>();
        public static IEnumerable<Type> CrosslinkModels() => Concrete<CrosslinkFragmentIntensityModel>();
        public static IEnumerable<Type> DetectabilityModels() => Concrete<DetectabilityModel>();

        private static object Instantiate(Type t)
        {
            var ctor = t.GetConstructors().First(c => c.GetParameters().All(p => p.HasDefaultValue));
            return ctor.Invoke(ctor.GetParameters().Select(p => p.DefaultValue).ToArray());
        }

        private static void AssertBuildsBatches(object model, object inputList)
        {
            var method = model.GetType().GetMethod("ToBatchedRequests", BindingFlags.NonPublic | BindingFlags.Instance)!;
            var batches = (IEnumerable)method.Invoke(model, new[] { inputList })!;

            int count = 0;
            foreach (var batch in batches)
            {
                var dict = (IDictionary<string, object>)batch;
                Assert.That(dict.ContainsKey("id"), Is.True);
                Assert.That(dict["inputs"], Is.InstanceOf<IEnumerable>());
                count++;
            }
            Assert.That(count, Is.GreaterThanOrEqualTo(1), $"{model.GetType().Name}.ToBatchedRequests produced no batches");
        }

        [TestCaseSource(nameof(FragmentModels))]
        public void FragmentIntensity_ToBatchedRequests_BuildsBatches(Type modelType)
        {
            var inputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 30, "QE", "HCD") { ValidatedFullSequence = "PEPTIDEK" }
            };
            AssertBuildsBatches(Instantiate(modelType), inputs);
        }

        [TestCaseSource(nameof(RtModels))]
        public void RetentionTime_ToBatchedRequests_BuildsBatches(Type modelType)
        {
            var inputs = new List<RetentionTimePredictionInput>
            {
                new("PEPTIDEK") { ValidatedFullSequence = "PEPTIDEK" }
            };
            AssertBuildsBatches(Instantiate(modelType), inputs);
        }

        [TestCaseSource(nameof(CcsModels))]
        public void Ccs_ToBatchedRequests_BuildsBatches(Type modelType)
        {
            var inputs = new List<CCSPredictionInput>
            {
                new("PEPTIDEK", 2) { ValidatedFullSequence = "PEPTIDEK" }
            };
            AssertBuildsBatches(Instantiate(modelType), inputs);
        }

        [TestCaseSource(nameof(CrosslinkModels))]
        public void Crosslink_ToBatchedRequests_BuildsBatches(Type modelType)
        {
            var inputs = new List<CrosslinkIntensityPredictionInput>
            {
                new("PEPTIDEK[UNIMOD:1896]", "ACDEK[UNIMOD:1896]", 2, 30)
                {
                    ValidatedAlphaSequence = "PEPTIDEK[UNIMOD:1896]",
                    ValidatedBetaSequence = "ACDEK[UNIMOD:1896]"
                }
            };
            AssertBuildsBatches(Instantiate(modelType), inputs);
        }

        [TestCaseSource(nameof(DetectabilityModels))]
        public void Detectability_ToBatchedRequests_BuildsBatches(Type modelType)
        {
            var inputs = new List<DetectabilityPredictionInput>
            {
                new("PEPTIDEK") { ValidatedFullSequence = "PEPTIDEK" }
            };
            AssertBuildsBatches(Instantiate(modelType), inputs);
        }

        // ── Model-specific request-building overrides ───────────────────────────────

        [Test]
        public void Tmt_TryCleanSequence_RejectsSequenceWithoutNTerminalLabel()
        {
            var model = new TmtProbe();
            // Return value null = rejected; that's what the prediction pipeline keys off.
            var result = model.Clean("PEPTIDEK", out _, out var warning);

            Assert.That(result, Is.Null);
            Assert.That(warning, Is.Not.Null);
            Assert.That(warning!.Message, Does.Contain("N-terminal"));
        }

        [Test]
        public void UniSpec_ValidateModelSpecificInputs_RejectsChargeOutsideInstrumentRange()
        {
            // VELOS supports charges 2-4; charge 5 must be rejected client-side (no network).
            var model = new UniSpec();
            var predictions = model.Predict(new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 5, 30, "VELOS", null)
            });

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].FragmentIntensities, Is.Null);
            Assert.That(predictions[0].Warning, Is.Not.Null);
        }

        [Test]
        public void Lac_ValidateModelSpecificInputs_RejectsUnknownInstrument()
        {
            // "orbitrap" is not in {ECLIPSE, ASTRAL, LUMOS}; the override upper-cases then rejects.
            var model = new Prosit2025IntensityLac();
            var predictions = model.Predict(new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 30, "orbitrap", "HCD")
            });

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].FragmentIntensities, Is.Null);
            Assert.That(predictions[0].Warning, Is.Not.Null);
        }

        private sealed class TmtProbe : Prosit2020IntensityTMT
        {
            public string? Clean(string sequence, out string? api, out WarningException? warning)
                => TryCleanSequence(sequence, out api, out warning);
        }
    }
}
