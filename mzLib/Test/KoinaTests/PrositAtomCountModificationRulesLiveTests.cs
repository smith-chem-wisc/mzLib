using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Client;
using PredictionClients.Koina.SupportedModels;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using Readers.ProForma;

namespace Test.KoinaTests
{
    /// <summary>
    /// Live canary for PrositAtomCountModificationRules: one batch per Prosit 2025 22PTM/40PTM and 2024 PTMs_gl model, with
    /// good peptides, N-terminal pyro-Glu from E and from Q (sent on the residue), ids only zero-shot encodes, and
    /// E[UNIMOD:954], which Koina answers with HTTP 400 for the whole batch. The batch must succeed, with E[UNIMOD:954] and a
    /// C-terminal modification refused client-side.
    /// </summary>
    [TestFixture]
    [Category("ExternalService")]
    [Category("Koina")]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class PrositAtomCountModificationRulesLiveTests : KoinaLiveTestFixture
    {
        private static readonly string[] Batch =
        {
            "PEPTIDEK",
            "PEPM[UNIMOD:35]TIDEK",
            "PEPC[UNIMOD:21]TIDEK",  // zero-shot
            "PEPTIDEK-[UNIMOD:34]",  // C-terminal, refused before sending
            "[UNIMOD:27]-EPEPTIDEK", // sent as E[UNIMOD:27]PEPTIDEK
            "[UNIMOD:28]-QPEPTIDEK", // sent as Q[UNIMOD:28]PEPTIDEK
            "PEPK[UNIMOD:64]TIDEK",  // zero-shot
            "PEPK[UNIMOD:259]TIDEK", // zero-shot
            "[UNIMOD:5]-PEPTIDEK",   // N-terminal, zero-shot
            "PEPE[UNIMOD:954]TIDEK", // refused before sending
        };

        private static readonly bool[] ExpectedMask = { true, true, true, false, true, true, true, true, true, false };

        private static IEnumerable<TestCaseData> RetentionTimeModels()
        {
            yield return new TestCaseData(new Prosit2025iRT40PTM()).SetArgDisplayNames("Prosit_2025_irt_40PTM");
            yield return new TestCaseData(new Prosit2025iRT22PTM()).SetArgDisplayNames("Prosit_2025_irt_22PTM");
            yield return new TestCaseData(new Prosit2024iRTPTMsGl()).SetArgDisplayNames("Prosit_2024_irt_PTMs_gl");
        }

        private static IEnumerable<TestCaseData> IntensityModels()
        {
            yield return new TestCaseData(new Prosit2025Intensity40PTM()).SetArgDisplayNames("Prosit_2025_intensity_40PTM");
            yield return new TestCaseData(new Prosit2025Intensity22PTM()).SetArgDisplayNames("Prosit_2025_intensity_22PTM");
            yield return new TestCaseData(new Prosit2024IntensityPTMsGl()).SetArgDisplayNames("Prosit_2024_intensity_PTMs_gl");
        }

        [TestCaseSource(nameof(RetentionTimeModels))]
        public void RetentionTime_BatchWithARefusedPairSucceeds(RetentionTimeModel model)
        {
            var predictions = model.Predict(Batch.Select(s => new RetentionTimePredictionInput(s) { SequenceParser = ProFormaSequenceParser.Instance }).ToList());

            Assert.That(model.ValidInputsMask, Is.EqualTo(ExpectedMask));
            Assert.That(predictions.Where((_, i) => ExpectedMask[i]).All(p => p.PredictedRetentionTime != null), Is.True);
            Assert.That(predictions[4].ValidatedFullSequence, Is.EqualTo("E[UNIMOD:27]PEPTIDEK"));
            Assert.That(predictions[5].ValidatedFullSequence, Is.EqualTo("Q[UNIMOD:28]PEPTIDEK"));
            Assert.That(predictions[3].PredictedRetentionTime, Is.Null);
            Assert.That(predictions[3].Warning!.Message, Does.Contain("C-terminus"));
            Assert.That(predictions[9].PredictedRetentionTime, Is.Null);
            Assert.That(predictions[9].Warning!.Message, Does.Contain("E[UNIMOD:954]"));
        }

        [TestCaseSource(nameof(IntensityModels))]
        public void Intensity_BatchWithARefusedPairSucceeds(FragmentIntensityModel model)
        {
            var predictions = model.Predict(Batch.Select(s => new FragmentIntensityPredictionInput(s, 2, 30, null, "HCD") { SequenceParser = ProFormaSequenceParser.Instance }).ToList());

            Assert.That(model.ValidInputsMask, Is.EqualTo(ExpectedMask));
            Assert.That(predictions.Where((_, i) => ExpectedMask[i]).All(p => p.FragmentIntensities is { Count: > 0 }), Is.True);
            Assert.That(predictions[4].ValidatedFullSequence, Is.EqualTo("E[UNIMOD:27]PEPTIDEK"));
            Assert.That(predictions[5].ValidatedFullSequence, Is.EqualTo("Q[UNIMOD:28]PEPTIDEK"));
            Assert.That(predictions[3].FragmentIntensities, Is.Null);
            Assert.That(predictions[3].Warning!.Message, Does.Contain("C-terminus"));
            Assert.That(predictions[9].FragmentIntensities, Is.Null);
            Assert.That(predictions[9].Warning!.Message, Does.Contain("E[UNIMOD:954]"));
        }
    }

    /// <summary>
    /// Sends one modification the rule refuses, for each kind of refusal, straight to each model that uses the rule. Koina
    /// must still answer HTTP 400; if it starts predicting one, this fails and the rule is out of date for that model.
    /// </summary>
    [TestFixture]
    [Category("ExternalService")]
    [Category("Koina")]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class PrositAtomCountRefusalLiveTests : KoinaLiveTestFixture
    {
        private static readonly (string Token, string Peptide, string Reason)[] Refusals =
        {
            ("K[UNIMOD:2073]", "PEPK[UNIMOD:2073]TIDEK", "is not in Koina's UNIMOD table"),
            ("A[UNIMOD:35]", "PEPA[UNIMOD:35]TIDEK", "has no attachment group"),
            ("K[UNIMOD:954]", "PEPK[UNIMOD:954]TIDEK", "Koina encodes only H, C, N, O, P and S"),
            ("S[UNIMOD:295]", "PEPS[UNIMOD:295]TIDEK", "cannot read the composition"),
            ("[UNIMOD:27]-", "[UNIMOD:27]-EPEPTIDEK", "negative"),
        };

        private static IEnumerable<TestCaseData> Cases()
        {
            var models = new (string Name, IRawSender Model)[]
            {
                ("Prosit_2025_irt_22PTM", new Irt22()),
                ("Prosit_2025_irt_40PTM", new Irt40()),
                ("Prosit_2024_irt_PTMs_gl", new IrtPtmsGl()),
                ("Prosit_2025_intensity_22PTM", new Intensity22()),
                ("Prosit_2025_intensity_40PTM", new Intensity40()),
                ("Prosit_2024_intensity_PTMs_gl", new IntensityPtmsGl()),
            };
            foreach (var (name, model) in models)
                foreach (var (token, peptide, reason) in Refusals)
                    yield return new TestCaseData(model, token, peptide, reason).SetArgDisplayNames(name, token);
        }

        [TestCaseSource(nameof(Cases))]
        public void RefusedModification_IsRefusedByKoina(IRawSender model, string token, string peptide, string reason)
        {
            Assert.That(PrositAtomCountModificationRules.RefusalReason(token), Does.Contain(reason));

            Assert.ThrowsAsync<KoinaInputRejectedException>(() => model.Send(peptide), $"Koina predicted {peptide}");
        }

        public interface IRawSender
        {
            /// <summary>Sends one peptide in Koina's notation to the model, without any check.</summary>
            Task<string> Send(string koinaSequence);
        }

        private static InputField[] IrtFields(string koinaSequence) =>
            new[] { new InputField("peptide_sequences", "BYTES", new[] { koinaSequence }) };

        private static InputField[] IntensityFields(string koinaSequence) => new[]
        {
            new InputField("peptide_sequences", "BYTES", new[] { koinaSequence }),
            new InputField("precursor_charges", "INT32", new[] { 2 }),
            new InputField("collision_energies", "FP32", new[] { 30f }),
            new InputField("fragmentation_types", "BYTES", new[] { "HCD" }),
        };

        private sealed class Irt22 : Prosit2025iRT22PTM, IRawSender
        {
            public Task<string> Send(string koinaSequence) => SendInferenceRequestAsync(ModelName, BuildBatchedRequest(0, IrtFields(koinaSequence)), CancellationToken.None);
        }

        private sealed class Irt40 : Prosit2025iRT40PTM, IRawSender
        {
            public Task<string> Send(string koinaSequence) => SendInferenceRequestAsync(ModelName, BuildBatchedRequest(0, IrtFields(koinaSequence)), CancellationToken.None);
        }

        private sealed class IrtPtmsGl : Prosit2024iRTPTMsGl, IRawSender
        {
            public Task<string> Send(string koinaSequence) => SendInferenceRequestAsync(ModelName, BuildBatchedRequest(0, IrtFields(koinaSequence)), CancellationToken.None);
        }

        private sealed class Intensity22 : Prosit2025Intensity22PTM, IRawSender
        {
            public Task<string> Send(string koinaSequence) => SendInferenceRequestAsync(ModelName, BuildBatchedRequest(0, IntensityFields(koinaSequence)), CancellationToken.None);
        }

        private sealed class Intensity40 : Prosit2025Intensity40PTM, IRawSender
        {
            public Task<string> Send(string koinaSequence) => SendInferenceRequestAsync(ModelName, BuildBatchedRequest(0, IntensityFields(koinaSequence)), CancellationToken.None);
        }

        private sealed class IntensityPtmsGl : Prosit2024IntensityPTMsGl, IRawSender
        {
            public Task<string> Send(string koinaSequence) => SendInferenceRequestAsync(ModelName, BuildBatchedRequest(0, IntensityFields(koinaSequence)), CancellationToken.None);
        }
    }
}
