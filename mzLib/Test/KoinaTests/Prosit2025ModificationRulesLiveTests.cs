using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using Readers.ProForma;

namespace Test.KoinaTests
{
    /// <summary>
    /// Live canary for Prosit2025ModificationRules: one batch per Prosit 2025 22PTM/40PTM model, with good peptides,
    /// N-terminal pyro-Glu from E and from Q (sent on the residue), and E[UNIMOD:954], which Koina answers with
    /// HTTP 400 for the whole batch. The batch must succeed, with E[UNIMOD:954] refused client-side.
    /// </summary>
    [TestFixture]
    [Category("ExternalService")]
    [Category("Koina")]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class Prosit2025ModificationRulesLiveTests : KoinaLiveTestFixture
    {
        private static readonly string[] Batch =
        {
            "PEPTIDEK",
            "PEPM[UNIMOD:35]TIDEK",
            "PEPC[UNIMOD:21]TIDEK",  // zero-shot
            "PEPTIDEK-[UNIMOD:34]",  // C-terminal, zero-shot
            "[UNIMOD:27]-EPEPTIDEK", // sent as E[UNIMOD:27]PEPTIDEK
            "[UNIMOD:28]-QPEPTIDEK", // sent as Q[UNIMOD:28]PEPTIDEK
            "PEPE[UNIMOD:954]TIDEK", // refused before sending
        };

        private static readonly bool[] ExpectedMask = { true, true, true, true, true, true, false };

        private static IEnumerable<TestCaseData> RetentionTimeModels()
        {
            yield return new TestCaseData(new Prosit2025iRT40PTM()).SetArgDisplayNames("Prosit_2025_irt_40PTM");
            yield return new TestCaseData(new Prosit2025iRT22PTM()).SetArgDisplayNames("Prosit_2025_irt_22PTM");
        }

        private static IEnumerable<TestCaseData> IntensityModels()
        {
            yield return new TestCaseData(new Prosit2025Intensity40PTM()).SetArgDisplayNames("Prosit_2025_intensity_40PTM");
            yield return new TestCaseData(new Prosit2025Intensity22PTM()).SetArgDisplayNames("Prosit_2025_intensity_22PTM");
        }

        [TestCaseSource(nameof(RetentionTimeModels))]
        public void RetentionTime_BatchWithARefusedPairSucceeds(RetentionTimeModel model)
        {
            var predictions = model.Predict(Batch.Select(s => new RetentionTimePredictionInput(s) { SequenceParser = ProFormaSequenceParser.Instance }).ToList());

            Assert.That(model.ValidInputsMask, Is.EqualTo(ExpectedMask));
            Assert.That(predictions.Take(6).All(p => p.PredictedRetentionTime != null), Is.True);
            Assert.That(predictions[4].ValidatedFullSequence, Is.EqualTo("E[UNIMOD:27]PEPTIDEK"));
            Assert.That(predictions[5].ValidatedFullSequence, Is.EqualTo("Q[UNIMOD:28]PEPTIDEK"));
            Assert.That(predictions[6].PredictedRetentionTime, Is.Null);
            Assert.That(predictions[6].Warning!.Message, Does.Contain("E[UNIMOD:954]"));
        }

        [TestCaseSource(nameof(IntensityModels))]
        public void Intensity_BatchWithARefusedPairSucceeds(FragmentIntensityModel model)
        {
            var predictions = model.Predict(Batch.Select(s => new FragmentIntensityPredictionInput(s, 2, 30, null, "HCD") { SequenceParser = ProFormaSequenceParser.Instance }).ToList());

            Assert.That(model.ValidInputsMask, Is.EqualTo(ExpectedMask));
            Assert.That(predictions.Take(6).All(p => p.FragmentIntensities is { Count: > 0 }), Is.True);
            Assert.That(predictions[4].ValidatedFullSequence, Is.EqualTo("E[UNIMOD:27]PEPTIDEK"));
            Assert.That(predictions[5].ValidatedFullSequence, Is.EqualTo("Q[UNIMOD:28]PEPTIDEK"));
            Assert.That(predictions[6].FragmentIntensities, Is.Null);
            Assert.That(predictions[6].Warning!.Message, Does.Contain("E[UNIMOD:954]"));
        }
    }
}
