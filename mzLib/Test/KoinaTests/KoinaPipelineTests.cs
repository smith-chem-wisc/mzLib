using System.Linq;
using Chemistry;
using MassSpectrometry;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.CCSModels;
using PredictionClients.Koina.SupportedModels.CrosslinkIntensityModels;
using PredictionClients.Koina.SupportedModels.FlyabilityModels;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;

namespace Test.KoinaTests
{
    /// <summary>
    /// Drives the full batched prediction pipeline (validation -> ToBatchedRequests -> transport
    /// -> ResponseToPredictions -> realignment) for each model family without network, by
    /// overriding the SendInferenceRequestAsync seam with a canned response.
    /// </summary>
    [TestFixture]
    public class KoinaPipelineTests
    {
        // One peptide, two fragments — shared by the intensity-style families.
        private const string IntensityJson =
            "{\"outputs\":[" +
            "{\"name\":\"annotation\",\"datatype\":\"BYTES\",\"shape\":[2],\"data\":[\"b1+1\",\"y1+1\"]}," +
            "{\"name\":\"mz\",\"datatype\":\"FP32\",\"shape\":[2],\"data\":[100.0,200.0]}," +
            "{\"name\":\"intensities\",\"datatype\":\"FP32\",\"shape\":[2],\"data\":[0.5,0.6]}]}";

        [Test]
        public void FragmentIntensity_Predict_RunsPipelineWithFakeTransport()
        {
            var predictions = new FakeFragmentModel().Predict(new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 30, null, null)
            });

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].FragmentAnnotations, Is.EquivalentTo(new[] { "b1+1", "y1+1" }));
            Assert.That(predictions[0].FragmentIntensities, Is.EquivalentTo(new[] { 0.5, 0.6 }));
        }

        /// <summary>
        /// A modified peptide goes through prediction and into a library spectrum in the default mapping mode. The label
        /// keeps its mzLib notation, and the precursor and the C-containing y3 carry the carbamidomethyl mass.
        /// </summary>
        [Test]
        public void FragmentIntensity_ModifiedPeptide_BuildsALibrarySpectrum()
        {
            const string sequence = "PEPTIDEC[Common Fixed:Carbamidomethyl on C]K";
            var model = new FakeModifiedFragmentModel();
            model.Predict(new List<FragmentIntensityPredictionInput> { new(sequence, 2, 30, null, null) });

            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 12.0 }, out _);

            Assert.That(spectra.Count, Is.EqualTo(1));
            var peptide = new PeptideWithSetModifications(sequence);
            Assert.That(spectra[0].Sequence, Is.EqualTo(sequence));
            Assert.That(spectra[0].PrecursorMz, Is.EqualTo(peptide.MonoisotopicMass.ToMz(2)).Within(1e-4));
            var y3 = spectra[0].MatchedFragmentIons.Single(ion => ion.NeutralTheoreticalProduct.ProductType == ProductType.y
                && ion.NeutralTheoreticalProduct.FragmentNumber == 3);
            var products = new List<Product>();
            new PeptideWithSetModifications("PEPTIDECK").Fragment(DissociationType.HCD, FragmentationTerminus.C, products);
            double unmodifiedY3 = products.Single(p => p.ProductType == ProductType.y && p.FragmentNumber == 3).NeutralMass.ToMz(1);
            Assert.That(y3.Mz - unmodifiedY3, Is.EqualTo(57.02146).Within(1e-3));
        }

        [Test]
        public void RetentionTime_Predict_RunsPipelineWithFakeTransport()
        {
            var predictions = new FakeRtModel().Predict(new List<RetentionTimePredictionInput> { new("PEPTIDEK") });

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].PredictedRetentionTime, Is.EqualTo(42.5));
        }

        [Test]
        public void Ccs_Predict_RunsPipelineWithFakeTransport()
        {
            var predictions = new FakeCcsModel().Predict(new List<CCSPredictionInput> { new("PEPTIDEK", 2) });

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].PredictedCCS, Is.EqualTo(321.5));
        }

        [Test]
        public void Crosslink_Predict_RunsPipelineWithFakeTransport()
        {
            var predictions = new FakeCrosslinkModel().Predict(new List<CrosslinkIntensityPredictionInput>
            {
                new("PEPTIDEK[UNIMOD:1896]", "ACDEK[UNIMOD:1896]", 2, 30)
            });

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].FragmentIntensities, Is.EquivalentTo(new[] { 0.5, 0.6 }));
        }

        [Test]
        public void Detectability_Predict_RunsPipelineWithFakeTransport()
        {
            var predictions = new FakeDetectabilityModel().Predict(new List<DetectabilityPredictionInput> { new("PEPTIDEK") });

            Assert.That(predictions.Count, Is.EqualTo(1));
            Assert.That(predictions[0].DetectabilityProbabilities, Is.Not.Null);
            Assert.That(predictions[0].DetectabilityProbabilities!.Value.HighDetectability, Is.EqualTo(0.4));
        }

        // ── canned-transport subclasses of real models ─────────────────────────────

        private sealed class FakeFragmentModel : Prosit2019Intensity
        {
            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
                => Task.FromResult(IntensityJson);
        }

        private sealed class FakeModifiedFragmentModel : Prosit2020IntensityHCD
        {
            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
                => Task.FromResult(
                    "{\"outputs\":[" +
                    "{\"name\":\"annotation\",\"datatype\":\"BYTES\",\"shape\":[2],\"data\":[\"b2+1\",\"y3+1\"]}," +
                    "{\"name\":\"mz\",\"datatype\":\"FP32\",\"shape\":[2],\"data\":[227.1,450.2]}," +
                    "{\"name\":\"intensities\",\"datatype\":\"FP32\",\"shape\":[2],\"data\":[0.5,1.0]}]}");
        }

        private sealed class FakeRtModel : Prosit2019iRT
        {
            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
                => Task.FromResult("{\"outputs\":[{\"name\":\"irt\",\"datatype\":\"FP32\",\"shape\":[1],\"data\":[42.5]}]}");
        }

        private sealed class FakeCcsModel : IM2Deep
        {
            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
                => Task.FromResult("{\"outputs\":[{\"name\":\"ccs\",\"datatype\":\"FP32\",\"shape\":[1],\"data\":[321.5]}]}");
        }

        private sealed class FakeCrosslinkModel : Prosit2023IntensityXLCMS2
        {
            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
                => Task.FromResult(IntensityJson);
        }

        private sealed class FakeDetectabilityModel : PFly2024FineTuned
        {
            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
                => Task.FromResult("{\"outputs\":[{\"name\":\"detectability\",\"datatype\":\"FP32\",\"shape\":[4],\"data\":[0.1,0.2,0.3,0.4]}]}");
        }
    }
}
