using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Linq;
using System.Net;
using System.Net.Http;
using System.Threading;
using System.Threading.Tasks;
using Chromatography;
using Chromatography.RetentionTimePrediction;
using Newtonsoft.Json.Linq;
using NUnit.Framework;
using Omics.SequenceConversion;
using Readers.ProForma;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Client;
using PredictionClients.Koina.SupportedModels.CCSModels;
using PredictionClients.Koina.SupportedModels.CrosslinkIntensityModels;
using PredictionClients.Koina.SupportedModels.FlyabilityModels;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using PredictionClients.Koina.Util;
using WarningException = System.ComponentModel.WarningException;

namespace Test.KoinaTests
{
    /// <summary>
    /// Koina answers a whole request with 400 when the model raises on one of its inputs. These tests
    /// pin what the models now do about it (KoinaModelBase.RequestPredictionsAsync): split the refused
    /// request until the one input is found, report it, and keep every other prediction. And what they
    /// still must not do: split an outage, a malformed request, or a model refusing everything.
    ///
    /// No network. The transport is a fake Koina behind the SendInferenceRequestAsync seam, and the
    /// refusal it sends is the real body, kept verbatim.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class KoinaRefusedInputIsolationTests
    {
        private const string KoinasReason = "E_954: no attachment group for residue 'E'";

        // ── one refused input, every family ───────────────────────────────────────

        [TestCase("RetentionTime")]
        [TestCase("FragmentIntensity")]
        [TestCase("Ccs")]
        [TestCase("Detectability")]
        [TestCase("Crosslink")]
        public void OneRefusedInputAmongMany_IsIsolated_AndEveryOtherInputIsPredicted(string family)
        {
            // Ten inputs, four to a request: requests of 4, 4 and 2, and input 5 is in the second.
            var koina = FakeKoina.Refusing(5);

            var outcome = Predict(family, koina, count: 10, maxBatchSize: 4);

            AssertRefusedExactly(outcome, 5);
            // The three requests, then the refused one halved twice: [4,5] [6,7], then [4] [5].
            Assert.That(koina.Requests, Is.EqualTo(7));
        }

        [TestCase("RetentionTime")]
        [TestCase("FragmentIntensity")]
        public void TwoRefusedInputsInOneRequest_AreBothIsolated(string family)
        {
            var koina = FakeKoina.Refusing(1, 6);

            var outcome = Predict(family, koina, count: 8, maxBatchSize: 8);

            AssertRefusedExactly(outcome, 1, 6);
            // 1 request, then [0-3] [4-7], then [0,1] [2,3] [4,5] [6,7], then [0] [1] [6] [7].
            Assert.That(koina.Requests, Is.EqualTo(1 + 2 + 4 + 4));
        }

        [Test]
        public void TheOnlyInputRefused_IsReportedWithoutAnyFurtherRequest()
        {
            var koina = FakeKoina.Refusing(0);

            var outcome = Predict("RetentionTime", koina, count: 1);

            AssertRefusedExactly(outcome, 0);
            Assert.That(koina.Requests, Is.EqualTo(1));
        }

        /// <summary>
        /// An input refused on the client side never reaches Koina, so the inputs Koina sees are not
        /// numbered as the caller's are. Koina's refusal still has to land on the right input.
        /// </summary>
        [Test]
        public void ARefusalFromKoina_AfterAnInputRefusedOnTheClientSide_LandsOnTheRightInput()
        {
            var koina = FakeKoina.Refusing(3);
            var model = new RtModel(koina, maxBatchSize: 4);
            var inputs = Inputs(6);
            inputs[1] = new RetentionTimePredictionInput("PEPTIDEBK"); // B is not a residue the model reads

            var predictions = model.Predict(inputs);

            Assert.That(model.ValidInputsMask, Is.EqualTo(new[] { true, false, true, false, true, true }));
            Assert.That(predictions[1].Warning!.Message, Does.Contain("Invalid base sequence"));
            Assert.That(predictions[1].Warning!.InnerException, Is.Null, "refused here, not by Koina");
            Assert.That(predictions[3].Warning!.Message, Does.Contain(KoinasReason));
            Assert.That(predictions[3].PredictedRetentionTime, Is.Null);
            foreach (int i in new[] { 0, 2, 4, 5 })
                Assert.That(predictions[i].PredictedRetentionTime, Is.EqualTo(Value(i)));
            // Inputs 0 2 3 4 then 5; then [0,2] [3,4]; then [3] [4].
            Assert.That(koina.Requests, Is.EqualTo(2 + 2 + 2));
        }

        /// <summary>
        /// The happy path is the old path: one parse over every response, not one per request.
        /// </summary>
        [Test]
        public void NothingRefused_TheParserSeesEveryResponseAtOnce_AsBefore()
        {
            var koina = new FakeKoina(_ => null);
            var model = new RtModel(koina, maxBatchSize: 4);

            var predictions = model.Predict(Inputs(10));

            Assert.That(predictions.Select(p => p.PredictedRetentionTime), Is.EqualTo(Enumerable.Range(0, 10).Select(Value)));
            Assert.That(model.ParserCalls, Is.EqualTo(new[] { 3 }), "one call carrying all three responses");
            Assert.That(koina.Requests, Is.EqualTo(3));
        }

        // ── what is NOT split ─────────────────────────────────────────────────────

        /// <summary>
        /// An outage keeps today's behaviour: it throws, and nothing is sent again. Only the request
        /// carrying input 5 fails, which is the case splitting would otherwise chase.
        /// </summary>
        [TestCase("ServerError", typeof(HttpRequestException))]
        [TestCase("ServiceFault", typeof(KoinaServiceException))]
        [TestCase("MalformedRequest", typeof(HttpRequestException))]
        [TestCase("Timeout", typeof(OperationCanceledException))]
        public void AnOutageOrAMalformedRequest_IsNotSplit_AndThrowsAsBefore(string failure, Type expected)
        {
            Func<Exception> fail = failure switch
            {
                // The very body of a refusal, but from a 5xx: an outage, whatever it says.
                "ServerError" => () => KoinaServiceException.ForFailedResponse(503, "Service Unavailable", RealRefusal),
                "ServiceFault" => () => KoinaServiceException.ForFailedResponse(400, "Bad Request",
                    TestKoinaServiceFaultClassification.RealCudaFault),
                // Triton refuses this before any model runs; every request would carry the same fault.
                "MalformedRequest" => () => KoinaServiceException.ForFailedResponse(400, "Bad Request",
                    "{\"error\":\"unexpected inference input 'peptide_sequenes' for model 'Prosit_2019_irt'\"}"),
                _ => () => new TaskCanceledException("The request was canceled due to the configured HttpClient.Timeout elapsing."),
            };
            var koina = new FakeKoina(sequences => sequences.Any(s => IndexOf(s) == 5) ? fail() : null);
            var model = new RtModel(koina, maxBatchSize: 4);

            var thrown = Assert.Catch(() => model.Predict(Inputs(10)));

            Assert.That(thrown, Is.InstanceOf(expected));
            Assert.That(thrown, Is.Not.InstanceOf<KoinaInputRejectedException>());
            Assert.That(koina.Requests, Is.EqualTo(3), "the three requests the session always sends, and no more");
        }

        [Test]
        public void AnOutageAlongsideARefusal_StillThrows_AndNothingIsSplit()
        {
            // Input 1's request is refused; input 9's meets an outage in the same group of requests.
            var koina = new FakeKoina(sequences =>
                sequences.Any(s => IndexOf(s) == 9) ? KoinaServiceException.ForFailedResponse(503, "Service Unavailable", "")
                : sequences.Any(s => IndexOf(s) == 1) ? FakeKoina.Refusal()
                : null);
            var model = new RtModel(koina, maxBatchSize: 4);

            Assert.Throws<HttpRequestException>(() => model.Predict(Inputs(10)));
            Assert.That(koina.Requests, Is.EqualTo(3));
        }

        // ── handling modes ────────────────────────────────────────────────────────

        [Test]
        public void ThrowExceptionMode_ARefusedInputThrows_AsAClientSideRefusalDoes()
        {
            var koina = FakeKoina.Refusing(5);
            var model = new RtModel(koina, maxBatchSize: 4, mode: SequenceConversionHandlingMode.ThrowException);

            var thrown = Assert.Throws<ArgumentException>(() => model.Predict(Inputs(10)));

            Assert.That(thrown!.Message, Does.Contain(Sequence(5)).And.Contain(KoinasReason));
            Assert.That(thrown.InnerException, Is.InstanceOf<KoinaInputRejectedException>());
        }

        /// <summary>
        /// RemoveIncompatibleElements has nothing it could remove: Koina does not say which part of the
        /// input it refused. So it leaves the input unpredicted, as ReturnNull does.
        /// </summary>
        [TestCase(SequenceConversionHandlingMode.ReturnNull)]
        [TestCase(SequenceConversionHandlingMode.RemoveIncompatibleElements)]
        [TestCase(SequenceConversionHandlingMode.UsePrimarySequence)]
        public void EveryOtherMode_LeavesTheRefusedInputUnpredicted(SequenceConversionHandlingMode mode)
        {
            var koina = FakeKoina.Refusing(2);

            var outcome = Predict("RetentionTime", koina, count: 4, maxBatchSize: 4, mode: mode);

            AssertRefusedExactly(outcome, 2);
        }

        // ── the real case, and PredictRetentionTimeEquivalents ────────────────────

        /// <summary>
        /// The real case, offline: Prosit_2025_irt_40PTM refuses zinc on E, E[UNIMOD:954], sent as
        /// ProForma the way the MSL pipeline sends it. Prosit2025ModificationRules now refuses it before
        /// sending, so the one request carries the other three and there is nothing to split. The fake still
        /// refuses E[UNIMOD:954]: if the rule let it through, it would be sent and take five requests.
        /// </summary>
        [Test]
        public void Prosit2025iRT40PTM_AnE954_IsRefusedBeforeSending_AndFailsAlone()
        {
            var koina = new FakeKoina(sequences => sequences.Any(s => s.Contains("E[UNIMOD:954]")) ? FakeKoina.Refusal() : null);
            var model = new Rt40PtmModel(koina);
            var inputs = new[] { Sequence(0), "PEPE[UNIMOD:954]TIDEK", Sequence(2), Sequence(3) }
                .Select(s => new RetentionTimePredictionInput(s) { SequenceParser = ProFormaSequenceParser.Instance }).ToList();

            var predictions = model.Predict(inputs);

            Assert.That(predictions[1].ValidatedFullSequence, Is.Null, "refused before sending");
            Assert.That(predictions[1].PredictedRetentionTime, Is.Null);
            Assert.That(predictions[1].Warning!.Message, Does.Contain("Prosit_2025_irt_40PTM cannot place E[UNIMOD:954]"));
            Assert.That(predictions[1].Warning.InnerException, Is.Null, "mzLib refused it, not Koina");
            Assert.That(model.ValidInputsMask, Is.EqualTo(new[] { true, false, true, true }));
            foreach (int i in new[] { 0, 2, 3 })
                Assert.That(predictions[i].PredictedRetentionTime, Is.EqualTo(Value(i)));
            Assert.That(koina.Requests, Is.EqualTo(1));
        }

        /// <summary>
        /// Before, the catch in PredictRetentionTimeEquivalents turned one refusal into PredictionError
        /// for every peptide in the call.
        /// </summary>
        [Test]
        public void PredictRetentionTimeEquivalents_OneRefusedPeptide_FailsAlone()
        {
            var koina = FakeKoina.Refusing(1);
            var model = new RtModel(koina);
            var peptides = Enumerable.Range(0, 4).Select(i => new Peptide(Sequence(i))).ToList();

            var results = model.PredictRetentionTimeEquivalents(peptides);

            Assert.That(results.Count, Is.EqualTo(4));
            Assert.That(results[1].PredictedValue, Is.Null);
            Assert.That(results[1].FailureReason, Is.EqualTo(RetentionTimeFailureReason.IncompatibleModifications),
                "the model could not handle this peptide; nothing went wrong with the call");
            Assert.That(model.Predictions[1].Warning!.Message, Does.Contain(KoinasReason));
            foreach (int i in new[] { 0, 2, 3 })
            {
                Assert.That(results[i].PredictedValue, Is.EqualTo(Value(i)));
                Assert.That(results[i].FailureReason, Is.Null);
            }
        }

        [Test]
        public void PredictRetentionTimeEquivalents_AnOutage_StillFailsEveryPeptide()
        {
            var koina = new FakeKoina(_ => KoinaServiceException.ForFailedResponse(503, "Service Unavailable", ""));
            var model = new RtModel(koina);

            var results = model.PredictRetentionTimeEquivalents(new[] { Sequence(0), Sequence(1) }.Select(s => new Peptide(s)).ToList());

            Assert.That(results.All(r => r.PredictedValue is null && r.FailureReason == RetentionTimeFailureReason.PredictionError), Is.True);
            Assert.That(koina.Requests, Is.EqualTo(1));
        }

        // ── bounds ────────────────────────────────────────────────────────────────

        [Test]
        public void OneRefusedInputAmongAThousand_CostsTwoRequestsPerHalving()
        {
            var koina = FakeKoina.Refusing(617);

            var outcome = Predict("RetentionTime", koina, count: 1000, maxBatchSize: 1000);

            AssertRefusedExactly(outcome, 617);
            // 1000 -> 500 -> ... -> 1 is ceil(log2 1000) = 10 halvings, two requests each.
            Assert.That(koina.Requests, Is.EqualTo(1 + 2 * 10));
        }

        /// <summary>
        /// A model that refuses everything looks like a broken model, or a notation it no longer reads,
        /// not like bad inputs. Isolating all 1000 would take 1998 more requests for nothing.
        /// </summary>
        [Test]
        public void AModelThatRefusesEverything_StopsSplittingAndThrows_AfterABoundedNumberOfRequests()
        {
            var koina = new FakeKoina(_ => FakeKoina.Refusal());
            var model = new RtModel(koina, maxBatchSize: 1000);

            Assert.Throws<KoinaInputRejectedException>(() => model.Predict(Inputs(1000)));

            // At most 10 refused inputs are allowed here (1% of 1000, and never fewer than 10). Rounds of
            // 2, 4, 8 and 16 requests leave 16 still refused, so at least 16 bad inputs: it stops there.
            Assert.That(koina.Requests, Is.EqualTo(1 + 2 + 4 + 8 + 16));
        }

        [Test]
        public void SplitRequests_GoOutNoMoreThanMaxNumberOfBatchesPerRequestAtATime()
        {
            // Sixteen inputs, four requests, one refused input in each.
            var koina = FakeKoina.Refusing(1, 6, 9, 14);
            koina.Latency = TimeSpan.FromMilliseconds(20);

            var outcome = Predict("RetentionTime", koina, count: 16, maxBatchSize: 4, batchesPerRequest: 2);

            AssertRefusedExactly(outcome, 1, 6, 9, 14);
            Assert.That(koina.MaxInFlight, Is.EqualTo(2));
            // 4 requests, then 8 halves of two, then the 4 refused halves split into 8 singles.
            Assert.That(koina.Requests, Is.EqualTo(4 + 8 + 8));
        }

        // ── models that cannot be split ───────────────────────────────────────────

        [Test]
        public void AModelThatDoesNotBatchInOrder_IsNotSplit()
        {
            // All ten inputs in one request, whatever MaxBatchSize says: which inputs a request carried is unknown.
            var koina = FakeKoina.Refusing(5);
            var model = new RtModel(koina, maxBatchSize: 4) { Batching = BatchingRule.AllInOneRequest };

            Assert.Throws<KoinaInputRejectedException>(() => model.Predict(Inputs(10)));
            Assert.That(koina.Requests, Is.EqualTo(1));
        }

        [Test]
        public void AModelWhoseHalvesDoNotFitOneRequest_IsNotSplitFurther()
        {
            var koina = FakeKoina.Refusing(1);
            var model = new RtModel(koina, maxBatchSize: 4) { Batching = BatchingRule.OneRequestPerInputBelowMaxBatchSize };

            Assert.Throws<KoinaInputRejectedException>(() => model.Predict(Inputs(4)));
            Assert.That(koina.Requests, Is.EqualTo(1));
        }

        // ── helpers ───────────────────────────────────────────────────────────────

        /// <summary>
        /// What Prosit_2025_irt_40PTM answered on 2026-11-09 for a request holding PEPE[UNIMOD:954]TIDEK,
        /// verbatim, with HTTP 400.
        /// </summary>
        internal const string RealRefusal =
            """{"error":"in ensemble 'Prosit_2025_irt_40PTM', Failed to process the request(s) for model instance 'Prosit_Preprocess_ac_loss_0_1', message: KeyError: \"E_954: no attachment group for residue 'E'\"\n\nAt:\n  /models/repo/Prosit_Preprocess_ac_loss/1/zero_shot.py(77): zero_shot_counts\n  /models/repo/Prosit_Preprocess_ac_loss/1/model.py(119): get_ac\n  /models/repo/Prosit_Preprocess_ac_loss/1/model.py(127): <listcomp>\n  /models/repo/Prosit_Preprocess_ac_loss/1/model.py(127): get_ac_all\n  /models/repo/Prosit_Preprocess_ac_loss/1/model.py(148): execute\n"}""";

        private static void AssertRefusedExactly(Outcome outcome, params int[] refused)
        {
            Assert.Multiple(() =>
            {
                for (int i = 0; i < outcome.Values.Length; i++)
                {
                    if (refused.Contains(i))
                    {
                        Assert.That(outcome.Values[i], Is.Null, $"input {i} was refused");
                        Assert.That(outcome.Mask[i], Is.False, $"input {i} was refused");
                        Assert.That(outcome.Warnings[i]?.Message, Does.Contain(KoinasReason), $"input {i} carries Koina's reason");
                        Assert.That(outcome.Warnings[i]?.Message, Does.Not.Contain("zero_shot.py"), "without the Python traceback");
                        Assert.That(outcome.Warnings[i]?.InnerException, Is.InstanceOf<KoinaInputRejectedException>());
                    }
                    else
                    {
                        Assert.That(outcome.Values[i], Is.EqualTo(Value(i)), $"input {i} keeps its own prediction");
                        Assert.That(outcome.Mask[i], Is.True, $"input {i} was predicted");
                    }
                }
            });
        }

        private sealed record Outcome(double?[] Values, WarningException?[] Warnings, bool[] Mask);

        /// <summary>
        /// Runs one model family's Predict over inputs 0..count-1 against <paramref name="koina"/>.
        /// </summary>
        private static Outcome Predict(string family, FakeKoina koina, int count, int maxBatchSize = 1000,
            SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ReturnNull, int batchesPerRequest = 500)
        {
            var indices = Enumerable.Range(0, count).ToList();
            switch (family)
            {
                case "RetentionTime":
                {
                    var model = new RtModel(koina, maxBatchSize, mode, batchesPerRequest);
                    var p = model.Predict(Inputs(count));
                    return new(p.Select(x => x.PredictedRetentionTime).ToArray(), p.Select(x => x.Warning).ToArray(), model.ValidInputsMask);
                }
                case "FragmentIntensity":
                {
                    var model = new FragmentModel(koina, maxBatchSize, mode, batchesPerRequest);
                    var p = model.Predict(indices.Select(i => new FragmentIntensityPredictionInput(Sequence(i), 2, 30, null, null)).ToList());
                    return new(p.Select(x => x.FragmentIntensities?[0]).ToArray(), p.Select(x => x.Warning).ToArray(), model.ValidInputsMask);
                }
                case "Ccs":
                {
                    var model = new CcsModel(koina, maxBatchSize, mode, batchesPerRequest);
                    var p = model.Predict(indices.Select(i => new CCSPredictionInput(Sequence(i), 2)).ToList());
                    return new(p.Select(x => x.PredictedCCS).ToArray(), p.Select(x => x.Warning).ToArray(), model.ValidInputsMask);
                }
                case "Detectability":
                {
                    var model = new DetectabilityTestModel(koina, maxBatchSize, mode, batchesPerRequest);
                    var p = model.Predict(indices.Select(i => new DetectabilityPredictionInput(Sequence(i))).ToList());
                    return new(p.Select(x => x.DetectabilityProbabilities?.NotDetectable).ToArray(), p.Select(x => x.Warning).ToArray(), model.ValidInputsMask);
                }
                case "Crosslink":
                {
                    var model = new CrosslinkModel(koina, maxBatchSize, mode, batchesPerRequest);
                    var p = model.Predict(indices.Select(i => new CrosslinkIntensityPredictionInput(
                        Sequence(i).TrimEnd('K') + "K[UNIMOD:1896]", "ACDEK[UNIMOD:1896]", 2, 30)).ToList());
                    return new(p.Select(x => x.FragmentIntensities?[0]).ToArray(), p.Select(x => x.Warning).ToArray(), model.ValidInputsMask);
                }
                default:
                    throw new ArgumentException($"unknown family {family}");
            }
        }

        private static List<RetentionTimePredictionInput> Inputs(int count) =>
            Enumerable.Range(0, count).Select(i => new RetentionTimePredictionInput(Sequence(i))).ToList();

        // Each input's sequence spells its index in three residues every model here accepts unmodified,
        // and the fake Koina answers each sequence with its index plus a half. So a prediction that lands
        // on the wrong input after splitting shows up as the wrong number.
        private const string Alphabet = "ADEFGHILNQRSTVWY";

        internal static string Sequence(int index) =>
            "PEP" + Alphabet[index / 256 % 16] + Alphabet[index / 16 % 16] + Alphabet[index % 16] + "K";

        internal static int IndexOf(string sequence)
        {
            if (sequence.Length < 6 || !sequence.StartsWith("PEP"))
                return -1;
            int a = Alphabet.IndexOf(sequence[3]), b = Alphabet.IndexOf(sequence[4]), c = Alphabet.IndexOf(sequence[5]);
            return a < 0 || b < 0 || c < 0 ? -1 : a * 256 + b * 16 + c;
        }

        private static double? Value(int index) => index + 0.5;

        private static string Number(double value) => value.ToString("0.0###", CultureInfo.InvariantCulture);

        private static string OneValuePerInput(string outputName, IReadOnlyList<string> sequences, int valuesPerInput = 1)
        {
            var data = sequences.SelectMany(s => Enumerable.Range(0, valuesPerInput).Select(k => Number(IndexOf(s) + 0.5 + k)));
            return $"{{\"outputs\":[{{\"name\":\"{outputName}\",\"datatype\":\"FP32\",\"shape\":[{sequences.Count},{valuesPerInput}],\"data\":[{string.Join(",", data)}]}}]}}";
        }

        private static string TwoFragmentsPerInput(IReadOnlyList<string> sequences)
        {
            int n = 2 * sequences.Count;
            var annotations = sequences.SelectMany(_ => new[] { "\"b1+1\"", "\"y1+1\"" });
            var mz = sequences.SelectMany(_ => new[] { "100.5", "200.5" });
            var intensities = sequences.SelectMany(s => new[] { Number(IndexOf(s) + 0.5), "0.25" });
            return "{\"outputs\":["
                + $"{{\"name\":\"annotation\",\"datatype\":\"BYTES\",\"shape\":[{n}],\"data\":[{string.Join(",", annotations)}]}},"
                + $"{{\"name\":\"mz\",\"datatype\":\"FP32\",\"shape\":[{n}],\"data\":[{string.Join(",", mz)}]}},"
                + $"{{\"name\":\"intensities\",\"datatype\":\"FP32\",\"shape\":[{n}],\"data\":[{string.Join(",", intensities)}]}}]}}";
        }

        /// <summary>
        /// A stand-in for Koina: counts the requests, and either throws what <c>fail</c> picks for a
        /// request's sequences or answers it.
        /// </summary>
        private sealed class FakeKoina
        {
            private readonly Func<IReadOnlyList<string>, Exception?> _fail;
            private int _requests;
            private int _inFlight;
            private int _maxInFlight;

            public FakeKoina(Func<IReadOnlyList<string>, Exception?> fail) => _fail = fail;

            public static FakeKoina Refusing(params int[] indices) =>
                new(sequences => sequences.Any(s => indices.Contains(IndexOf(s))) ? Refusal() : null);

            /// <summary>The exception HTTP.InferenceRequest throws for the real refusal.</summary>
            public static Exception Refusal() => KoinaServiceException.ForFailedResponse(400, "Bad Request", RealRefusal);

            public int Requests => _requests;
            public int MaxInFlight => _maxInFlight;
            public TimeSpan Latency { get; set; } = TimeSpan.Zero;

            public async Task<string> Answer(Dictionary<string, object> request, Func<IReadOnlyList<string>, string> respond)
            {
                Interlocked.Increment(ref _requests);
                int inFlight = Interlocked.Increment(ref _inFlight);
                int seen;
                while ((seen = _maxInFlight) < inFlight && Interlocked.CompareExchange(ref _maxInFlight, inFlight, seen) != seen) { }
                try
                {
                    await Task.Delay(Latency);
                    var sequences = JObject.FromObject(request)["inputs"]!
                        .First(input => ((string)input["name"]!).StartsWith("peptide_sequences"))["data"]!
                        .Select(token => (string)token!).ToList();
                    if (_fail(sequences) is { } failure)
                        throw failure;
                    return respond(sequences);
                }
                finally
                {
                    Interlocked.Decrement(ref _inFlight);
                }
            }
        }

        private sealed class Peptide : IRetentionPredictable
        {
            public Peptide(string fullSequence) => FullSequence = fullSequence;
            public string FullSequence { get; }
            public string BaseSequence => FullSequence;
            public string FullSequenceWithMassShifts => FullSequence;
            public double MonoisotopicMass => 1000.0;
        }

        // ── the five families, each a real model with the fake Koina behind it ────

        private enum BatchingRule { InOrder, AllInOneRequest, OneRequestPerInputBelowMaxBatchSize }

        private sealed class RtModel : Prosit2019iRT
        {
            private readonly FakeKoina _koina;
            private readonly int _maxBatchSize;

            public RtModel(FakeKoina koina, int maxBatchSize = 1000,
                SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ReturnNull, int batchesPerRequest = 500)
                : base(mode, batchesPerRequest, 0)
            {
                _koina = koina;
                _maxBatchSize = maxBatchSize;
            }

            public override int MaxBatchSize => _maxBatchSize;
            public BatchingRule Batching { get; init; }
            public List<int> ParserCalls { get; } = new();

            protected override List<Dictionary<string, object>> ToBatchedRequests(List<RetentionTimePredictionInput> validInputs)
            {
                InputField Sequences(IEnumerable<RetentionTimePredictionInput> inputs) =>
                    new("peptide_sequences", "BYTES", inputs.Select(p => p.ValidatedFullSequence!).ToArray());

                return Batching switch
                {
                    BatchingRule.AllInOneRequest => new() { BuildBatchedRequest(0, Sequences(validInputs)) },
                    BatchingRule.OneRequestPerInputBelowMaxBatchSize when validInputs.Count < MaxBatchSize =>
                        validInputs.Select((p, i) => BuildBatchedRequest(i, Sequences(new[] { p }))).ToList(),
                    _ => base.ToBatchedRequests(validInputs),
                };
            }

            protected override List<PeptideRTPrediction> ResponseToPredictions(IReadOnlyList<string> responses, List<RetentionTimePredictionInput> requestInputs)
            {
                ParserCalls.Add(responses.Count);
                return base.ResponseToPredictions(responses, requestInputs);
            }

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken cancellationToken)
                => _koina.Answer(request, sequences => OneValuePerInput("irt", sequences));
        }

        private sealed class Rt40PtmModel : Prosit2025iRT40PTM
        {
            private readonly FakeKoina _koina;
            public Rt40PtmModel(FakeKoina koina) : base(SequenceConversionHandlingMode.ReturnNull, 500, 0) => _koina = koina;

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken cancellationToken)
                => _koina.Answer(request, sequences => OneValuePerInput("irt", sequences));
        }

        private sealed class FragmentModel : Prosit2019Intensity
        {
            private readonly FakeKoina _koina;
            private readonly int _maxBatchSize;

            public FragmentModel(FakeKoina koina, int maxBatchSize, SequenceConversionHandlingMode mode, int batchesPerRequest)
                : base(mode, IncompatibleParameterHandlingMode.ReturnNull, PredictionClients.Koina.Util.FragmentIonMappingMode.MapToValidatedFullSequence, batchesPerRequest, 0)
            {
                _koina = koina;
                _maxBatchSize = maxBatchSize;
            }

            public override int MaxBatchSize => _maxBatchSize;

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken cancellationToken)
                => _koina.Answer(request, TwoFragmentsPerInput);
        }

        private sealed class CcsModel : IM2Deep
        {
            private readonly FakeKoina _koina;
            private readonly int _maxBatchSize;

            public CcsModel(FakeKoina koina, int maxBatchSize, SequenceConversionHandlingMode mode, int batchesPerRequest)
                : base(mode, IncompatibleParameterHandlingMode.ReturnNull, batchesPerRequest, 0)
            {
                _koina = koina;
                _maxBatchSize = maxBatchSize;
            }

            public override int MaxBatchSize => _maxBatchSize;

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken cancellationToken)
                => _koina.Answer(request, sequences => OneValuePerInput("ccs", sequences));
        }

        private sealed class DetectabilityTestModel : PFly2024FineTuned
        {
            private readonly FakeKoina _koina;
            private readonly int _maxBatchSize;

            public DetectabilityTestModel(FakeKoina koina, int maxBatchSize, SequenceConversionHandlingMode mode, int batchesPerRequest)
                : base(mode, batchesPerRequest, 0)
            {
                _koina = koina;
                _maxBatchSize = maxBatchSize;
            }

            public override int MaxBatchSize => _maxBatchSize;

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken cancellationToken)
                => _koina.Answer(request, sequences => OneValuePerInput("detectability", sequences, valuesPerInput: 4));
        }

        private sealed class CrosslinkModel : Prosit2023IntensityXLCMS2
        {
            private readonly FakeKoina _koina;
            private readonly int _maxBatchSize;

            public CrosslinkModel(FakeKoina koina, int maxBatchSize, SequenceConversionHandlingMode mode, int batchesPerRequest)
                : base(mode, IncompatibleParameterHandlingMode.ReturnNull, batchesPerRequest, 0)
            {
                _koina = koina;
                _maxBatchSize = maxBatchSize;
            }

            public override int MaxBatchSize => _maxBatchSize;

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken cancellationToken)
                => _koina.Answer(request, TwoFragmentsPerInput);
        }
    }

    /// <summary>
    /// The classification behind the splitting: which 400 means the model refused an input, as
    /// opposed to a malformed request or a server that failed. No network; the bodies are real.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public class TestKoinaInputRejectedClassification
    {
        [Test]
        public void TheModelRefusingAnInputBecomesAKoinaInputRejectedException()
        {
            var thrown = KoinaServiceException.ForFailedResponse(400, "Bad Request", KoinaRefusedInputIsolationTests.RealRefusal);

            Assert.That(thrown, Is.InstanceOf<KoinaInputRejectedException>());
            Assert.That(thrown, Is.Not.InstanceOf<KoinaServiceException>());
            Assert.That(((KoinaInputRejectedException)thrown).ServerError, Does.Contain("E_954: no attachment group for residue 'E'"));
            Assert.That(thrown.StatusCode, Is.EqualTo(HttpStatusCode.BadRequest));
            Assert.That(thrown.Message, Does.StartWith("Request failed with status 400 Bad Request: "),
                "the same message as every other failed response");
        }

        [Test]
        public void TheSameBodyFromA5xxIsAnOutageNotARefusal()
        {
            var thrown = KoinaServiceException.ForFailedResponse(503, "Service Unavailable", KoinaRefusedInputIsolationTests.RealRefusal);

            Assert.That(thrown.GetType(), Is.EqualTo(typeof(HttpRequestException)));
            Assert.That(KoinaInputRejectedException.IsInputRejection(503, KoinaRefusedInputIsolationTests.RealRefusal), Is.False);
        }

        [TestCase("{\"error\":\"unexpected inference input 'peptide_sequenes' for model 'Prosit_2025_irt_40PTM'\"}")]
        [TestCase("{\"error\":\"expected 2 inputs but got 1\"}")]
        [TestCase("")]
        public void AMalformedRequestIsNotARefusal(string body)
        {
            Assert.That(KoinaInputRejectedException.IsInputRejection(400, body), Is.False);
            Assert.That(KoinaServiceException.ForFailedResponse(400, "Bad Request", body).GetType(), Is.EqualTo(typeof(HttpRequestException)));
        }

        /// <summary>
        /// A model's own code can also fail on the GPU. That is the server's fault, not the input's,
        /// and it stays a service fault even though the refusal marker is there too.
        /// </summary>
        [Test]
        public void AServiceFaultRaisedInsideTheModelIsStillAServiceFault()
        {
            const string body = "{\"error\":\"in ensemble 'X', Failed to process the request(s) for model instance 'X_0', message: RuntimeError: CUDA error: device-side assert triggered\"}";

            Assert.That(KoinaInputRejectedException.IsInputRejection(400, body), Is.False);
            Assert.That(KoinaServiceException.ForFailedResponse(400, "Bad Request", body), Is.InstanceOf<KoinaServiceException>());
        }

        [Test]
        public void ExistingCatchClausesStillSeeItAsAnHttpRequestException()
        {
            var e = new KoinaInputRejectedException("Request failed with status 400", "KeyError");

            Assert.That(e, Is.InstanceOf<HttpRequestException>(), "callers already catching HttpRequestException must keep working unchanged");
            Assert.That(e.ServerError, Is.EqualTo("KeyError"));
        }
    }

    /// <summary>
    /// The live canary for the case that started this: Prosit_2025_irt_40PTM refuses zinc on E,
    /// E[UNIMOD:954]. Prosit2025ModificationRules now refuses it before sending, so this sends it through
    /// a copy of the model without that rule. If Koina ever accepts it, this fails and says the case has
    /// changed, and the rule is out of date; if Koina is down or faulting, KoinaLiveTestFixture skips it.
    /// </summary>
    [TestFixture]
    [Category("ExternalService")]
    [Category("Koina")]
    [ExcludeFromCodeCoverage]
    public class KoinaRefusedInputLiveTests : KoinaLiveTestFixture
    {
        [Test]
        public void Prosit2025iRT40PTM_RefusesZincOnE_AndPredictsTheOthers()
        {
            // ProForma, the way the MSL pipeline sends modified peptides to Koina.
            var sequences = new[] { "PEPTIDEK", "LGGNEQVTR", "PEPE[UNIMOD:954]TIDEK", "AGGSSEPVTGLADK", "VEATFGVDESNAK" };
            var model = new Prosit2025iRT40PTMWithoutModificationRules();

            var predictions = model.Predict(sequences
                .Select(s => new RetentionTimePredictionInput(s) { SequenceParser = ProFormaSequenceParser.Instance }).ToList());

            Assert.That(predictions.Count, Is.EqualTo(sequences.Length));
            Assert.That(predictions[2].ValidatedFullSequence, Is.EqualTo("PEPE[UNIMOD:954]TIDEK"),
                "the canary only means something if the refused peptide reached Koina");
            Assert.That(predictions[2].PredictedRetentionTime, Is.Null);
            Assert.That(model.ValidInputsMask[2], Is.False);
            Assert.That(predictions[2].Warning?.Message, Does.Contain("no attachment group for residue 'E'"));
            Assert.That(predictions[2].Warning?.InnerException, Is.InstanceOf<KoinaInputRejectedException>());
            foreach (int i in new[] { 0, 1, 3, 4 })
            {
                Assert.That(predictions[i].PredictedRetentionTime, Is.Not.Null, $"{sequences[i]} is predicted");
                Assert.That(model.ValidInputsMask[i], Is.True);
            }
        }

        /// <summary>
        /// Prosit2025iRT40PTM as it was before Prosit2025ModificationRules: the same converter and request,
        /// with base's cleaning only, so E[UNIMOD:954] reaches Koina.
        /// </summary>
        private sealed class Prosit2025iRT40PTMWithoutModificationRules : RetentionTimeModel
        {
            private static readonly ISequenceConverter Converter =
                CreateUnimodConverterAcceptAll(new UnimodSequenceFormatSchema(UnimodLabelStyle.UpperCase, '[', ']', "-", "-"));

            public Prosit2025iRT40PTMWithoutModificationRules() : base(Converter) { }

            public override string ModelName => "Prosit_2025_irt_40PTM";
            public override int MaxBatchSize => 1000;
            public override int MaxNumberOfBatchesPerRequest { get; init; } = 500;
            public override int ThrottlingDelayInMilliseconds { get; init; } = 100;
            public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1500;
            public override int MaxPeptideLength => 30;
            public override int MinPeptideLength => 1;
            public override bool IsIndexedRetentionTimeModel => true;
            public override bool AcceptsAllUnimodModifications => true;
            public override SequenceConversionHandlingMode ModHandlingMode { get; init; } = SequenceConversionHandlingMode.ReturnNull;

            protected override List<Dictionary<string, object>> ToBatchedRequests(List<RetentionTimePredictionInput> validInputs) =>
                validInputs.Select(p => p.ValidatedFullSequence!).Chunk(MaxBatchSize)
                    .Select((batch, i) => BuildBatchedRequest(i, new InputField("peptide_sequences", "BYTES", batch)))
                    .ToList();
        }
    }
}
