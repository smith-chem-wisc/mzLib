using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;
using NUnit.Framework;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using Readers.ProForma;

namespace Test.KoinaTests
{
    /// <summary>
    /// Offline tests for the (residue, UNIMOD) pairs Koina's Prosit 2025 22PTM and 40PTM models refuse, and for how the
    /// four models refuse them before sending. Expected outcomes were checked against the live server on 2026-11-09:
    /// every "refused" case below returned HTTP 400 there, and every "accepted" case returned 200.
    /// </summary>
    [TestFixture]
    public class Prosit2025ModificationRulesTests
    {
        private const string ZincOnE = "PEPE[UNIMOD:954]TIDEK";
        private const string OxidizedM = "PEPM[Common Variable:Oxidation on M]TIDEK";
        private static readonly ISequenceParser ProForma = ProFormaSequenceParser.Instance;

        #region The rule, token by token

        /// <summary>
        /// Trained pairs, including those on residues with no attachment group (P, N, E) and on either terminus.
        /// K[UNIMOD:5634] is trained and Koina skips its UNIMOD lookup, so it passes although mzLib has no such id.
        /// </summary>
        [TestCase("M[UNIMOD:35]")]
        [TestCase("C[UNIMOD:4]")]
        [TestCase("P[UNIMOD:35]")]
        [TestCase("N[UNIMOD:7]")]
        [TestCase("E[UNIMOD:27]")]
        [TestCase("Q[UNIMOD:28]")]
        [TestCase("R[UNIMOD:267]")]
        [TestCase("K[UNIMOD:5634]")]
        [TestCase("[UNIMOD:737]-")]
        [TestCase("-[UNIMOD:1]")]
        public void TrainedTokens_AreAccepted(string token)
        {
            Assert.That(Prosit2025ModificationRules.RefusalReason(token), Is.Null);
        }

        /// <summary>
        /// Untrained pairs Koina encodes zero-shot: an attachment group, a readable H C N O P S composition, and no
        /// negative count. K[UNIMOD:259] checks that isotopes count as their element (13C(6) 15N(2) minus C(6) N(2)).
        /// </summary>
        [TestCase("C[UNIMOD:21]")]
        [TestCase("K[UNIMOD:64]")]
        [TestCase("K[UNIMOD:259]")]
        [TestCase("K[UNIMOD:122]")]
        [TestCase("M[UNIMOD:5]")]
        [TestCase("Y[UNIMOD:35]")]
        [TestCase("[UNIMOD:36]-")]
        [TestCase("-[UNIMOD:34]")]
        public void ZeroShotTokens_AreAccepted(string token)
        {
            Assert.That(Prosit2025ModificationRules.RefusalReason(token), Is.Null);
        }

        [TestCase("E[UNIMOD:954]")]
        [TestCase("N[UNIMOD:35]")]
        [TestCase("N[UNIMOD:385]")]
        [TestCase("A[UNIMOD:35]")]
        [TestCase("E[UNIMOD:299]")]
        [TestCase("G[UNIMOD:21]")]
        public void ResidueWithoutAttachmentGroup_IsRefused(string token)
        {
            Assert.That(Prosit2025ModificationRules.RefusalReason(token),
                Does.Contain($"residue {token[0]} has no attachment group").And.Contain(token));
        }

        [TestCase("K[UNIMOD:954]", "Zn")]
        [TestCase("D[UNIMOD:954]", "Zn")]
        [TestCase("D[UNIMOD:530]", "K")]
        [TestCase("D[UNIMOD:956]", "Mg")]
        [TestCase("D[UNIMOD:30]", "Na")]
        [TestCase("Y[UNIMOD:340]", "Br")]
        public void ElementOutsideHCNOPS_IsRefused(string token, string element)
        {
            Assert.That(Prosit2025ModificationRules.RefusalReason(token),
                Does.Contain($"contains {element}").And.Contain("only H, C, N, O, P and S"));
        }

        [TestCase("S[UNIMOD:295]")]   // dHex
        [TestCase("S[UNIMOD:1436]")]  // Me dHex Hex HexNAc
        [TestCase("-[UNIMOD:1948]")]  // Sulf
        public void CompositionKoinaCannotRead_IsRefused(string token)
        {
            Assert.That(Prosit2025ModificationRules.RefusalReason(token), Does.Contain("cannot read the composition"));
        }

        [TestCase("[UNIMOD:27]-", "the N-terminus")]
        [TestCase("[UNIMOD:28]-", "the N-terminus")]
        [TestCase("[UNIMOD:385]-", "the N-terminus")]
        [TestCase("-[UNIMOD:2]", "the C-terminus")]
        [TestCase("S[UNIMOD:23]", "residue S")]
        [TestCase("W[UNIMOD:351]", "residue W")]
        public void NegativeAtomCount_IsRefused(string token, string site)
        {
            Assert.That(Prosit2025ModificationRules.RefusalReason(token),
                Does.Contain($"the attachment group of {site}").And.Contain("negative"));
        }

        /// <summary>
        /// Ids newer than Koina's unimod.obo (2072), ids in its trained table that the file lacks (12118), and ids no
        /// catalog has all fail the server's UNIMOD lookup, on any residue.
        /// </summary>
        [TestCase("K[UNIMOD:2073]")]
        [TestCase("K[UNIMOD:12118]")]
        [TestCase("E[UNIMOD:99999]")]
        [TestCase("K[UNIMOD:15]")]    // a gap in UNIMOD's numbering
        public void IdKoinaDoesNotKnow_IsRefused(string token)
        {
            Assert.That(Prosit2025ModificationRules.RefusalReason(token), Does.Contain("is not in Koina's UNIMOD table"));
        }

        [TestCase("E954")]
        [TestCase("[UNIMOD:1]")]
        [TestCase("EE[UNIMOD:954]")]
        public void RefusalReason_RejectsAStringThatIsNotAToken(string notAToken)
        {
            Assert.Throws<ArgumentException>(() => Prosit2025ModificationRules.RefusalReason(notAToken));
        }

        #endregion

        #region Sequences

        [Test]
        public void FindRefusedModifications_NamesEachRefusedTokenWhereItSits()
        {
            var refused = Prosit2025ModificationRules.FindRefusedModifications("[UNIMOD:27]-EPEM[UNIMOD:35]E[UNIMOD:954]K-[UNIMOD:2]");

            Assert.That(refused.Select(r => r.Token), Is.EqualTo(new[] { "[UNIMOD:27]-", "E[UNIMOD:954]", "-[UNIMOD:2]" }));
            Assert.That(refused.Select(r => r.Position), Is.EqualTo(new[]
            {
                ModificationPositionType.NTerminus, ModificationPositionType.Residue, ModificationPositionType.CTerminus
            }));
            Assert.That(refused.Select(r => r.ResidueIndex), Is.EqualTo(new int?[] { null, 4, null }));
        }

        [TestCase("PEPTIDEK")]
        [TestCase("[UNIMOD:1]-PEPC[UNIMOD:4]TM[UNIMOD:35]IDEK-[UNIMOD:34]")]
        public void FindRefusedModifications_IsEmptyWhenKoinaAcceptsEverything(string koinaSequence)
        {
            Assert.That(Prosit2025ModificationRules.FindRefusedModifications(koinaSequence), Is.Empty);
        }

        /// <summary>
        /// Sequences the MSL rat pilot sent to Prosit_2025_irt_40PTM, each of which Koina answered with HTTP 400, and the
        /// pair it named. Each is refused here, so it would not have been sent.
        /// </summary>
        [TestCase("GGDDLD[UNIMOD:530]PNYVLSSR", "D[UNIMOD:530]")]
        [TestCase("LISWYD[UNIMOD:954]NEYGYSNR", "D[UNIMOD:954]")]
        [TestCase("[UNIMOD:5]-HIADLAGNPE[UNIMOD:30]VILPVPAFNVINGGSHAGNK", "E[UNIMOD:30]")]
        [TestCase("IGEHTPSSLAIM[UNIMOD:5]ENAN[UNIMOD:385]VLAR", "N[UNIMOD:385]")]
        [TestCase("SQC[UNIMOD:275]LAEIEHDN[UNIMOD:35]IPADLPSIAADFVEDKEVC[UNIMOD:4]KNYAEAK", "N[UNIMOD:35]")]
        [TestCase("WGDAGAEYVVE[UNIMOD:299]STGIFTTM[UNIMOD:5]EK", "E[UNIMOD:299]")]
        [TestCase("VDGDIC[UNIMOD:4]SLNE[UNIMOD:952]FSITGSTYAPEGEVLK", "E[UNIMOD:952]")]
        [TestCase("LPADTE[UNIMOD:1870]VVC[UNIMOD:4]APPTAYIDFAR", "E[UNIMOD:1870]")]
        [TestCase("AGHPFMWNE[UNIMOD:531]HLGY[UNIMOD:21]VLTC[UNIMOD:4]PSNLGTGLR", "E[UNIMOD:531]")]
        [TestCase("LSVE[UNIMOD:951]ALNSLTGEFK", "E[UNIMOD:951]")]
        [TestCase("K[UNIMOD:122]AGHPFMWNE[UNIMOD:956]HLGYVLTC[UNIMOD:4]PSNLGTGLR", "E[UNIMOD:956]")]
        public void PilotSequencesKoinaRejected_AreRefused(string sent, string pairKoinaNamed)
        {
            Assert.That(Prosit2025ModificationRules.FindRefusedModifications(sent).Select(r => r.Token), Is.EqualTo(new[] { pairKoinaNamed }));
        }

        [TestCase("[UNIMOD:27]-EPEPTIDEK", "E[UNIMOD:27]PEPTIDEK")]
        [TestCase("[UNIMOD:28]-QPEPTIDEK", "Q[UNIMOD:28]PEPTIDEK")]
        [TestCase("[UNIMOD:27]-E", "E[UNIMOD:27]")]
        [TestCase("[UNIMOD:28]-QPEPTIDEK-[UNIMOD:34]", "Q[UNIMOD:28]PEPTIDEK-[UNIMOD:34]")]
        public void PlacePyroGluOnItsResidue_MovesItWhereKoinaTrainedIt(string nTerminal, string onResidue)
        {
            Assert.That(Prosit2025ModificationRules.PlacePyroGluOnItsResidue(nTerminal), Is.EqualTo(onResidue));
            Assert.That(Prosit2025ModificationRules.FindRefusedModifications(onResidue), Is.Empty);
        }

        /// <summary>
        /// Only Glu->pyro-Glu (27) on E and Gln->pyro-Glu (28) on Q move. Anything else, or a residue that already
        /// carries a modification, goes to the rule as written.
        /// </summary>
        [TestCase("[UNIMOD:27]-QPEPTIDEK")]
        [TestCase("[UNIMOD:28]-EPEPTIDEK")]
        [TestCase("[UNIMOD:27]-E[UNIMOD:34]PEPTIDEK")]
        [TestCase("[UNIMOD:385]-CPEPTIDEK")]
        [TestCase("PEPTIDEK")]
        public void PlacePyroGluOnItsResidue_LeavesEverythingElseAlone(string koinaSequence)
        {
            Assert.That(Prosit2025ModificationRules.PlacePyroGluOnItsResidue(koinaSequence), Is.EqualTo(koinaSequence));
        }

        #endregion

        #region The four models

        private static IEnumerable<TestCaseData> FourModels()
        {
            yield return new TestCaseData((Func<SequenceConversionHandlingMode, IProbe>)(mode => new Irt40Probe(mode))).SetArgDisplayNames("Prosit_2025_irt_40PTM");
            yield return new TestCaseData((Func<SequenceConversionHandlingMode, IProbe>)(mode => new Irt22Probe(mode))).SetArgDisplayNames("Prosit_2025_irt_22PTM");
            yield return new TestCaseData((Func<SequenceConversionHandlingMode, IProbe>)(mode => new Intensity40Probe(mode))).SetArgDisplayNames("Prosit_2025_intensity_40PTM");
            yield return new TestCaseData((Func<SequenceConversionHandlingMode, IProbe>)(mode => new Intensity22Probe(mode))).SetArgDisplayNames("Prosit_2025_intensity_22PTM");
        }

        /// <summary>
        /// The original failure: one E[UNIMOD:954] in a batch. It is masked out and never sent, and the rest of the batch
        /// is predicted. Every UNIMOD id is still accepted (K[UNIMOD:64] is trained nowhere); only the pair Koina cannot
        /// place is refused. ProForma input, as the MSL pipeline sends it: on master, base resolves "Metal:Zinc on E"
        /// in mzLib notation to no UNIMOD id and refuses it before this rule runs.
        /// </summary>
        [TestCaseSource(nameof(FourModels))]
        public void ReturnNull_RefusedInputIsMaskedAndNotSent(Func<SequenceConversionHandlingMode, IProbe> create)
        {
            var probe = create(SequenceConversionHandlingMode.ReturnNull);

            var warnings = probe.Predict(ProForma, "PEPTIDEK", ZincOnE, "PEPM[UNIMOD:35]TIDEK", "PEPK[UNIMOD:64]TIDEK");

            Assert.That(probe.Mask, Is.EqualTo(new[] { true, false, true, true }));
            Assert.That(probe.Sent, Is.EqualTo(new[] { "PEPTIDEK", "PEPM[UNIMOD:35]TIDEK", "PEPK[UNIMOD:64]TIDEK" }));
            Assert.That(warnings[1], Does.Contain(probe.ModelName).And.Contain("E[UNIMOD:954]").And.Contain("residue E has no attachment group"));
        }

        [TestCaseSource(nameof(FourModels))]
        public void ReturnNull_ReturnsNullWithAWarning(Func<SequenceConversionHandlingMode, IProbe> create)
        {
            var probe = create(SequenceConversionHandlingMode.ReturnNull);

            var cleaned = probe.Clean(ZincOnE, ProForma, out var api, out var warning);

            Assert.That(cleaned, Is.Null);
            Assert.That(api, Is.Null);
            Assert.That(warning!.Message, Does.Contain("E[UNIMOD:954]"));
        }

        [TestCaseSource(nameof(FourModels))]
        public void ThrowException_Throws(Func<SequenceConversionHandlingMode, IProbe> create)
        {
            var probe = create(SequenceConversionHandlingMode.ThrowException);

            var ex = Assert.Throws<ArgumentException>(() => probe.Clean(ZincOnE, ProForma, out _, out _));
            Assert.That(ex!.Message, Does.Contain("E[UNIMOD:954]"));
            Assert.That(probe.Clean(OxidizedM, null, out var api, out _), Is.EqualTo(OxidizedM));
            Assert.That(api, Is.EqualTo("PEPM[UNIMOD:35]TIDEK"));
        }

        /// <summary>
        /// Each refused modification is dropped and the rest kept, in the sequence sent and in the mzLib sequence
        /// returned, which keep describing the same peptide. The pyro-Glu still moves to its residue first.
        /// </summary>
        [TestCaseSource(nameof(FourModels))]
        public void RemoveIncompatibleElements_DropsOnlyTheRefusedModifications(Func<SequenceConversionHandlingMode, IProbe> create)
        {
            var probe = create(SequenceConversionHandlingMode.RemoveIncompatibleElements);
            const string sequence = "[UNIMOD:27]-EPEC[UNIMOD:4]TE[UNIMOD:954]DEM[UNIMOD:35]K-[UNIMOD:2]";

            // The third also has a mass shift base removes itself; both warnings are kept.
            var warnings = probe.Predict(ProForma, sequence, "[UNIMOD:385]-CPEM[UNIMOD:35]TIDEK", "PEPE[UNIMOD:954]TIDE[+12.3456]K");
            var cleaned = probe.Clean(sequence, ProForma, out var api, out var warning);

            Assert.That(probe.Mask, Is.EqualTo(new[] { true, true, true }));
            Assert.That(probe.Sent, Is.EqualTo(new[] { "E[UNIMOD:27]PEC[UNIMOD:4]TEDEM[UNIMOD:35]K", "CPEM[UNIMOD:35]TIDEK", "PEPETIDEK" }));
            Assert.That(api, Is.EqualTo("E[UNIMOD:27]PEC[UNIMOD:4]TEDEM[UNIMOD:35]K"));
            Assert.That(cleaned, Is.EqualTo("[UNIMOD:27]EPEC[UNIMOD:4]TEDEM[UNIMOD:35]K"));
            Assert.That(warning!.Message, Does.StartWith($"Removing modifications {probe.ModelName} cannot place")
                .And.Contain("E[UNIMOD:954]").And.Contain("-[UNIMOD:2]"));
            Assert.That(warnings[1], Does.Contain("[UNIMOD:385]-"));
            Assert.That(warnings[2], Does.Contain("E[UNIMOD:954]").And.Contain("+12.3456"));
        }

        /// <summary>
        /// mzLib's N-terminal pyro-Glu is sent on its residue, where Koina trained it. The mzLib sequence, which is what a
        /// library entry is built from, keeps it N-terminal.
        /// </summary>
        [TestCaseSource(nameof(FourModels))]
        public void NTerminalPyroGlu_IsSentOnItsResidue(Func<SequenceConversionHandlingMode, IProbe> create)
        {
            var probe = create(SequenceConversionHandlingMode.ReturnNull);

            var warnings = probe.Predict(ProForma, "[UNIMOD:27]-EPEPTIDEK", "[UNIMOD:28]-QPEPTIDEK", "[UNIMOD:28]-EPEPTIDEK");

            Assert.That(probe.Mask, Is.EqualTo(new[] { true, true, false }));
            Assert.That(probe.Sent, Is.EqualTo(new[] { "E[UNIMOD:27]PEPTIDEK", "Q[UNIMOD:28]PEPTIDEK" }));
            Assert.That(warnings[2], Does.Contain("[UNIMOD:28]-").And.Contain("the N-terminus"));
            Assert.That(probe.Clean("[UNIMOD:27]-EPEPTIDEK", ProForma, out _, out _), Is.EqualTo("[UNIMOD:27]EPEPTIDEK"));
        }

        /// <summary>
        /// Koina accepts a C-terminal modification it can encode and refuses one it cannot.
        /// </summary>
        [TestCaseSource(nameof(FourModels))]
        public void CTerminalModifications_FollowTheRule(Func<SequenceConversionHandlingMode, IProbe> create)
        {
            var probe = create(SequenceConversionHandlingMode.ReturnNull);

            probe.Predict(ProForma, "PEPTIDEK-[UNIMOD:34]", "PEPTIDEK-[UNIMOD:2]");

            Assert.That(probe.Mask, Is.EqualTo(new[] { true, false }));
            Assert.That(probe.Sent, Is.EqualTo(new[] { "PEPTIDEK-[UNIMOD:34]" }));
        }

        /// <summary>
        /// Base's own refusals are untouched: an invalid residue is still refused by base, before this rule runs.
        /// </summary>
        [TestCaseSource(nameof(FourModels))]
        public void BaseRefusals_AreLeftToBase(Func<SequenceConversionHandlingMode, IProbe> create)
        {
            var probe = create(SequenceConversionHandlingMode.ReturnNull);

            Assert.That(probe.Clean("PEPTIDEBK", null, out var api, out var warning), Is.Null);
            Assert.That(api, Is.Null);
            Assert.That(warning!.Message, Does.Contain("Invalid base sequence"));
        }

        #endregion

        #region Probes: the four models with a captured transport

        public interface IProbe
        {
            string ModelName { get; }
            List<string> Sent { get; }
            bool[] Mask { get; }

            /// <summary>Runs Predict, reading each input with the given parser, and returns each input's warning ("" when none).</summary>
            List<string> Predict(ISequenceParser? parser, params string[] sequences);

            string? Clean(string sequence, ISequenceParser? parser, out string? api, out WarningException? warning);
        }

        /// <summary>
        /// Records the peptide_sequences of each request and answers with one value per peptide, in the output layout
        /// of the model family.
        /// </summary>
        private static Task<string> Capture(Dictionary<string, object> request, List<string> sent, bool intensity)
        {
            var inputs = (JArray)JObject.Parse(JsonConvert.SerializeObject(request))["inputs"]!;
            var peptides = inputs.First(i => (string)i["name"]! == "peptide_sequences")["data"]!.Select(t => (string)t!).ToList();
            sent.AddRange(peptides);
            int n = peptides.Count;
            string json = intensity
                ? "{\"outputs\":[" +
                  $"{{\"name\":\"annotation\",\"datatype\":\"BYTES\",\"shape\":[{n}],\"data\":[{string.Join(",", Enumerable.Repeat("\"b1+1\"", n))}]}}," +
                  $"{{\"name\":\"mz\",\"datatype\":\"FP32\",\"shape\":[{n}],\"data\":[{string.Join(",", Enumerable.Repeat("100.0", n))}]}}," +
                  $"{{\"name\":\"intensities\",\"datatype\":\"FP32\",\"shape\":[{n}],\"data\":[{string.Join(",", Enumerable.Repeat("0.5", n))}]}}]}}"
                : $"{{\"outputs\":[{{\"name\":\"irt\",\"datatype\":\"FP32\",\"shape\":[{n},1],\"data\":[{string.Join(",", Enumerable.Repeat("42.5", n))}]}}]}}";
            return Task.FromResult(json);
        }

        private static List<string> Warnings(IEnumerable<WarningException?> warnings) => warnings.Select(w => w?.Message ?? "").ToList();

        private sealed class Irt40Probe : Prosit2025iRT40PTM, IProbe
        {
            public Irt40Probe(SequenceConversionHandlingMode mode) : base(mode) { }
            public List<string> Sent { get; } = new();
            public bool[] Mask => ValidInputsMask;
            public List<string> Predict(ISequenceParser? parser, params string[] sequences) =>
                Warnings(Predict(sequences.Select(s => new RetentionTimePredictionInput(s) { SequenceParser = parser }).ToList()).Select(p => p.Warning));
            public string? Clean(string sequence, ISequenceParser? parser, out string? api, out WarningException? warning) => TryCleanSequence(sequence, parser, out api, out warning);
            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct) => Capture(request, Sent, false);
        }

        private sealed class Irt22Probe : Prosit2025iRT22PTM, IProbe
        {
            public Irt22Probe(SequenceConversionHandlingMode mode) : base(mode) { }
            public List<string> Sent { get; } = new();
            public bool[] Mask => ValidInputsMask;
            public List<string> Predict(ISequenceParser? parser, params string[] sequences) =>
                Warnings(Predict(sequences.Select(s => new RetentionTimePredictionInput(s) { SequenceParser = parser }).ToList()).Select(p => p.Warning));
            public string? Clean(string sequence, ISequenceParser? parser, out string? api, out WarningException? warning) => TryCleanSequence(sequence, parser, out api, out warning);
            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct) => Capture(request, Sent, false);
        }

        private sealed class Intensity40Probe : Prosit2025Intensity40PTM, IProbe
        {
            public Intensity40Probe(SequenceConversionHandlingMode mode) : base(mode) { }
            public List<string> Sent { get; } = new();
            public bool[] Mask => ValidInputsMask;
            public List<string> Predict(ISequenceParser? parser, params string[] sequences) =>
                Warnings(Predict(sequences.Select(s => new FragmentIntensityPredictionInput(s, 2, 30, null, "HCD") { SequenceParser = parser }).ToList()).Select(p => p.Warning));
            public string? Clean(string sequence, ISequenceParser? parser, out string? api, out WarningException? warning) => TryCleanSequence(sequence, parser, out api, out warning);
            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct) => Capture(request, Sent, true);
        }

        private sealed class Intensity22Probe : Prosit2025Intensity22PTM, IProbe
        {
            public Intensity22Probe(SequenceConversionHandlingMode mode) : base(mode) { }
            public List<string> Sent { get; } = new();
            public bool[] Mask => ValidInputsMask;
            public List<string> Predict(ISequenceParser? parser, params string[] sequences) =>
                Warnings(Predict(sequences.Select(s => new FragmentIntensityPredictionInput(s, 2, 30, null, "HCD") { SequenceParser = parser }).ToList()).Select(p => p.Warning));
            public string? Clean(string sequence, ISequenceParser? parser, out string? api, out WarningException? warning) => TryCleanSequence(sequence, parser, out api, out warning);
            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct) => Capture(request, Sent, true);
        }

        #endregion
    }
}
