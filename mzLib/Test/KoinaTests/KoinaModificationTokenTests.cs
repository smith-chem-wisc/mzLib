using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text.RegularExpressions;
using NUnit.Framework;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using Readers.ProForma;
using WarningException = System.ComponentModel.WarningException;

namespace Test.KoinaTests
{
    /// <summary>
    /// Builds, for any model, the peptides that exercise its declared modification tokens. The peptides are
    /// written with UNIMOD ids, which is both ProForma and what the crosslink models take as they are.
    /// </summary>
    internal static class ModificationTokenCases
    {
        private static readonly Assembly KoinaAssembly = typeof(FragmentIntensityModel).Assembly;

        public static IEnumerable<Type> AllModels() =>
            KoinaAssembly.GetTypes().Where(t => !t.IsAbstract && t.GetProperty("AllowedModificationTokens") != null);

        public static IEnumerable<Type> TokenDeclaringModels() => AllModels().Where(t => Tokens(Instantiate(t)) != null);

        public static object Instantiate(Type modelType, SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ReturnNull)
        {
            var ctor = modelType.GetConstructors().First(c => c.GetParameters().All(p => p.HasDefaultValue));
            return ctor.Invoke(ctor.GetParameters()
                .Select(p => p.ParameterType == typeof(SequenceConversionHandlingMode) ? mode : p.DefaultValue)
                .ToArray());
        }

        public static IReadOnlySet<string>? Tokens(object model) =>
            (IReadOnlySet<string>?)model.GetType().GetProperty("AllowedModificationTokens")!.GetValue(model);

        public static IReadOnlySet<int> AllowedIds(object model) =>
            (IReadOnlySet<int>)model.GetType().GetProperty("AllowedUnimodIds")!.GetValue(model)!;

        public static bool AcceptsAllIds(object model) =>
            (bool)model.GetType().GetProperty("AcceptsAllUnimodModifications")!.GetValue(model)!;

        public static int IdOf(string token) => int.Parse(Regex.Match(token, @"\d+").Value);

        /// <summary>
        /// A peptide carrying the given tokens, plus the N-terminal label the model requires when none of them is one.
        /// </summary>
        public static string PeptideWith(object model, params string[] tokens)
        {
            var nTerminal = tokens.SingleOrDefault(t => t.EndsWith('-'));
            if (nTerminal == null
                && model.GetType().GetProperty("RequiredNTerminalUnimodIds")!.GetValue(model) is IReadOnlySet<int> { Count: > 0 } required)
            {
                nTerminal = $"[UNIMOD:{required.First()}]-";
            }

            return $"{nTerminal}AEPT{string.Concat(tokens.Where(t => !t.EndsWith('-')))}IDER";
        }

        /// <summary>
        /// An id the model allows, on a residue it has no token for.
        /// </summary>
        public static string UndeclaredToken(object model)
        {
            var tokens = Tokens(model)!;
            int id = IdOf(tokens.First(t => !t.EndsWith('-')));
            return "ACDEFGHIKLMNPQRSTVWY".Select(residue => $"{residue}[UNIMOD:{id}]").First(t => !tokens.Contains(t));
        }

        public static string? Clean(object model, string sequence, ISequenceParser? parser, out string? apiSequence, out WarningException? warning)
        {
            var method = model.GetType().GetMethod("TryCleanSequence", BindingFlags.NonPublic | BindingFlags.Instance)!;
            var args = new object?[] { sequence, parser, null, null };
            var result = (string?)method.Invoke(model, args);
            warning = (WarningException?)args[3];

            // The crosslink models send the validated sequence as it is.
            if (model is CrosslinkFragmentIntensityModel)
                apiSequence = result;
            else if (args[2] is CanonicalSequence cleaned)
                apiSequence = (string?)model.GetType()
                    .GetMethod("SerializeKoinaSequence", BindingFlags.NonPublic | BindingFlags.Instance)!
                    .Invoke(model, new object?[] { cleaned, null });
            else
                apiSequence = null;
            return result;
        }
    }

    [TestFixture]
    public class KoinaModificationTokenTests
    {
        // UniSpec's preprocessing keys on the id alone, so it has no residue tokens to declare.
        [Test]
        public void EveryModelThatListsItsIds_DeclaresTokens_ExceptUniSpec()
        {
            var listingIdsWithoutTokens = ModificationTokenCases.AllModels()
                .Where(t =>
                {
                    var model = ModificationTokenCases.Instantiate(t);
                    return !ModificationTokenCases.AcceptsAllIds(model)
                        && ModificationTokenCases.AllowedIds(model).Count > 0
                        && ModificationTokenCases.Tokens(model) == null;
                })
                .Select(t => t.Name);

            Assert.That(listingIdsWithoutTokens, Is.EquivalentTo(new[] { nameof(UniSpec) }));
        }

        [TestCaseSource(typeof(ModificationTokenCases), nameof(ModificationTokenCases.TokenDeclaringModels))]
        public void AllowedUnimodIds_AreTheIdsOfTheDeclaredTokens(Type modelType)
        {
            var model = ModificationTokenCases.Instantiate(modelType);
            var tokens = ModificationTokenCases.Tokens(model)!;

            // A residue or the N-terminus with one id: the only forms a single modification can take.
            Assert.That(tokens, Has.All.Match(@"^([A-Z]\[UNIMOD:\d+\]|\[UNIMOD:\d+\]-)$"));
            Assert.That(ModificationTokenCases.AllowedIds(model), Is.EquivalentTo(tokens.Select(ModificationTokenCases.IdOf).Distinct()));
            Assert.That(ModificationTokenCases.AcceptsAllIds(model), Is.False);
        }

        [TestCaseSource(typeof(ModificationTokenCases), nameof(ModificationTokenCases.TokenDeclaringModels))]
        public void EveryDeclaredToken_PassesValidation_AndTheSameIdOnAnUndeclaredResidueDoesNot(Type modelType)
        {
            var model = ModificationTokenCases.Instantiate(modelType);

            Assert.Multiple(() =>
            {
                foreach (var token in ModificationTokenCases.Tokens(model)!)
                {
                    var sequence = ModificationTokenCases.PeptideWith(model, token);
                    var result = ModificationTokenCases.Clean(model, sequence, ProFormaSequenceParser.Instance, out var api, out var warning);

                    Assert.That(result, Is.Not.Null, $"{sequence}: {warning?.Message}");
                    Assert.That(api, Does.Contain(token), sequence);
                }
            });

            var undeclared = ModificationTokenCases.PeptideWith(model, ModificationTokenCases.UndeclaredToken(model));
            var rejected = ModificationTokenCases.Clean(model, undeclared, ProFormaSequenceParser.Instance, out var rejectedApi, out var rejectedWarning);

            Assert.That(rejected, Is.Null, undeclared);
            Assert.That(rejectedApi, Is.Null, undeclared);
            Assert.That(rejectedWarning?.Message, Does.Contain("unsupported modification"), undeclared);
        }

        [Test]
        public void Hcd_AllowedIdOnDeclaredResidue_PassesFromMzLibAndProFormaSources()
        {
            var model = new HcdProbe();

            var fromMzLib = model.Clean("PEPTM[Common Variable:Oxidation on M]IDEK", null, out var mzLibApi, out var mzLibWarning);
            var fromProForma = model.Clean("PEPTM[UNIMOD:35]IDEK", ProFormaSequenceParser.Instance, out var proFormaApi, out var proFormaWarning);

            Assert.That(fromMzLib, Is.EqualTo("PEPTM[Common Variable:Oxidation on M]IDEK"));
            Assert.That(fromProForma, Is.EqualTo("PEPTM[UNIMOD:35]IDEK"));
            Assert.That(mzLibApi, Is.EqualTo("PEPTM[UNIMOD:35]IDEK"));
            Assert.That(proFormaApi, Is.EqualTo(mzLibApi));
            Assert.That(mzLibWarning, Is.Null);
            Assert.That(proFormaWarning, Is.Null);
        }

        // Oxidation (UNIMOD:35) is allowed, but the model only has a token for it on M.
        [TestCase(SequenceConversionHandlingMode.ReturnNull, null)]
        [TestCase(SequenceConversionHandlingMode.RemoveIncompatibleElements, "PEPTIDEK")]
        [TestCase(SequenceConversionHandlingMode.UsePrimarySequence, "PEPTIDEK")]
        public void Hcd_AllowedIdOnUndeclaredResidue_SameOutcomeFromMzLibAndProFormaSources(SequenceConversionHandlingMode mode, string? expected)
        {
            var model = new HcdProbe(mode);

            var fromMzLib = model.Clean("PEP[Unimod:Oxidation on P]TIDEK", null, out var mzLibApi, out var mzLibWarning);
            var fromProForma = model.Clean("PEP[UNIMOD:35]TIDEK", ProFormaSequenceParser.Instance, out var proFormaApi, out var proFormaWarning);

            Assert.That(fromMzLib, Is.EqualTo(expected));
            Assert.That(fromProForma, Is.EqualTo(expected));
            Assert.That(mzLibApi, Is.EqualTo(expected));
            Assert.That(proFormaApi, Is.EqualTo(mzLibApi));
            Assert.That(mzLibWarning, Is.Not.Null);
            Assert.That(proFormaWarning, Is.Not.Null);
        }

        [Test]
        public void Hcd_AllowedIdOnUndeclaredResidueThrowMode_ThrowsFromMzLibAndProFormaSources()
        {
            var model = new HcdProbe(SequenceConversionHandlingMode.ThrowException);

            Assert.That(() => model.Clean("PEP[Unimod:Oxidation on P]TIDEK", null, out _, out _),
                Throws.ArgumentException.With.Message.Contains("unsupported modifications"));
            Assert.That(() => model.Clean("PEP[UNIMOD:35]TIDEK", ProFormaSequenceParser.Instance, out _, out _),
                Throws.ArgumentException.With.Message.Contains("unsupported modifications"));
        }

        [Test]
        public void Hcd_RemoveIncompatibleElements_KeepsTheSameIdOnItsDeclaredResidue()
        {
            var model = new HcdProbe(SequenceConversionHandlingMode.RemoveIncompatibleElements);

            var fromMzLib = model.Clean("PEP[Unimod:Oxidation on P]TM[Common Variable:Oxidation on M]IDEK", null, out var mzLibApi, out var mzLibWarning);
            var fromProForma = model.Clean("PEP[UNIMOD:35]TM[UNIMOD:35]IDEK", ProFormaSequenceParser.Instance, out var proFormaApi, out var proFormaWarning);

            Assert.That(fromMzLib, Is.EqualTo("PEPTM[Common Variable:Oxidation on M]IDEK"));
            Assert.That(fromProForma, Is.EqualTo("PEPTM[UNIMOD:35]IDEK"));
            Assert.That(mzLibApi, Is.EqualTo("PEPTM[UNIMOD:35]IDEK"));
            Assert.That(proFormaApi, Is.EqualTo(mzLibApi));
            Assert.That(mzLibWarning?.Message, Does.Contain("Oxidation on P"));
            Assert.That(proFormaWarning?.Message, Does.Contain("UNIMOD:35"));
        }

        // The TMT label is allowed and required at the N-terminus, and has a token on K but none on S.
        [TestCase(SequenceConversionHandlingMode.ReturnNull, null, null)]
        [TestCase(SequenceConversionHandlingMode.RemoveIncompatibleElements, "[Multiplex Label:TMT6-plex on X]PEPSIDEK[Multiplex Label:TMT6-plex on K]", "[UNIMOD:737]-PEPSIDEK[UNIMOD:737]")]
        public void Tmt_LabelOnUndeclaredResidue_LeavesTheRequiredNTerminalLabelAlone(SequenceConversionHandlingMode mode, string? expectedMzLib, string? expected)
        {
            var model = new TmtProbe(mode);

            var fromMzLib = model.Clean("[Multiplex Label:TMT6-plex on X]PEPS[Unimod:TMT6plex on S]IDEK[Multiplex Label:TMT6-plex on K]", null, out var mzLibApi, out _);
            var fromProForma = model.Clean("[UNIMOD:737]-PEPS[UNIMOD:737]IDEK[UNIMOD:737]", ProFormaSequenceParser.Instance, out var proFormaApi, out _);

            Assert.That(fromMzLib, Is.EqualTo(expectedMzLib));
            Assert.That(fromProForma, Is.EqualTo(expected));
            Assert.That(mzLibApi, Is.EqualTo(expected));
            Assert.That(proFormaApi, Is.EqualTo(mzLibApi));
        }

        private sealed class HcdProbe : Prosit2020IntensityHCD
        {
            public HcdProbe(SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ReturnNull)
                : base(modHandlingMode: mode) { }

            public string? Clean(string sequence, ISequenceParser? sourceParser, out string? api, out WarningException? warning)
            {
                var validated = TryCleanSequence(sequence, sourceParser, out var koinaSequence, out warning);
                api = koinaSequence is { } cleaned ? SerializeKoinaSequence(cleaned, out _) : null;
                return validated;
            }
        }

        private sealed class TmtProbe : Prosit2020IntensityTMT
        {
            public TmtProbe(SequenceConversionHandlingMode mode)
                : base(modHandlingMode: mode) { }

            public string? Clean(string sequence, ISequenceParser? sourceParser, out string? api, out WarningException? warning)
            {
                var validated = TryCleanSequence(sequence, sourceParser, out var koinaSequence, out warning);
                api = koinaSequence is { } cleaned ? SerializeKoinaSequence(cleaned, out _) : null;
                return validated;
            }
        }
    }

    /// <summary>
    /// Sends every declared token of every token-declaring model to Koina in one request per model, together with
    /// one peptide carrying an allowed id on a residue the model has no token for. Every declared token must come
    /// back predicted, and that last peptide unpredicted with mzLib's own unsupported-modification warning.
    /// </summary>
    [TestFixture]
    [Category("ExternalService")]
    [Category("Koina")]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class KoinaModificationTokenLiveTests : KoinaLiveTestFixture
    {
        [TestCaseSource(typeof(ModificationTokenCases), nameof(ModificationTokenCases.TokenDeclaringModels))]
        public void EveryDeclaredTokenIsPredicted_AndAnUndeclaredOneIsHeldBack(Type modelType)
        {
            var model = ModificationTokenCases.Instantiate(modelType);
            var tokens = ModificationTokenCases.Tokens(model)!.ToList();
            var undeclared = ModificationTokenCases.PeptideWith(model, ModificationTokenCases.UndeclaredToken(model));
            var proForma = ProFormaSequenceParser.Instance;

            List<string> sequences;
            List<(bool Predicted, WarningException? Warning)> outcomes;
            switch (model)
            {
                case FragmentIntensityModel fragmentModel:
                    sequences = tokens.Select(t => ModificationTokenCases.PeptideWith(model, t)).Append(undeclared).ToList();
                    outcomes = fragmentModel
                        .Predict(sequences.Select(s => new FragmentIntensityPredictionInput(s, 2, 30, "LUMOS", "HCD") { SequenceParser = proForma }).ToList())
                        .Select(p => (p.FragmentIntensities is { Count: > 0 }, p.Warning)).ToList();
                    break;
                case RetentionTimeModel retentionTimeModel:
                    sequences = tokens.Select(t => ModificationTokenCases.PeptideWith(model, t)).Append(undeclared).ToList();
                    outcomes = retentionTimeModel
                        .Predict(sequences.Select(s => new RetentionTimePredictionInput(s) { SequenceParser = proForma }).ToList())
                        .Select(p => (p.PredictedRetentionTime.HasValue, p.Warning)).ToList();
                    break;
                case CrosslinkFragmentIntensityModel crosslinkModel:
                    // A crosslinked peptide carries its crosslinker exactly once, so one peptide holds every token.
                    var crosslinked = ModificationTokenCases.PeptideWith(model, tokens.ToArray());
                    sequences = new List<string> { crosslinked, undeclared };
                    outcomes = crosslinkModel
                        .Predict(sequences.Select(s => new CrosslinkIntensityPredictionInput(s, crosslinked, 2, 30)).ToList())
                        .Select(p => (p.FragmentIntensities is { Count: > 0 }, p.Warning)).ToList();
                    break;
                default:
                    Assert.Fail($"{modelType.Name} declares modification tokens but this test cannot build its inputs.");
                    return;
            }

            Assert.That(outcomes, Has.Count.EqualTo(sequences.Count));
            Assert.Multiple(() =>
            {
                for (int i = 0; i < sequences.Count - 1; i++)
                    Assert.That(outcomes[i].Predicted, Is.True, $"{sequences[i]}: {outcomes[i].Warning?.Message}");

                Assert.That(outcomes[^1].Predicted, Is.False, undeclared);
                Assert.That(outcomes[^1].Warning?.Message, Does.Contain("unsupported modification"), undeclared);
            });
        }
    }
}
