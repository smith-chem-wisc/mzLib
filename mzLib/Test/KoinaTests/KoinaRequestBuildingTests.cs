using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using Omics.Modifications;
using Omics.SequenceConversion;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.CrosslinkIntensityModels;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers.ProForma;

namespace Test.KoinaTests
{
    /// <summary>
    /// Drives each concrete model's request-building code (ToBatchedRequests and the
    /// model-specific TryCleanSequence/Validate overrides) without network access, so the Koina
    /// request-construction paths count toward coverage. ToBatchedRequests does no validation —
    /// it only reads the Koina payload — so one fully-populated input per family is enough.
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

        private static void AssertBuildsBatches(object model, object inputList, string? koinaSequence = null)
        {
            var method = model.GetType().GetMethod("ToBatchedRequests", BindingFlags.NonPublic | BindingFlags.Instance)!;
            var batches = (IEnumerable)method.Invoke(model, new[] { inputList })!;

            int count = 0;
            foreach (var batch in batches)
            {
                var dict = (IDictionary<string, object>)batch;
                Assert.That(dict.ContainsKey("id"), Is.True);
                Assert.That(dict["inputs"], Is.InstanceOf<IEnumerable>());
                if (koinaSequence != null)
                {
                    var sent = ((IEnumerable<object>)dict["inputs"])
                        .SelectMany(input => ((Array)input.GetType().GetProperty("data")!.GetValue(input)!).Cast<object>())
                        .ToList();
                    Assert.That(sent, Does.Contain(koinaSequence), $"{model.GetType().Name} must send the Koina payload");
                    Assert.That(sent, Does.Not.Contain(SourceFormatSequence), $"{model.GetType().Name} must not send ValidatedFullSequence");
                }
                count++;
            }
            Assert.That(count, Is.GreaterThanOrEqualTo(1), $"{model.GetType().Name}.ToBatchedRequests produced no batches");
        }

        private const string SourceFormatSequence = "PEPTIDEK[source format]";

        [TestCaseSource(nameof(FragmentModels))]
        public void FragmentIntensity_ToBatchedRequests_BuildsBatches(Type modelType)
        {
            var inputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 30, "QE", "HCD") { ValidatedFullSequence = SourceFormatSequence, KoinaSequence = "PEPTIDEK" }
            };
            AssertBuildsBatches(Instantiate(modelType), inputs, "PEPTIDEK");
        }

        [TestCaseSource(nameof(RtModels))]
        public void RetentionTime_ToBatchedRequests_BuildsBatches(Type modelType)
        {
            var inputs = new List<RetentionTimePredictionInput>
            {
                new("PEPTIDEK") { ValidatedFullSequence = SourceFormatSequence, KoinaSequence = "PEPTIDEK" }
            };
            AssertBuildsBatches(Instantiate(modelType), inputs, "PEPTIDEK");
        }

        [TestCaseSource(nameof(CcsModels))]
        public void Ccs_ToBatchedRequests_BuildsBatches(Type modelType)
        {
            var inputs = new List<CCSPredictionInput>
            {
                new("PEPTIDEK", 2) { ValidatedFullSequence = SourceFormatSequence, KoinaSequence = "PEPTIDEK" }
            };
            AssertBuildsBatches(Instantiate(modelType), inputs, "PEPTIDEK");
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
                new("PEPTIDEK") { ValidatedFullSequence = SourceFormatSequence, KoinaSequence = "PEPTIDEK" }
            };
            AssertBuildsBatches(Instantiate(modelType), inputs, "PEPTIDEK");
        }

        // ── Model-specific request-building overrides ───────────────────────────────

        [Test]
        public void Tmt_TryCleanSequence_RejectsSequenceWithoutNTerminalLabel(
            [Values(SequenceConversionHandlingMode.ReturnNull, SequenceConversionHandlingMode.RemoveIncompatibleElements)] SequenceConversionHandlingMode mode)
        {
            // Return value null = rejected; that's what the prediction pipeline keys off.
            var intensity = new TmtProbe(mode).Clean("PEPTIDEK", out _, out var intensityWarning);
            var irt = new IrtTmtProbe(mode).Clean("PEPTIDEK", out _, out var irtWarning);

            Assert.That(intensity, Is.Null);
            Assert.That(intensityWarning?.Message, Does.Contain("N-terminal"));
            Assert.That(irt, Is.Null);
            Assert.That(irtWarning?.Message, Does.Contain("N-terminal"));
        }

        [Test]
        public void Tmt_TryCleanSequence_ThrowExceptionModeThrowsWithoutNTerminalLabel()
        {
            const SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException;

            Assert.That(() => new TmtProbe(mode).Clean("PEPTIDEK", out _, out _),
                Throws.ArgumentException.With.Message.Contains("N-terminal"));
            Assert.That(() => new IrtTmtProbe(mode).Clean("PEPTIDEK", out _, out _),
                Throws.ArgumentException.With.Message.Contains("N-terminal"));
        }

        [Test]
        public void Tmt_IrtTryCleanSequence_AcceptsSupportedNTerminalLabel()
        {
            var result = new IrtTmtProbe().Clean("[Multiplex Label:TMT18 on X]PEPTIDEK", out var api, out var warning);

            Assert.That(result, Is.Not.Null);
            Assert.That(warning, Is.Null);
            Assert.That(api, Does.StartWith("[UNIMOD:2016]-"), "TMT18 on X should serialize to UNIMOD:2016.");
        }

        [Test]
        public void Tmt_Construction_AcceptsModesThatKeepTheLabel(
            [Values(SequenceConversionHandlingMode.ThrowException, SequenceConversionHandlingMode.ReturnNull, SequenceConversionHandlingMode.RemoveIncompatibleElements)] SequenceConversionHandlingMode mode)
        {
            Assert.That(new Prosit2020IntensityTMT(mode).ModHandlingMode, Is.EqualTo(mode));
            Assert.That(new Prosit2020iRTTMT(mode).ModHandlingMode, Is.EqualTo(mode));
        }

        [Test]
        public void Tmt_TryCleanSequence_AcceptsSupportedNTerminalLabel()
        {
            // Positive branch: a supported N-terminal TMT label must survive cleaning (offline).
            var model = new TmtProbe();
            var result = model.Clean("[Multiplex Label:TMT6-plex on X]PEPTIDEK", out var api, out var warning);

            Assert.That(result, Is.Not.Null, "A supported N-terminal TMT label must survive cleaning.");
            Assert.That(warning, Is.Null);
            Assert.That(api, Does.StartWith("[UNIMOD:737]-"), "TMT6-plex on X should serialize to UNIMOD:737.");
        }

        [Test]
        public void Tmt_TryCleanSequence_ProFormaSourceReachesSameApiSequenceAsMzLib()
        {
            var model = new TmtProbe();

            var mzLibResult = model.Clean("[Multiplex Label:TMT6-plex on X]PEPTIDEK", out var mzLibApi, out var mzLibWarning);
            var proFormaResult = model.CleanWithParser("[UNIMOD:737]-PEPTIDEK", ProFormaSequenceParser.Instance, out var proFormaApi, out var proFormaWarning);

            Assert.That(mzLibResult, Is.Not.Null);
            Assert.That(proFormaResult, Is.Not.Null, "A ProForma-sourced N-terminal TMT label must survive cleaning.");
            Assert.That(proFormaWarning, Is.Null);
            Assert.That(proFormaApi, Does.StartWith("[UNIMOD:737]-"));
            Assert.That(proFormaApi, Is.EqualTo(mzLibApi), "mzLib and ProForma sources of the same peptide must serialize byte-identically.");
        }

        [Test]
        public void Tmt_TryCleanSequence_ProFormaSourceWithoutNTerminalLabel_IsRejected()
        {
            var model = new TmtProbe();

            var result = model.CleanWithParser("PEPTIDEK", ProFormaSequenceParser.Instance, out var api, out var warning);

            Assert.That(result, Is.Null);
            Assert.That(api, Is.Null);
            Assert.That(warning, Is.Not.Null);
            Assert.That(warning!.Message, Does.Contain("N-terminal"));
        }

        [Test]
        public void IrtTmt_TryCleanSequence_RequiresNTerminalLabelFromEitherSource()
        {
            var model = new IrtTmtProbe();

            var mzLib = model.Clean("[Multiplex Label:TMT18 on X]PEPTIDEK", out var mzLibApi, out _);
            var proForma = model.CleanWithParser("[UNIMOD:2016]-PEPTIDEK", ProFormaSequenceParser.Instance, out var proFormaApi, out _);
            var unlabeled = model.CleanWithParser("PEPTIDEK", ProFormaSequenceParser.Instance, out _, out var unlabeledWarning);

            Assert.That(mzLib, Is.Not.Null);
            Assert.That(proForma, Is.Not.Null);
            Assert.That(proFormaApi, Is.EqualTo(mzLibApi).And.StartWith("[UNIMOD:2016]-"));
            Assert.That(unlabeled, Is.Null);
            Assert.That(unlabeledWarning?.Message, Does.Contain("N-terminal"));
        }

        [Test]
        public void EveryModel_RequiredNTerminalModifications_AreAllowedByThatModel()
        {
            var models = FragmentModels().Concat(RtModels()).Concat(CcsModels()).Concat(CrosslinkModels()).Concat(DetectabilityModels());
            var checkedModels = 0;
            Assert.Multiple(() =>
            {
                foreach (var modelType in models)
                {
                    var model = Instantiate(modelType);
                    var required = (IReadOnlySet<int>?)modelType.GetProperty("RequiredNTerminalUnimodIds")!.GetValue(model);
                    var acceptsAll = (bool)modelType.GetProperty("AcceptsAllUnimodModifications")!.GetValue(model)!;
                    if (required == null || acceptsAll)
                        continue;

                    var allowed = (IReadOnlySet<int>)modelType.GetProperty("AllowedUnimodIds")!.GetValue(model)!;
                    Assert.That(required, Is.SubsetOf(allowed), modelType.Name);
                    checkedModels++;
                }
            });
            Assert.That(checkedModels, Is.GreaterThanOrEqualTo(2), "Expected at least the two Prosit 2020 TMT models to declare required N-terminal mods.");
        }

        [Test]
        public void EveryModel_RequiringAModification_RejectsUsePrimarySequenceAtConstruction()
        {
            var models = FragmentModels().Concat(RtModels()).Concat(CcsModels()).Concat(CrosslinkModels()).Concat(DetectabilityModels());
            var checkedModels = 0;
            Assert.Multiple(() =>
            {
                foreach (var modelType in models)
                {
                    var model = Instantiate(modelType);
                    if (modelType.GetProperty("RequiredNTerminalUnimodIds")!.GetValue(model) == null)
                        continue;

                    var ctor = modelType.GetConstructors().First(c => c.GetParameters().All(p => p.HasDefaultValue)
                        && c.GetParameters().Any(p => p.ParameterType == typeof(SequenceConversionHandlingMode)));
                    var args = ctor.GetParameters()
                        .Select(p => p.ParameterType == typeof(SequenceConversionHandlingMode) ? SequenceConversionHandlingMode.UsePrimarySequence : p.DefaultValue)
                        .ToArray();
                    Assert.That(() => ctor.Invoke(args),
                        Throws.TypeOf<TargetInvocationException>().With.InnerException.TypeOf<ArgumentException>(), modelType.Name);
                    Assert.That(() => modelType.GetProperty("ModHandlingMode")!.SetValue(Instantiate(modelType), SequenceConversionHandlingMode.UsePrimarySequence),
                        Throws.TypeOf<TargetInvocationException>().With.InnerException.TypeOf<ArgumentException>(), modelType.Name + " (init)");
                    checkedModels++;
                }
            });
            Assert.That(checkedModels, Is.GreaterThanOrEqualTo(2), "Expected at least the two Prosit 2020 TMT models to declare required N-terminal mods.");
        }

        [Test]
        public void EveryModel_ModificationItsConverterResolves_IsAllowedByThatModel()
        {
            var oxidation = CanonicalModification.AtResidue(3, 'M', "Common Variable:Oxidation on M", mzLibId: "Common Variable:Oxidation on M");
            var models = FragmentModels().Concat(RtModels()).Concat(CcsModels()).Concat(CrosslinkModels()).Concat(DetectabilityModels());
            var checkedAcceptAllModels = 0;
            Assert.Multiple(() =>
            {
                foreach (var modelType in models)
                {
                    var model = Instantiate(modelType);
                    var converter = (ISequenceConverter)modelType.GetProperty("SequenceConverter", BindingFlags.NonPublic | BindingFlags.Instance)!.GetValue(model)!;
                    if (converter.Serializer.ModificationLookup?.TryResolve(oxidation)?.UnimodId is not int id)
                        continue;

                    var acceptsAll = (bool)modelType.GetProperty("AcceptsAllUnimodModifications")!.GetValue(model)!;
                    var allowed = (IReadOnlySet<int>)modelType.GetProperty("AllowedUnimodIds")!.GetValue(model)!;
                    Assert.That(acceptsAll || allowed.Contains(id), Is.True,
                        $"{modelType.Name}'s converter resolves UNIMOD:{id}, which the model doesn't allow.");
                    if (acceptsAll)
                        checkedAcceptAllModels++;
                }
            });
            Assert.That(checkedAcceptAllModels, Is.GreaterThanOrEqualTo(16), "Expected all 16 accept-all models to resolve oxidation and be checked.");
        }

        [Test]
        public void Tmt_TryCleanSequence_ProFormaSourceWithOutOfSetUnimodId_ReturnNullRejectsBeforeSerialization()
        {
            var model = new TmtProbe();

            var result = model.CleanWithParser("[UNIMOD:737]-PEPS[UNIMOD:21]IDEK", ProFormaSequenceParser.Instance, out var api, out var warning);

            Assert.That(result, Is.Null);
            Assert.That(api, Is.Null);
            Assert.That(warning, Is.Not.Null);
            Assert.That(warning!.Message, Does.Contain("UNIMOD:21"));
        }

        [TestCase(SequenceConversionHandlingMode.ReturnNull, null)]
        [TestCase(SequenceConversionHandlingMode.RemoveIncompatibleElements, "[UNIMOD:737]-PEPSIDEK")]
        public void Tmt_TryCleanSequence_OutOfSetModification_SameOutcomeFromMzLibAndProFormaSources(SequenceConversionHandlingMode mode, string? expected)
        {
            var model = new TmtProbe(mode);

            var fromMzLib = model.Clean("[Multiplex Label:TMT6-plex on X]PEPS[Common Biological:Phosphorylation on S]IDEK", out var mzLibApi, out _);
            var fromProForma = model.CleanWithParser("[UNIMOD:737]-PEPS[UNIMOD:21]IDEK", ProFormaSequenceParser.Instance, out var proFormaApi, out var proFormaWarning);

            Assert.That(fromMzLib is null, Is.EqualTo(expected is null));
            Assert.That(fromProForma is null, Is.EqualTo(expected is null));
            Assert.That(mzLibApi, Is.EqualTo(expected));
            Assert.That(proFormaApi, Is.EqualTo(expected));
            Assert.That(proFormaWarning?.Message, Does.Contain("UNIMOD:21"));
        }

        [Test]
        public void Tmt_TryCleanSequence_OutOfSetModificationThrowMode_ThrowsFromMzLibAndProFormaSources()
        {
            var model = new TmtProbe(SequenceConversionHandlingMode.ThrowException);

            Assert.Throws<ArgumentException>(() => model.Clean("[Multiplex Label:TMT6-plex on X]PEPS[Common Biological:Phosphorylation on S]IDEK", out _, out _));
            Assert.Throws<ArgumentException>(() => model.CleanWithParser("[UNIMOD:737]-PEPS[UNIMOD:21]IDEK", ProFormaSequenceParser.Instance, out _, out _));
        }

        [Test]
        public void Tmt_ToBatchedRequests_SendsFragmentationType([Values("HCD", "CID")] string fragType)
        {
            // Both supported fragmentation types must flow through to the Koina request offline.
            var model = new TmtProbe();
            var inputs = new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, 30, null, fragType) { KoinaSequence = "PEPTIDEK" }
            };

            var batches = model.Build(inputs);

            Assert.That(batches, Has.Count.GreaterThanOrEqualTo(1));
            var fragData = DataForInput(batches[0], "fragmentation_types");
            Assert.That(fragData, Is.Not.Null, "TMT request must include a fragmentation_types input.");
            Assert.That(fragData!.Cast<string>(), Does.Contain(fragType));
        }

        [Test]
        public void XlNms2_TryCleanSequence_AcceptsModelSpecificCrosslinkerUnimod()
        {
            // XLNMS2 accepts only {4, 35, 1898}; the 1898 crosslinker must pass its real validation.
            var model = new XlNms2Probe();
            var result = model.Clean("PEPTIDEK[UNIMOD:1898]", out var warning);

            Assert.That(result, Is.Not.Null);
            Assert.That(warning, Is.Null);
        }

        [Test]
        public void XlNms2_TryCleanSequence_RejectsForeignCrosslinkerUnimod()
        {
            // 1896 is the CMS2 crosslinker, not accepted by XLNMS2 ({4, 35, 1898}).
            // The shared smoke test bypasses this by pre-seeding Validated* fields, so assert it here.
            var model = new XlNms2Probe();
            var result = model.Clean("PEPTIDEK[UNIMOD:1896]", out var warning);

            Assert.That(result, Is.Null);
            Assert.That(warning, Is.Not.Null);
            Assert.That(warning!.Message, Does.Contain("1896"));
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

        private static Modification CatalogMod(string type, string id) =>
            Mods.AllProteinModsList.Single(m => m.ModificationType == type && m.IdWithMotif == id);

        private static string Digested(string peptide)
        {
            var oxidation = CatalogMod("Common Variable", "Oxidation on M");
            var phospho = CatalogMod("Common Biological", "Phosphorylation on S");
            var tmt = CatalogMod("Multiplex Label", "TMT6-plex on X");
            var tmtK = CatalogMod("Multiplex Label", "TMT6-plex on K");
            var localized = new Dictionary<int, List<Modification>>
            {
                [1] = new() { CatalogMod("UniProt", "N-acetylalanine on A") },
                [12] = new() { CatalogMod("UniProt", "Arginine amide on R") }
            };
            var digestionParams = new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, minPeptideLength: 1, maxModsForPeptides: 2);
            var peptides = new Protein("PEPMSIDEK", "P").Digest(digestionParams, new List<Modification>(), new List<Modification> { oxidation, phospho })
                .Concat(new Protein("PEPMSIDEK", "P").Digest(digestionParams, new List<Modification> { tmt, tmtK }, new List<Modification> { oxidation, phospho }))
                .Concat(new Protein("APEPTIDEKAAR", "P", oneBasedModifications: localized).Digest(digestionParams, new List<Modification>(), new List<Modification>()))
                .Select(p => p.FullSequence);
            return peptides.Single(s => s == peptide);
        }

        private const string OxidizedPhosphopeptide = "PEPM[Common Variable:Oxidation on M]S[Common Biological:Phosphorylation on S]IDEK";
        private const string AcetylatedNTerminus = "[UniProt:N-acetylalanine on A]APEPTIDEK";
        private const string AmidatedCTerminus = "AAR-[UniProt:Arginine amide on R]";
        private const string TmtLabeled = "[Multiplex Label:TMT6-plex on X]PEPMSIDEK[Multiplex Label:TMT6-plex on K]";
        private const string TmtLabeledOxidizedPhosphopeptide =
            "[Multiplex Label:TMT6-plex on X]PEPM[Common Variable:Oxidation on M]S[Common Biological:Phosphorylation on S]IDEK[Multiplex Label:TMT6-plex on K]";

        [TestCase(OxidizedPhosphopeptide, "PEPM[UNIMOD:35]S[UNIMOD:21]IDEK")]
        [TestCase(AcetylatedNTerminus, "[UNIMOD:1]APEPTIDEK")]
        [TestCase(AmidatedCTerminus, "AAR-[UNIMOD:2]")]
        [TestCase(TmtLabeled, "[UNIMOD:737]PEPMSIDEK[UNIMOD:737]")]
        public void AcceptAllModel_MetaMorpheusFullSequence_SendsTheCatalogUnimodIds(string fullSequence, string expected)
        {
            var result = new Ms2PipProbe(SequenceConversionHandlingMode.ThrowException).Clean(Digested(fullSequence), out var api, out var warning);

            Assert.That(result, Is.EqualTo(fullSequence));
            Assert.That(api, Is.EqualTo(expected));
            Assert.That(warning, Is.Null);
        }

        [Test]
        public void AcceptAllModel_UniProtModificationCitingAUnimodRecordOfAnotherMass_IsNotSentAsThatId()
        {
            // UniProt's N,N-dimethylproline (+28.031) cites UNIMOD:529 (+29.039).
            var dimethylproline = CatalogMod("UniProt", "N,N-dimethylproline on P");
            var record529 = Mods.UnimodModifications.First(m => m.ModificationType == "Unimod" && CanonicalModification.GetUnimodId(m) == 529);
            Assert.That(CanonicalModification.GetUnimodId(dimethylproline) == 529
                        && Math.Abs(record529.MonoisotopicMass!.Value - dimethylproline.MonoisotopicMass!.Value) > 0.01,
                "N,N-dimethylproline no longer cites a UNIMOD record of another mass, which is what this test needs.");
            var localized = new Dictionary<int, List<Modification>> { [1] = new() { dimethylproline } };
            var fullSequence = new Protein("PEPTIDEK", "P", oneBasedModifications: localized)
                .Digest(new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, minPeptideLength: 1), new List<Modification>(), new List<Modification>())
                .Single(p => p.AllModsOneIsNterminus.Count == 1).FullSequence;

            var result = new Ms2PipProbe(SequenceConversionHandlingMode.ReturnNull).Clean(fullSequence, out _, out var warning);

            Assert.That(fullSequence, Is.EqualTo("[UniProt:N,N-dimethylproline on P]PEPTIDEK"));
            Assert.That(result, Is.Null);
            Assert.That(warning?.Message, Does.Contain("UniProt:N,N-dimethylproline on P"));
        }

        [TestCase(SequenceConversionHandlingMode.ReturnNull, null, null)]
        [TestCase(SequenceConversionHandlingMode.RemoveIncompatibleElements, "PEPM[UNIMOD:35]SIDEK", "PEPM[Common Variable:Oxidation on M]SIDEK")]
        [TestCase(SequenceConversionHandlingMode.UsePrimarySequence, "PEPMSIDEK", "PEPMSIDEK")]
        public void AllowListModel_MetaMorpheusPhosphopeptide_KeepsOxidationAndHandlesPhosphoPerMode(SequenceConversionHandlingMode mode, string? expected, string? expectedValidated)
        {
            var result = new PrositHcdProbe(mode).Clean(Digested(OxidizedPhosphopeptide), out var api, out var warning);

            Assert.That(result, Is.EqualTo(expectedValidated));
            Assert.That(api, Is.EqualTo(expected));
            if (mode != SequenceConversionHandlingMode.UsePrimarySequence)
                Assert.That(warning?.Message, Does.Contain("Common Biological:Phosphorylation on S"));
        }

        [TestCase(OxidizedPhosphopeptide)]
        [TestCase(AcetylatedNTerminus)]
        [TestCase(AmidatedCTerminus)]
        [TestCase(TmtLabeled)]
        public void AllowListModel_MetaMorpheusModificationItDoesNotAllow_ThrowsInThrowMode(string fullSequence)
        {
            Assert.That(() => new PrositHcdProbe(SequenceConversionHandlingMode.ThrowException).Clean(Digested(fullSequence), out _, out _),
                Throws.ArgumentException.With.Message.Contains("unsupported modifications"));
        }

        [TestCase(SequenceConversionHandlingMode.ThrowException, "[UNIMOD:737]-PEPMSIDEK[UNIMOD:737]")]
        [TestCase(SequenceConversionHandlingMode.ReturnNull, "[UNIMOD:737]-PEPMSIDEK[UNIMOD:737]")]
        public void TmtModel_MetaMorpheusTmtLabels_AreTheRequiredNTerminalLabel(SequenceConversionHandlingMode mode, string expected)
        {
            var result = new TmtProbe(mode).Clean(Digested(TmtLabeled), out var api, out var warning);

            Assert.That(result, Is.EqualTo(TmtLabeled));
            Assert.That(api, Is.EqualTo(expected));
            Assert.That(warning, Is.Null);
        }

        [TestCase(SequenceConversionHandlingMode.ReturnNull, null)]
        [TestCase(SequenceConversionHandlingMode.RemoveIncompatibleElements, "[UNIMOD:737]-PEPM[UNIMOD:35]SIDEK[UNIMOD:737]")]
        public void TmtModel_MetaMorpheusTmtPhosphopeptide_HandlesPhosphoPerMode(SequenceConversionHandlingMode mode, string? expected)
        {
            var result = new TmtProbe(mode).Clean(Digested(TmtLabeledOxidizedPhosphopeptide), out var api, out var warning);

            Assert.That(result is null, Is.EqualTo(expected is null));
            Assert.That(api, Is.EqualTo(expected));
            Assert.That(warning?.Message, Does.Contain("Common Biological:Phosphorylation on S"));
        }

        private sealed class Ms2PipProbe : Ms2PipHCD2021
        {
            public Ms2PipProbe(SequenceConversionHandlingMode mode) : base(modHandlingMode: mode) { }

            public string? Clean(string sequence, out string? api, out WarningException? warning)
            {
                var validated = TryCleanSequence(sequence, null, out var koinaSequence, out warning);
                api = koinaSequence is { } cleaned ? SerializeKoinaSequence(cleaned, out _) : null;
                return validated;
            }
        }

        private sealed class PrositHcdProbe : Prosit2020IntensityHCD
        {
            public PrositHcdProbe(SequenceConversionHandlingMode mode) : base(modHandlingMode: mode) { }

            public string? Clean(string sequence, out string? api, out WarningException? warning)
            {
                var validated = TryCleanSequence(sequence, null, out var koinaSequence, out warning);
                api = koinaSequence is { } cleaned ? SerializeKoinaSequence(cleaned, out _) : null;
                return validated;
            }
        }

        private sealed class TmtProbe : Prosit2020IntensityTMT
        {
            public TmtProbe(SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ReturnNull)
                : base(modHandlingMode: mode) { }

            public string? Clean(string sequence, out string? api, out WarningException? warning)
                => CleanWithParser(sequence, null, out api, out warning);

            public string? CleanWithParser(string sequence, ISequenceParser? sourceParser, out string? api, out WarningException? warning)
            {
                var validated = TryCleanSequence(sequence, sourceParser, out var koinaSequence, out warning);
                api = koinaSequence is { } cleaned ? SerializeKoinaSequence(cleaned, out _) : null;
                return validated;
            }

            public List<Dictionary<string, object>> Build(List<FragmentIntensityPredictionInput> inputs)
                => ToBatchedRequests(inputs);
        }

        private sealed class IrtTmtProbe : Prosit2020iRTTMT
        {
            public IrtTmtProbe(SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ReturnNull)
                : base(modHandlingMode: mode) { }

            public string? Clean(string sequence, out string? api, out WarningException? warning)
                => CleanWithParser(sequence, null, out api, out warning);

            public string? CleanWithParser(string sequence, ISequenceParser? sourceParser, out string? api, out WarningException? warning)
            {
                var validated = TryCleanSequence(sequence, sourceParser, out var koinaSequence, out warning);
                api = koinaSequence is { } cleaned ? SerializeKoinaSequence(cleaned, out _) : null;
                return validated;
            }
        }

        /// <summary>
        /// Reads the flat data array for a named input field out of a built batch request.
        /// BuildBatchedRequest stores each input as an anonymous { name, shape, datatype, data }
        /// object, so this reflects over those properties to find the requested field.
        /// </summary>
        private static Array? DataForInput(Dictionary<string, object> batch, string inputName)
        {
            foreach (var input in (IEnumerable<object>)batch["inputs"])
            {
                var t = input.GetType();
                var name = (string)t.GetProperty("name")!.GetValue(input)!;
                if (name == inputName)
                    return (Array)t.GetProperty("data")!.GetValue(input)!;
            }
            return null;
        }

        private sealed class XlNms2Probe : Prosit2024IntensityXLNMS2
        {
            public string? Clean(string sequence, out WarningException? warning)
                => TryCleanSequence(sequence, null, out _, out warning);
        }
    }
}
