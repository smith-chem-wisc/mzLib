using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.CCSModels;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;

namespace Test.KoinaTests
{
    [TestFixture]
    [Category("ExternalService")]
    [Category("Koina")]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class KoinaDashedNTerminalLiveTests : KoinaLiveTestFixture
    {
        private const string Bare = "PEPTIDEK";
        private const string Acetylated = "[Unimod:Acetyl on X]PEPTIDEK";
        private const string TmtLabelled = "[Multiplex Label:TMT6-plex on X]PEPTIDEK";
        private const double AcetylMass = 42.0106;

        public static IEnumerable<Type> RetentionTimeModels() => new[] { typeof(AlphaPeptDeepRTGeneric), typeof(ChronologerRT) };

        public static IEnumerable<Type> CcsModels() => new[] { typeof(AlphaPeptDeepCCSGeneric) };

        public static IEnumerable<TestCaseData> FragmentCases()
        {
            yield return new TestCaseData(typeof(AlphaPeptDeepMs2Generic), "QE");
            foreach (var t in new[] { typeof(Ms2PipHCD2021), typeof(Ms2PipImmunoHCD), typeof(Ms2PipTimsTOF2023), typeof(Ms2PipTimsTOF2024), typeof(Ms2PipTTOF5600) })
                yield return new TestCaseData(t, null);
        }

        private static T Instantiate<T>(Type t)
        {
            var ctor = t.GetConstructors().First(c => c.GetParameters().All(p => p.HasDefaultValue));
            return (T)ctor.Invoke(ctor.GetParameters().Select(p => p.DefaultValue).ToArray());
        }

        // Acetyl must differ from both the bare peptide and a second N-terminal id, so a model that reads the id as residues
        // or ignores it can't pass.
        private static void AssertUsesTheNTerminalModification(double? bare, double? acetylated, double? tmtLabelled, string? warning)
        {
            Assert.That(bare, Is.Not.Null);
            Assert.That(acetylated, Is.Not.Null, warning);
            Assert.That(tmtLabelled, Is.Not.Null, warning);
            Assert.That(acetylated, Is.Not.EqualTo(bare).Within(0.5));
            Assert.That(acetylated, Is.Not.EqualTo(tmtLabelled).Within(0.5));
        }

        [TestCaseSource(nameof(RetentionTimeModels))]
        public void RetentionTime_DashedNTerminalModification_IsUsed(Type modelType)
        {
            var rt = Instantiate<RetentionTimeModel>(modelType)
                .Predict(new List<RetentionTimePredictionInput> { new(Bare), new(Acetylated), new(TmtLabelled) });

            AssertUsesTheNTerminalModification(rt[0].PredictedRetentionTime, rt[1].PredictedRetentionTime, rt[2].PredictedRetentionTime,
                rt[1].Warning?.Message ?? rt[2].Warning?.Message);
        }

        [TestCaseSource(nameof(CcsModels))]
        public void Ccs_DashedNTerminalModification_IsUsed(Type modelType)
        {
            var ccs = Instantiate<CollisionalCrossSectionModel>(modelType)
                .Predict(new List<CCSPredictionInput> { new(Bare, 2), new(Acetylated, 2), new(TmtLabelled, 2) });

            AssertUsesTheNTerminalModification(ccs[0].PredictedCCS, ccs[1].PredictedCCS, ccs[2].PredictedCCS,
                ccs[1].Warning?.Message ?? ccs[2].Warning?.Message);
        }

        [TestCaseSource(nameof(FragmentCases))]
        public void Fragment_DashedNTerminalAcetyl_ShiftsBIons(Type modelType, string? instrument)
        {
            var spectra = Instantiate<FragmentIntensityModel>(modelType).Predict(new List<FragmentIntensityPredictionInput>
            {
                new(Bare, 2, 35, instrument, null), new(Acetylated, 2, 35, instrument, null)
            });

            Assert.That(spectra[1].FragmentMZs, Is.Not.Null, spectra[1].Warning?.Message);
            double B2(PeptideFragmentIntensityPrediction p) => p.FragmentMZs![p.FragmentAnnotations!.IndexOf("b2+1")];
            Assert.That(B2(spectra[1]) - B2(spectra[0]), Is.EqualTo(AcetylMass).Within(0.01));
        }
    }
}
