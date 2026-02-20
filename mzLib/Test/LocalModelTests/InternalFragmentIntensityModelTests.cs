using NUnit.Framework;
using Omics.Fragmentation;
using PredictionClients.LocalModels;
using Readers.SpectralLibrary;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using CategoryAttribute = NUnit.Framework.CategoryAttribute;

namespace Test.LocalModelTests
{
    /// <summary>
    /// Tests for InternalFragmentIntensityModel.
    ///
    /// TEST ORGANISATION
    /// -----------------
    /// 1. Constructor / input validation   — matches Prosit test patterns exactly
    /// 2. Inference (ONNX)                — requires the .onnx file on disk
    /// 3. Spectral library generation     — in-memory and file round-trips
    /// 4. Utility methods                 — BareSequence, feature vector
    ///
    /// Tests that require the ONNX model file are decorated with [Category("RequiresOnnxModel")]
    /// so they can be excluded in CI environments where the model is not present.
    /// All other tests are pure unit tests that run without any files.
    ///
    /// The ONNX model is resolved via InternalFragmentIntensityModel.DefaultOnnxModelPath,
    /// which checks embedded resources, then Resources folder, then LocalModels folder.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class InternalFragmentIntensityModelTests
    {
        // ── Shared valid inputs ──────────────────────────────────────────────────
        private static readonly List<string> ValidPeptides = new() { "PEPTIDEK", "ELVISLIVESK", "SAMPLER" };
        private static readonly List<int> ValidCharges = new() { 2, 2, 2 };
        private static readonly List<double?> ValidRetentionTimes = new() { 100.0, 200.0, 300.0 };

        // ════════════════════════════════════════════════════════════════════════
        // 1. CONSTRUCTOR / INPUT VALIDATION
        // ════════════════════════════════════════════════════════════════════════

        [Test]
        public static void Constructor_ValidInputs_NoWarning()
        {
            var model = new InternalFragmentIntensityModel(
                ValidPeptides, ValidCharges, ValidRetentionTimes,
                out var warnings);

            Assert.That(warnings, Is.Null,
                "Valid inputs should produce no warning");
            Assert.That(model.PeptideSequences.Count, Is.EqualTo(3),
                "All three valid peptides should be accepted");
        }

        [Test]
        public static void Constructor_MismatchedListLengths_ThrowsArgumentException()
        {
            var peptides = new List<string> { "PEPTIDEK", "ELVISLIVESK" };
            var charges = new List<int> { 2 };          // wrong length
            var rts = new List<double?> { 100.0, null };

            Assert.Throws<ArgumentException>(() =>
                new InternalFragmentIntensityModel(
                    peptides, charges, rts,
                    out _));
        }

        [Test]
        public static void Constructor_EmptyInput_ReturnsWarning()
        {
            var model = new InternalFragmentIntensityModel(
                new List<string>(), new List<int>(), new List<double?>(),
                out var warnings);

            Assert.That(warnings, Is.Not.Null,
                "Empty input should produce a warning");
            Assert.That(model.PeptideSequences.Count, Is.EqualTo(0));
        }

        [Test]
        public static void Constructor_AllValidCharges_NoneFiltered()
        {
            foreach (int charge in new[] { 1, 2, 3, 4, 5, 6 })
            {
                var model = new InternalFragmentIntensityModel(
                    new List<string> { "PEPTIDEK" },
                    new List<int> { charge },
                    new List<double?> { null },
                    out var warnings);

                Assert.That(model.PeptideSequences.Count, Is.EqualTo(1),
                    $"Charge {charge} should be accepted");
                Assert.That(warnings, Is.Null,
                    $"Charge {charge} should not produce a warning");
            }
        }

        [Test]
        public static void Constructor_InvalidCharge_Filtered_WithWarning()
        {
            foreach (int badCharge in new[] { 0, -1, 7 })
            {
                var model = new InternalFragmentIntensityModel(
                    new List<string> { "PEPTIDEK" },
                    new List<int> { badCharge },
                    new List<double?> { null },
                    out var warnings);

                Assert.That(model.PeptideSequences.Count, Is.EqualTo(0),
                    $"Charge {badCharge} should be filtered out");
                Assert.That(warnings, Is.Not.Null,
                    $"Charge {badCharge} should produce a warning");
            }
        }

        [Test]
        public static void Constructor_PeptideTooShort_Filtered()
        {
            // MinPeptideLength = 5; shorter sequences cannot yield any internal fragment of length 3
            var model = new InternalFragmentIntensityModel(
                new List<string> { "AAAA" },
                new List<int> { 2 },
                new List<double?> { null },
                out var warnings);

            Assert.That(model.PeptideSequences.Count, Is.EqualTo(0),
                "Peptide shorter than MinPeptideLength should be filtered");
            Assert.That(warnings, Is.Not.Null);
        }

        [Test]
        public static void Constructor_PeptideAtMinLength_Accepted()
        {
            // MinPeptideLength = 5 — "AAAAA" is exactly at the boundary
            var model = new InternalFragmentIntensityModel(
                new List<string> { "AAAAA" },
                new List<int> { 2 },
                new List<double?> { null },
                out var warnings);

            Assert.That(model.PeptideSequences.Count, Is.EqualTo(1));
            Assert.That(warnings, Is.Null);
        }

        [Test]
        public static void Constructor_InvalidAminoAcid_Filtered()
        {
            // 'X' (unknown AA), 'B' (Asx), 'Z' (Glx) are not in AaMass
            foreach (var badSeq in new[] { "PEPTXIDEK", "PEPTBIDEK", "PEPTZIDEK" })
            {
                var model = new InternalFragmentIntensityModel(
                    new List<string> { badSeq },
                    new List<int> { 2 },
                    new List<double?> { null },
                    out var warnings);

                Assert.That(model.PeptideSequences.Count, Is.EqualTo(0),
                    $"'{badSeq}' with non-canonical AA should be filtered");
                Assert.That(warnings, Is.Not.Null);
            }
        }

        [Test]
        public static void Constructor_MixedValidAndInvalid_OnlyValidKept()
        {
            var peptides = new List<string> { "PEPTIDEK", "PEPTXIDEK", "SAMPLER" };
            var charges = new List<int> { 2, 2, 2 };
            var rts = new List<double?> { null, null, null };

            var model = new InternalFragmentIntensityModel(
                peptides, charges, rts,
                out var warnings);

            Assert.That(model.PeptideSequences.Count, Is.EqualTo(2),
                "Only the two valid peptides should be retained");
            Assert.That(warnings, Is.Not.Null,
                "Filtered entry should produce a warning");
        }

        [Test]
        public static void Constructor_ModifiedSequences_BareSequenceUsedForValidation()
        {
            // mzLib modification brackets — stripped before length / AA validation
            var modifiedSeq = "PEPTM[Common Variable:Oxidation on M]IDER";
            var model = new InternalFragmentIntensityModel(
                new List<string> { modifiedSeq },
                new List<int> { 2 },
                new List<double?> { null },
                out var warnings);

            // Bare sequence = "PEPTMIDER" (9 AA) — valid
            Assert.That(model.PeptideSequences.Count, Is.EqualTo(1),
                "Modified sequence with valid bare AA sequence should be accepted");
            Assert.That(warnings, Is.Null);
        }

        [Test]
        public static void Constructor_NullRetentionTimes_Accepted()
        {
            var rts = new List<double?> { null, null, null };
            var model = new InternalFragmentIntensityModel(
                ValidPeptides, ValidCharges, rts,
                out var warnings);

            Assert.That(model.PeptideSequences.Count, Is.EqualTo(3));
            Assert.That(warnings, Is.Null,
                "Null retention times should not cause filtering");
        }

        [Test]
        public static void DefaultOnnxModelPath_ReturnsValidPath()
        {
            var defaultPath = InternalFragmentIntensityModel.DefaultOnnxModelPath;

            Assert.That(defaultPath, Is.Not.Null.And.Not.Empty,
                "DefaultOnnxModelPath should return a non-empty path");
            Assert.That(defaultPath, Does.EndWith(InternalFragmentIntensityModel.DefaultModelFileName),
                "Default path should end with the expected model filename");
        }

        // ════════════════════════════════════════════════════════════════════════
        // 2. UTILITY METHODS
        // ════════════════════════════════════════════════════════════════════════

        [Test]
        public static void BareSequence_NoMods_ReturnsIdentical()
        {
            Assert.That(InternalFragmentIntensityModel.BareSequence("PEPTIDEK"),
                Is.EqualTo("PEPTIDEK"));
        }

        [Test]
        public static void BareSequence_WithMods_StripsBrackets()
        {
            var full = "PEPTM[Common Variable:Oxidation on M]IDER";
            Assert.That(InternalFragmentIntensityModel.BareSequence(full),
                Is.EqualTo("PEPTMIDER"));
        }

        [Test]
        public static void BareSequence_MultipleMods_AllStripped()
        {
            var full = "M[Common Variable:Oxidation on M]PEPTC[Common Fixed:Carbamidomethyl on C]IDE";
            Assert.That(InternalFragmentIntensityModel.BareSequence(full),
                Is.EqualTo("MPEPTCIDE"));
        }

        [Test]
        public static void BareSequence_NTerminalMod_Stripped()
        {
            var full = "[Common Biological:Acetylation on X]SAMPLER";
            Assert.That(InternalFragmentIntensityModel.BareSequence(full),
                Is.EqualTo("SAMPLER"));
        }

        // ════════════════════════════════════════════════════════════════════════
        // 3. ONNX INFERENCE  (requires model file — tagged for selective CI runs)
        // ════════════════════════════════════════════════════════════════════════

        [Test, Category("RequiresOnnxModel")]
        public static async Task RunInferenceAsync_ValidPeptides_PopulatesPredictions()
        {
            var model = new InternalFragmentIntensityModel(
                ValidPeptides, ValidCharges, ValidRetentionTimes,
                out _);

            await model.RunInferenceAsync();

            Assert.That(model.Predictions.Count, Is.EqualTo(3),
                "Should have one prediction per peptide");

            foreach (var pred in model.Predictions)
            {
                Assert.That(pred.FragmentAnnotations.Count, Is.GreaterThan(0),
                    "Each peptide should have at least one candidate fragment");
                Assert.That(pred.FragmentAnnotations.Count, Is.EqualTo(pred.FragmentMZs.Count));
                Assert.That(pred.FragmentAnnotations.Count, Is.EqualTo(pred.FragmentIntensities.Count));
            }
        }

        [Test, Category("RequiresOnnxModel")]
        public static async Task RunInferenceAsync_PredictionValues_AreChemicallyReasonable()
        {
            var model = new InternalFragmentIntensityModel(
                ValidPeptides, ValidCharges, ValidRetentionTimes,
                out _);

            await model.RunInferenceAsync();

            foreach (var pred in model.Predictions)
            {
                foreach (var mz in pred.FragmentMZs)
                    Assert.That(mz, Is.GreaterThan(0).And.LessThan(3000),
                        "Fragment m/z should be positive and within peptide mass range");

                foreach (var intensity in pred.FragmentIntensities)
                    Assert.That(intensity, Is.GreaterThanOrEqualTo(0),
                        "Predicted intensities should be non-negative");
            }
        }

        [Test, Category("RequiresOnnxModel")]
        public static async Task RunInferenceAsync_AnnotationsMatchMzLibFormat()
        {
            // All annotation strings must follow "bIb[start-end]+1"
            var annotRegex = new System.Text.RegularExpressions.Regex(
                @"^[a-zA-Z]+I[a-zA-Z]+\[\d+-\d+\]\+\d+$");

            var model = new InternalFragmentIntensityModel(
                ValidPeptides, ValidCharges, ValidRetentionTimes,
                out _);

            await model.RunInferenceAsync();

            foreach (var pred in model.Predictions)
                foreach (var ann in pred.FragmentAnnotations)
                    Assert.That(annotRegex.IsMatch(ann), Is.True,
                        $"Annotation '{ann}' does not match expected bIb[start-end]+1 format");
        }

        [Test, Category("RequiresOnnxModel")]
        public static async Task RunInferenceAsync_EmptyInput_ReturnsNoError()
        {
            var model = new InternalFragmentIntensityModel(
                new List<string>(), new List<int>(), new List<double?>(),
                out _);

            Assert.DoesNotThrowAsync(async () => await model.RunInferenceAsync());
            Assert.That(model.PredictedSpectra.Count, Is.EqualTo(0));
        }

        [Test, Category("RequiresOnnxModel")]
        public static async Task RunInferenceAsync_DisposesAfterInference()
        {
            var model = new InternalFragmentIntensityModel(
                ValidPeptides, ValidCharges, ValidRetentionTimes,
                out _);

            await model.RunInferenceAsync();

            // Second call on disposed model should throw
            Assert.ThrowsAsync<ObjectDisposedException>(
                async () => await model.RunInferenceAsync(),
                "Second RunInferenceAsync call should throw ObjectDisposedException");

            // Results remain accessible after disposal
            Assert.That(model.PredictedSpectra.Count, Is.GreaterThan(0),
                "PredictedSpectra should remain accessible after disposal");
        }

        [Test, Category("RequiresOnnxModel")]
        public static void RunInferenceAsync_DisposeBeforeInference_ThrowsObjectDisposedException()
        {
            var model = new InternalFragmentIntensityModel(
                ValidPeptides, ValidCharges, ValidRetentionTimes,
                out _);

            model.Dispose();

            Assert.ThrowsAsync<ObjectDisposedException>(
                async () => await model.RunInferenceAsync());
        }

        // ════════════════════════════════════════════════════════════════════════
        // 4. SPECTRAL LIBRARY GENERATION
        // ════════════════════════════════════════════════════════════════════════

        [Test, Category("RequiresOnnxModel")]
        public static async Task PredictedSpectra_InternalIonsHaveIsInternalFragmentTrue()
        {
            var model = new InternalFragmentIntensityModel(
                ValidPeptides, ValidCharges, ValidRetentionTimes,
                out _);

            await model.RunInferenceAsync();

            Assert.That(model.PredictedSpectra.Count, Is.GreaterThan(0),
                "At least one peptide should produce passing internal ions");

            foreach (var spectrum in model.PredictedSpectra)
            {
                foreach (var ion in spectrum.MatchedFragmentIons)
                {
                    Assert.That(ion.IsInternalFragment, Is.True,
                        "All ions in this model's spectra are internal fragments");
                    Assert.That(ion.NeutralTheoreticalProduct.SecondaryProductType, Is.Not.Null,
                        "Internal fragment Product must have a SecondaryProductType");
                }
            }
        }

        [Test, Category("RequiresOnnxModel")]
        public static async Task PredictedSpectra_IntensitiesNormalizedToOne()
        {
            var model = new InternalFragmentIntensityModel(
                ValidPeptides, ValidCharges, ValidRetentionTimes,
                out _);

            await model.RunInferenceAsync();

            foreach (var spectrum in model.PredictedSpectra)
            {
                if (spectrum.MatchedFragmentIons.Count == 0) continue;
                double maxIntensity = spectrum.MatchedFragmentIons.Max(ion => ion.Intensity);
                Assert.That(maxIntensity, Is.EqualTo(1.0).Within(1e-9),
                    $"Max intensity in spectrum {spectrum.Name} should be normalized to 1.0");
            }
        }

        [Test, Category("RequiresOnnxModel")]
        public static async Task PredictedSpectra_MaxInternalIonsLimitRespected()
        {
            const int maxIons = 5;
            var model = new InternalFragmentIntensityModel(
                ValidPeptides, ValidCharges, ValidRetentionTimes,
                out _,
                maxInternalIonsPerPeptide: maxIons);

            await model.RunInferenceAsync();

            foreach (var spectrum in model.PredictedSpectra)
                Assert.That(spectrum.MatchedFragmentIons.Count, Is.LessThanOrEqualTo(maxIons),
                    $"Spectrum {spectrum.Name} exceeds MaxInternalIonsPerPeptide = {maxIons}");
        }

        [Test, Category("RequiresOnnxModel")]
        public static async Task PredictedSpectra_MinIntensityFilterRespected()
        {
            const double minTicNI = 0.005;
            var model = new InternalFragmentIntensityModel(
                ValidPeptides, ValidCharges, ValidRetentionTimes,
                out _,
                minIntensityFilter: minTicNI);

            await model.RunInferenceAsync();

            // After normalization, intensities are relative within each spectrum.
            // What we can verify is that the *raw* predicted TicNIs above the threshold
            // were the ones included — proxy check: all spectra should have fewer ions
            // than when using the default 0.002 threshold.
            var modelDefault = new InternalFragmentIntensityModel(
                ValidPeptides, ValidCharges, ValidRetentionTimes,
                out _);
            await modelDefault.RunInferenceAsync();

            int totalHighThreshold = model.PredictedSpectra.Sum(s => s.MatchedFragmentIons.Count);
            int totalDefaultThreshold = modelDefault.PredictedSpectra.Sum(s => s.MatchedFragmentIons.Count);

            Assert.That(totalHighThreshold, Is.LessThanOrEqualTo(totalDefaultThreshold),
                "Higher MinIntensityFilter should produce fewer ions than default");
        }

        [Test, Category("RequiresOnnxModel")]
        public static async Task PredictedSpectra_DuplicatePeptides_DeduplicatedWithWarning()
        {
            var duplicatePeptides = new List<string> { "PEPTIDEK", "PEPTIDEK", "SAMPLER" };
            var duplicateCharges = new List<int> { 2, 2, 2 };
            var duplicateRts = new List<double?> { 100.0, 100.0, 200.0 };

            var model = new InternalFragmentIntensityModel(
                duplicatePeptides, duplicateCharges, duplicateRts,
                out _);

            var warning = await model.RunInferenceAsync();

            Assert.That(model.PredictedSpectra.Count, Is.LessThanOrEqualTo(2),
                "Duplicate spectra should be deduplicated");
            Assert.That(warning, Is.Not.Null,
                "Deduplication should produce a warning");
        }

        [Test, Category("RequiresOnnxModel")]
        public static async Task PredictedSpectra_ProductTypeIsB_SecondaryIsB()
        {
            // Our model uses b-type internal fragments (double b-type cleavage)
            var model = new InternalFragmentIntensityModel(
                ValidPeptides, ValidCharges, ValidRetentionTimes,
                out _);

            await model.RunInferenceAsync();

            foreach (var spectrum in model.PredictedSpectra)
                foreach (var ion in spectrum.MatchedFragmentIons)
                {
                    Assert.That(ion.NeutralTheoreticalProduct.ProductType, Is.EqualTo(ProductType.b),
                        "Primary product type should be 'b'");
                    Assert.That(ion.NeutralTheoreticalProduct.SecondaryProductType, Is.EqualTo(ProductType.b),
                        "Secondary product type should be 'b'");
                }
        }

        // ════════════════════════════════════════════════════════════════════════
        // 5. MSP ROUND-TRIP (write library → read back → compare)
        // ════════════════════════════════════════════════════════════════════════

        [Test, Category("RequiresOnnxModel")]
        public static async Task MspRoundTrip_WrittenLibraryMatchesInMemorySpectra()
        {
            var outPath = Path.Combine(
                TestContext.CurrentContext.TestDirectory,
                "internalFragmentLibraryTest.msp");

            SpectralLibrary? savedLib = null;
            try
            {
                var model = new InternalFragmentIntensityModel(
                    ValidPeptides, ValidCharges, ValidRetentionTimes,
                    out _,
                    spectralLibrarySavePath: outPath);

                await model.RunInferenceAsync();

                Assert.That(File.Exists(outPath), Is.True, "MSP file should have been written");

                savedLib = new SpectralLibrary(new List<string> { outPath });
                var savedSpectra = savedLib.GetAllLibrarySpectra().ToList();

                // Same number of spectra
                Assert.That(savedSpectra.Count, Is.EqualTo(model.PredictedSpectra.Count),
                    "Saved library should have same number of spectra as in-memory");

                var sortedInMemory = model.PredictedSpectra
                    .OrderBy(s => s.Sequence).ThenBy(s => s.ChargeState).ToList();
                var sortedSaved = savedSpectra
                    .OrderBy(s => s.Sequence).ThenBy(s => s.ChargeState).ToList();

                for (int i = 0; i < sortedInMemory.Count; i++)
                {
                    var mem = sortedInMemory[i];
                    var saved = sortedSaved[i];

                    Assert.That(mem.Sequence, Is.EqualTo(saved.Sequence));
                    Assert.That(mem.ChargeState, Is.EqualTo(saved.ChargeState));
                    Assert.That(mem.MatchedFragmentIons.Count,
                        Is.EqualTo(saved.MatchedFragmentIons.Count),
                        $"Ion count mismatch for {mem.Name}");

                    var memFrags = mem.MatchedFragmentIons.OrderBy(f => f.Mz).ToList();
                    var savedFrags = saved.MatchedFragmentIons.OrderBy(f => f.Mz).ToList();

                    for (int j = 0; j < memFrags.Count; j++)
                    {
                        Assert.That(memFrags[j].Mz,
                            Is.EqualTo(savedFrags[j].Mz).Within(1e-4),
                            "m/z should round-trip within 0.1 mDa");
                        Assert.That(memFrags[j].Intensity,
                            Is.EqualTo(savedFrags[j].Intensity).Within(1e-6),
                            "Normalized intensity should round-trip");
                        Assert.That(memFrags[j].NeutralTheoreticalProduct.ProductType,
                            Is.EqualTo(savedFrags[j].NeutralTheoreticalProduct.ProductType),
                            "ProductType should round-trip through MSP");
                        Assert.That(memFrags[j].NeutralTheoreticalProduct.FragmentNumber,
                            Is.EqualTo(savedFrags[j].NeutralTheoreticalProduct.FragmentNumber),
                            "Start residue (FragmentNumber) should round-trip through MSP");
                        Assert.That(memFrags[j].NeutralTheoreticalProduct.SecondaryFragmentNumber,
                            Is.EqualTo(savedFrags[j].NeutralTheoreticalProduct.SecondaryFragmentNumber),
                            "End residue (SecondaryFragmentNumber) should round-trip through MSP");
                    }
                }
            }
            finally
            {
                savedLib?.CloseConnections();
                if (File.Exists(outPath)) File.Delete(outPath);
            }
        }

        [Test, Category("RequiresOnnxModel")]
        public static async Task MspRoundTrip_ParsedIonsAreInternalFragments()
        {
            var outPath = Path.Combine(
                TestContext.CurrentContext.TestDirectory,
                "internalFragmentRoundTripTest.msp");

            SpectralLibrary? savedLib = null;
            try
            {
                var model = new InternalFragmentIntensityModel(
                    ValidPeptides, ValidCharges, ValidRetentionTimes,
                    out _,
                    spectralLibrarySavePath: outPath);

                await model.RunInferenceAsync();

                savedLib = new SpectralLibrary(new List<string> { outPath });
                foreach (var spectrum in savedLib.GetAllLibrarySpectra())
                    foreach (var ion in spectrum.MatchedFragmentIons)
                        Assert.That(ion.IsInternalFragment, Is.True,
                            $"Ion in saved/re-read library should have IsInternalFragment=true; " +
                            $"annotation was '{ion.NeutralTheoreticalProduct.Annotation}'");
            }
            finally
            {
                savedLib?.CloseConnections();
                if (File.Exists(outPath)) File.Delete(outPath);
            }
        }

        // ════════════════════════════════════════════════════════════════════════
        // 6. SPECTRUM ANGLE SANITY CHECK
        // ════════════════════════════════════════════════════════════════════════

        [Test, Category("RequiresOnnxModel")]
        public static async Task SpectralAngle_SameSpectrumWithItself_IsOne()
        {
            // A library spectrum compared against its own fragment ions should give
            // spectral angle = 1.0 (perfect match).
            var model = new InternalFragmentIntensityModel(
                new List<string> { "PEPTIDEK" },
                new List<int> { 2 },
                new List<double?> { null },
                out _);

            await model.RunInferenceAsync();

            Assume.That(model.PredictedSpectra.Count, Is.EqualTo(1),
                "Need at least one spectrum for this test");

            var spectrum = model.PredictedSpectra[0];
            var result = spectrum.CalculateSpectralAngleOnTheFly(spectrum.MatchedFragmentIons);

            // SpectralContrastAngle of identical spectra = 1.0
            Assert.That(double.Parse(result, System.Globalization.CultureInfo.InvariantCulture),
                Is.EqualTo(1.0).Within(1e-4),
                "Self-comparison spectral angle should be 1.0");
        }
    }
}