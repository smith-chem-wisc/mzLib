using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using Chemistry;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SequenceConversion;
using Omics.SpectrumMatch;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Util;
using Proteomics.ProteolyticDigestion;
using Readers.ProForma;
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

        private const string CarbamidomethylPeptide = "PEPTIDEC[Common Fixed:Carbamidomethyl on C]K";

        /// <summary>
        /// A modified peptide goes through prediction into a library spectrum whichever format it came in and whichever
        /// sequence the fragments are mapped onto. Koina is sent the UNIMOD string, the validated sequence stays in the
        /// input's format, and the precursor and the C-containing y3 carry the carbamidomethyl mass.
        /// </summary>
        [TestCase(CarbamidomethylPeptide, false, CarbamidomethylPeptide, FragmentIonMappingMode.MapToValidatedFullSequence)]
        [TestCase(CarbamidomethylPeptide, false, CarbamidomethylPeptide, FragmentIonMappingMode.MapToInputFullSequence)]
        [TestCase("PEPTIDEC[UNIMOD:4]K", true, "PEPTIDEC[Unimod:Carbamidomethyl on C]K", FragmentIonMappingMode.MapToValidatedFullSequence)]
        [TestCase("PEPTIDEC[UNIMOD:4]K", true, "PEPTIDEC[Unimod:Carbamidomethyl on C]K", FragmentIonMappingMode.MapToInputFullSequence)]
        public void FragmentIntensity_ModifiedPeptide_BuildsALibrarySpectrum(string sequence, bool proForma, string expectedLabel, FragmentIonMappingMode mappingMode)
        {
            var model = new CannedHcdModel(mappingMode);
            var input = new FragmentIntensityPredictionInput(sequence, 2, 30, null, null) { SequenceParser = proForma ? ProFormaSequenceParser.Instance : null };

            var predictions = model.Predict(new List<FragmentIntensityPredictionInput> { input });
            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 12.0 }, out _);

            Assert.That(SentSequences(model.Requests), Is.EqualTo(new[] { "PEPTIDEC[UNIMOD:4]K" }));
            Assert.That(predictions[0].ValidatedFullSequence, Is.EqualTo(sequence));
            Assert.That(predictions[0].SequenceParser, Is.SameAs(input.SequenceParser));
            Assert.That(spectra.Count, Is.EqualTo(1));
            Assert.That(spectra[0].Sequence, Is.EqualTo(expectedLabel));
            Assert.That(spectra[0].PrecursorMz, Is.EqualTo(new PeptideWithSetModifications(CarbamidomethylPeptide).MonoisotopicMass.ToMz(2)).Within(1e-4));
            Assert.That(spectra[0].PrecursorMz - new PeptideWithSetModifications("PEPTIDECK").MonoisotopicMass.ToMz(2), Is.EqualTo(57.02146 / 2).Within(1e-3));
            Assert.That(FragmentMz(spectra[0], ProductType.y, 3) - UnmodifiedFragmentMz("PEPTIDECK", ProductType.y, 3), Is.EqualTo(57.02146).Within(1e-3));
        }

        /// <summary>
        /// The Prosit TMT model takes its required label under mzLib's catalog name and builds the labeled peptide in
        /// both mapping modes.
        /// </summary>
        [Test]
        public void FragmentIntensity_TmtLabeledPeptide_BuildsALibrarySpectrum(
            [Values(FragmentIonMappingMode.MapToValidatedFullSequence, FragmentIonMappingMode.MapToInputFullSequence)] FragmentIonMappingMode mappingMode)
        {
            const string sequence = "[Multiplex Label:TMT6-plex on X]PEPTIDEK";
            var model = new CannedTmtModel(mappingMode);

            var predictions = model.Predict(new List<FragmentIntensityPredictionInput> { new(sequence, 2, 30, null, "HCD") });
            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 12.0 }, out _);

            Assert.That(SentSequences(model.Requests), Is.EqualTo(new[] { "[UNIMOD:737]-PEPTIDEK" }));
            Assert.That(predictions[0].ValidatedFullSequence, Is.EqualTo("[Multiplex Label:TMT6-plex on X]PEPTIDEK"));
            Assert.That(spectra.Count, Is.EqualTo(1));
            var labeled = new PeptideWithSetModifications("[Multiplex Label:TMT6-plex on X]PEPTIDEK");
            Assert.That(spectra[0].PrecursorMz, Is.EqualTo(labeled.MonoisotopicMass.ToMz(2)).Within(1e-4));
            Assert.That(FragmentMz(spectra[0], ProductType.b, 2) - UnmodifiedFragmentMz("PEPTIDEK", ProductType.b, 2), Is.EqualTo(229.16293).Within(1e-3));
            Assert.That(spectra[0].Sequence, Is.EqualTo(sequence), "a name mzLib knows keeps its label");
        }

        /// <summary>
        /// A terminal modification stays at the terminus it was written at: an N-terminal methyl, by id or by mass, is on
        /// the b ions and not the y ions, and a C-terminal methyl or ethyl ester is on the y ions and not the b ions.
        /// </summary>
        [Test]
        public void FragmentIntensity_TerminalModification_StaysOnItsTerminus(
            [Values("[UNIMOD:34]-PEPTIDEK", "[+14.0157]-PEPTIDEK", "PEPTIDEK-[Methyl]", "PEPTIDEK-[Ethyl]")] string sequence,
            [Values(FragmentIonMappingMode.MapToValidatedFullSequence, FragmentIonMappingMode.MapToInputFullSequence)] FragmentIonMappingMode mappingMode)
        {
            var model = new CannedMs2PipModel(mappingMode);
            var input = new FragmentIntensityPredictionInput(sequence, 2, null, null, null) { SequenceParser = ProFormaSequenceParser.Instance };

            model.Predict(new List<FragmentIntensityPredictionInput> { input });
            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 1.0 }, out _);

            Assert.That(spectra.Count, Is.EqualTo(1));
            var shift = sequence.Contains("Ethyl") ? 28.0313 : 14.01565;
            var nTerminal = sequence.StartsWith("[");
            Assert.That(FragmentMz(spectra[0], ProductType.b, 2) - UnmodifiedFragmentMz("PEPTIDEK", ProductType.b, 2), Is.EqualTo(nTerminal ? shift : 0).Within(1e-3));
            Assert.That(FragmentMz(spectra[0], ProductType.y, 3) - UnmodifiedFragmentMz("PEPTIDEK", ProductType.y, 3), Is.EqualTo(nTerminal ? 0 : shift).Within(1e-3));
        }

        /// <summary>
        /// ProteaseGuru's configuration: an allow-list model, RemoveIncompatibleElements, fragments mapped onto the input.
        /// Koina is sent the peptide without the modifications the model doesn't allow, but the spectrum is built from
        /// the input, which still carries them, so they're found outside the model's allow-list.
        /// </summary>
        [TestCase("PEPS[UNIMOD:21]IDEC[UNIMOD:4]K", true, "PEPSIDEC[UNIMOD:4]K", "PEPS[Common Biological:Phosphorylation on S]IDEC[Common Fixed:Carbamidomethyl on C]K")]
        [TestCase("[Common Fixed:TMT6plex on N-terminus]PEPTIDEK", false, "PEPTIDEK", "[Multiplex Label:TMT6-plex on X]PEPTIDEK")]
        public void FragmentIntensity_AllowListModelRemovingAModification_BuildsTheInputWithIt(string sequence, bool proForma, string sent, string equivalent)
        {
            var model = new CannedHcdModel(FragmentIonMappingMode.MapToInputFullSequence, SequenceConversionHandlingMode.RemoveIncompatibleElements);
            var input = new FragmentIntensityPredictionInput(sequence, 2, 30, null, null) { SequenceParser = proForma ? ProFormaSequenceParser.Instance : null };

            var predictions = model.Predict(new List<FragmentIntensityPredictionInput> { input });
            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 1.0 }, out _);

            Assert.That(SentSequences(model.Requests), Is.EqualTo(new[] { sent }));
            Assert.That(model.ValidInputsMask, Is.EqualTo(new[] { true }), predictions[0].Warning?.Message);
            Assert.That(spectra.Count, Is.EqualTo(1));
            Assert.That(spectra[0].PrecursorMz, Is.EqualTo(new PeptideWithSetModifications(equivalent).MonoisotopicMass.ToMz(2)).Within(1e-4));
        }

        /// <summary>
        /// A mass only a UniProt modification matches (3'-geranyl-2',N2-cyclotryptophan on W), removed for Koina by an
        /// allow-list model, is found in mzLib's UniProt catalog when the spectrum is built from the input.
        /// </summary>
        [Test]
        public void FragmentIntensity_AllowListModelRemovingAUniProtOnlyMass_BuildsTheInputWithTheUniProtModification()
        {
            var model = new CannedHcdModel(FragmentIonMappingMode.MapToInputFullSequence, SequenceConversionHandlingMode.RemoveIncompatibleElements);
            var input = new FragmentIntensityPredictionInput("GLSW[+136.1252]DEFK", 2, 30, null, null) { SequenceParser = ProFormaSequenceParser.Instance };

            var predictions = model.Predict(new List<FragmentIntensityPredictionInput> { input });
            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 1.0 }, out _);

            Assert.That(SentSequences(model.Requests), Is.EqualTo(new[] { "GLSWDEFK" }));
            Assert.That(model.ValidInputsMask, Is.EqualTo(new[] { true }), predictions[0].Warning?.Message);
            Assert.That(spectra.Count, Is.EqualTo(1));
            Assert.That(spectra[0].Sequence, Does.StartWith("GLSW[UniProt:"));
            Assert.That(spectra[0].PrecursorMz - new PeptideWithSetModifications("GLSWDEFK").MonoisotopicMass.ToMz(2), Is.EqualTo(136.1252 / 2).Within(1e-3));
        }

        /// <summary>
        /// mzLib names are read the way mzLib reads them, except that a modification at a terminus is one of that
        /// terminus's: the N-terminal <c>Methyl on X</c> stays on the b ions (mzLib's dictionary holds the C-terminal
        /// one), and a C-terminal modification written on the last residue moves to the C-terminus, as reading the name
        /// moves it. A known id under a type mzLib doesn't pair it with is labeled with mzLib's own type. Koina's lookup
        /// can't take the N-terminal methyl, so that one is removed for Koina and mapped onto the input.
        /// </summary>
        [TestCase("[Unimod:Methyl on X]AGLSDEFK", "[Unimod:Methyl on X]AGLSDEFK", "AGLSDEFK", 14.01565, 0, FragmentIonMappingMode.MapToInputFullSequence, SequenceConversionHandlingMode.RemoveIncompatibleElements)]
        [TestCase("GLSDEFA[Unimod:Amidated on X]", "GLSDEFA-[Unimod:Amidated on X]", "GLSDEFA", 0, -0.98402, FragmentIonMappingMode.MapToValidatedFullSequence, SequenceConversionHandlingMode.ReturnNull)]
        [TestCase("GLSDEFA[Unimod:Amidated on X]", "GLSDEFA-[Unimod:Amidated on X]", "GLSDEFA", 0, -0.98402, FragmentIonMappingMode.MapToInputFullSequence, SequenceConversionHandlingMode.ReturnNull)]
        [TestCase("PEPTIDEK[Common Fixed:TMT6plex on K]", "PEPTIDEK[Unimod:TMT6plex on K]", "PEPTIDEK", 0, 229.16293, FragmentIonMappingMode.MapToValidatedFullSequence, SequenceConversionHandlingMode.ReturnNull)]
        [TestCase("PEPTIDEK[Common Fixed:TMT6plex on K]", "PEPTIDEK[Unimod:TMT6plex on K]", "PEPTIDEK", 0, 229.16293, FragmentIonMappingMode.MapToInputFullSequence, SequenceConversionHandlingMode.ReturnNull)]
        public void FragmentIntensity_MzLibName_IsBuiltAndLabeledTheWayMzLibReadsIt(string sequence, string expectedLabel, string baseSequence,
            double bShift, double yShift, FragmentIonMappingMode mappingMode, SequenceConversionHandlingMode modHandlingMode)
        {
            var model = new CannedMs2PipModel(mappingMode, modHandlingMode);

            var predictions = model.Predict(new List<FragmentIntensityPredictionInput> { new(sequence, 2, null, null, null) });
            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 1.0 }, out var warning);

            Assert.That(spectra.Count, Is.EqualTo(1), predictions[0].Warning?.Message ?? warning?.Message);
            Assert.That(spectra[0].Sequence, Is.EqualTo(expectedLabel));
            Assert.That(FragmentMz(spectra[0], ProductType.b, 2) - UnmodifiedFragmentMz(baseSequence, ProductType.b, 2), Is.EqualTo(bShift).Within(1e-3));
            Assert.That(FragmentMz(spectra[0], ProductType.y, 3) - UnmodifiedFragmentMz(baseSequence, ProductType.y, 3), Is.EqualTo(yShift).Within(1e-3));
        }

        /// <summary>
        /// Koina resolves modifications to protein modifications only. Masses that only RNA modifications of mzLib's
        /// match (terminal dephosphorylation, the 3' cyclic phosphate's water loss, dimethyladenosine), on an
        /// allow-list model that removes them for Koina and maps the fragments onto the input: the peptide is built
        /// with a protein modification of that mass where there is one, and otherwise fails with a warning.
        /// </summary>
        [TestCase("[-79.9663]-PEPTIDEK", null)]
        [TestCase("PEPTIDEK-[-18.0106]", -18.01056)]
        [TestCase("PEPTIDEA[+28.0313]K", 28.0313)]
        public void FragmentIntensity_MassOnlyAnRnaModificationMatches_GetsAProteinModificationOrFails(string sequence, double? proteinModificationMass)
        {
            var model = new CannedHcdModel(FragmentIonMappingMode.MapToInputFullSequence, SequenceConversionHandlingMode.RemoveIncompatibleElements);
            var input = new FragmentIntensityPredictionInput(sequence, 2, 30, null, null) { SequenceParser = ProFormaSequenceParser.Instance };

            var predictions = model.Predict(new List<FragmentIntensityPredictionInput> { input });
            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 1.0 }, out _);

            var rnaNames = Mods.AllRnaModsList.Select(m => $"[{m.ModificationType}:{m.IdWithMotif}]").ToHashSet();
            if (proteinModificationMass == null)
            {
                Assert.That(model.ValidInputsMask, Is.EqualTo(new[] { false }));
                Assert.That(predictions[0].Warning?.Message, Does.Contain(sequence));
                Assert.That(spectra, Is.Empty);
                return;
            }

            Assert.That(spectra.Count, Is.EqualTo(1), predictions[0].Warning?.Message);
            Assert.That(rnaNames.Any(spectra[0].Sequence.Contains), Is.False, spectra[0].Sequence);
            var baseSequence = sequence.Contains("PEPTIDEA") ? "PEPTIDEAK" : "PEPTIDEK";
            Assert.That(spectra[0].PrecursorMz - new PeptideWithSetModifications(baseSequence).MonoisotopicMass.ToMz(2),
                Is.EqualTo(proteinModificationMass.Value / 2).Within(1e-3));
        }

        /// <summary>
        /// A dropped modification is gone from the request sent to Koina and from the validated sequence, and the
        /// spectrum is built without it.
        /// </summary>
        [TestCase("PEPS[Common Biological:Phosphorylation on S]IDEC[Common Fixed:Carbamidomethyl on C]K", false, "PEPSIDEC[Common Fixed:Carbamidomethyl on C]K")]
        [TestCase("PEPS[UNIMOD:21]IDEC[UNIMOD:4]K", true, "PEPSIDEC[UNIMOD:4]K")]
        public void FragmentIntensity_RemoveIncompatibleElements_DropsTheModificationFromWhatIsSentAndValidated(string sequence, bool proForma, string expectedValidated)
        {
            var model = new CannedHcdModel(FragmentIonMappingMode.MapToValidatedFullSequence, SequenceConversionHandlingMode.RemoveIncompatibleElements);
            var input = new FragmentIntensityPredictionInput(sequence, 2, 30, null, null) { SequenceParser = proForma ? ProFormaSequenceParser.Instance : null };

            var predictions = model.Predict(new List<FragmentIntensityPredictionInput> { input });
            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 12.0 }, out _);

            Assert.That(SentSequences(model.Requests), Is.EqualTo(new[] { "PEPSIDEC[UNIMOD:4]K" }));
            Assert.That(predictions[0].ValidatedFullSequence, Is.EqualTo(expectedValidated));
            Assert.That(spectra[0].PrecursorMz, Is.EqualTo(new PeptideWithSetModifications("PEPSIDEC[Common Fixed:Carbamidomethyl on C]K").MonoisotopicMass.ToMz(2)).Within(1e-4));
        }

        /// <summary>
        /// Under MapToInputFullSequence the fragments are mapped onto the input when Koina answers. A peptide mzLib can't
        /// build (UNIMOD:2114, which an accept-all model sends but mzLib has no definition for) fails alone, with a
        /// warning, like a rejected input; the rest of the batch still gets its spectrum.
        /// </summary>
        [Test]
        public void FragmentIntensity_MapToInputFullSequence_PeptideMzLibCannotBuildFailsAlone(
            [Values(SequenceConversionHandlingMode.ReturnNull, SequenceConversionHandlingMode.RemoveIncompatibleElements)] SequenceConversionHandlingMode mode)
        {
            var model = new CannedMs2PipModel(FragmentIonMappingMode.MapToInputFullSequence, mode);
            var unbuildable = new FragmentIntensityPredictionInput("PEPTIDEK[UNIMOD:2114]", 2, null, null, null) { SequenceParser = ProFormaSequenceParser.Instance };
            var buildable = new FragmentIntensityPredictionInput("PEPTIDECK", 2, null, null, null);

            var predictions = model.Predict(new List<FragmentIntensityPredictionInput> { unbuildable, buildable });
            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 1.0, 2.0 }, out _);

            Assert.That(SentSequences(model.Requests), Is.EqualTo(new[] { "PEPTIDEK[UNIMOD:2114]", "PEPTIDECK" }));
            Assert.That(predictions[0].FragmentAnnotations, Is.Null);
            Assert.That(predictions[0].Warning?.Message, Does.Contain("PEPTIDEK[UNIMOD:2114]"));
            Assert.That(model.ValidInputsMask, Is.EqualTo(new[] { false, true }));
            Assert.That(predictions[1].FragmentAnnotations, Is.Not.Empty);
            Assert.That(spectra.Select(s => s.Sequence), Is.EqualTo(new[] { "PEPTIDECK" }));
        }

        /// <summary>
        /// A ProForma mass-only modification is resolved once, for Koina, and ValidatedFullSequence names that same
        /// modification, so the library spectrum built from it has the mass of what was predicted: Thr->Val for -1.9793,
        /// and phospho, not sulfo, for +79.9568.
        /// </summary>
        [TestCase("GLST[-1.9793]DEFK", "GLST[UNIMOD:1210]DEFK", "GLSTDEFK", -1.97926)]
        [TestCase("GLSS[+79.9568]DEFK", "GLSS[UNIMOD:21]DEFK", "GLSSDEFK", 79.96633)]
        public void FragmentIntensity_ProFormaMassOnlyModification_LibrarySpectrumHasThePredictedMass(string sequence, string expectedValidated, string baseSequence, double modificationMass)
        {
            var model = new CannedMs2PipModel(FragmentIonMappingMode.MapToValidatedFullSequence);
            var input = new FragmentIntensityPredictionInput(sequence, 2, null, null, null) { SequenceParser = ProFormaSequenceParser.Instance };

            var predictions = model.Predict(new List<FragmentIntensityPredictionInput> { input });
            var spectra = model.GenerateLibrarySpectraFromPredictions(new double?[] { 1.0 }, out _);

            Assert.That(SentSequences(model.Requests), Is.EqualTo(new[] { expectedValidated }));
            Assert.That(predictions[0].ValidatedFullSequence, Is.EqualTo(expectedValidated));
            Assert.That(spectra.Count, Is.EqualTo(1));
            Assert.That(spectra[0].PrecursorMz,
                Is.EqualTo((new PeptideWithSetModifications(baseSequence).MonoisotopicMass + modificationMass).ToMz(2)).Within(1e-3));
        }

        [Test]
        public void FragmentIntensity_MapToInputFullSequence_PeptideMzLibCannotBuildThrowsInThrowExceptionMode()
        {
            var model = new CannedMs2PipModel(FragmentIonMappingMode.MapToInputFullSequence, SequenceConversionHandlingMode.ThrowException);
            var unbuildable = new FragmentIntensityPredictionInput("PEPTIDEK[UNIMOD:2114]", 2, null, null, null) { SequenceParser = ProFormaSequenceParser.Instance };

            Assert.That(() => model.Predict(new List<FragmentIntensityPredictionInput> { unbuildable, new("PEPTIDECK", 2, null, null, null) }),
                Throws.ArgumentException.With.Message.Contains("PEPTIDEK[UNIMOD:2114]"));
            Assert.That(SentSequences(model.Requests), Is.EqualTo(new[] { "PEPTIDEK[UNIMOD:2114]", "PEPTIDECK" }), "it fails after Koina answered");
        }

        /// <summary>
        /// The Koina payload is serialized when the requests are built. An input the model's serializer can't write is
        /// dropped there with a warning and marked invalid, and the rest of the batch is still sent; in ThrowException
        /// mode it throws instead. The same for every family.
        /// </summary>
        [Test]
        public void FragmentIntensity_InputThatCannotBeSerializedForKoina_IsDroppedAlone(
            [Values(SequenceConversionHandlingMode.ReturnNull, SequenceConversionHandlingMode.RemoveIncompatibleElements)] SequenceConversionHandlingMode mode)
        {
            var model = new UnserializableFragmentModel("PEPTIDEK") { ModHandlingMode = mode };

            var predictions = model.Predict(new List<FragmentIntensityPredictionInput>
            {
                new("PEPTIDEK", 2, null, null, null),
                new("ELVISK", 2, null, null, null)
            });

            Assert.That(SentSequences(model.Requests), Is.EqualTo(new[] { "ELVISK" }));
            Assert.That(model.ValidInputsMask, Is.EqualTo(new[] { false, true }));
            Assert.That(model.ModelInputs[0].ValidatedFullSequence, Is.Null);
            Assert.That(predictions[0].FragmentAnnotations, Is.Null);
            Assert.That(predictions[0].Warning?.Message, Does.Contain(UnserializableMessage));
            Assert.That(predictions[1].FragmentAnnotations, Is.Not.Empty);
        }

        [Test]
        public void RetentionTime_InputThatCannotBeSerializedForKoina_IsDroppedAlone(
            [Values(SequenceConversionHandlingMode.ReturnNull, SequenceConversionHandlingMode.RemoveIncompatibleElements)] SequenceConversionHandlingMode mode)
        {
            var model = new UnserializableRtModel("PEPTIDEK") { ModHandlingMode = mode };

            var predictions = model.Predict(new List<RetentionTimePredictionInput> { new("PEPTIDEK"), new("ELVISK") });

            Assert.That(SentSequences(model.Requests), Is.EqualTo(new[] { "ELVISK" }));
            Assert.That(model.ValidInputsMask, Is.EqualTo(new[] { false, true }));
            Assert.That(model.ModelInputs[0].ValidatedFullSequence, Is.Null);
            Assert.That(predictions[0].PredictedRetentionTime, Is.Null);
            Assert.That(predictions[0].Warning?.Message, Does.Contain(UnserializableMessage));
            Assert.That(predictions[1].PredictedRetentionTime, Is.EqualTo(1.0));
        }

        [Test]
        public void Ccs_InputThatCannotBeSerializedForKoina_IsDroppedAlone(
            [Values(SequenceConversionHandlingMode.ReturnNull, SequenceConversionHandlingMode.RemoveIncompatibleElements)] SequenceConversionHandlingMode mode)
        {
            var model = new UnserializableCcsModel("PEPTIDEK") { ModHandlingMode = mode };

            var predictions = model.Predict(new List<CCSPredictionInput> { new("PEPTIDEK", 2), new("ELVISK", 2) });

            Assert.That(SentSequences(model.Requests), Is.EqualTo(new[] { "ELVISK" }));
            Assert.That(model.ValidInputsMask, Is.EqualTo(new[] { false, true }));
            Assert.That(model.ModelInputs[0].ValidatedFullSequence, Is.Null);
            Assert.That(predictions[0].PredictedCCS, Is.Null);
            Assert.That(predictions[0].Warning?.Message, Does.Contain(UnserializableMessage));
            Assert.That(predictions[1].PredictedCCS, Is.EqualTo(1.0));
        }

        [Test]
        public void Detectability_InputThatCannotBeSerializedForKoina_IsDroppedAlone(
            [Values(SequenceConversionHandlingMode.ReturnNull, SequenceConversionHandlingMode.RemoveIncompatibleElements)] SequenceConversionHandlingMode mode)
        {
            var model = new UnserializableDetectabilityModel("PEPTIDEK") { ModHandlingMode = mode };

            var predictions = model.Predict(new List<DetectabilityPredictionInput> { new("PEPTIDEK"), new("ELVISK") });

            Assert.That(SentSequences(model.Requests), Is.EqualTo(new[] { "ELVISK" }));
            Assert.That(model.ValidInputsMask, Is.EqualTo(new[] { false, true }));
            Assert.That(model.ModelInputs[0].ValidatedFullSequence, Is.Null);
            Assert.That(predictions[0].DetectabilityProbabilities, Is.Null);
            Assert.That(predictions[0].Warning?.Message, Does.Contain(UnserializableMessage));
            Assert.That(predictions[1].DetectabilityProbabilities, Is.Not.Null);
        }

        [Test]
        public void EveryFamily_InputThatCannotBeSerializedForKoina_ThrowsInThrowExceptionMode()
        {
            const SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException;
            var fragment = new UnserializableFragmentModel("PEPTIDEK") { ModHandlingMode = mode };
            var rt = new UnserializableRtModel("PEPTIDEK") { ModHandlingMode = mode };
            var ccs = new UnserializableCcsModel("PEPTIDEK") { ModHandlingMode = mode };
            var detectability = new UnserializableDetectabilityModel("PEPTIDEK") { ModHandlingMode = mode };

            Assert.Multiple(() =>
            {
                Assert.That(() => fragment.Predict(new List<FragmentIntensityPredictionInput> { new("ELVISK", 2, null, null, null), new("PEPTIDEK", 2, null, null, null) }),
                    Throws.ArgumentException.With.Message.Contains(UnserializableMessage));
                Assert.That(() => rt.Predict(new List<RetentionTimePredictionInput> { new("ELVISK"), new("PEPTIDEK") }),
                    Throws.ArgumentException.With.Message.Contains(UnserializableMessage));
                Assert.That(() => ccs.Predict(new List<CCSPredictionInput> { new("ELVISK", 2), new("PEPTIDEK", 2) }),
                    Throws.ArgumentException.With.Message.Contains(UnserializableMessage));
                Assert.That(() => detectability.Predict(new List<DetectabilityPredictionInput> { new("ELVISK"), new("PEPTIDEK") }),
                    Throws.ArgumentException.With.Message.Contains(UnserializableMessage));
                Assert.That(fragment.Requests.Concat(rt.Requests).Concat(ccs.Requests).Concat(detectability.Requests), Is.Empty,
                    "nothing is sent once the batch has thrown");
            });
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

        [Test]
        public void EveryFamily_Predict_ParsesWithTheInputsOwnSequenceParser()
        {
            // Each family's validation loop must hand the input record's parser to TryCleanSequence, not null.
            var fragmentParser = new CountingParser();
            var rtParser = new CountingParser();
            var ccsParser = new CountingParser();
            var detectabilityParser = new CountingParser();

            new FakeFragmentModel().Predict(new List<FragmentIntensityPredictionInput> { new("PEPTIDEK", 2, 30, null, null) { SequenceParser = fragmentParser } });
            new FakeRtModel().Predict(new List<RetentionTimePredictionInput> { new("PEPTIDEK") { SequenceParser = rtParser } });
            new FakeCcsModel().Predict(new List<CCSPredictionInput> { new("PEPTIDEK", 2) { SequenceParser = ccsParser } });
            new FakeDetectabilityModel().Predict(new List<DetectabilityPredictionInput> { new("PEPTIDEK") { SequenceParser = detectabilityParser } });

            Assert.Multiple(() =>
            {
                Assert.That(fragmentParser.Calls, Is.EqualTo(1), "fragment intensity");
                Assert.That(rtParser.Calls, Is.EqualTo(1), "retention time");
                Assert.That(ccsParser.Calls, Is.EqualTo(1), "CCS");
                Assert.That(detectabilityParser.Calls, Is.EqualTo(1), "detectability");
            });
        }

        private static double FragmentMz(LibrarySpectrum spectrum, ProductType type, int number) =>
            spectrum.MatchedFragmentIons.Single(ion => ion.NeutralTheoreticalProduct.ProductType == type
                && ion.NeutralTheoreticalProduct.FragmentNumber == number).Mz;

        private static double UnmodifiedFragmentMz(string baseSequence, ProductType type, int number)
        {
            var products = new List<Product>();
            new PeptideWithSetModifications(baseSequence).Fragment(DissociationType.HCD, FragmentationTerminus.Both, products);
            return products.Single(p => p.ProductType == type && p.FragmentNumber == number).NeutralMass.ToMz(1);
        }

        private static List<string> SentSequences(List<Dictionary<string, object>> requests) =>
            requests.SelectMany(request => PeptideSequences(request)).ToList();

        private static string[] PeptideSequences(Dictionary<string, object> request)
        {
            var input = ((IEnumerable<object>)request["inputs"]).Single(i => (string)i.GetType().GetProperty("name")!.GetValue(i)! == "peptide_sequences");
            return ((Array)input.GetType().GetProperty("data")!.GetValue(input)!).Cast<string>().ToArray();
        }

        // Answers b2+1 and y3+1 for every peptide in the request, and keeps the request.
        private static Task<string> CannedFragments(List<Dictionary<string, object>> requests, Dictionary<string, object> request)
        {
            requests.Add(request);
            int peptides = PeptideSequences(request).Length;
            string Repeat(string values) => string.Join(",", Enumerable.Repeat(values, peptides));
            string Output(string name, string datatype, string data) =>
                "{\"name\":\"" + name + "\",\"datatype\":\"" + datatype + "\",\"shape\":[" + 2 * peptides + "],\"data\":[" + data + "]}";
            return Task.FromResult("{\"outputs\":[" +
                Output("annotation", "BYTES", Repeat("\"b2+1\",\"y3+1\"")) + "," +
                Output("mz", "FP32", Repeat("0.0,0.0")) + "," +
                Output("intensities", "FP32", Repeat("0.5,1.0")) + "]}");
        }

        private sealed class CannedHcdModel : Prosit2020IntensityHCD
        {
            public CannedHcdModel(FragmentIonMappingMode mappingMode, SequenceConversionHandlingMode modHandlingMode = SequenceConversionHandlingMode.ReturnNull)
                : base(modHandlingMode, fragmentIonMappingMode: mappingMode) { }

            public List<Dictionary<string, object>> Requests { get; } = new();

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
                => CannedFragments(Requests, request);
        }

        private sealed class CannedTmtModel : Prosit2020IntensityTMT
        {
            public CannedTmtModel(FragmentIonMappingMode mappingMode) : base(fragmentIonMappingMode: mappingMode) { }

            public List<Dictionary<string, object>> Requests { get; } = new();

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
                => CannedFragments(Requests, request);
        }

        private sealed class CannedMs2PipModel : Ms2PipHCD2021
        {
            public CannedMs2PipModel(FragmentIonMappingMode mappingMode, SequenceConversionHandlingMode modHandlingMode = SequenceConversionHandlingMode.ReturnNull)
                : base(modHandlingMode, fragmentIonMappingMode: mappingMode) { }

            public List<Dictionary<string, object>> Requests { get; } = new();

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
                => CannedFragments(Requests, request);
        }

        private const string UnserializableMessage = "This sequence is fine for mzLib but not for Koina.";

        // Answers 1.0 for every peptide in the request (four values per peptide for detectability), and keeps the request.
        private static Task<string> CannedValues(List<Dictionary<string, object>> requests, Dictionary<string, object> request, string name, int valuesPerPeptide = 1)
        {
            requests.Add(request);
            int values = PeptideSequences(request).Length * valuesPerPeptide;
            return Task.FromResult("{\"outputs\":[{\"name\":\"" + name + "\",\"datatype\":\"FP32\",\"shape\":[" + values + "],\"data\":["
                + string.Join(",", Enumerable.Repeat("1.0", values)) + "]}]}");
        }

        // Models whose serializer can't write one base sequence for Koina, though it cleans fine.
        private sealed class UnserializableFragmentModel : FragmentIntensityModel
        {
            public UnserializableFragmentModel(string unserializable)
                : base(new UnserializableConverter(CreateUnimodConverter(UnimodSequenceFormatSchema.Instance, new HashSet<int>()), unserializable)) { }

            public List<Dictionary<string, object>> Requests { get; } = new();
            public override string ModelName => "Unserializable";
            public override int MaxBatchSize => 1000;
            public override int MaxNumberOfBatchesPerRequest { get; init; } = 1;
            public override int ThrottlingDelayInMilliseconds { get; init; } = 0;
            public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1000;
            public override int MaxPeptideLength => 30;
            public override int MinPeptideLength => 1;
            public override HashSet<int> AllowedPrecursorCharges => new() { 2 };
            public override SequenceConversionHandlingMode ModHandlingMode { get; init; } = SequenceConversionHandlingMode.ReturnNull;
            public override IncompatibleParameterHandlingMode ParameterHandlingMode { get; init; } = IncompatibleParameterHandlingMode.ReturnNull;
            public override FragmentIonMappingMode FragmentIonMappingMode { get; init; } = FragmentIonMappingMode.MapToValidatedFullSequence;

            protected override List<Dictionary<string, object>> ToBatchedRequests(List<FragmentIntensityPredictionInput> validInputs)
                => new() { BuildBatchedRequest(0, new InputField("peptide_sequences", "BYTES", validInputs.Select(p => GetKoinaSequence(p)).ToArray())) };

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
                => CannedFragments(Requests, request);
        }

        private sealed class UnserializableRtModel : RetentionTimeModel
        {
            public UnserializableRtModel(string unserializable)
                : base(new UnserializableConverter(CreateUnimodConverter(UnimodSequenceFormatSchema.Instance, new HashSet<int>()), unserializable)) { }

            public List<Dictionary<string, object>> Requests { get; } = new();
            public override string ModelName => "UnserializableRt";
            public override int MaxBatchSize => 1000;
            public override int MaxNumberOfBatchesPerRequest { get; init; } = 1;
            public override int ThrottlingDelayInMilliseconds { get; init; } = 0;
            public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1000;
            public override int MaxPeptideLength => 30;
            public override int MinPeptideLength => 1;
            public override bool IsIndexedRetentionTimeModel => true;
            public override SequenceConversionHandlingMode ModHandlingMode { get; init; } = SequenceConversionHandlingMode.ReturnNull;

            protected override List<Dictionary<string, object>> ToBatchedRequests(List<RetentionTimePredictionInput> validInputs)
                => new() { BuildBatchedRequest(0, new InputField("peptide_sequences", "BYTES", validInputs.Select(p => GetKoinaSequence(p)).ToArray())) };

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
                => CannedValues(Requests, request, "irt");
        }

        private sealed class UnserializableCcsModel : CollisionalCrossSectionModel
        {
            public UnserializableCcsModel(string unserializable)
                : base(new UnserializableConverter(CreateUnimodConverter(UnimodSequenceFormatSchema.Instance, new HashSet<int>()), unserializable)) { }

            public List<Dictionary<string, object>> Requests { get; } = new();
            public override string ModelName => "UnserializableCcs";
            public override int MaxBatchSize => 1000;
            public override int MaxNumberOfBatchesPerRequest { get; init; } = 1;
            public override int ThrottlingDelayInMilliseconds { get; init; } = 0;
            public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1000;
            public override int MaxPeptideLength => 30;
            public override int MinPeptideLength => 1;
            public override SequenceConversionHandlingMode ModHandlingMode { get; init; } = SequenceConversionHandlingMode.ReturnNull;

            protected override List<Dictionary<string, object>> ToBatchedRequests(List<CCSPredictionInput> validInputs)
                => new() { BuildBatchedRequest(0, new InputField("peptide_sequences", "BYTES", validInputs.Select(p => GetKoinaSequence(p)).ToArray())) };

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
                => CannedValues(Requests, request, "ccs");
        }

        private sealed class UnserializableDetectabilityModel : DetectabilityModel
        {
            public UnserializableDetectabilityModel(string unserializable)
                : base(new UnserializableConverter(CreateUnimodConverter(UnimodSequenceFormatSchema.Instance, new HashSet<int>()), unserializable)) { }

            public List<Dictionary<string, object>> Requests { get; } = new();
            public override string ModelName => "UnserializableDetectability";
            public override int MaxBatchSize => 1000;
            public override int MaxNumberOfBatchesPerRequest { get; init; } = 1;
            public override int ThrottlingDelayInMilliseconds { get; init; } = 0;
            public override int BenchmarkedTimeForOneMaxBatchSizeInMilliseconds => 1000;
            public override int MaxPeptideLength => 30;
            public override int MinPeptideLength => 1;
            public override int NumberOfDetectabilityClasses => 4;
            public override List<string> DetectabilityClasses => new() { "Not Detectable", "Low", "Intermediate", "High" };
            public override SequenceConversionHandlingMode ModHandlingMode { get; init; } = SequenceConversionHandlingMode.ReturnNull;

            protected override List<Dictionary<string, object>> ToBatchedRequests(List<DetectabilityPredictionInput> validInputs)
                => new() { BuildBatchedRequest(0, new InputField("peptide_sequences", "BYTES", validInputs.Select(p => GetKoinaSequence(p)).ToArray())) };

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
                => CannedValues(Requests, request, "detectability", 4);
        }

        private sealed class UnserializableConverter(ISequenceConverter inner, string unserializable) : ISequenceConverter
        {
            public string FormatName => inner.FormatName;
            public string SourceFormatName => inner.SourceFormatName;
            public string TargetFormatName => inner.TargetFormatName;
            public ISequenceParser Parser => inner.Parser;
            public ISequenceSerializer Serializer => inner.Serializer;

            public CanonicalSequence? Parse(string input, ConversionWarnings? warnings = null, SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
                => inner.Parse(input, warnings, mode);

            public string? Serialize(CanonicalSequence sequence, ConversionWarnings? warnings = null, SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
                => sequence.BaseSequence == unserializable
                    ? throw new SequenceConversionException(UnserializableMessage, ConversionFailureReason.IncompatibleModifications)
                    : inner.Serialize(sequence, warnings, mode);

            public string? Convert(string input, ConversionWarnings? warnings = null, SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
                => inner.Convert(input, warnings, mode);
        }

        private sealed class CountingParser : ISequenceParser
        {
            public int Calls { get; private set; }
            public string FormatName => MzLibSequenceParser.Instance.FormatName;
            public SequenceFormatSchema Schema => MzLibSequenceParser.Instance.Schema;
            public bool CanParse(string input) => MzLibSequenceParser.Instance.CanParse(input);

            public CanonicalSequence? Parse(string input, ConversionWarnings? warnings = null, SequenceConversionHandlingMode mode = SequenceConversionHandlingMode.ThrowException)
            {
                Calls++;
                return MzLibSequenceParser.Instance.Parse(input, warnings, mode);
            }
        }

        // ── canned-transport subclasses of real models ─────────────────────────────

        private sealed class FakeFragmentModel : Prosit2019Intensity
        {
            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
                => Task.FromResult(IntensityJson);
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
