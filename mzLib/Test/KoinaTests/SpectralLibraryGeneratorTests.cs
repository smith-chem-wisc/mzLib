using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using Chemistry;
using Chromatography;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.CZE;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.SupportedModels.FlyabilityModels;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using PredictionClients.Koina.Util;
using PredictionClients.SpectralLibraryGeneration;
using SpectralLibraryReader = Readers.SpectralLibrary.SpectralLibrary;

namespace Test.KoinaTests
{
    /// <summary>
    /// Runs the generator against real Koina model classes whose transport returns canned responses,
    /// so nothing here needs the network.
    /// </summary>
    [TestFixture]
    public class SpectralLibraryGeneratorTests
    {
        private const string Oxidized = "PEPTIDEM[Common Variable:Oxidation on M]K";

        [Test]
        public void SpectraCarryTheirRetentionTimesAndTheModelNames()
        {
            var fragmentModel = new FakeFragmentModel();
            var generator = new SpectralLibraryGenerator(fragmentModel, new FakeKoinaRetentionTimeModel());

            List<LibrarySpectrum> spectra = generator.GenerateLibrarySpectra(new List<FragmentIntensityPredictionInput>
            {
                Input("PEPTIDEK", 2),
                Input("PEPTIDEK", 3),
                Input("ELVISK", 2)
            });

            Assert.That(spectra.Select(s => s.Name), Is.EqualTo(new[] { "PEPTIDEK/2", "PEPTIDEK/3", "ELVISK/2" }));
            Assert.That(spectra.Select(s => s.RetentionTime), Is.EqualTo(new double?[] { 30.5, 30.5, 12.25 }));
            Assert.That(spectra.Select(s => s.FragmentModelName), Is.All.EqualTo("Prosit_2020_intensity_HCD"));
            Assert.That(spectra.Select(s => s.RetentionTimeModelName), Is.All.EqualTo("Prosit_2019_irt"));
        }

        [Test]
        public void WithoutARetentionTimePredictorSpectraHaveNoRetentionTime()
        {
            var generator = new SpectralLibraryGenerator(new FakeFragmentModel());

            LibrarySpectrum spectrum = generator.GenerateLibrarySpectra(new List<FragmentIntensityPredictionInput> { Input("PEPTIDEK", 2) }).Single();

            Assert.That(spectrum.RetentionTime, Is.Null);
            Assert.That(spectrum.RetentionTimeModelName, Is.Null);
            Assert.That(spectrum.FragmentModelName, Is.EqualTo("Prosit_2020_intensity_HCD"));
        }

        [Test]
        public void APeptideWithNoUsableRetentionTimeIsGivenNoneAndNoRetentionTimeModel()
        {
            var predictor = new FakeRetentionTimePredictor(new Dictionary<string, double?>
            {
                ["PEPTIDEK"] = 30.5,
                ["ELVISK"] = null,
                ["PEPTIDER"] = double.NaN
            });
            var generator = new SpectralLibraryGenerator(new FakeFragmentModel(), predictor);

            List<LibrarySpectrum> spectra = generator.GenerateLibrarySpectra(new List<FragmentIntensityPredictionInput>
            {
                Input("PEPTIDEK", 2),
                Input("ELVISK", 2),
                Input("PEPTIDER", 2)
            });

            Assert.That(spectra.Select(s => s.RetentionTime), Is.EqualTo(new double?[] { 30.5, null, null }));
            Assert.That(spectra.Select(s => s.RetentionTimeModelName), Is.EqualTo(new[] { "FakePredictor", null, null }));
        }

        [Test]
        public void EachPeptideIsSentToTheRetentionTimePredictorOnceWithItsBaseSequence()
        {
            var predictor = new FakeRetentionTimePredictor(new Dictionary<string, double?> { ["PEPTIDEK"] = 1, [Oxidized] = 2 });
            var generator = new SpectralLibraryGenerator(new FakeFragmentModel(), predictor);

            generator.GenerateLibrarySpectra(new List<FragmentIntensityPredictionInput>
            {
                Input("PEPTIDEK", 2),
                Input("PEPTIDEK", 3),
                Input(Oxidized, 2)
            });

            Assert.That(predictor.Received.Select(p => p.FullSequence), Is.EqualTo(new[] { "PEPTIDEK", Oxidized }));
            Assert.That(predictor.Received.Select(p => p.BaseSequence), Is.EqualTo(new[] { "PEPTIDEK", "PEPTIDEMK" }));
        }

        [Test]
        public void UndetectablePeptidesAreLeftOutAndUnjudgedOnesAreKept()
        {
            // pFly takes no modifications, so it cannot judge the oxidized peptide.
            var detectability = new FakeDetectabilityModel(new Dictionary<string, double>
            {
                ["PEPTIDEK"] = 0.5,
                ["ELVISK"] = 0.51
            });
            var fragmentModel = new FakeFragmentModel();
            var generator = new SpectralLibraryGenerator(fragmentModel, detectabilityModel: detectability, detectabilityThreshold: 0.5);

            List<LibrarySpectrum> spectra = generator.GenerateLibrarySpectra(new List<FragmentIntensityPredictionInput>
            {
                Input("PEPTIDEK", 2),
                Input("ELVISK", 2),
                Input("ELVISK", 3),
                Input(Oxidized, 2)
            });

            Assert.That(detectability.Received, Is.EqualTo(new[] { "PEPTIDEK", "ELVISK" }));
            Assert.That(fragmentModel.Received, Is.EqualTo(new[] { "PEPTIDEK", "PEPTIDEM[UNIMOD:35]K" }));
            Assert.That(spectra.Select(s => s.Name), Is.EqualTo(new[] { "PEPTIDEK/2", Oxidized + "/2" }));
        }

        [Test]
        public void WhenEveryPeptideIsUndetectableNothingIsPredicted()
        {
            var detectability = new FakeDetectabilityModel(new Dictionary<string, double> { ["PEPTIDEK"] = 0.9 });
            var fragmentModel = new FakeFragmentModel();
            var predictor = new FakeRetentionTimePredictor(new Dictionary<string, double?>());
            var generator = new SpectralLibraryGenerator(fragmentModel, predictor, detectability);

            List<LibrarySpectrum> spectra = generator.GenerateLibrarySpectra(new List<FragmentIntensityPredictionInput> { Input("PEPTIDEK", 2) });

            Assert.That(spectra, Is.Empty);
            Assert.That(fragmentModel.Received, Is.Empty);
            Assert.That(predictor.Received, Is.Empty);
        }

        [Test]
        public void FragmentsAreKeptDownToTheModelsOwnIntensityFloor()
        {
            var generator = new SpectralLibraryGenerator(new FakeFragmentModel());

            LibrarySpectrum spectrum = generator.GenerateLibrarySpectra(new List<FragmentIntensityPredictionInput> { Input("PEPTIDEK", 2) }).Single();

            Assert.That(spectrum.MatchedFragmentIons.Select(ion => ion.Intensity), Is.EquivalentTo(FakeFragmentModel.Intensities));
        }

        [Test]
        public void ACapillaryElectrophoresisPredictorIsRefused()
        {
            using var cze = new CZERetentionTimePredictor();

            Assert.Throws<NotSupportedException>(() => new SpectralLibraryGenerator(new FakeFragmentModel(), cze));
        }

        [Test]
        public void TheWrittenMspHoldsWhatWasReturned()
        {
            string path = TemporaryPath(".msp");
            try
            {
                var generator = new SpectralLibraryGenerator(new FakeFragmentModel(), new FakeKoinaRetentionTimeModel());
                List<LibrarySpectrum> written = generator.GenerateAndWrite(
                    new List<FragmentIntensityPredictionInput> { Input("PEPTIDEK", 2), Input("ELVISK", 2) },
                    path, SpectralLibraryFormat.Msp);

                Assert.That(File.ReadAllText(path), Does.Contain("RT=30.5 FragmentModel=Prosit_2020_intensity_HCD RTModel=Prosit_2019_irt"));
                Assert.That(ReadBack(path), Is.EqualTo(written));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void TheWrittenMslHoldsWhatWasReturned()
        {
            string path = TemporaryPath(".msl");
            try
            {
                List<LibrarySpectrum> written = SpectralLibraryGenerator.WriteSpectralLibrary(
                    new List<LibrarySpectrum> { Spectrum("PEPTIDEK", 30.5, isDecoy: true) }, path, SpectralLibraryFormat.Msl);

                List<LibrarySpectrum> read = ReadBack(path);
                Assert.That(read.Select(s => s.Name), Is.EqualTo(written.Select(s => s.Name)));
                Assert.That(read.Single().MatchedFragmentIons.Count, Is.EqualTo(written.Single().MatchedFragmentIons.Count));
                Assert.That(read.Single().IsDecoy, Is.True);
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void FragmentFiltersApplyInOrderAndLeaveTheGivenSpectraAlone()
        {
            // Ion m/z 100..500 with intensities 0.1, 1.0, 0.4, 0.05, 0.6.
            LibrarySpectrum spectrum = Spectrum("PEPTIDEK", 30.5);
            spectrum.FragmentModelName = "Fragments";
            spectrum.RetentionTimeModelName = "Times";
            string path = TemporaryPath(".msp");
            try
            {
                List<LibrarySpectrum> written = SpectralLibraryGenerator.WriteSpectralLibrary(
                    new List<LibrarySpectrum> { spectrum }, path, SpectralLibraryFormat.Msp,
                    minMz: 150, maxMz: 450, minRelativeIntensity: 0.1, maxFragmentCount: 1);

                // m/z keeps 200, 300, 400; 0.1 of the strongest (1.0) drops 400 (0.05); the top one is 200.
                LibrarySpectrum filtered = written.Single();
                Assert.That(filtered.MatchedFragmentIons.Select(ion => ion.Mz), Is.EqualTo(new[] { 200.0 }));
                Assert.That(filtered.RetentionTime, Is.EqualTo(30.5));
                Assert.That(filtered.FragmentModelName, Is.EqualTo("Fragments"));
                Assert.That(filtered.RetentionTimeModelName, Is.EqualTo("Times"));
                Assert.That(spectrum.MatchedFragmentIons.Count, Is.EqualTo(5));
                Assert.That(ReadBack(path), Is.EqualTo(written));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [TestCase(200, null, null, null, new[] { 200.0, 300, 400, 500 })]
        [TestCase(null, 300, null, null, new[] { 100.0, 200, 300 })]
        [TestCase(null, null, 0.4, null, new[] { 200.0, 300, 500 })]
        [TestCase(null, null, null, 3, new[] { 200.0, 300, 500 })]
        [TestCase(250, null, 0.5, null, new[] { 500.0 })]
        public void EachFragmentFilterWorksAlone(double? minMz, double? maxMz, double? minRelativeIntensity, int? maxFragmentCount, double[] keptMz)
        {
            string path = TemporaryPath(".msp");
            try
            {
                List<LibrarySpectrum> written = SpectralLibraryGenerator.WriteSpectralLibrary(
                    new List<LibrarySpectrum> { Spectrum("PEPTIDEK", null, intensityScale: 2) }, path, SpectralLibraryFormat.Msp,
                    minMz, maxMz, minRelativeIntensity, maxFragmentCount);

                Assert.That(written.Single().MatchedFragmentIons.Select(ion => ion.Mz), Is.EqualTo(keptMz));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [Test]
        public void ASpectrumLeftWithNoFragmentsIsNotWritten()
        {
            string path = TemporaryPath(".msp");
            try
            {
                List<LibrarySpectrum> written = SpectralLibraryGenerator.WriteSpectralLibrary(
                    new List<LibrarySpectrum>
                    {
                        Spectrum("PEPTIDEK", null),
                        new LibrarySpectrum("ELVISK", 350.2, 2, new List<MatchedFragmentIon>(), null)
                    },
                    path, SpectralLibraryFormat.Msp);

                Assert.That(written.Select(s => s.Name), Is.EqualTo(new[] { "PEPTIDEK/2" }));
                Assert.That(ReadBack(path).Select(s => s.Name), Is.EqualTo(new[] { "PEPTIDEK/2" }));
            }
            finally
            {
                File.Delete(path);
            }
        }

        [TestCase(SpectralLibraryFormat.Msp, ".msp")]
        [TestCase(SpectralLibraryFormat.Msl, ".msl")]
        public void NoInputsWriteAnEmptyLibrary(SpectralLibraryFormat format, string extension)
        {
            string path = TemporaryPath(extension);
            try
            {
                var generator = new SpectralLibraryGenerator(new FakeFragmentModel(), new FakeKoinaRetentionTimeModel());

                List<LibrarySpectrum> written = generator.GenerateAndWrite(new List<FragmentIntensityPredictionInput>(), path, format);

                Assert.That(written, Is.Empty);
                Assert.That(ReadBack(path), Is.Empty);
            }
            finally
            {
                File.Delete(path);
            }
        }

        private static FragmentIntensityPredictionInput Input(string sequence, int charge) =>
            new(sequence, charge, 30, null, null);

        private static LibrarySpectrum Spectrum(string sequence, double? retentionTime, bool isDecoy = false, double intensityScale = 1)
        {
            double[] intensities = { 0.1, 1.0, 0.4, 0.05, 0.6 };
            var ions = intensities.Select((intensity, i) => new MatchedFragmentIon(
                new Product(ProductType.y, FragmentationTerminus.C, (100.0 * (i + 1)).ToMass(1), i + 1, sequence.Length - i, 0),
                100 * (i + 1), intensity * intensityScale, 1)).ToList();
            return new LibrarySpectrum(sequence, 450.7, 2, ions, retentionTime, isDecoy);
        }

        private static string TemporaryPath(string extension) =>
            Path.Combine(Path.GetTempPath(), "SpectralLibraryGeneratorTests_" + Guid.NewGuid().ToString("N") + extension);

        private static List<LibrarySpectrum> ReadBack(string path)
        {
            var library = new SpectralLibraryReader(new List<string> { path });
            try
            {
                return library.GetAllLibrarySpectra().ToList();
            }
            finally
            {
                library.CloseConnections();
            }
        }

        private static List<string> RequestedSequences(Dictionary<string, object> request) =>
            JObject.FromObject(request)["inputs"]!
                .First(field => (string)field["name"]! == "peptide_sequences")["data"]!
                .Select(sequence => (string)sequence!)
                .ToList();

        private static string Output(string name, int rows, IEnumerable<object> data) =>
            JsonConvert.SerializeObject(new { name, datatype = name == "annotation" ? "BYTES" : "FP32", shape = new[] { rows, data.Count() / rows }, data });

        private static string Response(params string[] outputs) => "{\"outputs\":[" + string.Join(",", outputs) + "]}";

        private sealed class FakeFragmentModel : Prosit2020IntensityHCD
        {
            public static readonly double[] Intensities = { 0.2, 1.0, 0.5, 0.00005 };
            private static readonly string[] Annotations = { "b2+1", "y1+1", "y2+1", "y3+1" };

            public List<string> Received { get; } = new();

            public FakeFragmentModel() : base(fragmentIonMappingMode: FragmentIonMappingMode.MapToInputFullSequence) { }

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
            {
                List<string> sequences = RequestedSequences(request);
                Received.AddRange(sequences);
                int n = sequences.Count;
                return Task.FromResult(Response(
                    Output("annotation", n, Enumerable.Repeat(Annotations, n).SelectMany(a => a)),
                    Output("mz", n, Enumerable.Repeat(new object[] { 100.0, 200.0, 300.0, 400.0 }, n).SelectMany(a => a)),
                    Output("intensities", n, Enumerable.Repeat(Intensities, n).SelectMany(a => a).Cast<object>())));
            }
        }

        private sealed class FakeKoinaRetentionTimeModel : Prosit2019iRT
        {
            private static readonly Dictionary<string, double> Values = new() { ["PEPTIDEK"] = 30.5, ["ELVISK"] = 12.25 };

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
            {
                List<string> sequences = RequestedSequences(request);
                return Task.FromResult(Response(Output("irt", sequences.Count, sequences.Select(s => (object)Values[s]))));
            }
        }

        private sealed class FakeDetectabilityModel : PFly2024FineTuned
        {
            private readonly Dictionary<string, double> _notDetectable;

            public List<string> Received { get; } = new();

            public FakeDetectabilityModel(Dictionary<string, double> notDetectable) => _notDetectable = notDetectable;

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken ct)
            {
                List<string> sequences = RequestedSequences(request);
                Received.AddRange(sequences);
                var data = sequences.SelectMany(s => new object[] { _notDetectable[s], 0.2 * (1 - _notDetectable[s]), 0.3 * (1 - _notDetectable[s]), 0.5 * (1 - _notDetectable[s]) });
                return Task.FromResult(Response(Output("detectability", sequences.Count, data)));
            }
        }

        private sealed class FakeRetentionTimePredictor : IRetentionTimePredictor
        {
            private readonly Dictionary<string, double?> _values;

            public List<IRetentionPredictable> Received { get; } = new();

            public FakeRetentionTimePredictor(Dictionary<string, double?> values) => _values = values;

            public string PredictorName => "FakePredictor";
            public SeparationType SeparationType => SeparationType.HPLC;

            public IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)>
                PredictRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
            {
                var list = peptides.ToList();
                Received.AddRange(list);
                return list.Select(p => (_values[p.FullSequence], p, (RetentionTimeFailureReason?)null)).ToList();
            }

            public double? PredictRetentionTimeEquivalent(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason) =>
                throw new NotImplementedException();

            public string? GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason) =>
                throw new NotImplementedException();

            public void Dispose() { }
        }
    }
}
