using Chemistry;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using MassSpectrometry.Proteomics;
using MassSpectrometry.Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using Easy.Common.Extensions;
using MassSpectrometry.Deconvolution;
using MassSpectrometry.Deconvolution.Algorithms;
using MassSpectrometry.Deconvolution.Parameters;
using TopDownProteomics.MassSpectrometry;
using IsotopicDistribution = Chemistry.IsotopicDistribution;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public sealed class TestDeconvolution
    {

        #region Old Deconvolution

        [Test]
        [TestCase(586.2143122, 24, 41983672, 586.2)]//This is a lesser abundant charge state envelope at the low mz end
        [TestCase(740.372202090153, 19, 108419280, 740.37)]//This is the most abundant charge state envelope
        [TestCase(1081.385183, 13, 35454636, 1081.385)]//This is a lesser abundant charge state envelope at the high mz end
        public void TestDeconvolutionProteoformMultiChargeState(double selectedIonMz, int selectedIonChargeStateGuess, double selectedIonIntensity, double isolationMz)
        {
            MsDataScan[] Scans = new MsDataScan[1];

            //txt file, not mgf, because it's an MS1. Most intense proteoform has mass of ~14037.9 Da
            string Ms1SpectrumPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\14kDaProteoformMzIntensityMs1.txt");

            string[] spectrumLines = File.ReadAllLines(Ms1SpectrumPath);

            int mzIntensityPairsCount = spectrumLines.Length;
            double[] ms1mzs = new double[mzIntensityPairsCount];
            double[] ms1intensities = new double[mzIntensityPairsCount];

            for (int i = 0; i < mzIntensityPairsCount; i++)
            {
                string[] pair = spectrumLines[i].Split('\t');
                ms1mzs[i] = Convert.ToDouble(pair[0], CultureInfo.InvariantCulture);
                ms1intensities[i] = Convert.ToDouble(pair[1], CultureInfo.InvariantCulture);
            }

            MzSpectrum spectrum = new MzSpectrum(ms1mzs, ms1intensities, false);

            Scans[0] = new MsDataScan(spectrum, 1, 1, false, Polarity.Positive, 1.0, new MzRange(495, 1617), "first spectrum", MZAnalyzerType.Unknown, spectrum.SumOfAllY, null, null, null, selectedIonMz, selectedIonChargeStateGuess, selectedIonIntensity, isolationMz, 4);

            var myMsDataFile = new FakeMsDataFile(Scans);

            MsDataScan scan = myMsDataFile.GetAllScansList()[0];

            List<IsotopicEnvelope> isolatedMasses = scan.GetIsolatedMassesAndCharges(spectrum, 1, 60, 4, 3).ToList();

            List<double> monoIsotopicMasses = isolatedMasses.Select(m => m.MonoisotopicMass).ToList();

            //The primary monoisotopic mass should be the same regardless of which peak in which charge state was selected for isolation.
            //this case is interesting because other monoisotopic mass may have a sodium adduct. The unit test could be expanded to consider this.
            Assert.That(monoIsotopicMasses[0], Is.EqualTo(14037.926829).Within(.0005));
        }

        [Test]
        [TestCase("APSGGKK", "12-18-17_frac7_calib_ms1_663_665.mzML", 2)]
        [TestCase("PKRKAEGDAKGDKAKVKDEPQRRSARLSAKPAPPKPEPKPKKAPAKKGEKVPKGKKGKADAGKEGNNPAENGDAKTDQAQKAEGAGDAK", "FXN11_tr1_032017-calib_ms1_scans716_718.mzML", 8)]
        [TestCase("PKRKVSSAEGAAKEEPKRRSARLSAKPPAKVEAKPKKAAAKDKSSDKKVQTKGKRGAKGKQAEVANQETKEDLPAENGETKTEESPASDEAGEKEAKSD", "FXN11_tr1_032017-calib_ms1_scans781_783.mzML", 16)]
        public static void CheckGetMostAbundantObservedIsotopicMass(string peptide, string file, int charge)
        {
            Protein test1 = new Protein(peptide, "Accession");
            DigestionParams d = new DigestionParams();
            PeptideWithSetModifications pw = new PeptideWithSetModifications(test1, d, 1, test1.Length, CleavageSpecificity.None, "", 0, new Dictionary<int, Modification>(), 0);
            double mostAbundantMz = pw.MostAbundantMass.ToMz(charge);

            string singleScan = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", file);
            Mzml singleMZML = Mzml.LoadAllStaticData(singleScan);

            List<MsDataScan> singlescan = singleMZML.GetAllScansList();

            MzSpectrum singlespec = singlescan[0].MassSpectrum;
            MzRange singleRange = new MzRange(singlespec.XArray.Min(), singlespec.XArray.Max());
            int minAssumedChargeState = 1;
            int maxAssumedChargeState = 60;
            double deconvolutionTolerancePpm = 20;
            double intensityRatioLimit = 3;

            //check assigned correctly
            List<IsotopicEnvelope> lie2 = singlespec.Deconvolute(singleRange, minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatioLimit).ToList();
            List<IsotopicEnvelope> lie2_charge = lie2.Where(p => p.Charge == charge).ToList();
            Assert.That(lie2_charge[0].MostAbundantObservedIsotopicMz, Is.EqualTo(mostAbundantMz).Within(0.1));

            //check that if already assigned, skips assignment and just recalls same value
            List<IsotopicEnvelope> lie3 = singlespec.Deconvolute(singleRange, minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatioLimit).ToList();
            Assert.AreEqual(lie2.Select(p => p.MostAbundantObservedIsotopicMass), lie3.Select(p => p.MostAbundantObservedIsotopicMass));
        }

        #endregion

        #region Classic Deconvolution

        [Test]
        [TestCase(586.2143122, 24, 41983672, 586.2)]//This is a lesser abundant charge state envelope at the low mz end
        [TestCase(740.372202090153, 19, 108419280, 740.37)]//This is the most abundant charge state envelope
        [TestCase(1081.385183, 13, 35454636, 1081.385)]//This is a lesser abundant charge state envelope at the high mz end
        public void TestClassicDeconvolutionProteoformMultiChargeState(double selectedIonMz, int selectedIonChargeStateGuess, double selectedIonIntensity, double isolationMz)
        {
            MsDataScan[] Scans = new MsDataScan[1];

            //txt file, not mgf, because it's an MS1. Most intense proteoform has mass of ~14037.9 Da
            string Ms1SpectrumPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\14kDaProteoformMzIntensityMs1.txt");

            string[] spectrumLines = File.ReadAllLines(Ms1SpectrumPath);

            int mzIntensityPairsCount = spectrumLines.Length;
            double[] ms1mzs = new double[mzIntensityPairsCount];
            double[] ms1intensities = new double[mzIntensityPairsCount];

            for (int i = 0; i < mzIntensityPairsCount; i++)
            {
                string[] pair = spectrumLines[i].Split('\t');
                ms1mzs[i] = Convert.ToDouble(pair[0], CultureInfo.InvariantCulture);
                ms1intensities[i] = Convert.ToDouble(pair[1], CultureInfo.InvariantCulture);
            }

            MzSpectrum spectrum = new MzSpectrum(ms1mzs, ms1intensities, false);

            Scans[0] = new MsDataScan(spectrum, 1, 1, false, Polarity.Positive, 1.0, new MzRange(495, 1617), "first spectrum", MZAnalyzerType.Unknown, spectrum.SumOfAllY, null, null, null, selectedIonMz, selectedIonChargeStateGuess, selectedIonIntensity, isolationMz, 4);

            var myMsDataFile = new FakeMsDataFile(Scans);

            MsDataScan scan = myMsDataFile.GetAllScansList()[0];

            // The ones marked 2 are for checking an overload method

            DeconvolutionParameters deconParameters = new ClassicDeconvolutionParameters(1, 60, 4, 3);
            Deconvoluter deconvoluter = new Deconvoluter(DeconvolutionTypes.ClassicDeconvolution, deconParameters);

            List<IsotopicEnvelope> isolatedMasses = scan.GetIsolatedMassesAndCharges(DeconvolutionTypes.ClassicDeconvolution, deconParameters).ToList();
            List<IsotopicEnvelope> isolatedMasses2 = scan.GetIsolatedMassesAndCharges(deconvoluter).ToList();

            List<double> monoIsotopicMasses = isolatedMasses.Select(m => m.MonoisotopicMass).ToList();
            List<double> monoIsotopicMasses2 = isolatedMasses2.Select(m => m.MonoisotopicMass).ToList();

            //The primary monoisotopic mass should be the same regardless of which peak in which charge state was selected for isolation.
            //this case is interesting because other monoisotopic mass may have a sodium adduct. The unit test could be expanded to consider this.
            Assert.That(monoIsotopicMasses[0], Is.EqualTo(14037.926829).Within(.0005));
            Assert.That(monoIsotopicMasses2[0], Is.EqualTo(14037.926829).Within(.0005));
        }

        [Test]
        [TestCase("APSGGKK", "12-18-17_frac7_calib_ms1_663_665.mzML", 2)]
        [TestCase("PKRKAEGDAKGDKAKVKDEPQRRSARLSAKPAPPKPEPKPKKAPAKKGEKVPKGKKGKADAGKEGNNPAENGDAKTDQAQKAEGAGDAK", "FXN11_tr1_032017-calib_ms1_scans716_718.mzML", 8)]
        [TestCase("PKRKVSSAEGAAKEEPKRRSARLSAKPPAKVEAKPKKAAAKDKSSDKKVQTKGKRGAKGKQAEVANQETKEDLPAENGETKTEESPASDEAGEKEAKSD", "FXN11_tr1_032017-calib_ms1_scans781_783.mzML", 16)]
        public static void CheckClassicGetMostAbundantObservedIsotopicMass(string peptide, string file, int charge)
        {
            Protein test1 = new Protein(peptide, "Accession");
            DigestionParams d = new DigestionParams();
            PeptideWithSetModifications pw = new PeptideWithSetModifications(test1, d, 1, test1.Length, CleavageSpecificity.None, "", 0, new Dictionary<int, Modification>(), 0);
            double pwsmMonoisotopicMass = pw.MostAbundantMass;

            string singleScan = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", file);
            Mzml singleMZML = Mzml.LoadAllStaticData(singleScan);

            List<MsDataScan> singlescan = singleMZML.GetAllScansList();

            MzSpectrum singlespec = singlescan[0].MassSpectrum;
            MzRange singleRange = new MzRange(singlespec.XArray.Min(), singlespec.XArray.Max());
            int minAssumedChargeState = 1;
            int maxAssumedChargeState = 60;
            double deconvolutionTolerancePpm = 20;
            double intensityRatioLimit = 3;

            DeconvolutionParameters deconParameters =
                new ClassicDeconvolutionParameters(minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm,
                    intensityRatioLimit);
            Deconvoluter deconvoluter = new Deconvoluter(DeconvolutionTypes.ClassicDeconvolution, deconParameters);

            //check assigned correctly

            List<IsotopicEnvelope> lie2 = deconvoluter.ClassicDeconvoluteMzSpectra(singlespec, singleRange).ToList();

            List<IsotopicEnvelope> lie2_charge = lie2.Where(p => p.Charge == charge).ToList();
            Assert.That(lie2_charge[0].MostAbundantObservedIsotopicMass, Is.EqualTo(pwsmMonoisotopicMass).Within(0.05));

            //check that if already assigned, skips assignment and just recalls same value
            List<IsotopicEnvelope> lie3 = deconvoluter.ClassicDeconvoluteMzSpectra(singlespec, singleRange).ToList();
            Assert.AreEqual(lie2.Select(p => p.MostAbundantObservedIsotopicMass), lie3.Select(p => p.MostAbundantObservedIsotopicMass));
        }

        #endregion

        #region SpectralDecon

        [Test]
        public static void TestGetSecondMostAbundantSpecies()
        {
            Protein testProtein = new Protein("PEPTIDEFPEPTIDEK", "Accession");
            DigestionParams digestionParams = new DigestionParams();
            PeptideWithSetModifications pwsm = new PeptideWithSetModifications(testProtein, digestionParams, 1, testProtein.Length, CleavageSpecificity.None, "", 0, new Dictionary<int, Modification>(), 0);
            IsotopicDistribution isotopicDistribution = IsotopicDistribution.GetDistribution(pwsm.FullChemicalFormula, 
                fineResolution: 0.125, minProbability: 0.001);

            int charge = 1;
            IsotopicEnvelope testEnvelope = new IsotopicEnvelope(isotopicDistribution, charge);
            Assert.That(testEnvelope.MostAbundantObservedIsotopicMz, Is.EqualTo(pwsm.MonoisotopicMass.ToMz(charge)).Within(0.1));
            Assert.That(testEnvelope.SecondMostAbundantObservedIsotopicMz, Is.EqualTo(isotopicDistribution.Masses[1].ToMz(charge)).Within(0.1));
        }

        [Test]

        public void TestIndexingForSpectralDecon()
        {
            //PEPTIDEK vs PEPTIDEFPEPTIDEK (the longer peptide has ~ twice the mass of the shorter, enabling a test of the indexing system
            Protein myProtein = new Protein("PEPTIDEKPEPTIDEFPEPTIDEK", "accession");
            DigestionParams digest1 = new DigestionParams(protease: "trypsin", maxMissedCleavages: 0, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain);

            List<PeptideWithSetModifications> pep = myProtein.Digest(digest1, new List<Modification>(), new List<Modification>()).ToList();

            int minAssumedChargeState = 1;
            int maxAssumedChargeState = 60;
            double deconvolutionTolerancePpm = 20;
            int binsPerDalton = 1;
            int scanMinimum = 460;

            SpectralDeconvolutionParameters spectralDeconParams = new SpectralDeconvolutionParameters(
                minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm,
                new List<Protein>() { myProtein },
                new List<Modification>(), new List<Modification>(), digest1,
                new List<SilacLabel>(), false, scanMinimumMz: scanMinimum, scanMaximumMz: 2000,
                ambiguityThresholdForIsotopicDistribution: 0.9, binsPerDalton: binsPerDalton);

            SpectralDeconvolutionAlgorithm spectralDecon = new SpectralDeconvolutionAlgorithm(spectralDeconParams);

            PeptideWithSetModifications peptidek = pep.Where(p => p.FullSequence.Equals("PEPTIDEK")).First();
            PeptideWithSetModifications doublePeptidek = pep.Where(p => p.FullSequence.Equals("PEPTIDEFPEPTIDEK")).First();
            Assert.That(spectralDecon.EnvelopeDictionary.ContainsKey(peptidek) & spectralDecon.EnvelopeDictionary.ContainsKey(doublePeptidek));
            Assert.That(spectralDecon.EnvelopeDictionary[peptidek].Count == 2 & spectralDecon.EnvelopeDictionary[doublePeptidek].Count == 4);

            List<List<MinimalSpectrum>> indexedSpectra = new();
            //Iterate over 2d Array
            foreach (var listOfSpectra in spectralDecon.IndexedLibrarySpectra)
            {
                if (listOfSpectra.IsNotNullOrEmpty()) indexedSpectra.Add(listOfSpectra);
            }
            
            // For the longer peptide, the first and second isotopes have extremely similar abundances,
            // so they should be stored in different bins for the +1, +2, and +4 charge states (the +3 charge state masses fall within the same bin [619 Thompsons])
            // The shorter peptide is approximately 1/2 the mass of the longer peptide. With bin sizes of one dalton,
            // every spectra for the shorter peptide should share a bin with a peptide from a longer spectra.
            // This assertion does a lot of heavy lifting in testing the indexing engine. 
            // DO NOT CHANGE unless you understand what is being tested here
            //Assert.That(indexedSpectra.Count == 7);

            int peptideCharge = 2;
            int binIndex = (int)Math.Floor(binsPerDalton * (peptidek.MonoisotopicMass.ToMz(charge: peptideCharge) - scanMinimum));
            int chargeIndex = peptideCharge - minAssumedChargeState;
            Assert.That(spectralDecon.SpectrumIndexToPwsmMap.TryGetValue((binIndex, chargeIndex, 0), out var peptidek2Charge));
            Assert.That(spectralDecon.SpectrumIndexToPwsmMap.TryGetValue((binIndex, 3, 0), out var doublePeptideK4Charge));
            Assert.That(peptidek2Charge.BaseSequence.Equals(peptidek.BaseSequence));
            Assert.That(doublePeptideK4Charge.BaseSequence.Equals(doublePeptidek.BaseSequence));

            peptideCharge = 1;
            chargeIndex = peptideCharge - minAssumedChargeState;
            binIndex = (int)Math.Floor(binsPerDalton * (doublePeptidek.MonoisotopicMass.ToMz(charge: 1) - scanMinimum));
            Assert.That(spectralDecon.SpectrumIndexToPwsmMap.TryGetValue((binIndex, chargeIndex, 0), out var doublePeptideK1Charge) && 
                        !spectralDecon.SpectrumIndexToPwsmMap.TryGetValue((binIndex, chargeIndex, 1), out var doesNotExist));
            Assert.That(doublePeptideK1Charge == doublePeptideK4Charge);

        }

        [Test]
        [TestCase("APSGGKK", "12-18-17_frac7_calib_ms1_663_665.mzML", 2)]
        [TestCase("PKRKAEGDAKGDKAKVKDEPQRRSARLSAKPAPPKPEPKPKKAPAKKGEKVPKGKKGKADAGKEGNNPAENGDAKTDQAQKAEGAGDAK", "FXN11_tr1_032017-calib_ms1_scans716_718.mzML", 8)]
        [TestCase("PKRKVSSAEGAAKEEPKRRSARLSAKPPAKVEAKPKKAAAKDKSSDKKVQTKGKRGAKGKQAEVANQETKEDLPAENGETKTEESPASDEAGEKEAKSD", "FXN11_tr1_032017-calib_ms1_scans781_783.mzML", 16)]
        public static void CheckSpectralGetMostAbundantObservedIsotopicMass(string peptide, string file, int charge)
        {
            Protein myProtein = new Protein(peptide, "Accession");
            DigestionParams digest1 = new DigestionParams("top-down");
            PeptideWithSetModifications pw = new PeptideWithSetModifications(myProtein, digest1, 1, myProtein.Length, CleavageSpecificity.None, "", 0, new Dictionary<int, Modification>(), 0);
            double pwsmMonoisotopicMass = pw.MostAbundantMass;

            string singleScan = Path.Combine(TestContext.CurrentContext.TestDirectory, "DataFiles", file);
            Mzml singleMZML = Mzml.LoadAllStaticData(singleScan);

            List<MsDataScan> singlescan = singleMZML.GetAllScansList();

            MzSpectrum singlespec = singlescan[0].MassSpectrum;
            MzRange singleRange = new MzRange(singlespec.XArray.Min(), singlespec.XArray.Max());


            int minAssumedChargeState = 1;
            int maxAssumedChargeState = 60;
            double deconvolutionTolerancePpm = 20;
            int binsPerDalton = 1;
            int scanMinimum = 460;

            DeconvolutionParameters deconParams = new SpectralDeconvolutionParameters(
                minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm,
                new List<Protein>() { myProtein },
                new List<Modification>(), new List<Modification>(), digest1,
                new List<SilacLabel>(), false, scanMinimumMz: singleRange.Minimum, scanMaximumMz: singleRange.Maximum,
                ambiguityThresholdForIsotopicDistribution: 0.9, binsPerDalton: binsPerDalton);

            Deconvoluter deconvoluter = new Deconvoluter(DeconvolutionTypes.SpectralDeconvolution, deconParams);

            //check assigned correctly

            List<IsotopicEnvelope> lie2 = deconvoluter.SpectralDeconvoluteMzSpectra(singlespec, singleRange).ToList();

            List<IsotopicEnvelope> lie2_charge = lie2.Where(p => p.Charge == charge).ToList();
            Assert.That(lie2_charge[0].MostAbundantObservedIsotopicMass, Is.EqualTo(pwsmMonoisotopicMass).Within(0.05));

            //check that if already assigned, skips assignment and just recalls same value
            List<IsotopicEnvelope> lie3 = deconvoluter.ClassicDeconvoluteMzSpectra(singlespec, singleRange).ToList();
            Assert.AreEqual(lie2.Select(p => p.MostAbundantObservedIsotopicMass), lie3.Select(p => p.MostAbundantObservedIsotopicMass));
        }
        #endregion
    }
}