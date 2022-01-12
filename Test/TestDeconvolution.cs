using Chemistry;
using FlashLFQ;
using IO.ThermoRawFileReader;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;

namespace Test
{
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public sealed class TestDeconvolution
    {

        [Test]
        public static void CheckGetMostAbundantObservedIsotopicMass()
        {
            
            string fullFilePathWithExtension = @"D:\TDBU\Jurkat\TD-Projects-JurkatTopDownSeanDaiPaper\FXN11_tr1_032017.raw"; //@"D:\Other\averageMass\for_shortreed_protein_mix_sid 15_topdown_standard-qb.raw";
            ThermoRawFileReader staticRaw = ThermoRawFileReader.LoadAllStaticData(fullFilePathWithExtension);
            List<MsDataScan> scan = staticRaw.GetAllScansList(); //.Where(p => p.RetentionTime == 35.88835).ToList();
            scan = scan.Where(p => p.MsnOrder == 1).ToList();
            scan = scan.Where(p => p.OneBasedScanNumber == 587).ToList();

            //scan = scan.Where(p => p.RetentionTime < 35.888).Where(p => p.RetentionTime > 35.8).ToList();
            MzSpectrum spec = scan[0].MassSpectrum;
            MzRange theRange = new MzRange(spec.XArray.Min(), spec.XArray.Max());
            int minAssumedChargeState = 1;
            int maxAssumedChargeState = 60;
            double deconvolutionTolerancePpm = 20;
            double intensityRatioLimit = 3;

            

            List<MassSpectrometry.IsotopicEnvelope> lie = spec.Deconvolute(theRange, minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatioLimit).ToList();

            //List<MassSpectrometry.IsotopicEnvelope> lie = scan[3000].GetIsolatedMassesAndCharges(spec, 1, 60, 4, 3).ToList();
            lie = lie.Where(p => p.MostAbundantObservedIsotopicMass < 8246.5+1).OrderByDescending(p => p.MostAbundantObservedIsotopicMass).ToList();
            var compare = lie.Select(p => p.MostAbundantObservedIsotopicMass).ToList()[0];

            Assert.That(compare, Is.EqualTo(8246.5).Within(1));


            //List<IsotopicEnvelope>

            //List<MassSpectrometry.IsotopicEnvelope> lie = spec.Deconvolute(theRange, minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatioLimit).ToList();

            //Deconvolute(int ? minScan, int ? maxScan, int minAssumedChargeState, int maxAssumedChargeState, double deconvolutionTolerancePpm, double intensityRatioLimit, double aggregationTolerancePpm, Func<MsDataScan, bool> scanFilterFunc, int maxThreads = -1)
            /*
            ProteinGroup pg = new ProteinGroup("Accession", "Gene", "Organism");
            Peptide pep1 = new Peptide("AAAAAA", "AAAAAA", true, new HashSet<ProteinGroup> { pg });
            Peptide pep2 = new Peptide("AAA[H]AAA", "AAAAAA", true, new HashSet<ProteinGroup> { pg });

            var dist1 = IsotopicDistribution.GetDistribution(pep1.GetChemicalFormula(), 0.1, 0.01);

            var dist2 = IsotopicDistribution.GetDistribution(pep2.GetChemicalFormula(), 0.1, 0.01);

            MsDataScan[] Scans = new MsDataScan[2];
            double[] ms1intensities = new double[] { 0.8, 0.8, 0.2, 0.02, 0.2, 0.02 };
            double[] ms1mzs = dist1.Masses.Concat(dist2.Masses).OrderBy(b => b).Select(b => b.ToMz(1)).ToArray();

            double selectedIonMz = ms1mzs[1];

            MzSpectrum MS1 = new MzSpectrum(ms1mzs, ms1intensities, false);

            Scans[0] = new MsDataScan(MS1, 1, 1, false, Polarity.Positive, 1.0, new MzRange(300, 2000), "first spectrum", MZAnalyzerType.Unknown, MS1.SumOfAllY, null, null, null);

            // Horrible fragmentation, but we don't care about this!
            double[] ms2intensities = new double[] { 1000 };
            double[] ms2mzs = new double[] { 1000 };
            MzSpectrum MS2 = new MzSpectrum(ms2mzs, ms2intensities, false);
            double isolationMZ = selectedIonMz;
            Scans[1] = new MsDataScan(MS2, 2, 2, false, Polarity.Positive, 2.0, new MzRange(100, 1500), "second spectrum", MZAnalyzerType.Unknown, MS2.SumOfAllY, null, null, null, selectedIonMz, null, null, isolationMZ, 2.5, DissociationType.HCD, 1, null);

            var myMsDataFile = new FakeMsDataFile(Scans);

            var cool = myMsDataFile.GetAllScansList().Last();

            int maxAssumedChargeState = 1;
            Tolerance massTolerance = Tolerance.ParseToleranceString("10 PPM");

            var isolatedMasses = cool.GetIsolatedMassesAndCharges(myMsDataFile.GetOneBasedScan(cool.OneBasedPrecursorScanNumber.Value).MassSpectrum, 1, maxAssumedChargeState, 10, 5).ToList();
            */

            //Assert.AreEqual(lie[0].MostAbundantObservedIsotopicMass,);

            //Assert.That(lie[0].MostAbundantObservedIsotopicMass, Is.EqualTo(9112.40295).Within(0.1));
        }



        [Test]
        [TestCase(586.2143122,24, 41983672, 586.2)]//This is a lesser abundant charge state envelope at the low mz end
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

            List<MassSpectrometry.IsotopicEnvelope> isolatedMasses = scan.GetIsolatedMassesAndCharges(spectrum, 1, 60, 4, 3).ToList();

            List<double> monoIsotopicMasses = isolatedMasses.Select(m => m.MonoisotopicMass).ToList();

            //The primary monoisotopic mass should be the same regardless of which peak in which charge state was selected for isolation.
            //this case is interesting because other monoisotopic mass may have a sodium adduct. The unit test could be expanded to consider this.
            Assert.That(monoIsotopicMasses[0], Is.EqualTo(14037.926829).Within(.0005));
        }
    }
}