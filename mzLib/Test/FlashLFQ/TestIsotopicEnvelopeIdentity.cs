using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using IsotopicEnvelope = FlashLFQ.IsotopicEnvelope;

namespace Test.FlashLFQ
{
    /// <summary>
    /// GetIsotopicEnvelopes and CheckIsotopicEnvelopeCorrelation are FlashLFQ's innermost loop. They were rewritten to
    /// stop allocating (LINQ projections, a list that was filled and never read, MathNet's IEnumerable Pearson) and must
    /// return exactly what they did before. These tests hold them against verbatim copies of the previous code.
    /// </summary>
    [TestFixture]
    public static class TestIsotopicEnvelopeIdentity
    {
        private static string DataDirectory => Path.Combine(TestContext.CurrentContext.TestDirectory, "FlashLFQ", "TestData");

        [Test]
        public static void PearsonMatchesMathNetBitForBit()
        {
            var random = new Random(11);
            for (int trial = 0; trial < 200_000; trial++)
            {
                int length = random.Next(0, 10);
                var points = new List<(double expIntensity, double theorIntensity, double theorMass)>();
                for (int i = 0; i < length; i++)
                {
                    double exp = random.Next(5) == 0 ? 0 : Math.Exp(random.NextDouble() * 25);
                    double theor = random.Next(6) == 0 ? 0 : Math.Exp(random.NextDouble() * 25);
                    if (random.Next(20) == 0) theor = exp; // perfectly correlated runs
                    points.Add((exp, theor, 1000 + i));
                }

                double expected = Correlation.Pearson(points.Select(p => p.expIntensity), points.Select(p => p.theorIntensity));
                double actual = FlashLfqEngine.PearsonCorrelationOfIntensities(points);
                Assert.That(BitConverter.DoubleToInt64Bits(actual), Is.EqualTo(BitConverter.DoubleToInt64Bits(expected)),
                    $"trial {trial}: {string.Join(", ", points)}");
            }
        }

        [Test]
        public static void MinTheoreticalMassMatchesEnumerableMin()
        {
            var random = new Random(5);
            for (int trial = 0; trial < 20_000; trial++)
            {
                var points = Enumerable.Range(0, random.Next(1, 8))
                    .Select(_ => (0.0, 0.0, random.Next(30) == 0 ? double.NaN : 500 + random.NextDouble() * 3000))
                    .ToList();

                double expected = points.Min(p => p.Item3);
                double actual = FlashLfqEngine.MinTheoreticalMass(points);
                Assert.That(BitConverter.DoubleToInt64Bits(actual), Is.EqualTo(BitConverter.DoubleToInt64Bits(expected)));
            }
        }

        /// <summary>
        /// Every identification in a real file, at every charge state, on its own XIC and on XICs shifted by one isotope
        /// either way (which exercise the off-by-one branches), must give the same envelopes as the previous code: the same
        /// peak objects, charge states, and bitwise-identical intensities and correlations.
        /// </summary>
        [Test]
        public static void GetIsotopicEnvelopesMatchesThePreviousCode()
        {
            var file = new SpectraFileInfo(Path.Combine(DataDirectory, "20100614_Velos1_TaGe_SA_K562_3.mzML"), "a", 0, 0, 0);
            var ids = new List<Identification>();
            foreach (var psm in SpectrumMatchTsvReader.ReadPsmTsv(Path.Combine(DataDirectory, "AllPSMs.psmtsv"), out _))
            {
                if (psm.FileNameWithoutExtension != file.FilenameWithoutExtension || psm.FullSequence.Contains('|')) continue;
                ids.Add(new Identification(file, psm.BaseSeq, psm.FullSequence, (double)psm.MonoisotopicMass, (double)psm.RetentionTime,
                    psm.PrecursorCharge, new List<ProteinGroup> { new ProteinGroup(psm.ProteinAccession, "", "") }));
            }
            Assert.That(ids.Count, Is.GreaterThan(100));

            var engine = new FlashLfqEngine(ids, silent: true, maxThreads: 1);
            engine.CalculateTheoreticalIsotopeDistributions();
            var index = PeakIndexingEngine.InitializeIndexingEngine(file);
            engine.IndexingEngineDictionary[file] = index;

            var peakfindingTolerance = new PpmTolerance(FlashLfqEngine.PeakfindingPpmTolerance);
            int envelopes = 0;
            foreach (var id in ids)
            {
                foreach (int z in new[] { 1, 2, 3, 4 })
                {
                    foreach (double isotopeShift in new[] { 0.0, -Constants.C13MinusC12, Constants.C13MinusC12 })
                    {
                        var xic = index.GetXic((id.PeakfindingMass + isotopeShift).ToMz(z), id.Ms2RetentionTimeInMinutes,
                            peakfindingTolerance, FlashLfqEngine.MissedScansAllowed);

                        var expected = PreviousGetIsotopicEnvelopes(engine, xic, id, z, file);
                        var actual = engine.GetIsotopicEnvelopes(xic, id, z, file);

                        Assert.That(actual.Count, Is.EqualTo(expected.Count));
                        for (int i = 0; i < expected.Count; i++)
                        {
                            Assert.That(actual[i].IndexedPeak, Is.SameAs(expected[i].IndexedPeak));
                            Assert.That(actual[i].ChargeState, Is.EqualTo(expected[i].ChargeState));
                            Assert.That(BitConverter.DoubleToInt64Bits(actual[i].Intensity), Is.EqualTo(BitConverter.DoubleToInt64Bits(expected[i].Intensity)));
                            Assert.That(BitConverter.DoubleToInt64Bits(actual[i].PearsonCorrelation), Is.EqualTo(BitConverter.DoubleToInt64Bits(expected[i].PearsonCorrelation)));
                        }
                        envelopes += expected.Count;
                    }
                }
            }
            Assert.That(envelopes, Is.GreaterThan(1000), "the comparison must cover a real number of envelopes");
        }

        #region Verbatim copy of GetIsotopicEnvelopes and CheckIsotopicEnvelopeCorrelation before they stopped allocating

        private static List<IsotopicEnvelope> PreviousGetIsotopicEnvelopes(
            FlashLfqEngine engine,
            List<IIndexedPeak> xic,
            Identification identification,
            int chargeState,
            SpectraFileInfo spectraFile)
        {
            var isotopicEnvelopes = new List<IsotopicEnvelope>();
            var isotopeMassShifts = engine.ModifiedSequenceToIsotopicDistribution[identification.ModifiedSequence];

            if (isotopeMassShifts.Count < engine.FlashParams.NumIsotopesRequired)
            {
                return isotopicEnvelopes;
            }

            PpmTolerance isotopeTolerance = new PpmTolerance(engine.FlashParams.IsotopePpmTolerance);

            double[] experimentalIsotopeIntensities = new double[isotopeMassShifts.Count];
            double[] theoreticalIsotopeMassShifts = isotopeMassShifts.Select(p => p.Item1).ToArray();
            double[] theoreticalIsotopeAbundances = isotopeMassShifts.Select(p => p.Item2).ToArray();
            int peakfindingMassIndex = (int)Math.Round(identification.PeakfindingMass - identification.MonoisotopicMass, 0);
            List<IIndexedPeak> isotopologuePeaks = new List<IIndexedPeak>();

            var massShiftToIsotopePeaks = new Dictionary<int, List<(double expIntensity, double theorIntensity, double theorMass)>>
            {
                { -1, new List<(double, double, double)>() },
                { 0, new List<(double, double, double)>() },
                { 1, new List<(double, double, double)>() },
            };

            List<int> directions = new List<int> { -1, 1 };

            foreach (IIndexedPeak peak in xic)
            {
                Array.Clear(experimentalIsotopeIntensities, 0, experimentalIsotopeIntensities.Length);
                foreach (var kvp in massShiftToIsotopePeaks)
                {
                    kvp.Value.Clear();
                }

                double observedMass = peak.M.ToMass(chargeState);
                double observedMassError = observedMass - identification.PeakfindingMass;

                foreach (var shift in massShiftToIsotopePeaks)
                {
                    foreach (int direction in directions)
                    {
                        int start = direction == -1
                            ? peakfindingMassIndex - 1
                            : peakfindingMassIndex;

                        for (int i = start; i < theoreticalIsotopeAbundances.Length && i >= 0; i += direction)
                        {
                            double isotopeMass = identification.MonoisotopicMass + observedMassError +
                                                 theoreticalIsotopeMassShifts[i] + shift.Key * Constants.C13MinusC12;
                            double theoreticalIsotopeIntensity = theoreticalIsotopeAbundances[i] * peak.Intensity;

                            IIndexedPeak isotopePeak = engine.IndexingEngineDictionary[spectraFile]
                                .GetIndexedPeak(isotopeMass.ToMz(chargeState), peak.ZeroBasedScanIndex, isotopeTolerance);

                            if (isotopePeak == null
                                || isotopePeak.Intensity < theoreticalIsotopeIntensity / 4.0
                                || isotopePeak.Intensity > theoreticalIsotopeIntensity * 4.0)
                            {
                                break;
                            }

                            shift.Value.Add((isotopePeak.Intensity, theoreticalIsotopeIntensity, isotopeMass));
                            if (shift.Key == 0)
                            {
                                experimentalIsotopeIntensities[i] = isotopePeak.Intensity;
                                isotopologuePeaks.Add(isotopePeak);
                            }
                        }
                    }
                }

                if (massShiftToIsotopePeaks[0].Count < engine.FlashParams.NumIsotopesRequired)
                {
                    continue;
                }

                if (PreviousCheckIsotopicEnvelopeCorrelation(engine, massShiftToIsotopePeaks, peak, chargeState, isotopeTolerance, spectraFile, out var pearsonCorr))
                {
                    for (int i = 0; i < experimentalIsotopeIntensities.Length; i++)
                    {
                        if (experimentalIsotopeIntensities[i] == 0)
                        {
                            experimentalIsotopeIntensities[i] = theoreticalIsotopeAbundances[i] * experimentalIsotopeIntensities[peakfindingMassIndex];
                        }
                    }

                    isotopicEnvelopes.Add(new IsotopicEnvelope(peak, chargeState, experimentalIsotopeIntensities.Sum(), pearsonCorr));
                }
            }

            return isotopicEnvelopes;
        }

        private static bool PreviousCheckIsotopicEnvelopeCorrelation(
            FlashLfqEngine engine,
            Dictionary<int, List<(double expIntensity, double theorIntensity, double theorMass)>> massShiftToIsotopePeaks,
            IIndexedPeak peak,
            int chargeState,
            PpmTolerance isotopeTolerance,
            SpectraFileInfo spectraFile,
            out double pearsonCorrelation)
        {
            pearsonCorrelation = Correlation.Pearson(
                massShiftToIsotopePeaks[0].Select(p => p.expIntensity),
                massShiftToIsotopePeaks[0].Select(p => p.theorIntensity));

            foreach (var shift in massShiftToIsotopePeaks)
            {
                if (!shift.Value.Any())
                {
                    continue;
                }

                double unexpectedMass = shift.Value.Min(p => p.theorMass) - Constants.C13MinusC12;
                IIndexedPeak unexpectedPeak = engine.IndexingEngineDictionary[spectraFile]
                    .GetIndexedPeak(unexpectedMass.ToMz(chargeState), peak.ZeroBasedScanIndex, isotopeTolerance);

                if (unexpectedPeak == null)
                {
                    shift.Value.Add((0, 0, unexpectedMass));
                }
                else
                {
                    shift.Value.Add((unexpectedPeak.Intensity, 0, unexpectedMass));
                }
            }

            double corrWithPadding = Correlation.Pearson(
                massShiftToIsotopePeaks[0].Select(p => p.expIntensity),
                massShiftToIsotopePeaks[0].Select(p => p.theorIntensity));
            double corrShiftedLeft = Correlation.Pearson(
                massShiftToIsotopePeaks[-1].Select(p => p.expIntensity),
                massShiftToIsotopePeaks[-1].Select(p => p.theorIntensity));
            double corrShiftedRight = Correlation.Pearson(
                massShiftToIsotopePeaks[1].Select(p => p.expIntensity),
                massShiftToIsotopePeaks[1].Select(p => p.theorIntensity));

            if (double.IsNaN(corrShiftedLeft))
            {
                corrShiftedLeft = -1;
            }
            if (double.IsNaN(corrShiftedRight))
            {
                corrShiftedRight = -1;
            }

            return pearsonCorrelation > 0.7 && corrShiftedLeft - corrWithPadding < 0.1 && corrShiftedRight - corrWithPadding < 0.1;
        }

        #endregion
    }
}
