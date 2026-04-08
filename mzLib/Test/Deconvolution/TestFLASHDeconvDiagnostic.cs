using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;

namespace Test
{
    /// <summary>
    /// Hard-assertion diagnostic tests for FLASHDeconv Step 3–5 mass calculation.
    /// Each test checks one specific assumption. Run all of these when
    /// Deconvolute_SyntheticCytoC_FindsMassWithin20Ppm fails — the first failing
    /// assertion names the broken assumption precisely.
    /// Delete or [Ignore] once Steps 3–5 produce correct masses.
    /// </summary>
    [TestFixture]
    public sealed class TestFLASHDeconvDiagnostic
    {
        private static FLASHDeconvolutionParameters DefaultParams() =>
            new FLASHDeconvolutionParameters(
                minCharge: 1, maxCharge: 60,
                deconvolutionTolerancePpm: 10.0,
                minIsotopicPeakCount: 3,
                minCosineScore: 0.4,
                minMassRange: 50.0,
                maxMassRange: 100_000.0);

        // ── Constants from the first diagnostic run (now asserted) ────────────
        private const double InputMass = 12_223.2;
        private const int ExpectedIdx = 219;
        private const double ExpectedApexMass = 12_236.1822;
        private const double ExpectedDiffToMono = 7.0180;
        // formula.MonoisotopicMass = apex - DiffToMono = 12229.1642
        private const double ExpectedFormulaMono = ExpectedApexMass - ExpectedDiffToMono;

        // ══════════════════════════════════════════════════════════════════════
        // A. Averagine field assumptions
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void A1_GetMostIntenseMassIndex_ReturnsExpectedIndex()
        {
            var p = DefaultParams();
            Assert.That(p.AverageResidueModel.GetMostIntenseMassIndex(InputMass),
                Is.EqualTo(ExpectedIdx),
                $"GetMostIntenseMassIndex({InputMass}) should return {ExpectedIdx}");
        }

        [Test]
        public void A2_AllMasses_Index0_IsApex_IntensityDescending()
        {
            // AllMasses[idx] is sorted intensity-descending, so index 0 = most abundant = apex.
            var p = DefaultParams();
            int idx = p.AverageResidueModel.GetMostIntenseMassIndex(InputMass);
            double[] masses = p.AverageResidueModel.GetAllTheoreticalMasses(idx);
            double[] intens = p.AverageResidueModel.GetAllTheoreticalIntensities(idx);

            // Index 0 should have maximum intensity
            Assert.That(intens[0], Is.EqualTo(intens.Max()),
                "AllIntensities[idx][0] should be the maximum intensity (array is intensity-descending)");
            Assert.That(masses[0], Is.EqualTo(ExpectedApexMass).Within(0.001),
                $"AllMasses[idx][0] (apex) should be ≈ {ExpectedApexMass}");
        }

        [Test]
        public void A3_DiffToMonoisotopic_IsApexMinusFormulaMono()
        {
            // DiffToMonoisotopic[i] = masses[0] - chemicalFormula.MonoisotopicMass
            var p = DefaultParams();
            int idx = p.AverageResidueModel.GetMostIntenseMassIndex(InputMass);
            double diff = p.AverageResidueModel.GetDiffToMonoisotopic(idx);
            double apex = p.AverageResidueModel.GetAllTheoreticalMasses(idx)[0];

            Assert.That(diff, Is.EqualTo(ExpectedDiffToMono).Within(0.001),
                $"DiffToMonoisotopic should be ≈ {ExpectedDiffToMono}");
            // formula.Mono = apex - diff should match the expected value
            Assert.That(apex - diff, Is.EqualTo(ExpectedFormulaMono).Within(0.001),
                $"apex ({apex:F4}) - diff ({diff:F4}) = formula.Mono should be ≈ {ExpectedFormulaMono}");
        }

        [Test]
        public void A4_SortedMasses_Index0_IsFormulaMono()
        {
            // After mass-ascending sort, the lightest peak = formula.MonoisotopicMass.
            var p = DefaultParams();
            int idx = p.AverageResidueModel.GetMostIntenseMassIndex(InputMass);
            double[] masses = p.AverageResidueModel.GetAllTheoreticalMasses(idx);
            double[] intens = p.AverageResidueModel.GetAllTheoreticalIntensities(idx);
            double sortedMono = masses.Zip(intens).OrderBy(x => x.First).First().First;

            Assert.That(sortedMono, Is.EqualTo(ExpectedFormulaMono).Within(0.001),
                $"sortedMasses[0] should be formula.Mono ≈ {ExpectedFormulaMono}. " +
                $"Got {sortedMono:F4}.");
        }

        [Test]
        public void A5_FormulaMono_DiffersFromInputMass_ByMoreThan20Ppm()
        {
            // Critical: the formula's monoisotopic mass (12229.16) is NOT the same as
            // the input mass (12223.2). They differ by ~6 Da = ~490 ppm.
            // Therefore the algorithm should recover 12229.16, not 12223.2.
            double errorPpm = Math.Abs(ExpectedFormulaMono - InputMass) / InputMass * 1e6;
            Assert.That(errorPpm, Is.GreaterThan(20.0),
                $"formula.Mono ({ExpectedFormulaMono:F4}) should differ from input ({InputMass}) " +
                $"by > 20 ppm. Actual difference: {errorPpm:F1} ppm. " +
                $"The test must assert against formula.Mono, not the input mass.");
        }

        [Test]
        public void A6_DiffToMono_IsNearIntegerMultiple_OfC13MinusC12()
        {
            // apex - mono = nApex × C13MinusC12 within ~0.01 Da.
            // If this fails, the Averagine apex doesn't fall on a 1-Da isotope boundary.
            var p = DefaultParams();
            int idx = p.AverageResidueModel.GetMostIntenseMassIndex(InputMass);
            double diff = p.AverageResidueModel.GetDiffToMonoisotopic(idx);
            double nApexExact = diff / Constants.C13MinusC12;
            int nApex = (int)Math.Round(nApexExact);

            Console.WriteLine($"DiffToMono={diff:F4}, C13={Constants.C13MinusC12:F6}");
            Console.WriteLine($"nApex = {nApexExact:F4} ≈ {nApex}");

            Assert.That(nApexExact, Is.EqualTo(nApex).Within(0.1),
                $"DiffToMono ({diff:F4}) / C13MinusC12 = {nApexExact:F4} should be near-integer. " +
                $"If it's not, the apex is NOT a simple n*C13 above the monoisotopic.");
        }

        // ══════════════════════════════════════════════════════════════════════
        // B. Step 2 candidate assumptions
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void B1_Step2_FindsCandidatesNearFormulaMono()
        {
            // At least one candidate should be within a few C13 steps of formula.Mono.
            var p = DefaultParams();
            var spectrum = BuildSpectrumFromFormulaMono(ExpectedFormulaMono, new[] { 9, 10, 11, 12, 13 }, p);
            var cands = RunStep2(spectrum, p);

            Console.WriteLine($"Step 2 found {cands.Count} candidates:");
            foreach (var c in cands.Take(8))
            {
                double offsetDa = c.Mass - ExpectedFormulaMono;
                double offsetSteps = offsetDa / Constants.C13MinusC12;
                Console.WriteLine($"  Mass={c.Mass:F4}, offset from formulaMono={offsetDa:F4} Da = {offsetSteps:F2} C13 steps");
            }

            Assert.That(cands.Count, Is.GreaterThan(0), "Step 2 must find at least one candidate");

            // Any candidate should be reasonably close to the isotope envelope
            bool anyClose = cands.Any(c => Math.Abs(c.Mass - ExpectedFormulaMono) < 15.0);
            Assert.That(anyClose, Is.True,
                $"At least one candidate should be within 15 Da of formulaMono ({ExpectedFormulaMono:F4})");
        }

        [Test]
        public void B2_CandidateMass_IsNotAtApex_ItIsAtSomeLowerIsotope()
        {
            // Step 2 bins on the mz peaks that trigger the continuous-charge run.
            // These are typically NOT the apex peak but lighter isotopes.
            // candidate.Mass should be BETWEEN formula.Mono and apex.
            var p = DefaultParams();
            var spectrum = BuildSpectrumFromFormulaMono(ExpectedFormulaMono, new[] { 9, 10, 11, 12, 13 }, p);
            var cands = RunStep2(spectrum, p);
            Assert.That(cands.Count, Is.GreaterThan(0));

            var cand = cands[0];
            Console.WriteLine($"candidate.Mass = {cand.Mass:F4}");
            Console.WriteLine($"formula.Mono   = {ExpectedFormulaMono:F4}");
            Console.WriteLine($"apex           = {ExpectedApexMass:F4}");
            Console.WriteLine($"candidate is {cand.Mass - ExpectedFormulaMono:F4} Da above mono, " +
                              $"{ExpectedApexMass - cand.Mass:F4} Da below apex");

            // candidate.Mass should be strictly less than the apex
            Assert.That(cand.Mass, Is.LessThan(ExpectedApexMass),
                $"candidate.Mass ({cand.Mass:F4}) should be below the apex ({ExpectedApexMass:F4})");
            // candidate.Mass should be at or above formula.Mono
            Assert.That(cand.Mass, Is.GreaterThanOrEqualTo(ExpectedFormulaMono - 1.0),
                $"candidate.Mass ({cand.Mass:F4}) should be at or above formula.Mono ({ExpectedFormulaMono:F4})");
        }

        // ══════════════════════════════════════════════════════════════════════
        // C. New mass calculation: apex-observed approach
        // ══════════════════════════════════════════════════════════════════════

        [Test]
        public void C1_MostIntensePeak_IsNearApex_InNeutralMassSpace()
        {
            // After recruiting peaks for a single representative charge,
            // the most intense peak should be near the apex neutral mass.
            var p = DefaultParams();
            var spectrum = BuildSpectrumFromFormulaMono(ExpectedFormulaMono, new[] { 9, 10, 11, 12, 13 }, p);
            var cands = RunStep2(spectrum, p);
            Assert.That(cands.Count, Is.GreaterThan(0));

            var cand = cands[0];
            int polSign = 1;
            var tolerance = new PpmTolerance(p.DeconvolutionTolerancePpm);
            int nIso = p.AverageResidueModel
                .GetAllTheoreticalMasses(p.AverageResidueModel.GetMostIntenseMassIndex(cand.Mass))
                .Length;

            // Recruit for all charges
            var recruited = new List<(double mz, double intensity, int z)>();
            for (int z = cand.MinAbsCharge; z <= cand.MaxAbsCharge; z++)
            {
                double baseMz = cand.Mass.ToMz(polSign * z);
                double isoDelta = Constants.C13MinusC12 / z;
                for (int n = -10; n < nIso; n++)
                {
                    double target = baseMz + n * isoDelta;
                    if (target <= 0) continue;
                    var hits = spectrum.GetPeakIndicesWithinTolerance(target, tolerance);
                    if (hits.Count == 0 && n >= 0) break;
                    if (hits.Count == 0) continue;
                    int best = hits.OrderByDescending(i => spectrum.YArray[i]).First();
                    recruited.Add((spectrum.XArray[best], spectrum.YArray[best], z));
                }
            }

            Assert.That(recruited.Count, Is.GreaterThan(0), "Must recruit at least one peak");
            var apexPeak = recruited.OrderByDescending(pk => pk.intensity).First();
            double apexObsMass = apexPeak.mz.ToMass(polSign * apexPeak.z);

            Console.WriteLine($"Most intense recruited: mz={apexPeak.mz:F4}, z={apexPeak.z}, I={apexPeak.intensity:G4}");
            Console.WriteLine($"apexObsMass = {apexObsMass:F4}");
            Console.WriteLine($"ExpectedApexMass = {ExpectedApexMass:F4}");
            Console.WriteLine($"error = {Math.Abs(apexObsMass - ExpectedApexMass):F4} Da = {Math.Abs(apexObsMass - ExpectedApexMass) / ExpectedApexMass * 1e6:F1} ppm");

            Assert.That(apexObsMass, Is.EqualTo(ExpectedApexMass).Within(ExpectedApexMass * 50e-6),
                $"Most intense recruited peak mass ({apexObsMass:F4}) should be within 50 ppm of apex ({ExpectedApexMass:F4})");
        }

        [Test]
        public void C2_ApexObsMass_Minus_DiffToMono_GivesFormulaMono_Within5Ppm()
        {
            // The proposed new formula: monoMass = apexObsMass - DiffToMonoisotopic
            // This should recover ExpectedFormulaMono within 5 ppm.
            var p = DefaultParams();
            var spectrum = BuildSpectrumFromFormulaMono(ExpectedFormulaMono, new[] { 9, 10, 11, 12, 13 }, p);
            var cands = RunStep2(spectrum, p);
            Assert.That(cands.Count, Is.GreaterThan(0));

            var cand = cands[0];
            int polSign = 1;
            var tolerance = new PpmTolerance(p.DeconvolutionTolerancePpm);
            int nIso = p.AverageResidueModel
                .GetAllTheoreticalMasses(p.AverageResidueModel.GetMostIntenseMassIndex(cand.Mass))
                .Length;

            var recruited = new List<(double mz, double intensity, int z)>();
            for (int z = cand.MinAbsCharge; z <= cand.MaxAbsCharge; z++)
            {
                double baseMz = cand.Mass.ToMz(polSign * z);
                double isoDelta = Constants.C13MinusC12 / z;
                for (int n = -10; n < nIso; n++)
                {
                    double target = baseMz + n * isoDelta;
                    if (target <= 0) continue;
                    var hits = spectrum.GetPeakIndicesWithinTolerance(target, tolerance);
                    if (hits.Count == 0 && n >= 0) break;
                    if (hits.Count == 0) continue;
                    int best = hits.OrderByDescending(i => spectrum.YArray[i]).First();
                    recruited.Add((spectrum.XArray[best], spectrum.YArray[best], z));
                }
            }

            var apexPeak = recruited.OrderByDescending(pk => pk.intensity).First();
            double apexObsMass = apexPeak.mz.ToMass(polSign * apexPeak.z);
            int avgIdxCand = p.AverageResidueModel.GetMostIntenseMassIndex(cand.Mass);
            double diffCand = p.AverageResidueModel.GetDiffToMonoisotopic(avgIdxCand);
            double monoComputed = apexObsMass - diffCand;

            Console.WriteLine($"apexObsMass       = {apexObsMass:F4}");
            Console.WriteLine($"DiffToMono(cand)  = {diffCand:F4}");
            Console.WriteLine($"monoComputed      = {monoComputed:F4}");
            Console.WriteLine($"ExpectedFormulaMono = {ExpectedFormulaMono:F4}");
            Console.WriteLine($"error = {Math.Abs(monoComputed - ExpectedFormulaMono):F4} Da = {Math.Abs(monoComputed - ExpectedFormulaMono) / ExpectedFormulaMono * 1e6:F1} ppm");

            Assert.That(monoComputed, Is.EqualTo(ExpectedFormulaMono).Within(ExpectedFormulaMono * 5e-6),
                $"apexObsMass ({apexObsMass:F4}) - DiffToMono ({diffCand:F4}) = {monoComputed:F4} " +
                $"should be within 5 ppm of ExpectedFormulaMono ({ExpectedFormulaMono:F4})");
        }

        // ══════════════════════════════════════════════════════════════════════
        // Helpers
        // ══════════════════════════════════════════════════════════════════════

        private static MzSpectrum BuildSpectrumFromFormulaMono(
            double formulaMono, int[] charges, FLASHDeconvolutionParameters p,
            double baseIntensity = 100_000.0)
        {
            // Build Averagine intensity profile from the apex lookup
            int idx = p.AverageResidueModel.GetMostIntenseMassIndex(formulaMono + 7.0); // look up near apex
            double[] m = p.AverageResidueModel.GetAllTheoreticalMasses(idx);
            double[] ins = p.AverageResidueModel.GetAllTheoreticalIntensities(idx);
            var sorted = m.Zip(ins).OrderBy(x => x.First).ToArray();
            double[] sortedIntens = sorted.Select(x => x.Second).ToArray();
            int nIso = sortedIntens.Length;

            var mzList = new List<double>();
            var itList = new List<double>();
            foreach (int z in charges)
                for (int n = 0; n < nIso; n++)
                {
                    double mz = (formulaMono + n * Constants.C13MinusC12).ToMz(z);
                    double it = baseIntensity * sortedIntens[n];
                    if (it < baseIntensity * 0.001) continue;
                    mzList.Add(mz);
                    itList.Add(it);
                }

            var pairs = mzList.Zip(itList).OrderBy(x => x.First).ToArray();
            return new MzSpectrum(pairs.Select(x => x.First).ToArray(),
                                  pairs.Select(x => x.Second).ToArray(), false);
        }

        private static List<FLASHDeconvolutionAlgorithm.CandidateMass> RunStep2(
            MzSpectrum spectrum, FLASHDeconvolutionParameters p)
        {
            var logPeaks = FLASHDeconvolutionAlgorithm.BuildLogMzPeaks(spectrum, spectrum.Range, p.Polarity);
            int minAbs = Math.Abs(p.MinAssumedChargeState);
            int cr = Math.Abs(p.MaxAssumedChargeState) - minAbs + 1;
            var up = FLASHDeconvolutionAlgorithm.BuildUniversalPattern(minAbs, cr);
            var hp = FLASHDeconvolutionAlgorithm.BuildHarmonicPatterns(minAbs, cr);
            double bmf = 1.0 / (p.DeconvolutionTolerancePpm * 1e-6);
            return FLASHDeconvolutionAlgorithm.FindCandidateMasses(logPeaks, up, hp, bmf, p);
        }
    }
}