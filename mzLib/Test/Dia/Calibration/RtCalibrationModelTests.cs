// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using MassSpectrometry.Dia;
using NUnit.Framework;

namespace Test.Dia.Calibration
{
    [TestFixture]
    public class RtCalibrationModelTests
    {
        // ── Construction & basic properties ──────────────────────────────────

        [Test]
        public void Constructor_ValidParameters_SetsProperties()
        {
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 2.0, sigmaMinutes: 0.3, rSquared: 0.98, anchorCount: 50);

            Assert.That(model.Slope, Is.EqualTo(0.5));
            Assert.That(model.Intercept, Is.EqualTo(2.0));
            Assert.That(model.SigmaMinutes, Is.EqualTo(0.3));
            Assert.That(model.RSquared, Is.EqualTo(0.98));
            Assert.That(model.AnchorCount, Is.EqualTo(50));
        }

        [Test]
        public void Constructor_ZeroSlope_Throws()
        {
            Assert.Throws<ArgumentException>(() =>
                new RtCalibrationModel(slope: 0.0, intercept: 1.0, sigmaMinutes: 0.1,
                    rSquared: 0.99, anchorCount: 10));
        }

        [Test]
        public void SigmaIrt_ComputedCorrectly()
        {
            // σ_iRT = σ_minutes / |slope|
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 0.0, sigmaMinutes: 1.0, rSquared: 0.95, anchorCount: 20);

            Assert.That(model.SigmaIrt, Is.EqualTo(2.0).Within(1e-10));
        }

        [Test]
        public void SigmaIrt_NegativeSlope_UsesAbsoluteValue()
        {
            var model = new RtCalibrationModel(
                slope: -0.25, intercept: 100.0, sigmaMinutes: 0.5, rSquared: 0.95, anchorCount: 20);

            Assert.That(model.SigmaIrt, Is.EqualTo(2.0).Within(1e-10));
        }

        // ── IsReliable ──────────────────────────────────────────────────────

        [Test]
        public void IsReliable_HighRSquaredEnoughAnchors_True()
        {
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 2.0, sigmaMinutes: 0.3, rSquared: 0.95, anchorCount: 50);

            Assert.That(model.IsReliable, Is.True);
        }

        [Test]
        public void IsReliable_LowRSquared_False()
        {
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 2.0, sigmaMinutes: 0.3, rSquared: 0.80, anchorCount: 50);

            Assert.That(model.IsReliable, Is.False);
        }

        [Test]
        public void IsReliable_TooFewAnchors_False()
        {
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 2.0, sigmaMinutes: 0.3, rSquared: 0.99, anchorCount: 5);

            Assert.That(model.IsReliable, Is.False);
        }

        // ── ToIrt / ToMinutes round-trip ─────────────────────────────────────

        [Test]
        public void ToIrt_ToMinutes_RoundTrip()
        {
            // observed_RT = 0.5 * iRT + 2.0
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 2.0, sigmaMinutes: 0.3, rSquared: 0.98, anchorCount: 50);

            double iRT = 40.0;
            double expectedMinutes = 0.5 * 40.0 + 2.0; // = 22.0

            Assert.That(model.ToMinutes(iRT), Is.EqualTo(expectedMinutes).Within(1e-10));
            Assert.That(model.ToIrt(expectedMinutes), Is.EqualTo(iRT).Within(1e-10));
        }

        [Test]
        public void ToIrt_KnownValues()
        {
            // observed_RT = 0.5 * iRT + 2.0
            // iRT = (observed_RT - 2.0) / 0.5
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 2.0, sigmaMinutes: 0.1, rSquared: 0.99, anchorCount: 100);

            Assert.That(model.ToIrt(2.0), Is.EqualTo(0.0).Within(1e-10));  // intercept maps to iRT=0
            Assert.That(model.ToIrt(12.0), Is.EqualTo(20.0).Within(1e-10));
            Assert.That(model.ToIrt(52.0), Is.EqualTo(100.0).Within(1e-10));
        }

        // ── Window half-width ───────────────────────────────────────────────

        [Test]
        public void GetIrtWindowHalfWidth_DefaultK3()
        {
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 0.0, sigmaMinutes: 1.0, rSquared: 0.95, anchorCount: 20);

            // σ_iRT = 1.0 / 0.5 = 2.0, window = 3 * 2.0 = 6.0
            Assert.That(model.GetIrtWindowHalfWidth(), Is.EqualTo(6.0).Within(1e-10));
        }

        [Test]
        public void GetMinutesWindowHalfWidth_DefaultK3()
        {
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 0.0, sigmaMinutes: 1.0, rSquared: 0.95, anchorCount: 20);

            // window = 3 * 1.0 = 3.0 minutes
            Assert.That(model.GetMinutesWindowHalfWidth(), Is.EqualTo(3.0).Within(1e-10));
        }

        [Test]
        public void GetIrtWindowHalfWidth_CustomK()
        {
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 0.0, sigmaMinutes: 1.0, rSquared: 0.95, anchorCount: 20);

            Assert.That(model.GetIrtWindowHalfWidth(k: 2.0), Is.EqualTo(4.0).Within(1e-10));
        }

        // ── RT scoring ──────────────────────────────────────────────────────

        [Test]
        public void ComputeRtScore_ZeroResidual_ReturnsZero()
        {
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 0.0, sigmaMinutes: 1.0, rSquared: 0.95, anchorCount: 20);

            Assert.That(model.ComputeRtScore(0.0), Is.EqualTo(0.0).Within(1e-10));
        }

        [Test]
        public void ComputeRtScore_OneSigma_ReturnsMinusHalf()
        {
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 0.0, sigmaMinutes: 1.0, rSquared: 0.95, anchorCount: 20);

            double sigmaIrt = model.SigmaIrt; // = 2.0
            // rtScore = -(σ^2) / (2 * σ^2) = -0.5
            Assert.That(model.ComputeRtScore(sigmaIrt), Is.EqualTo(-0.5).Within(1e-10));
        }

        [Test]
        public void ComputeRtScore_LargerResidual_MoreNegative()
        {
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 0.0, sigmaMinutes: 1.0, rSquared: 0.95, anchorCount: 20);

            double score1Sigma = model.ComputeRtScore(model.SigmaIrt);
            double score2Sigma = model.ComputeRtScore(2.0 * model.SigmaIrt);

            Assert.That(score2Sigma, Is.LessThan(score1Sigma));
        }

        [Test]
        public void ComputeRtScore_TwoArgOverload_MatchesSingleArg()
        {
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 2.0, sigmaMinutes: 0.5, rSquared: 0.95, anchorCount: 20);

            double libraryIrt = 50.0;
            double scanIrt = 53.0;
            double residual = scanIrt - libraryIrt; // = 3.0

            Assert.That(model.ComputeRtScore(libraryIrt, scanIrt),
                Is.EqualTo(model.ComputeRtScore(residual)).Within(1e-10));
        }

        // ── Provisional model ───────────────────────────────────────────────

        [Test]
        public void CreateProvisional_MapsRangeCorrectly()
        {
            // Run: 5–65 min, Library iRT: 0–100
            var model = RtCalibrationModel.CreateProvisional(
                runRtMin: 5.0, runRtMax: 65.0,
                libraryIrtMin: 0.0, libraryIrtMax: 100.0,
                initialWindowIrt: 20.0);

            // slope = 60 / 100 = 0.6 min/iRT
            Assert.That(model.Slope, Is.EqualTo(0.6).Within(1e-10));

            // intercept: 5 = 0.6 * 0 + b → b = 5
            Assert.That(model.Intercept, Is.EqualTo(5.0).Within(1e-10));

            // Verify endpoints map correctly
            Assert.That(model.ToMinutes(0.0), Is.EqualTo(5.0).Within(1e-10));
            Assert.That(model.ToMinutes(100.0), Is.EqualTo(65.0).Within(1e-10));
        }

        [Test]
        public void CreateProvisional_IsNotReliable()
        {
            var model = RtCalibrationModel.CreateProvisional(5.0, 65.0, 0.0, 100.0);
            Assert.That(model.IsReliable, Is.False);
            Assert.That(model.AnchorCount, Is.EqualTo(0));
        }

        [Test]
        public void CreateProvisional_WindowHalfWidth_MatchesInitialWindow()
        {
            var model = RtCalibrationModel.CreateProvisional(
                runRtMin: 5.0, runRtMax: 65.0,
                libraryIrtMin: 0.0, libraryIrtMax: 100.0,
                initialWindowIrt: 21.0);

            // k=3 window should equal initialWindowIrt = 21
            Assert.That(model.GetIrtWindowHalfWidth(k: 3.0), Is.EqualTo(21.0).Within(1e-6));
        }

        [Test]
        public void CreateProvisional_DegenerateRange_DoesNotThrow()
        {
            // Zero RT range — degenerate but should not crash
            var model = RtCalibrationModel.CreateProvisional(
                runRtMin: 30.0, runRtMax: 30.0,
                libraryIrtMin: 0.0, libraryIrtMax: 100.0);

            Assert.That(model, Is.Not.Null);
            Assert.That(model.IsReliable, Is.False);
        }

        // ── ToString ─────────────────────────────────────────────────────────

        [Test]
        public void ToString_ContainsKeyInfo()
        {
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 2.0, sigmaMinutes: 0.3, rSquared: 0.98, anchorCount: 50);

            string s = model.ToString();
            Assert.That(s, Does.Contain("0.5"));
            Assert.That(s, Does.Contain("R²"));
            Assert.That(s, Does.Contain("anchors=50"));
        }

        [Test]
        public void ToString_UnreliableModel_ShowsTag()
        {
            var model = new RtCalibrationModel(
                slope: 0.5, intercept: 2.0, sigmaMinutes: 0.3, rSquared: 0.50, anchorCount: 3);

            Assert.That(model.ToString(), Does.Contain("UNRELIABLE"));
        }
    }
}
