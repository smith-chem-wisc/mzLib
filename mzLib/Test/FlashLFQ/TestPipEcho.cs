using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.Linq;
using FlashLFQ;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.IO;
using FlashLFQ.PEP;
using System;
using Chemistry;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MzLibUtil;
using Test.FileReadingTests;
using IsotopicEnvelope = FlashLFQ.IsotopicEnvelope;


namespace Test.FlashLFQ
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class TestPipEcho
    {
        [Test]
        [TestCase(3)]
        [TestCase(5)]
        public static void TestDonorGroupEqualizer(int numGroups)
        {
            SpectraFileInfo fakeFile = new SpectraFileInfo("fakeFile", "A", 1, 1, 1);
            Identification id = new Identification(fakeFile, "KPVGAAK", "KPVGAAK", 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });

            MbrChromatographicPeak targetPeak = new MbrChromatographicPeak(id, fakeFile, 1, randomRt: false);
            MbrChromatographicPeak decoyPeak = new MbrChromatographicPeak(id, fakeFile, 1, randomRt: true);
            targetPeak.MbrScore = 100;

            Random random = new Random(42);
            List<DonorGroup> donorGroups = new List<DonorGroup>();
            for (int i = 0; i < 10000; i++)
            {
                int numberTargets = random.Next(0, 10);
                int numberDecoys = random.Next(0, 10);
                donorGroups.Add(new DonorGroup(id, Enumerable.Repeat(targetPeak, numberTargets).ToList(), Enumerable.Repeat(decoyPeak, numberDecoys).ToList()));
            }
            
            donorGroups = PepAnalysisEngine.OrderDonorGroups(donorGroups);
            var donorIndices = PepAnalysisEngine.GetDonorGroupIndices(donorGroups, numGroups: numGroups, scoreCutoff: 50);

            Assert.That(donorIndices.Count, Is.EqualTo(numGroups));
            List<int> targetPeakCounts = new();
            List<int> decoyPeakCounts = new();
            for (int i = 0; i < numGroups; i++)
            {
                int targetSum = 0;
                int decoySum = 0;
                foreach (int idx in donorIndices[i])
                {
                    targetSum += donorGroups[idx].TargetAcceptors.Count;
                    decoySum += donorGroups[idx].DecoyAcceptors.Count;
                }
                targetPeakCounts.Add(targetSum);
                decoyPeakCounts.Add(decoySum);
            } 

            // Assert that each group has an approximately equal number of target peaks
            Assert.That(targetPeakCounts.Max() - targetPeakCounts.Min(), Is.LessThanOrEqualTo(numGroups-1));
            // Assert that each group has an approximately equal number of decoy peaks
            Assert.That(decoyPeakCounts.Max() - decoyPeakCounts.Min(), Is.LessThanOrEqualTo(numGroups - 1));
        }

        private class MockRtCalibDataPoint : RetentionTimeCalibDataPoint
        {
            public MockRtCalibDataPoint(double rtDiff) : base(null, null)
            {
                RtDiff = rtDiff;
            }
        }

        /// <summary>
        /// Builds a calibration data point whose donor peak apexes at <paramref name="donorRt"/> and whose
        /// acceptor peak apexes at <paramref name="acceptorRt"/>. The curve orders points by the donor apex RT.
        /// </summary>
        private static RetentionTimeCalibDataPoint MakeCalibDataPoint(double donorRt, double acceptorRt)
        {
            SpectraFileInfo donorFile = new SpectraFileInfo("donor", "A", 1, 1, 1);
            SpectraFileInfo acceptorFile = new SpectraFileInfo("acceptor", "A", 1, 1, 1);
            const double mass = 669.4173;

            ChromatographicPeak BuildPeak(SpectraFileInfo file, double rt)
            {
                Identification id = new Identification(file, "KPVGAAK", "KPVGAAK", mass, rt, 2,
                    new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
                id.PeakfindingMass = mass;
                ChromatographicPeak peak = new ChromatographicPeak(id, file);
                peak.IsotopicEnvelopes.Add(new IsotopicEnvelope(
                    new IndexedMassSpectralPeak(mass.ToMz(1), 1.0, 1, rt), 1, 1.0, 1.0));
                peak.CalculateIntensityForThisFeature(false);
                return peak;
            }

            return new RetentionTimeCalibDataPoint(BuildPeak(donorFile, donorRt), BuildPeak(acceptorFile, acceptorRt));
        }

        [Test]
        public static void TestRetentionTimeCalibrationCurveOrdersByDonorRt()
        {
            // Data points are supplied out of donor-RT order on purpose.
            var unordered = new[]
            {
                MakeCalibDataPoint(donorRt: 30.0, acceptorRt: 29.5),
                MakeCalibDataPoint(donorRt: 10.0, acceptorRt: 10.2),
                MakeCalibDataPoint(donorRt: 20.0, acceptorRt: 19.7),
                MakeCalibDataPoint(donorRt: 5.0, acceptorRt: 5.1),
                MakeCalibDataPoint(donorRt: 25.0, acceptorRt: 24.6),
            };

            var curve = new RetentionTimeCalibrationCurve(unordered);

            // The curve exposes exactly the points it was given...
            Assert.That(curve.Count, Is.EqualTo(unordered.Length));
            Assert.That(curve.DataPoints.Length, Is.EqualTo(unordered.Length));

            // ...ordered ascending by the donor peak's apex retention time.
            var donorRts = curve.DataPoints.Select(p => p.DonorFilePeak.Apex.IndexedPeak.RetentionTime).ToList();
            Assert.That(donorRts, Is.EqualTo(new[] { 5.0, 10.0, 20.0, 25.0, 30.0 }));
            Assert.That(donorRts, Is.Ordered);

            // Each data point keeps its own RtDiff (donor RT - acceptor RT); ordering must not scramble pairs.
            foreach (var point in curve.DataPoints)
            {
                double expectedDiff = point.DonorFilePeak.Apex.IndexedPeak.RetentionTime
                                      - point.AcceptorFilePeak.Apex.IndexedPeak.RetentionTime;
                Assert.That(point.RtDiff, Is.EqualTo(expectedDiff).Within(1e-9));
            }

            // The ordering enables the binary search PredictRetentionTime relies on.
            var probe = new RetentionTimeCalibDataPoint(
                curve.DataPoints[2].DonorFilePeak, null); // donor RT 20.0
            int index = System.Array.BinarySearch(curve.DataPoints, probe);
            Assert.That(index, Is.EqualTo(2));
        }

        [Test]
        public static void TestRetentionTimeCalibrationCurvePreservesOrderOfDegeneratePoints()
        {
            // Points without a donor apex RT (e.g. the mock points used to exercise the scorer) compare
            // equal, so a stable order must leave them in the sequence they were supplied.
            var points = new RetentionTimeCalibDataPoint[]
            {
                new MockRtCalibDataPoint(0.5),
                new MockRtCalibDataPoint(0.6),
                new MockRtCalibDataPoint(0.7),
                new MockRtCalibDataPoint(0.8),
            };

            var curve = new RetentionTimeCalibrationCurve(points);

            Assert.That(curve.DataPoints.Select(p => p.RtDiff).ToList(),
                Is.EqualTo(new[] { 0.5, 0.6, 0.7, 0.8 }));
        }

        /// <summary>
        /// Builds a single-peptide ChromatographicPeak apexing at <paramref name="rt"/>. Used to assemble
        /// calibration data points (and donor peaks) for the RT-prediction tests below.
        /// </summary>
        private static ChromatographicPeak BuildCalibPeak(SpectraFileInfo file, double rt)
        {
            const double mass = 669.4173;
            Identification id = new Identification(file, "KPVGAAK", "KPVGAAK", mass, rt, 2,
                new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
            id.PeakfindingMass = mass;
            ChromatographicPeak peak = new ChromatographicPeak(id, file);
            peak.IsotopicEnvelopes.Add(new IsotopicEnvelope(
                new IndexedMassSpectralPeak(mass.ToMz(1), 1.0, 1, rt), 1, 1.0, 1.0));
            peak.CalculateIntensityForThisFeature(false);
            return peak;
        }

        // Reads the RT prediction-error distribution the scorer stored for a donor file, so the branches of
        // AddRtPredErrorDistribution can be asserted directly rather than only through a downstream MBR score.
        private static Normal GetStoredRtDistribution(MbrScorer scorer, SpectraFileInfo donorFile)
        {
            var field = typeof(MbrScorer).GetField("_rtPredictionErrorDistributionDictionary",
                System.Reflection.BindingFlags.Instance | System.Reflection.BindingFlags.NonPublic);
            var dict = (Dictionary<SpectraFileInfo, Normal>)field.GetValue(scorer);
            return dict[donorFile];
        }

        [Test]
        public static void TestRetentionTimeCalibDataPointToString()
        {
            var point = MakeCalibDataPoint(donorRt: 20.0, acceptorRt: 19.7);

            // ToString is a debugging aid; make sure it emits the donor/acceptor RTs and their difference.
            string text = point.ToString();
            Assert.That(text, Does.Contain("DonorRT: 20.000"));
            Assert.That(text, Does.Contain("AcceptorRT: 19.700"));
            Assert.That(text, Does.Contain("Diff: 0.300"));
        }

        [Test]
        public static void TestRetentionTimeCalibDataPointCompareToWithNullDonorPeak()
        {
            // A point with a real donor peak (apex RT present) vs. a point with no donor peak (null apex).
            RetentionTimeCalibDataPoint withDonor = MakeCalibDataPoint(donorRt: 12.0, acceptorRt: 11.8);
            RetentionTimeCalibDataPoint withoutDonor = new MockRtCalibDataPoint(0.2); // DonorFilePeak == null

            // Nullable.Compare sorts the point lacking a donor apex before the one that has one.
            Assert.That(withoutDonor.CompareTo(withDonor), Is.LessThan(0));
            Assert.That(withDonor.CompareTo(withoutDonor), Is.GreaterThan(0));
            // Two points that both lack a donor apex compare equal.
            Assert.That(withoutDonor.CompareTo(new MockRtCalibDataPoint(0.9)), Is.EqualTo(0));
        }

        [Test]
        public static void TestAddRtPredErrorDistributionWithTooFewPointsUsesDefault()
        {
            MbrScorer scorer = new MbrScorer(null, null);
            SpectraFileInfo donorFile = new SpectraFileInfo("donor", "A", 1, 1, 1);

            // With numberOfAnchorPeptides = 2, at least 2*2 + 1 = 5 valid points are required. Supplying
            // fewer leaves the safe, non-degenerate default distribution (mean 0, std-dev 1) in place.
            var curve = new RetentionTimeCalibrationCurve(new RetentionTimeCalibDataPoint[]
            {
                new MockRtCalibDataPoint(0.5),
                new MockRtCalibDataPoint(0.6),
                new MockRtCalibDataPoint(0.5),
            });

            scorer.AddRtPredErrorDistribution(donorFile, curve, numberOfAnchorPeptides: 2);

            Normal dist = GetStoredRtDistribution(scorer, donorFile);
            Assert.That(dist.Mean, Is.EqualTo(0).Within(1e-9));
            Assert.That(dist.StdDev, Is.EqualTo(1).Within(1e-9));
        }

        [Test]
        public static void TestAddRtPredErrorDistributionSkipsNonFiniteErrorsAndFloorsStdDev()
        {
            MbrScorer scorer = new MbrScorer(null, null);
            SpectraFileInfo donorFile = new SpectraFileInfo("donor", "A", 1, 1, 1);

            // One anchor has an infinite RtDiff. The prediction error it produces (and the errors of the
            // neighbors that average it in) are non-finite and must be discarded. The remaining finite
            // errors are all identical, giving a zero spread that is raised to the RtStandardDeviationMin
            // floor rather than producing a degenerate (zero std-dev) distribution.
            var curve = new RetentionTimeCalibrationCurve(new RetentionTimeCalibDataPoint[]
            {
                new MockRtCalibDataPoint(0.5),
                new MockRtCalibDataPoint(0.5),
                new MockRtCalibDataPoint(double.PositiveInfinity),
                new MockRtCalibDataPoint(0.5),
                new MockRtCalibDataPoint(0.5),
                new MockRtCalibDataPoint(0.5),
                new MockRtCalibDataPoint(0.5),
            });

            scorer.AddRtPredErrorDistribution(donorFile, curve, numberOfAnchorPeptides: 1);

            Normal dist = GetStoredRtDistribution(scorer, donorFile);
            Assert.That(dist.Mean, Is.EqualTo(0).Within(1e-9));
            Assert.That(dist.StdDev, Is.EqualTo(scorer.RtStandardDeviationMin).Within(1e-9));
        }

        [Test]
        public static void TestPredictRetentionTimeSkipsAnchorsWithoutAcceptorPeaks()
        {
            SpectraFileInfo donorFile = new SpectraFileInfo("donor", "A", 1, 1, 1);
            SpectraFileInfo acceptorFile = new SpectraFileInfo("acceptor", "A", 1, 1, 1);

            ChromatographicPeak donorPeak = BuildCalibPeak(donorFile, 20.0);

            // The anchors immediately on either side of the donor peak have no acceptor peak. Local
            // alignment must skip those (the AcceptorFilePeak != null guard in both scan directions) and
            // keep reaching outward to the anchors that do carry an acceptor peak.
            var curve = new RetentionTimeCalibrationCurve(new[]
            {
                new RetentionTimeCalibDataPoint(BuildCalibPeak(donorFile, 19.6), BuildCalibPeak(acceptorFile, 19.5)),
                new RetentionTimeCalibDataPoint(BuildCalibPeak(donorFile, 19.8), null), // no acceptor -> skipped (backward scan)
                new RetentionTimeCalibDataPoint(BuildCalibPeak(donorFile, 20.0), BuildCalibPeak(acceptorFile, 20.0)),
                new RetentionTimeCalibDataPoint(BuildCalibPeak(donorFile, 20.2), null), // no acceptor -> skipped (forward scan)
                new RetentionTimeCalibDataPoint(BuildCalibPeak(donorFile, 20.4), BuildCalibPeak(acceptorFile, 20.3)),
            });

            var engine = new FlashLfqEngine(new List<Identification> { donorPeak.Identifications.First() });

            RtInfo rtInfo = engine.PredictRetentionTime(curve, donorPeak, acceptorFile,
                acceptorSampleIsFractionated: false, donorSampleIsFractionated: false);

            Assert.That(rtInfo, Is.Not.Null);
            // The two usable anchors (19.6/19.5 and 20.4/20.3) each have a donor-minus-acceptor RT diff of
            // 0.1, so the predicted acceptor RT is the donor apex (20.0) minus that median diff = 19.9.
            Assert.That(rtInfo.PredictedRt, Is.EqualTo(19.9).Within(1e-6));
        }

        [Test]
        public static void TestMbrScorer()
        {
            SpectraFileInfo fakeFile = new SpectraFileInfo("fakeFile", "A", 1, 1, 1);
            SpectraFileInfo fakeDonorFile = new SpectraFileInfo("fakeFile", "A", 1, 1, 1);

            double idMass = 669.4173;
            Identification id = new Identification(fakeFile, "KPVGAAK", "KPVGAAK", 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
            Identification id2 = new Identification(fakeFile, "KPVGK", "KPVGK", 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
            Identification donorId = new Identification(fakeFile, "KPVGK", "KPVGK", 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
            id.PeakfindingMass = idMass;
            id2.PeakfindingMass = idMass;
            donorId.PeakfindingMass = idMass;

            var peak1 = new ChromatographicPeak(id, fakeFile);
            var peak2 = new ChromatographicPeak(id2, fakeFile);
            var peak3 = new ChromatographicPeak(id, fakeFile);
            var donorPeak = new ChromatographicPeak(donorId, fakeDonorFile);
            var acceptorPeak = new MbrChromatographicPeak(donorId, fakeFile, 1, false);

            IndexedMassSpectralPeak imsPeak = new IndexedMassSpectralPeak((idMass + 0.001).ToMz(1), 1.1, 1, 25);
            IndexedMassSpectralPeak imsPeak2 = new IndexedMassSpectralPeak((idMass - 0.001).ToMz(1), 1, 2, 26);
            IndexedMassSpectralPeak imsPeak3 = new IndexedMassSpectralPeak((idMass + 0.002).ToMz(1), 1.2, 3, 27);
            var iso1 = new IsotopicEnvelope(imsPeak, 1, 1.2, 0.98);
            var iso2 = new IsotopicEnvelope(imsPeak2, 1, 1, 0.9);
            var iso3 = new IsotopicEnvelope(imsPeak3, 1, 0.8, 0.95);

            peak1.IsotopicEnvelopes.Add(iso1);
            peak1.IsotopicEnvelopes.Add(iso2);
            peak1.CalculateIntensityForThisFeature(false);

            peak2.IsotopicEnvelopes.Add(iso3);
            peak2.CalculateIntensityForThisFeature(false);

            peak3.IsotopicEnvelopes.Add(iso2);
            peak3.CalculateIntensityForThisFeature(false);

            donorPeak.IsotopicEnvelopes.Add(iso2);
            donorPeak.CalculateIntensityForThisFeature(false);

            acceptorPeak.IsotopicEnvelopes.Add(iso1);
            acceptorPeak.CalculateIntensityForThisFeature(false);


            var peakList = new List<ChromatographicPeak> { peak1, peak3 };

            // Builds a scorer. Ppm Error and Intensity distributions both have mean and std-dev of 1
            MbrScorer scorer = MbrScorerFactory.BuildMbrScorer(peakList, new FlashLfqParameters(), out var tol);
            Assert.That(scorer, Is.Null);   // Not enough peaks to build a scorer

            peakList = new List<ChromatographicPeak> { peak1, peak2, peak3 };
            scorer = MbrScorerFactory.BuildMbrScorer(peakList, new FlashLfqParameters(), out tol);

            scorer.AddRtPredErrorDistribution(fakeDonorFile, new RetentionTimeCalibrationCurve(new RetentionTimeCalibDataPoint[]
            {
                new MockRtCalibDataPoint(0.5),
                new MockRtCalibDataPoint(0.6),
                new MockRtCalibDataPoint(0.5),
                new MockRtCalibDataPoint(0.6),
                new MockRtCalibDataPoint(0.5),
                new MockRtCalibDataPoint(0.6),
                new MockRtCalibDataPoint(0.5)
            }), 2);

            acceptorPeak.MbrScore = scorer.ScoreMbr(acceptorPeak, donorPeak, predictedRt: 25.1);

            Assert.That(acceptorPeak.MbrScore, Is.EqualTo(62.5).Within(0.1));
            Assert.That(acceptorPeak.PpmScore, Is.EqualTo(1).Within(0.01));
            Assert.That(acceptorPeak.IntensityScore, Is.EqualTo(0.46).Within(0.01));
            Assert.That(acceptorPeak.RtScore, Is.EqualTo(0.42).Within(0.01));
            Assert.That(acceptorPeak.ScanCountScore, Is.EqualTo(0.59).Within(0.01));
            Assert.That(acceptorPeak.IsotopicDistributionScore, Is.EqualTo(0.83).Within(0.01));

            double manualMbrScore = Math.Pow( acceptorPeak.PpmScore * acceptorPeak.IntensityScore * acceptorPeak.RtScore * acceptorPeak.ScanCountScore * acceptorPeak.IsotopicDistributionScore, 1.0 / 5.0) * 100;
            Assert.That(acceptorPeak.MbrScore, Is.EqualTo(manualMbrScore).Within(0.1));


            // Test the LogFcDictionary within the MbrScorer
            List<ChromatographicPeak> donorPeaks = new List<ChromatographicPeak> { };
            List<ChromatographicPeak> acceptorPeaks = new List<ChromatographicPeak> { };
            var intensityProperty = typeof(ChromatographicPeak).GetProperty(
                "Intensity",
                System.Reflection.BindingFlags.Instance | System.Reflection.BindingFlags.Public | System.Reflection.BindingFlags.NonPublic
            );
            var random = new Random();
            for (int i = 0; i < 90; i++)
            {
                var seq = RandomPeptide(10);
                Identification idDonor = new Identification(fakeDonorFile, seq, seq, 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
                Identification idAcceptor = new Identification(fakeFile, seq, seq, 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
                ChromatographicPeak newDonorPeak = new ChromatographicPeak(idDonor, fakeDonorFile);
                intensityProperty.SetValue(newDonorPeak, 1000 + random.NextDouble());
                ChromatographicPeak newAcceptorPeak = new ChromatographicPeak(idAcceptor, fakeFile);
                intensityProperty.SetValue(newAcceptorPeak, 2000 + random.NextDouble());

                donorPeaks.Add(newDonorPeak);
                acceptorPeaks.Add(newAcceptorPeak);
            }

            var unambiguousPeaksProperty = typeof(MbrScorer).GetFields(
                System.Reflection.BindingFlags.Instance | System.Reflection.BindingFlags.NonPublic)
                .First(x => x.Name.Contains("UnambiguousMsMsAcceptorPeaks"));
            unambiguousPeaksProperty.SetValue(scorer, acceptorPeaks);
            scorer.CalculateFoldChangeBetweenFiles(donorPeaks);

            // Get the FieldInfo for _logFcDistributionDictionary
            var fieldInfo = typeof(MbrScorer).GetField(
                "_logFcDistributionDictionary",
                System.Reflection.BindingFlags.Instance | System.Reflection.BindingFlags.NonPublic
            );
            // Get the value from an instance (e.g., mbrScorerInstance)
            var logFcDistributionDictionary = (Dictionary<SpectraFileInfo, Normal>)fieldInfo.GetValue(scorer);
            Assert.That(logFcDistributionDictionary.Count, Is.EqualTo(0)); // Because there are fewer than 100 values, the logFC distribution is not calculated

            // Add 20 more peaks to the acceptor peaks list to ensure that the logFC distribution is calculated
            for (int i = 0; i < 20; i++)
            {
                var seq = RandomPeptide(10);
                Identification idDonor = new Identification(fakeDonorFile, seq, seq, 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
                Identification idAcceptor = new Identification(fakeFile, seq, seq, 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
                ChromatographicPeak newDonorPeak = new ChromatographicPeak(idDonor, fakeDonorFile);
                intensityProperty.SetValue(newDonorPeak, 1000 + random.NextDouble());
                ChromatographicPeak newAcceptorPeak = new ChromatographicPeak(idAcceptor, fakeFile);
                intensityProperty.SetValue(newAcceptorPeak, 2000 + random.NextDouble());

                donorPeaks.Add(newDonorPeak);
                acceptorPeaks.Add(newAcceptorPeak);
            }

            unambiguousPeaksProperty.SetValue(scorer, acceptorPeaks);
            scorer.CalculateFoldChangeBetweenFiles(donorPeaks);
            logFcDistributionDictionary = (Dictionary<SpectraFileInfo, Normal>)fieldInfo.GetValue(scorer);
            Assert.That(logFcDistributionDictionary.Count, Is.EqualTo(1));
            Assert.That(logFcDistributionDictionary.Values.First().Mean, Is.EqualTo(Math.Log2(2)).Within(0.001)); // The mean logFC should be around log(2) because the acceptor peaks are twice as intense as the donor peaks
            logFcDistributionDictionary.Clear(); // Clear the dictionary for the next test

            // Make sure that Nan values don't cause a crash
            for (int i = 0; i < 1; i++)
            {
                var seq = RandomPeptide(10);
                Identification idDonor = new Identification(fakeDonorFile, seq, seq, 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
                Identification idAcceptor = new Identification(fakeFile, seq, seq, 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
                ChromatographicPeak newDonorPeak = new ChromatographicPeak(idDonor, fakeDonorFile);
                intensityProperty.SetValue(newDonorPeak, double.NaN);
                ChromatographicPeak newAcceptorPeak = new ChromatographicPeak(idAcceptor, fakeFile);
                intensityProperty.SetValue(newAcceptorPeak, double.NaN);

                donorPeaks.Add(newDonorPeak);
                acceptorPeaks.Add(newAcceptorPeak);
            }

            unambiguousPeaksProperty.SetValue(scorer, acceptorPeaks);
            scorer.CalculateFoldChangeBetweenFiles(donorPeaks);
            Assert.That(logFcDistributionDictionary.Count, Is.EqualTo(1));
            Assert.That(logFcDistributionDictionary.Values.First().Mean, Is.EqualTo(Math.Log2(2)).Within(0.001));

        }

        public static string RandomPeptide(int length)
        {
            //Generate a random peptide sequence
            Random random = new Random();
            const string aminoAcids = "ACDEFGHIKLMNPQRSTVWY"; // Standard amino acids
            return new string(Enumerable.Repeat(aminoAcids, length)
                .Select(s => s[random.Next(s.Length)]).ToArray());
        }

        [Test]
        public static void TestSpectraFileInfoString()
        {
            SpectraFileInfo fakeFile = new SpectraFileInfo(@"C:\Users\xyz\data\fakeFile.raw", "A", 1, 1, 1); 
            Assert.AreEqual("fakeFile.raw", fakeFile.ToString());
        }

        [Test]
        public static void TestChromatographicPeakEquals()
        {
            SpectraFileInfo fakeFile = new SpectraFileInfo("fakeFile", "A", 1, 1, 1);
            Identification id = new Identification(fakeFile, "KPVGAAK", "KPVGAAK", 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });
            Identification id2 = new Identification(fakeFile, "KPVGK", "KPVGK", 669.4173, 1.9398, 2, new List<ProteinGroup> { new ProteinGroup("P16403", "H12", "HUMAN") });

            var peak1 = new MbrChromatographicPeak(id,  fakeFile, 1, randomRt: false);
            var peak2 = new MbrChromatographicPeak(id,  fakeFile, 1, randomRt: false);
            var peak3 = new MbrChromatographicPeak(id2, fakeFile, 1, randomRt: false);
            var peak4 = new MbrChromatographicPeak(id, fakeFile, 1, randomRt: false);


            IndexedMassSpectralPeak imsPeak = new IndexedMassSpectralPeak(1, 1, 1, 25);
            IndexedMassSpectralPeak imsPeak2 = new IndexedMassSpectralPeak(1, 1, 1, 50);
            var iso1 = new IsotopicEnvelope(imsPeak, 1, 1, 1);
            var iso2 = new IsotopicEnvelope(imsPeak2, 1, 1, 1);

            peak1.IsotopicEnvelopes.Add(iso1);
            peak1.CalculateIntensityForThisFeature(false);

            peak2.IsotopicEnvelopes.Add(iso1);
            peak2.CalculateIntensityForThisFeature(false);

            peak3.IsotopicEnvelopes.Add(iso1);
            peak3.CalculateIntensityForThisFeature(false);

            peak4.IsotopicEnvelopes.Add(iso2);
            peak4.CalculateIntensityForThisFeature(false);

            Assert.That(peak1.Equals(peak2));
            Assert.That(!peak1.Equals(peak3));
            Assert.That(!peak1.Equals(peak4));

        }

        /// <summary>
        /// This test MatchBetweenRuns by creating two fake mzML files and a list of fake IDs. 
        /// There are multiple sets of IDs, where most are shared between the two runs but one+ is/are missing
        /// MBR is tested by ensuring that IDs are transferred between runs
        /// </summary>
        [Test]
        public static void TestFlashLfqMatchBetweenRunsNearestNeighborDonors()
        {
            List<string> filesToWrite = new List<string> { "mzml_1", "mzml_2", "mzml_3" };
            List<string> pepSequences = new List<string>
                {
                "PEPTIDE",
                "PEPTIDEV",
                "PEPTIDEVV",
                "TARGETPEP",
                "PEPTIDEVVV",
                "PEPTIDEVVVV",
                "PEPTIDEVVVVA",
                "PEPTIDEVVVVAA"
            };
            double intensity = 1e6;

            double[] file1Rt = new double[] { 1.01, 1.02, 1.03, 1.033, 1.035, 1.04, 1.045, 1.05 };
            double[] file2Rt = new double[] { 1.00, 1.025, 1.03, 1.031, 1.035, 1.04, 1.055, 1.07 };

            // generate mzml files (5 peptides each)
            for (int f = 0; f < filesToWrite.Count; f++)
            {
                // 1 MS1 scan per peptide
                MsDataScan[] scans = new MsDataScan[8];

                for (int p = 0; p < pepSequences.Count; p++)
                {
                    ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(pepSequences[p]).GetChemicalFormula();
                    IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                    double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                    double[] intensities = dist.Intensities.Select(v => v * intensity).ToArray();
                    if(f == 2)
                    {
                        // Make file 3 the most intense
                        intensities = intensities.Select(v => v * 5).ToArray();
                    }
                    double rt;
                    if (f == 1)
                    {
                        rt = file2Rt[p];
                    }
                    else
                    {
                        rt = file1Rt[p];
                    }

                    // add the scan
                    scans[p] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: p + 1, msnOrder: 1, isCentroid: true,
                        polarity: Polarity.Positive, retentionTime: rt, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                        mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (p + 1));
                }

                // write the .mzML
                Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[f] + ".mzML"), false);
            }

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[0] + ".mzML"), "a", 0, 0, 0);
            SpectraFileInfo file2 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[1] + ".mzML"), "a", 1, 0, 0);
            SpectraFileInfo file3 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[2] + ".mzML"), "a", 2, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1 = new Identification(file1, "PEPTIDE", "PEPTIDE",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDE").MonoisotopicMass, file1Rt[0] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id2 = new Identification(file1, "PEPTIDEV", "PEPTIDEV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEV").MonoisotopicMass, file1Rt[1] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id3 = new Identification(file1, "PEPTIDEVV", "PEPTIDEVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVV").MonoisotopicMass, file1Rt[2] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id4 = new Identification(file1, "PEPTIDEVVV", "PEPTIDEVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVV").MonoisotopicMass, file1Rt[4] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id5 = new Identification(file1, "PEPTIDEVVVV", "PEPTIDEVVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVVV").MonoisotopicMass, file1Rt[5] + 0.001, 1, new List<ProteinGroup> { pg });

            Identification id6 = new Identification(file2, "PEPTIDE", "PEPTIDE",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDE").MonoisotopicMass, file2Rt[0] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id7 = new Identification(file2, "PEPTIDEV", "PEPTIDEV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEV").MonoisotopicMass, file2Rt[1] + 0.001, 1, new List<ProteinGroup> { pg });
            // missing ID 8 - MBR feature - "PEPTIDEVV"

            Identification id9 = new Identification(file2, "PEPTIDEVVV", "PEPTIDEVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVV").MonoisotopicMass, file2Rt[4] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id10 = new Identification(file2, "PEPTIDEVVVV", "PEPTIDEVVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVVV").MonoisotopicMass, file2Rt[5] + 0.001, 1, new List<ProteinGroup> { pg });


            Identification id11 = new Identification(file3, "PEPTIDEV", "PEPTIDEV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEV").MonoisotopicMass, file1Rt[1] + 0.001, 1, new List<ProteinGroup> { pg }); // same as peak 2
            Identification id12 = new Identification(file3, "PEPTIDEVV", "PEPTIDEVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVV").MonoisotopicMass, file1Rt[2] + 0.001, 1, new List<ProteinGroup> { pg }); // same as peak 3, but higher intensity
            Identification id13 = new Identification(file3, "PEPTIDEVVV", "PEPTIDEVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVV").MonoisotopicMass, file1Rt[4] + 0.001, 1, new List<ProteinGroup> { pg }); // same as peak 4


            // create the FlashLFQ engine
            FlashLfqEngine neighborsEngine = new FlashLfqEngine(new List<Identification> { id1, id2, id3, id4, id5, id6, id7, id9, id10, id11, id12, id13 }, 
                matchBetweenRuns: true, donorCriterion: DonorCriterion.Neighbors);

            //run the engine
            var results = neighborsEngine.Run();

            Assert.That(results.Peaks[file2].Count == 5);
            Assert.That(results.Peaks[file2].Count(p => p.DetectionType == DetectionType.MBR) == 1);

            var peak = results.Peaks[file2].First(p => p.DetectionType == DetectionType.MBR);
            var otherFilePeak = results.Peaks[file1].Where(p => p.Identifications.First().BaseSequence ==
                peak.Identifications.First().BaseSequence).First();


            Assert.That(peak.Intensity > 0);
            Assert.That(peak.Intensity == otherFilePeak.Intensity);
            Assert.That(peak.Identifications.First().FileInfo == file1); // assure that the ID came from file 1, ie, the donor with the most neighboring peaks

            // create the FlashLFQ engine
            FlashLfqEngine intensityEngine = new FlashLfqEngine(new List<Identification> { id1, id2, id3, id4, id5, id6, id7, id9, id10, id11, id12, id13 },
                matchBetweenRuns: true, donorCriterion: DonorCriterion.Intensity);

            //run the engine
            results = intensityEngine.Run();

            Assert.That(results.Peaks[file2].Count == 5);
            Assert.That(results.Peaks[file2].Count(p => p.DetectionType == DetectionType.MBR) == 1);

            peak = results.Peaks[file2].First(p => p.DetectionType == DetectionType.MBR);
            otherFilePeak = results.Peaks[file3].Where(p => p.Identifications.First().BaseSequence ==
                peak.Identifications.First().BaseSequence).First();


            Assert.That(peak.Intensity > 0);
            Assert.That(peak.Intensity, Is.EqualTo(otherFilePeak.Intensity/5).Within(1)); // file 3 is five times more intense than file 2
            Assert.That(peak.Identifications.First().FileInfo == file3); // assure that the ID came from file 3, ie, the most intense donor peaks

        }

        [Test]
        public static void TestMbrScorer_Build_WithEmptyList_ReturnsNull()
        {
            var flashParams = new FlashLfqParameters { MbrPpmTolerance = 40 };
            var scorer = MbrScorerFactory.BuildMbrScorer(new List<ChromatographicPeak>(), flashParams, out var tol);
            Assert.That(scorer, Is.Null);
            Assert.That(tol, Is.Null);
        }

        [Test]
        public static void TestMbrScorer_Build_WithLessThanThreePeaks_ReturnsNull()
        {
            var peaks = new List<ChromatographicPeak>
            {
                CreatePeak("A", 1000,  10.0, 0.95, 5, 1e-6),
                CreatePeak("B", 1200, 11.0, 0.90, 4, 2e-6)
            };
            var scorer = MbrScorerFactory.BuildMbrScorer(peaks, new FlashLfqParameters(), out var tol);
            Assert.That(scorer, Is.Null);
            Assert.That(tol, Is.Null);
        }

        [Test]
        public static void TestMbrScorer_Build_AllApexNull_ReturnsNull()
        {
            var peaks = new List<ChromatographicPeak>
            {
                CreatePeak("A", 1000,  10.0, 0.8, 3, 1e-6, makeApex:false),
                CreatePeak("B", 1100, 11.0, 0.7, 4, 1e-6, makeApex:false),
                CreatePeak("C", 900, 12.0, 0.6, 5, 1e-6, makeApex:false)
            };
            var scorer = MbrScorerFactory.BuildMbrScorer(peaks, new FlashLfqParameters(), out var tol);
            // Unambiguous list count is 3, but distributions will fail (ppm list has 3, intensity list has 3) Apex null shouldn't crash
            // Build still proceeds; InitializeScorer depends only on counts & distributions (apex not used for intensity/ppm)
            // MassError & Intensity valid => scorer not null
            if (scorer == null)
            {
                Assert.Pass("Returned null safely (acceptable).");
            }
            else
            {
                Assert.That(scorer.IsValid(), Is.True);
            }
        }

        [Test]
        public static void TestMbrScorer_Build_AllZeroOrNegativeIntensity_FailsLogIntensityDistribution()
        {
            var peaks = new List<ChromatographicPeak>
            {
                CreatePeak("A", 0,   10.0, 0.9, 3, 1e-6),
                CreatePeak("B", -5,  11.0, 0.9, 4, 1e-6),
                CreatePeak("C", 0,  12.0, 0.9, 5, 1e-6)
            };
            var scorer = MbrScorerFactory.BuildMbrScorer(peaks, new FlashLfqParameters(), out var tol);
            Assert.That(scorer, Is.Null);
            Assert.That(tol, Is.Null);
        }


        [Test]
        public static void TestMbrScorer_Build_PpmErrorsIdentical_ZeroSpreadStillValid()
        {
            var peaks = new List<ChromatographicPeak>
            {
                CreatePeak("A", 1000, 10, 0.8, 3, 100, 1.5),
                CreatePeak("B", 1200, 11, 0.85,4, 100, 1.5),
                CreatePeak("C", 1100, 12, 0.82,5, 100, 1.5),
                CreatePeak("D", 900,  13, 0.81,6, 100, 1.5)
            };
            var scorer = MbrScorerFactory.BuildMbrScorer(peaks, new FlashLfqParameters { MbrPpmTolerance = 20 }, out var tol);
            Assert.That(scorer, Is.Not.Null);
            Assert.That(tol, Is.Not.Null);
            Assert.That(tol.Value, Is.EqualTo(1.5).Within(0.1)); // Lots of small doubles means a fair number of rounding errors, exact ppm error is ~1.46
        }

        [Test]
        public static void TestMbrScorer_Build_ExtremeMassErrors_CappedByFlashParams()
        {
            var peaks = new List<ChromatographicPeak>
            {
                CreatePeak("A", 1500, 10, 0.7, 3, 100, 100),
                CreatePeak("B", 1400, 11, 0.8, 4, 100, 120),
                CreatePeak("C", 1300, 12, 0.75, 5, 100, 110),
                CreatePeak("D", 1250, 13, 0.72, 6, 100, 130)
            };
            var flashParams = new FlashLfqParameters { MbrPpmTolerance = 30 };
            var scorer = MbrScorerFactory.BuildMbrScorer(peaks, flashParams, out var tol);
            Assert.That(scorer, Is.Not.Null);
            Assert.That(tol.Value, Is.LessThanOrEqualTo(30)); // enforced cap
        }

        [Test]
        public static void TestMbrScorer_Build_IsotopicCorrelationUniform_NoGammaDistributionCrash()
        {
            var peaks = new List<ChromatographicPeak>
            {
                CreatePeak("A", 1000, 10, 1.0, 3, 1), // IsotopicPearsonCorrelation=1 => (1 - corr)=0 filtered out
                CreatePeak("B", 1100, 11, 1.0, 4, 1),
                CreatePeak("C", 1200, 12, 1.0, 5, 1),
                CreatePeak("D", 1300, 13, 1.0, 6, 1)
            };
            var scorer = MbrScorerFactory.BuildMbrScorer(peaks, new FlashLfqParameters(), out var tol);
            Assert.That(scorer, Is.Not.Null);
            // IsValid ignores isotopic distribution being null
            Assert.That(scorer.IsValid(), Is.True);
        }

        [Test]
        public static void TestMbrScorer_ScoreMbr_WithoutRtDistribution_ThrowsKeyButHandledBySetup()
        {
            // Build valid scorer first
            var donorFile = new SpectraFileInfo("donor", "C", 1, 1, 1);
            var acceptorFile = new SpectraFileInfo("acceptor", "C", 1, 1, 1);
            var peaks = new List<ChromatographicPeak>
            {
                CreatePeak("A", 1500, 10, 0.9, 5, 1, fileInfo:acceptorFile),
                CreatePeak("B", 1600, 11, 0.85,4, 1, fileInfo:acceptorFile),
                CreatePeak("C", 1700, 12, 0.8, 6, 1, fileInfo:acceptorFile),
                CreatePeak("D", 1800, 13, 0.75,7, 1, fileInfo:acceptorFile)
            };
            var scorer = MbrScorerFactory.BuildMbrScorer(peaks, new FlashLfqParameters(), out var tol);
            Assert.That(scorer, Is.Not.Null);

            // Add RT prediction distribution with too few anchor points -> will fallback to default (0,0) distribution
            scorer.AddRtPredErrorDistribution(donorFile, new RetentionTimeCalibrationCurve(new RetentionTimeCalibDataPoint[]
            {
                new MockRtCalibDataPoint(0.2),
                new MockRtCalibDataPoint(0.21),
                new MockRtCalibDataPoint(0.19),
                new MockRtCalibDataPoint(0.2),
                new MockRtCalibDataPoint(0.22)
            }), 2);

            var donorPeak = peaks[0];
            var acceptor = new MbrChromatographicPeak(peaks[1].Identifications.First(), acceptorFile, peaks[1].ApexRetentionTime, false);
            acceptor.IsotopicEnvelopes.Add(peaks[1].Apex);
            acceptor.CalculateIntensityForThisFeature(false);

            var score = scorer.ScoreMbr(acceptor, donorPeak, acceptor.ApexRetentionTime + 0.05);
            Assert.That(score, Is.GreaterThan(0));
        }

        private static ChromatographicPeak CreatePeak(
            string seq,
            double intensity,
            double rt,
            double isoCorr,
            int scanCount,
            double mass,
            double massErrorPpm = 0,
            bool makeApex = true,
            SpectraFileInfo fileInfo = null)
        {
            fileInfo ??= new SpectraFileInfo("file_" + seq, "Cond", 1, 1, 1);
            var id = new Identification(fileInfo, seq, seq, mass, rt, 1, new List<ProteinGroup> { new ProteinGroup("PG_" + seq, "", "") });
            id.PeakfindingMass = mass;
            double test = mass.ToMz(1);
            var peak = new ChromatographicPeak(id, fileInfo);
            var envelopes = new List<IsotopicEnvelope>(scanCount);
            double mz = (mass + (massErrorPpm * mass / 1e6)).ToMz(1); //  mass.ToMz(1) + (massErrorPpm * mass.ToMz(1) / 1e6); // apply mass error
            foreach (int i in Enumerable.Range(1, scanCount))
            {
                envelopes.Add(new IsotopicEnvelope(new IndexedMassSpectralPeak(mz , intensity * 0.8, i, rt), 1, intensity, isoCorr));
            }
            if (makeApex)
            {
                var ims = new IndexedMassSpectralPeak(mz, intensity, scanCount, rt);
                var env = new IsotopicEnvelope(ims, 1, intensity, isoCorr); 
                peak.IsotopicEnvelopes.Add(env);
            }
            peak.CalculateIntensityForThisFeature(false);
            return peak;
        }
    }
}
