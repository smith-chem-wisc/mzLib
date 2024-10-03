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


namespace TestFlashLFQ
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

            ChromatographicPeak targetPeak = new ChromatographicPeak(id, isMbrPeak: true, fakeFile, randomRt: false);
            ChromatographicPeak decoyPeak = new ChromatographicPeak(id, isMbrPeak: true, fakeFile, randomRt: true);
            targetPeak.MbrScore = 100;

            Random random = new Random(42);
            List<DonorGroup> donorGroups = new List<DonorGroup>();
            for (int i = 0; i < 10000; i++)
            {
                int numberTargets = random.Next(0, 10);
                int numberDecoys = random.Next(0, 10);
                donorGroups.Add(new DonorGroup(id, Enumerable.Repeat(targetPeak, numberTargets).ToList(), Enumerable.Repeat(targetPeak, numberDecoys).ToList()));
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

            var peak1 = new ChromatographicPeak(id, isMbrPeak: false, fakeFile, randomRt: false);
            var peak2 = new ChromatographicPeak(id, isMbrPeak: false, fakeFile, randomRt: false);
            var peak3 = new ChromatographicPeak(id2, isMbrPeak: false, fakeFile, randomRt: false);
            var peak4 = new ChromatographicPeak(id, isMbrPeak: false, fakeFile, randomRt: false);
            var donorPeak = new ChromatographicPeak(donorId, isMbrPeak: false, fakeDonorFile, randomRt: false);
            var acceptorPeak = new ChromatographicPeak(donorId, isMbrPeak: true, fakeFile, randomRt: false);

            IndexedMassSpectralPeak imsPeak = new IndexedMassSpectralPeak((idMass + 0.001).ToMz(1), 1.1, 1, 25);
            IndexedMassSpectralPeak imsPeak2 = new IndexedMassSpectralPeak((idMass - 0.001).ToMz(1), 1, 2, 26);
            var iso1 = new IsotopicEnvelope(imsPeak, 1, 1, 0.98);
            var iso2 = new IsotopicEnvelope(imsPeak2, 1, 1, 0.9);

            peak1.IsotopicEnvelopes.Add(iso1);
            peak1.IsotopicEnvelopes.Add(iso2);
            peak1.CalculateIntensityForThisFeature(false);

            //peak2.IsotopicEnvelopes.Add(iso1);
            //peak2.CalculateIntensityForThisFeature(false);

            //peak3.IsotopicEnvelopes.Add(iso1);
            //peak3.CalculateIntensityForThisFeature(false);

            peak4.IsotopicEnvelopes.Add(iso2);
            peak4.CalculateIntensityForThisFeature(false);

            donorPeak.IsotopicEnvelopes.Add(iso2);
            donorPeak.CalculateIntensityForThisFeature(false);

            acceptorPeak.IsotopicEnvelopes.Add(iso1);
            acceptorPeak.CalculateIntensityForThisFeature(false);


            var peakList = new List<ChromatographicPeak> { peak1, peak4 };
            var peakDict = peakList.ToDictionary(keySelector: p => p.Apex.IndexedPeak, elementSelector: p => p);

            // Builds a scorer. Ppm Error and Intensity distributions both have mean and std-dev of 1
            MbrScorer scorer = new MbrScorer(peakDict, peakList, new MathNet.Numerics.Distributions.Normal(1, 1), new MathNet.Numerics.Distributions.Normal(1,1));
            
            scorer.AddRtPredErrorDistribution(fakeDonorFile, new List<double> { 0.5, 0.6, 0.5, 0.6, 0.5, 0.6, 0.5 }, 2);

            acceptorPeak.MbrScore = scorer.ScoreMbr(acceptorPeak, donorPeak, predictedRt: 25.1);

            Assert.That(acceptorPeak.MbrScore, Is.EqualTo(58.7).Within(0.1));
            Assert.That(acceptorPeak.PpmScore, Is.EqualTo(0.62).Within(0.01));
            Assert.That(acceptorPeak.IntensityScore, Is.EqualTo(0.32).Within(0.01));
            Assert.That(acceptorPeak.RtScore, Is.EqualTo(0.96).Within(0.01));
            Assert.That(acceptorPeak.ScanCountScore, Is.EqualTo(0.5).Within(0.01));
            Assert.That(acceptorPeak.IsotopicDistributionScore, Is.EqualTo(0.74).Within(0.01));
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

            var peak1 = new ChromatographicPeak(id, isMbrPeak: true, fakeFile, randomRt: false);
            var peak2 = new ChromatographicPeak(id, isMbrPeak: true, fakeFile, randomRt: false);
            var peak3 = new ChromatographicPeak(id2, isMbrPeak: true, fakeFile, randomRt: false);
            var peak4 = new ChromatographicPeak(id, isMbrPeak: true, fakeFile, randomRt: false);

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
    }
}
