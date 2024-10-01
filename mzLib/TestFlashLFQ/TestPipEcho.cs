using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.Linq;
using FlashLFQ;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using System.IO;
using FlashLFQ.PEP;
using System;


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

    }
}
