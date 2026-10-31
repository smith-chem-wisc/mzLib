using FlashLFQ;
using MassSpectrometry;
using NUnit.Framework;
using System.Collections.Generic;
using System.Linq;
using IsotopicEnvelope = FlashLFQ.IsotopicEnvelope;

namespace Test.FlashLFQ
{
    /// <summary>
    /// Issue #1484. Match-between-runs scores each transfer against a random-retention-time decoy, taken from
    /// the donor file's MS/MS peaks that lie 5 to 11 Da from the donor. The decoy used to be picked by an
    /// index computed modulo the number of candidates, so one more (or one fewer) donor peptide in that mass
    /// window re-drew the decoy of every donor whose window held it. On PXD009792, 25, 9 and 5 extra donor
    /// peptides per file flipped a fifth to a third of the accepted transfers.
    ///
    /// The decoy is now the candidate with the smallest stable hash of (donor file, donor sequence, acceptor
    /// file, candidate sequence): rendezvous hashing. A new candidate changes a donor's decoy only when it
    /// wins, and removing a candidate changes it only when it was the one drawn.
    /// </summary>
    [TestFixture]
    public static class TestMbrDecoyDraw
    {
        private const double DonorRetentionTime = 10;
        private const double MinimumRetentionTimeDifference = 5;

        // Every peptide here carries a modification, so a key read from the base sequence instead of the
        // modified one would change the pinned draws
        private const string Modified = "[Common Variable:Oxidation on M]";

        private static readonly SpectraFileInfo DonorFile = new SpectraFileInfo("donor.mzML", "A", 0, 0, 0);
        private static readonly SpectraFileInfo AcceptorFile = new SpectraFileInfo("acceptor.mzML", "A", 1, 0, 0);
        private static readonly SpectraFileInfo OtherAcceptorFile = new SpectraFileInfo("other.mzML", "A", 2, 0, 0);

        private static readonly FlashLfqEngine Engine = new FlashLfqEngine(
            new List<Identification> { Donor(0) }, matchBetweenRuns: true, silent: true);

        /// <summary>
        /// 200 donor identifications, all with candidates 6 to 8 Da heavier, so every candidate below sits in
        /// every donor's window.
        /// </summary>
        private static Identification Donor(int index)
        {
            // MS/MS retention times as a real search reports them (not whole minutes)
            var id = new Identification(DonorFile, "DONOR" + Letters(index), "DONOR" + Letters(index) + Modified,
                1000 + index * 0.0137, DonorRetentionTime + 0.123 + index * 0.0371, 2, new List<ProteinGroup>());
            id.PeakfindingMass = id.MonoisotopicMass;
            return id;
        }

        private static List<Identification> Donors() => Enumerable.Range(0, 200).Select(Donor).ToList();

        /// <summary>
        /// An MS/MS peak in the donor file, far enough from the donor in retention time to be a decoy.
        /// </summary>
        private static ChromatographicPeak Candidate(string sequence, double mass, double retentionTime)
        {
            var id = new Identification(DonorFile, sequence, sequence + Modified, mass, retentionTime, 2, new List<ProteinGroup>());
            id.PeakfindingMass = mass;
            var peak = new ChromatographicPeak(id, DonorFile);
            peak.IsotopicEnvelopes = new List<IsotopicEnvelope>
            {
                new IsotopicEnvelope(new IndexedMassSpectralPeak(mass / 2 + 1.007, 1e6, 1, retentionTime), 2, 1e6, 1)
            };
            peak.CalculateIntensityForThisFeature(false);
            return peak;
        }

        private static List<ChromatographicPeak> Candidates(int count) => Enumerable.Range(0, count)
            .Select(i => Candidate("CAND" + Letters(i), 1007 + i * 0.05, 50 + i))
            .ToList();

        private static string Letters(int index)
        {
            string letters = "";
            do
            {
                letters = (char)('A' + index % 26) + letters;
                index /= 26;
            } while (index > 0);
            return letters;
        }

        private static List<string> Draws(List<ChromatographicPeak> candidates, SpectraFileInfo acceptor = null) =>
            Donors().Select(donor => Engine.GetRandomPeak(candidates, DonorRetentionTime,
                    MinimumRetentionTimeDifference, donor, acceptor ?? AcceptorFile)
                ?.Identifications.First().BaseSequence)
            .ToList();

        [Test]
        public static void AnAddedCandidateChangesADonorsDecoyOnlyWhenTheNewCandidateIsDrawn()
        {
            List<ChromatographicPeak> candidates = Candidates(20);
            List<string> before = Draws(candidates);

            // An unrelated donor peptide that lands in every donor's mass window, inserted in mass order
            ChromatographicPeak added = Candidate("UNRELATED", 1007.525, 80);
            var withAdded = candidates.Take(10).Append(added).Concat(candidates.Skip(10)).ToList();
            List<string> after = Draws(withAdded);

            var changed = Enumerable.Range(0, before.Count).Where(i => before[i] != after[i]).ToList();
            Assert.That(changed.Select(i => after[i]), Is.All.EqualTo("UNRELATED"),
                "a donor whose decoy changed must now draw the added candidate, not another old one");
            Assert.That(changed.Count, Is.LessThan(before.Count / 4),
                "about 1 donor in 21 should move to the new candidate");
        }

        [Test]
        public static void RemovingACandidateChangesOnlyTheDecoysOfDonorsThatDrewIt()
        {
            List<ChromatographicPeak> candidates = Candidates(20);
            List<string> before = Draws(candidates);

            string removed = candidates[7].Identifications.First().BaseSequence;
            List<string> after = Draws(candidates.Where((_, i) => i != 7).ToList());

            for (int i = 0; i < before.Count; i++)
            {
                if (before[i] != removed)
                {
                    Assert.That(after[i], Is.EqualTo(before[i]), $"donor {i} did not draw {removed}, so its decoy must stay");
                }
                else
                {
                    Assert.That(after[i], Is.Not.EqualTo(removed));
                }
            }
        }

        [Test]
        public static void TheDecoyDoesNotDependOnTheOrderOfTheCandidates()
        {
            List<ChromatographicPeak> candidates = Candidates(20);
            List<string> inMassOrder = Draws(candidates);

            var reversed = Enumerable.Reverse(candidates).ToList();
            Assert.That(Draws(reversed), Is.EqualTo(inMassOrder));
        }

        [Test]
        public static void TheAcceptorFileIsPartOfTheDraw()
        {
            List<ChromatographicPeak> candidates = Candidates(20);
            List<string> intoAcceptor = Draws(candidates, AcceptorFile);
            List<string> intoOther = Draws(candidates, OtherAcceptorFile);

            Assert.That(Enumerable.Range(0, intoAcceptor.Count).Count(i => intoAcceptor[i] != intoOther[i]),
                Is.GreaterThan(intoAcceptor.Count / 2),
                "two acceptors should draw their decoys independently");
        }

        /// <summary>
        /// The draw must be a fixed function of its keys, the same in every process. A per-process hash such as
        /// string.GetHashCode would make this value change from run to run.
        /// </summary>
        [Test]
        public static void TheDecoyIsTheSameInEveryProcess()
        {
            List<string> draws = Draws(Candidates(20));
            Assert.That(draws.Take(8), Is.EqualTo(new[]
            {
                "CANDM", "CANDQ", "CANDB", "CANDO", "CANDC", "CANDS", "CANDI", "CANDQ"
            }), "computed independently from FNV-1a 64 + fmix64 over the UTF-8 keys (a Python reimplementation)");
        }

        [Test]
        public static void TheHashIsFnv1aOverUtf8()
        {
            // Published FNV-1a 64-bit test vectors
            Assert.That(StableHash.Fnv1a64(""), Is.EqualTo(0xcbf29ce484222325UL));
            Assert.That(StableHash.Fnv1a64("a"), Is.EqualTo(0xaf63dc4c8601ec8cUL));
            Assert.That(StableHash.Fnv1a64("foobar"), Is.EqualTo(0x85944171f73967e8UL));

            // Separated fields hash apart
            ulong ab_c = StableHash.Append(StableHash.AppendSeparator(StableHash.Append(StableHash.FnvOffsetBasis, "AB")), "C");
            ulong a_bc = StableHash.Append(StableHash.AppendSeparator(StableHash.Append(StableHash.FnvOffsetBasis, "A")), "BC");
            Assert.That(ab_c, Is.Not.EqualTo(a_bc));

            // A long, non-ASCII text goes through the pooled buffer and still matches its stack-buffer prefix
            string longText = new string('é', 400);
            Assert.That(StableHash.Fnv1a64(longText),
                Is.EqualTo(StableHash.Append(StableHash.Fnv1a64(new string('é', 100)), new string('é', 300))));
        }

        [Test]
        public static void TheWindowStillWidensWhenNoCandidateIsNearInMass()
        {
            // Every candidate is about 50 Da away: outside 5-11 Da, inside the first widened window
            var farCandidates = Enumerable.Range(0, 5)
                .Select(i => Candidate("FAR" + Letters(i), 1050 + i, 50 + i))
                .ToList();
            Identification donor = Donor(0);

            ChromatographicPeak drawn = Engine.GetRandomPeak(farCandidates, DonorRetentionTime,
                MinimumRetentionTimeDifference, donor, AcceptorFile);
            Assert.That(drawn, Is.Not.Null);
            Assert.That(farCandidates, Does.Contain(drawn));

            Assert.That(Engine.GetRandomPeak(new List<ChromatographicPeak>(), DonorRetentionTime,
                MinimumRetentionTimeDifference, donor, AcceptorFile), Is.Null);
        }

        [Test]
        public static void ACandidateWithTheDonorsBaseSequenceOrRetentionTimeIsNeverDrawn()
        {
            Identification donor = Donor(3);
            var candidates = new List<ChromatographicPeak>
            {
                Candidate(donor.BaseSequence, 1007, 50),        // same base sequence
                Candidate("TOOCLOSE", 1007.1, DonorRetentionTime + 1), // within the retention-time margin
                Candidate("ALLOWED", 1007.2, 60),
            };

            ChromatographicPeak drawn = Engine.GetRandomPeak(candidates, DonorRetentionTime,
                MinimumRetentionTimeDifference, donor, AcceptorFile);
            Assert.That(drawn.Identifications.First().BaseSequence, Is.EqualTo("ALLOWED"));
        }
    }
}
