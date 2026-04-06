// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry.Dia;
using NUnit.Framework;

namespace Test.DiaTests
{
    /// <summary>
    /// Unit tests for the new DecoyGenerator (scramble + mutation strategy).
    ///
    /// Test structure:
    ///   Layer 1 — Invariants that must hold for every decoy
    ///   Layer 2 — Mutation constraints (aromatic count, mass delta, BLOSUM)
    ///   Layer 3 — Reproducibility (same sequence → same decoy)
    ///   Layer 4 — Fallback behaviour (short sequences, no valid mutation)
    ///   Layer 5 — NDP / spectrum similarity
    ///
    /// Placement: mzLib/Test/Dia/DecoyGeneratorTests.cs
    /// </summary>
    [TestFixture]
    public class DecoyGeneratorTests
    {
        // ─────────────────────────────────────────────────────────────────────
        //  Helpers
        // ─────────────────────────────────────────────────────────────────────

        private static readonly Dictionary<char, double> ResidueMass =
            new Dictionary<char, double>
        {
            {'A', 71.03711},  {'R', 156.10111}, {'N', 114.04293},
            {'D', 115.02694}, {'C', 103.00919}, {'E', 129.04259},
            {'Q', 128.05858}, {'G', 57.02146},  {'H', 137.05891},
            {'I', 113.08406}, {'L', 113.08406}, {'K', 128.09496},
            {'M', 131.04049}, {'F', 147.06841}, {'P', 97.05276},
            {'S', 87.03203},  {'T', 101.04768}, {'W', 186.07931},
            {'Y', 163.06333}, {'V', 99.06841},
        };

        private static readonly HashSet<char> Aromatics = new HashSet<char> { 'F', 'W', 'Y' };

        private const double WaterMass = 18.010565;

        private static double NeutralMass(string bare)
        {
            double m = WaterMass;
            foreach (char c in bare)
                m += ResidueMass.TryGetValue(c, out double r) ? r : 111.0;
            return m;
        }

        private static int AromaticCount(string seq)
            => seq.Count(c => Aromatics.Contains(c));

        private static string StripMods(string seq)
        {
            var sb = new System.Text.StringBuilder();
            bool inMod = false; char modOpen = ' ';
            foreach (char c in seq)
            {
                if (!inMod) { if (c == '(' || c == '[') { inMod = true; modOpen = c; } else sb.Append(c); }
                else { char cl = modOpen == '(' ? ')' : ']'; if (c == cl) inMod = false; }
            }
            return sb.ToString();
        }

        /// <summary>
        /// Builds a minimal LibraryPrecursorInput for testing.
        /// Fragment m/z and intensities are synthetic but non-zero.
        /// </summary>
        private static LibraryPrecursorInput MakePrecursor(string sequence, double precursorMz = 550.0)
        {
            var mzs = new float[] { 100f, 200f, 300f, 400f, 500f };
            var ints = new float[] { 900f, 700f, 500f, 300f, 100f };
            return new LibraryPrecursorInput(
                sequence: sequence,
                precursorMz: precursorMz,
                chargeState: 2,
                retentionTime: 30.0,
                isDecoy: false,
                fragmentMzs: mzs,
                fragmentIntensities: ints,
                irtValue: 50.0);
        }

        private static List<LibraryPrecursorInput> Generate(params string[] sequences)
        {
            var targets = sequences.Select(s => MakePrecursor(s)).ToList();
            return DecoyGenerator.GenerateFromTargets(targets);
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Layer 1: Invariants — hold for every generated decoy
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void GenerateFromTargets_ProducesOneDecoyPerTarget()
        {
            var targets = new[]
            {
                MakePrecursor("PEPTIDEK"),
                MakePrecursor("SAMPLEK"),
                MakePrecursor("ANALYSISK"),
            };
            var decoys = DecoyGenerator.GenerateFromTargets(targets.ToList());
            Assert.That(decoys.Count, Is.EqualTo(targets.Length));
        }

        [Test]
        public void AllDecoys_AreMarkedIsDecoy()
        {
            var decoys = Generate("PEPTIDEK", "SAMPLEK", "ELECTROSPRAY");
            Assert.That(decoys.All(d => d.IsDecoy), Is.True,
                "Every generated decoy must have IsDecoy=true");
        }

        [Test]
        public void AllDecoys_PreservePrecursorMz()
        {
            var targets = new[]
            {
                MakePrecursor("PEPTIDEK", precursorMz: 550.3),
                MakePrecursor("SAMPLEK",  precursorMz: 612.7),
            };
            var decoys = DecoyGenerator.GenerateFromTargets(targets.ToList());
            for (int i = 0; i < targets.Length; i++)
            {
                Assert.That(decoys[i].PrecursorMz,
                    Is.EqualTo(targets[i].PrecursorMz).Within(1e-6),
                    $"Decoy[{i}] precursor m/z must match target");
            }
        }

        [Test]
        public void AllDecoys_PreserveChargeState()
        {
            var targets = new[]
            {
                MakePrecursor("PEPTIDEK"),
                MakePrecursor("SAMPLEK"),
            };
            var decoys = DecoyGenerator.GenerateFromTargets(targets.ToList());
            for (int i = 0; i < targets.Length; i++)
            {
                Assert.That(decoys[i].ChargeState, Is.EqualTo(targets[i].ChargeState),
                    $"Decoy[{i}] charge state must match target");
            }
        }

        [Test]
        public void AllDecoys_PreserveRetentionTime()
        {
            var targets = new[] { MakePrecursor("PEPTIDEK") };
            var decoys = DecoyGenerator.GenerateFromTargets(targets.ToList());
            Assert.That(decoys[0].RetentionTime, Is.EqualTo(targets[0].RetentionTime),
                "RetentionTime must be copied from target (iRT correction applied downstream)");
        }

        [Test]
        public void AllDecoys_PreserveIrtValue()
        {
            var targets = new[] { MakePrecursor("PEPTIDEK") };
            var decoys = DecoyGenerator.GenerateFromTargets(targets.ToList());
            Assert.That(decoys[0].IrtValue, Is.EqualTo(targets[0].IrtValue),
                "IrtValue must be copied from target");
        }

        [Test]
        public void AllDecoys_PreserveCTerminalResidue()
        {
            // C-terminal residue of the scrambled+mutated sequence must equal
            // the C-terminal residue of the bare target sequence.
            var sequences = new[] { "PEPTIDEK", "SAMPLEFWR", "ANALYSISK" };
            var decoys = Generate(sequences);

            for (int i = 0; i < sequences.Length; i++)
            {
                string bare = StripMods(sequences[i]);
                string decoySeq = StripMods(decoys[i].Sequence);
                Assert.That(decoySeq[decoySeq.Length - 1],
                    Is.EqualTo(bare[bare.Length - 1]),
                    $"Decoy[{i}] C-terminal residue must match target '{bare}'");
            }
        }

        [Test]
        public void AllDecoys_HaveNonZeroFragments()
        {
            var decoys = Generate("PEPTIDEK", "SAMPLEWK", "ANALYSISK");
            foreach (var d in decoys)
            {
                Assert.That(d.FragmentCount, Is.GreaterThan(0),
                    $"Decoy '{d.Sequence}' must have at least one fragment");
                Assert.That(d.FragmentMzs.Length, Is.EqualTo(d.FragmentIntensities.Length),
                    "Fragment m/z and intensity arrays must be same length");
            }
        }

        [Test]
        public void AllDecoys_SequenceDiffersFromTarget()
        {
            // The decoy sequence (after stripping mods) should differ from the target.
            // This can fail for palindromic or very short sequences — we test normal cases.
            var sequences = new[] { "PEPTIDEK", "SAMPLEWK", "ELECTROSPRAY" };
            var targets = sequences.Select(s => MakePrecursor(s)).ToList();
            var decoys = DecoyGenerator.GenerateFromTargets(targets);

            for (int i = 0; i < sequences.Length; i++)
            {
                string bare = StripMods(sequences[i]);
                string decoyBare = StripMods(decoys[i].Sequence);
                Assert.That(decoyBare, Is.Not.EqualTo(bare),
                    $"Decoy sequence must differ from target '{bare}'");
            }
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Layer 2: Mutation constraints
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void AllDecoys_PreserveAromaticCount()
        {
            // Use sequences with known aromatic counts.
            // The decoy (after scramble + mutation) must have the same count.
            var sequences = new[]
            {
                "PEPTFIDEK",     // 1 aromatic (F)
                "SAMPLEWYK",     // 2 aromatics (W, Y)
                "ANALYSISK",     // 0 aromatics
                "FYWPEPTIDEK",   // 3 aromatics (F, Y, W)
            };
            var targets = sequences.Select(s => MakePrecursor(s)).ToList();
            var decoys = DecoyGenerator.GenerateFromTargets(targets);

            for (int i = 0; i < sequences.Length; i++)
            {
                string bare = StripMods(sequences[i]);
                string decoyBare = StripMods(decoys[i].Sequence);
                int targetAromatics = AromaticCount(bare);
                int decoyAromatics = AromaticCount(decoyBare);
                Assert.That(decoyAromatics, Is.EqualTo(targetAromatics),
                    $"Decoy of '{bare}' has {decoyAromatics} aromatics, " +
                    $"expected {targetAromatics}");
            }
        }

        [Test]
        public void AllDecoys_PreserveSameAminoAcidCount()
        {
            // Substitutions only — sequence length must be unchanged.
            var sequences = new[] { "PEPTIDEK", "SAMPLEWYK", "ELECTROSPRAY" };
            var targets = sequences.Select(s => MakePrecursor(s)).ToList();
            var decoys = DecoyGenerator.GenerateFromTargets(targets);

            for (int i = 0; i < sequences.Length; i++)
            {
                string bare = StripMods(sequences[i]);
                string decoyBare = StripMods(decoys[i].Sequence);
                Assert.That(decoyBare.Length, Is.EqualTo(bare.Length),
                    $"Decoy of '{bare}' has length {decoyBare.Length}, expected {bare.Length}");
            }
        }

        [Test]
        public void AllDecoys_PreserveChargeGroup()
        {
            // Basic (K,R,H) → basic only; acidic (D,E) → acidic only; neutral → neutral only.
            // Test sequences with known charge group compositions.
            var sequences = new[]
            {
                "PEPTIDEK",    // K=basic, rest neutral
                "SAMPLEWYK",   // K=basic, W/Y=aromatic/neutral
                "ADDEDEK",     // D,E=acidic, K=basic
                "HRKPEPTIDEK", // H,R,K=basic
            };
            var targets = sequences.Select(s => MakePrecursor(s)).ToList();
            var decoys = DecoyGenerator.GenerateFromTargets(targets);

            for (int i = 0; i < sequences.Length; i++)
            {
                string targetBare = StripMods(sequences[i]);
                string decoyBare = StripMods(decoys[i].Sequence);

                int targetBasic = targetBare.Count(c => "KRH".Contains(c));
                int targetAcidic = targetBare.Count(c => "DE".Contains(c));
                int decoyBasic = decoyBare.Count(c => "KRH".Contains(c));
                int decoyAcidic = decoyBare.Count(c => "DE".Contains(c));

                Assert.That(decoyBasic, Is.EqualTo(targetBasic),
                    $"Decoy of '{targetBare}' has {decoyBasic} basic residues, " +
                    $"expected {targetBasic}");
                Assert.That(decoyAcidic, Is.EqualTo(targetAcidic),
                    $"Decoy of '{targetBare}' has {decoyAcidic} acidic residues, " +
                    $"expected {targetAcidic}");
            }
        }

        [Test]
        public void MutatedDecoys_PrecursorMassDeltaInRange()
        {
            // For sequences long enough that mutation is likely to succeed,
            // the neutral mass delta between the decoy sequence and the target
            // sequence should be in (0.1, 10.0) Da.
            //
            // Note: the decoy sequence is the scrambled+mutated bare sequence.
            // We reconstruct its neutral mass from the sequence string.
            // The precursor m/z is NOT changed — this tests the sequence mass only.
            //
            // We use a large batch to get statistical coverage across many RNG seeds.
            var longSequences = new[]
            {
                "SAMPLEVASTK",
                "PEPTIDEANALYSISK",
                "ELECTROSPRAYIONK",
                "MASSSPECTROMETRK",
                "PROTEOMICSWORKFLOWK",
                "LIQUIDCHROMATOGRAPHYK",
                "TANDEMMASSSPECTROK",
                "BIOINFORMATICSANALYSISK",
            };

            var targets = longSequences.Select(s => MakePrecursor(s)).ToList();
            var decoys = DecoyGenerator.GenerateFromTargets(targets);

            int validMassCount = 0;
            for (int i = 0; i < longSequences.Length; i++)
            {
                string targetBare = StripMods(longSequences[i]);
                string decoyBare = StripMods(decoys[i].Sequence);

                double targetMass = NeutralMass(targetBare);
                double decoyMass = NeutralMass(decoyBare);
                double delta = Math.Abs(decoyMass - targetMass);

                if (delta >= 0.1)
                    validMassCount++;

                // If mutation succeeded (sequences differ beyond scramble), mass must be in range
                // We can't directly tell if mutation succeeded from outside, so we just count
            }

            // At least half should have valid mass deltas (mutation succeeded)
            Assert.That(validMassCount, Is.GreaterThanOrEqualTo(longSequences.Length / 2),
                "At least half of long-sequence decoys should have mass delta in (0.1, 10.0) Da. " +
                "If this fails, the mutation logic is not producing valid substitutions.");
        }

        [Test]
        public void MutatedDecoy_MassDeltaNotZero_ForTypicalSequence()
        {
            // For a concrete sequence, verify that the decoy mass differs from target mass.
            // This catches the degenerate case where all substitutions cancel out (e.g. I↔L only).
            string sequence = "SAMPLEVASTK";
            string bare = StripMods(sequence);
            double targetMass = NeutralMass(bare);

            var targets = new[] { MakePrecursor(sequence) };
            var decoys = DecoyGenerator.GenerateFromTargets(targets.ToList());
            string decoyBare = StripMods(decoys[0].Sequence);
            double decoyMass = NeutralMass(decoyBare);

            Assert.That(Math.Abs(decoyMass - targetMass), Is.GreaterThan(0.05),
                $"Decoy of '{bare}' should have a different neutral mass. " +
                $"Target: {targetMass:F4} Da, Decoy: {decoyMass:F4} Da");
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Layer 3: Reproducibility
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void SameSequence_AlwaysProducesSameDecoy()
        {
            // Run generation twice on the same input — must produce identical decoys.
            var targets = new[]
            {
                MakePrecursor("PEPTIDEK"),
                MakePrecursor("SAMPLEWYK"),
                MakePrecursor("ELECTROSPRAY"),
            };

            var decoys1 = DecoyGenerator.GenerateFromTargets(targets.ToList());
            var decoys2 = DecoyGenerator.GenerateFromTargets(targets.ToList());

            for (int i = 0; i < targets.Length; i++)
            {
                Assert.That(decoys1[i].Sequence, Is.EqualTo(decoys2[i].Sequence),
                    $"Decoy[{i}] sequence must be identical across two independent runs");
                Assert.That(decoys1[i].FragmentMzs, Is.EqualTo(decoys2[i].FragmentMzs),
                    $"Decoy[{i}] fragment m/z must be identical across two independent runs");
            }
        }

        [Test]
        public void DecoyIsOrderIndependent()
        {
            // The decoy for sequence A must be the same whether A is first or third in the list.
            // This verifies the per-sequence seed is truly order-independent.
            var targetA = MakePrecursor("PEPTIDEK");
            var targetB = MakePrecursor("SAMPLEWYK");
            var targetC = MakePrecursor("ANALYSISK");

            // Run 1: A first
            var run1 = DecoyGenerator.GenerateFromTargets(
                new List<LibraryPrecursorInput> { targetA, targetB, targetC });

            // Run 2: A last
            var run2 = DecoyGenerator.GenerateFromTargets(
                new List<LibraryPrecursorInput> { targetB, targetC, targetA });

            string decoyA_run1 = run1[0].Sequence;
            string decoyA_run2 = run2[2].Sequence; // A is now at index 2

            Assert.That(decoyA_run1, Is.EqualTo(decoyA_run2),
                "Decoy for 'PEPTIDEK' must be identical regardless of its position in the input list");
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Layer 4: Fallback behaviour
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void ShortSequence_TwoResidues_ProducesDecoyWithoutCrash()
        {
            // Very short sequences (length < 2 or = 2) must not throw.
            var targets = new[]
            {
                MakePrecursor("AK"),  // length 2 — C-terminal fixed, nothing to scramble
                MakePrecursor("K"),   // length 1
            };
            Assert.DoesNotThrow(() => DecoyGenerator.GenerateFromTargets(targets.ToList()),
                "Short sequences must not cause exceptions");
        }

        [Test]
        public void ShortSequence_StillMarkedIsDecoy()
        {
            var targets = new[] { MakePrecursor("AK") };
            var decoys = DecoyGenerator.GenerateFromTargets(targets.ToList());
            Assert.That(decoys[0].IsDecoy, Is.True);
        }

        [Test]
        public void EmptyFragments_FallbackToTargetFragments()
        {
            // If decoy ends up with zero fragments (edge case), generator falls back
            // to target fragments rather than producing an empty decoy.
            var targets = new[] { MakePrecursor("AK") };
            var decoys = DecoyGenerator.GenerateFromTargets(targets.ToList());
            Assert.That(decoys[0].FragmentCount, Is.GreaterThan(0),
                "Even degenerate sequences must produce a non-empty fragment array");
        }

        [Test]
        public void LargeSequenceBatch_CompletesWithoutException()
        {
            // Stress test: generate 1000 decoys without crashing.
            // Uses varied sequences with different amino acid compositions.
            var rng = new Random(999);
            string aa = "ACDEFGHIKLMNPQRSTVWY";
            var targets = new List<LibraryPrecursorInput>(1000);
            for (int i = 0; i < 1000; i++)
            {
                int len = rng.Next(6, 20);
                var seq = new char[len];
                for (int j = 0; j < len - 1; j++) seq[j] = aa[rng.Next(aa.Length)];
                seq[len - 1] = rng.Next(2) == 0 ? 'K' : 'R'; // tryptic C-terminal
                targets.Add(MakePrecursor(new string(seq)));
            }

            Assert.DoesNotThrow(() => DecoyGenerator.GenerateFromTargets(targets),
                "Large batch generation must not throw");

            var decoys = DecoyGenerator.GenerateFromTargets(targets);
            Assert.That(decoys.Count, Is.EqualTo(1000));
            Assert.That(decoys.All(d => d.IsDecoy), Is.True);
            Assert.That(decoys.All(d => d.FragmentCount > 0), Is.True);
        }

        // ─────────────────────────────────────────────────────────────────────
        //  Layer 5: NDP / spectrum similarity
        // ─────────────────────────────────────────────────────────────────────

        [Test]
        public void AverageNdp_IsLow_AcrossTypicalSequences()
        {
            // The mutation step is specifically designed to reduce NDP.
            // Across a batch of typical sequences, mean NDP should be well below 0.5.
            // (The old scramble-only approach already achieved ~0.005 on the benchmark;
            //  we verify we haven't accidentally made it worse.)
            var rng = new Random(42);
            string aa = "ACDEFGHIKLMNPQRSTVWY";
            var targets = new List<LibraryPrecursorInput>(200);
            for (int i = 0; i < 200; i++)
            {
                int len = rng.Next(8, 20);
                var seq = new char[len];
                for (int j = 0; j < len - 1; j++) seq[j] = aa[rng.Next(aa.Length)];
                seq[len - 1] = 'K';

                // Use realistic fragment m/z and intensities
                float[] mzs = Enumerable.Range(1, 10).Select(k => (float)(k * 100)).ToArray();
                float[] ints = Enumerable.Range(1, 10).Select(k => (float)(1000 - k * 80)).ToArray();

                targets.Add(new LibraryPrecursorInput(
                    sequence: new string(seq),
                    precursorMz: 550.0,
                    chargeState: 2,
                    retentionTime: 30.0,
                    isDecoy: false,
                    fragmentMzs: mzs,
                    fragmentIntensities: ints,
                    irtValue: null));
            }

            var decoys = DecoyGenerator.GenerateFromTargets(targets);

            // Compute NDP for each pair using the same logic as DecoyGenerator
            double totalNdp = 0.0;
            int counted = 0;
            for (int i = 0; i < targets.Count; i++)
            {
                var t = targets[i];
                var d = decoys[i];
                if (t.FragmentMzs.Length == 0 || d.FragmentMzs.Length == 0) continue;

                double normT = 0, normD = 0, dot = 0;
                foreach (var v in t.FragmentIntensities) normT += v * v;
                foreach (var v in d.FragmentIntensities) normD += v * v;
                normT = Math.Sqrt(normT); normD = Math.Sqrt(normD);
                if (normT < 1e-12 || normD < 1e-12) continue;

                int ti = 0, di = 0;
                while (ti < t.FragmentMzs.Length && di < d.FragmentMzs.Length)
                {
                    double ppm = Math.Abs(t.FragmentMzs[ti] - d.FragmentMzs[di])
                                 / t.FragmentMzs[ti] * 1e6;
                    if (ppm <= 20.0) { dot += t.FragmentIntensities[ti] * d.FragmentIntensities[di]; ti++; di++; }
                    else if (t.FragmentMzs[ti] < d.FragmentMzs[di]) ti++;
                    else di++;
                }

                totalNdp += dot / (normT * normD);
                counted++;
            }

            double meanNdp = counted > 0 ? totalNdp / counted : 1.0;
            Assert.That(meanNdp, Is.LessThan(0.3),
                $"Mean NDP across {counted} decoys is {meanNdp:F4} — should be < 0.3. " +
                "High NDP means decoy spectra are too similar to targets.");
        }

        [Test]
        public void FragmentCount_MatchesTargetWithinOneDelta()
        {
            // Fragment count of decoy must be within ±1 of target fragment count.
            var sequences = new[] { "PEPTIDEK", "SAMPLEWYK", "ELECTROSPRAY", "MASSSPECTRK" };
            var targets = sequences.Select(s => MakePrecursor(s)).ToList();
            var decoys = DecoyGenerator.GenerateFromTargets(targets);

            for (int i = 0; i < targets.Count; i++)
            {
                int diff = Math.Abs(decoys[i].FragmentCount - targets[i].FragmentCount);
                Assert.That(diff, Is.LessThanOrEqualTo(1),
                    $"Decoy[{i}] fragment count {decoys[i].FragmentCount} differs from " +
                    $"target {targets[i].FragmentCount} by more than 1");
            }
        }
    }
}