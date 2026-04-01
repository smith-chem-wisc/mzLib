// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0
//
// Location: MassSpectrometry/Dia/DecoyGenerator.cs

using System;
using System.Collections.Generic;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Generates decoy LibraryPrecursorInputs from a target list.
    ///
    /// Algorithm per precursor:
    ///   1. Scramble the sequence, keeping the C-terminal residue fixed
    ///      (always K or R for tryptic peptides).
    ///   2. Compute theoretical singly-charged b and y fragment m/z values
    ///      for the scrambled sequence.
    ///   3. For each theoretical decoy fragment m/z, find the closest target
    ///      fragment m/z and assign its intensity. This gives the decoy a
    ///      realistic intensity distribution without matching real DIA signal.
    ///   4. Keep precursor m/z, charge, and RT identical to the target so
    ///      the decoy competes in the same isolation window and RT region.
    ///   5. Deduplicate scrambled sequences — if the scramble produces the
    ///      forward sequence (palindromes), re-scramble up to MaxScrambleAttempts.
    /// </summary>
    public static class DecoyGenerator
    {
        private const int MaxScrambleAttempts = 10;

        /// <summary>Maximum NDP retries when decoy spectrum is too similar to target.</summary>
        private const int NdpMaxRetries = 5;

        /// <summary>NDP threshold below which a scramble is accepted immediately.
        /// If all retries exceed this, the lowest-NDP attempt is kept.</summary>
        private const double NdpAcceptThreshold = 0.10;
        private const double ProtonMass = 1.007276;
        private const double WaterMass = 18.010565;

        // Monoisotopic residue masses (standard amino acids)
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

        /// <summary>
        /// Generates one decoy per target. Returns a list of the same length as
        /// <paramref name="targets"/> with all entries marked IsDecoy=true.
        /// </summary>
        public static List<LibraryPrecursorInput> GenerateFromTargets(
            IReadOnlyList<LibraryPrecursorInput> targets,
            int randomSeed = 42)
        {
            var rng = new Random(randomSeed);
            var decoys = new List<LibraryPrecursorInput>(targets.Count);
            int scrambleFallbacks = 0;
            int countMismatches = 0;

            for (int i = 0; i < targets.Count; i++)
            {
                var t = targets[i];
                int targetFragCount = t.FragmentMzs?.Length ?? 0;

                // ── Step 1: Scramble the target sequence ─────────────────────
                // Strip modifications to get bare sequence, keep C-terminal fixed.
                string bare = StripModifications(t.Sequence);

                // ── Steps 2-5 with NDP-guided retry ──────────────────────────
                // After each scramble attempt, compute NDP between the resulting
                // decoy spectrum and the target. If NDP is too high, try again.
                // Keep the attempt with the lowest NDP across all retries.
                string scrambled;
                float[] decoyMzs;
                float[] decoyInts;

                if (bare.Length < 2)
                {
                    // Too short to scramble meaningfully
                    scrambled = bare;
                    decoyMzs = ComputeFragmentMzs(scrambled);
                    decoyInts = AssignIntensities(decoyMzs, t.FragmentMzs, t.FragmentIntensities);
                    (decoyMzs, decoyInts) = FilterZeroIntensity(decoyMzs, decoyInts);
                }
                else
                {
                    // Best candidate across all attempts
                    string bestScrambled = null;
                    float[] bestMzs = null;
                    float[] bestInts = null;
                    double bestNdp = double.MaxValue;

                    for (int attempt = 0; attempt < NdpMaxRetries; attempt++)
                    {
                        // Step 2: scramble (guaranteed different from forward sequence)
                        string candidate = Scramble(bare, rng, out bool ok);
                        if (!ok) { scrambleFallbacks++; candidate = Reverse(bare); }

                        // Step 3: fragment
                        float[] cMzs = ComputeFragmentMzs(candidate);

                        // Step 4: assign intensities
                        float[] cInts = AssignIntensities(cMzs, t.FragmentMzs, t.FragmentIntensities);

                        // Step 5: remove zero-intensity fragments
                        (cMzs, cInts) = FilterZeroIntensity(cMzs, cInts);

                        // Compute NDP against target
                        double ndp = NormalizedDotProduct(
                            t.FragmentMzs, t.FragmentIntensities,
                            cMzs, cInts);

                        // Keep best (lowest NDP) candidate
                        if (ndp < bestNdp)
                        {
                            bestNdp = ndp;
                            bestScrambled = candidate;
                            bestMzs = cMzs;
                            bestInts = cInts;
                        }

                        // Accept immediately if NDP is below threshold
                        if (ndp <= NdpAcceptThreshold) break;
                    }

                    scrambled = bestScrambled;
                    decoyMzs = bestMzs;
                    decoyInts = bestInts;
                }

                // ── Step 6: Check that count of (mz, intensity) pairs is equal (+/-1)
                // If decoy has more fragments than target (from the larger theoretical
                // b/y set), trim to target count by keeping highest-intensity fragments.
                // If decoy has fewer (some theoreticals had no assignable intensity),
                // accept the smaller count — within 1 is fine.
                if (decoyMzs.Length > targetFragCount && targetFragCount > 0)
                {
                    // Sort by intensity descending, keep top targetFragCount
                    var idx = new int[decoyMzs.Length];
                    for (int j = 0; j < idx.Length; j++) idx[j] = j;
                    Array.Sort(idx, (a, b) => decoyInts[b].CompareTo(decoyInts[a]));

                    var topMzs = new float[targetFragCount];
                    var topInts = new float[targetFragCount];
                    for (int j = 0; j < targetFragCount; j++)
                    {
                        topMzs[j] = decoyMzs[idx[j]];
                        topInts[j] = decoyInts[idx[j]];
                    }
                    Array.Sort(topMzs, topInts);  // re-sort by m/z ascending
                    decoyMzs = topMzs;
                    decoyInts = topInts;
                }

                // Count check: flag if still off by more than 1
                int diff = Math.Abs(decoyMzs.Length - targetFragCount);
                if (diff > 1) countMismatches++;

                // Fallback: if decoy ended up empty, use target fragments
                if (decoyMzs.Length == 0)
                {
                    decoyMzs = t.FragmentMzs;
                    decoyInts = t.FragmentIntensities;
                }

                decoys.Add(MakeDecoy(t, scrambled, decoyMzs, decoyInts));
            }

            Console.WriteLine($"  Decoys generated: {decoys.Count:N0}" +
                              (scrambleFallbacks > 0 ? $"  ({scrambleFallbacks} used reverse fallback)" : "") +
                              (countMismatches > 0 ? $"  ({countMismatches} count mismatches >1)" : " (all counts matched within ±1)"));
            return decoys;
        }

        // ── Sequence manipulation ─────────────────────────────────────────────

        /// <summary>
        /// Scrambles the middle residues (everything except the last character),
        /// keeping the C-terminal residue fixed. Retries up to MaxScrambleAttempts
        /// to avoid producing the original sequence.
        /// </summary>
        private static string Scramble(string bare, Random rng, out bool success)
        {
            if (bare.Length <= 2)
            {
                // Only 1 middle residue — can't scramble further
                success = true;
                return bare;
            }

            char cTerm = bare[bare.Length - 1];
            char[] middle = bare.Substring(0, bare.Length - 1).ToCharArray();

            for (int attempt = 0; attempt < MaxScrambleAttempts; attempt++)
            {
                // Fisher-Yates shuffle
                for (int j = middle.Length - 1; j > 0; j--)
                {
                    int k = rng.Next(j + 1);
                    (middle[j], middle[k]) = (middle[k], middle[j]);
                }

                string candidate = new string(middle) + cTerm;
                if (candidate != bare)
                {
                    success = true;
                    return candidate;
                }
            }

            success = false;
            return new string(middle) + cTerm; // return last attempt even if same
        }

        private static string Reverse(string bare)
        {
            if (bare.Length <= 1) return bare;
            char cTerm = bare[bare.Length - 1];
            char[] middle = bare.Substring(0, bare.Length - 1).ToCharArray();
            Array.Reverse(middle);
            return new string(middle) + cTerm;
        }

        // ── Theoretical fragment m/z computation ──────────────────────────────

        /// <summary>
        /// Computes singly-charged b and y fragment m/z values for a bare sequence.
        /// Returns them sorted ascending.
        /// </summary>
        private static float[] ComputeFragmentMzs(string sequence)
        {
            int L = sequence.Length;
            if (L < 2) return Array.Empty<float>();

            var mzs = new List<float>((L - 1) * 2);

            // b ions: N-terminal fragments
            // b_i = sum(residues 0..i-1) + proton
            double bMass = ProtonMass;
            for (int i = 0; i < L - 1; i++)
            {
                if (!ResidueMass.TryGetValue(sequence[i], out double rm))
                    rm = 111.0; // unknown residue — use average
                bMass += rm;
                mzs.Add((float)bMass);
            }

            // y ions: C-terminal fragments
            // y_i = sum(residues L-i..L-1) + water + proton
            double yMass = WaterMass + ProtonMass;
            for (int i = L - 1; i >= 1; i--)
            {
                if (!ResidueMass.TryGetValue(sequence[i], out double rm))
                    rm = 111.0;
                yMass += rm;
                mzs.Add((float)yMass);
            }

            mzs.Sort();
            return mzs.ToArray();
        }

        // ── Intensity assignment ──────────────────────────────────────────────

        /// <summary>
        /// For each decoy fragment m/z, finds the nearest target fragment m/z
        /// and assigns its intensity. Uses binary search for efficiency.
        /// </summary>
        private static float[] AssignIntensities(
            float[] decoyMzs,
            float[] targetMzs,
            float[] targetInts)
        {
            if (targetMzs == null || targetMzs.Length == 0)
                return new float[decoyMzs.Length]; // all zero

            var result = new float[decoyMzs.Length];
            for (int i = 0; i < decoyMzs.Length; i++)
            {
                float dMz = decoyMzs[i];
                int idx = BinarySearchNearest(targetMzs, dMz);
                result[i] = targetInts[idx];
            }
            return result;
        }

        /// <summary>
        /// Binary search returning the index of the element nearest to target.
        /// Assumes arr is sorted ascending.
        /// </summary>
        private static int BinarySearchNearest(float[] arr, float target)
        {
            int lo = 0, hi = arr.Length - 1;
            while (lo < hi)
            {
                int mid = (lo + hi) / 2;
                if (arr[mid] < target) lo = mid + 1;
                else hi = mid;
            }
            // lo is now the first index >= target
            if (lo == 0) return 0;
            // Compare lo and lo-1
            float distLo = Math.Abs(arr[lo] - target);
            float distLo1 = Math.Abs(arr[lo - 1] - target);
            return distLo1 <= distLo ? lo - 1 : lo;
        }

        // ── Helpers ───────────────────────────────────────────────────────────

        private static (float[] mzs, float[] ints) FilterZeroIntensity(
            float[] mzs, float[] ints)
        {
            int count = 0;
            for (int i = 0; i < ints.Length; i++)
                if (ints[i] > 0f) count++;

            if (count == ints.Length) return (mzs, ints);

            var fMzs = new float[count];
            var fInts = new float[count];
            int j = 0;
            for (int i = 0; i < ints.Length; i++)
            {
                if (ints[i] > 0f)
                {
                    fMzs[j] = mzs[i];
                    fInts[j] = ints[i];
                    j++;
                }
            }
            return (fMzs, fInts);
        }

        /// <summary>
        /// Normalized dot product between two spectra (sorted ascending by m/z).
        /// Uses two-pointer matching within ppmTolerance. Range [0, 1].
        /// </summary>
        private static double NormalizedDotProduct(
            float[] tMzs, float[] tInts,
            float[] dMzs, float[] dInts,
            double ppmTolerance = 20.0)
        {
            if (tMzs == null || dMzs == null || tInts == null || dInts == null) return 0.0;
            double dot = 0.0, normT = 0.0, normD = 0.0;
            for (int i = 0; i < tInts.Length; i++) normT += tInts[i] * tInts[i];
            for (int j = 0; j < dInts.Length; j++) normD += dInts[j] * dInts[j];
            normT = Math.Sqrt(normT);
            normD = Math.Sqrt(normD);
            if (normT < 1e-12 || normD < 1e-12) return 0.0;
            int ti = 0, di = 0;
            while (ti < tMzs.Length && di < dMzs.Length)
            {
                double ppm = Math.Abs(tMzs[ti] - dMzs[di]) / tMzs[ti] * 1e6;
                if (ppm <= ppmTolerance) { dot += tInts[ti] * dInts[di]; ti++; di++; }
                else if (tMzs[ti] < dMzs[di]) ti++;
                else di++;
            }
            return dot / (normT * normD);
        }

        private static LibraryPrecursorInput MakeDecoy(
            LibraryPrecursorInput target,
            string sequence,
            float[] fragmentMzs,
            float[] fragmentIntensities)
        {
            return new LibraryPrecursorInput(
                sequence: sequence,
                precursorMz: target.PrecursorMz,
                chargeState: target.ChargeState,
                retentionTime: target.RetentionTime,
                isDecoy: true,
                fragmentMzs: fragmentMzs,
                fragmentIntensities: fragmentIntensities,
                irtValue: target.IrtValue);
        }

        /// <summary>
        /// Strips modification annotations like (UniMod:21) or [UniMod:21]
        /// to return the bare amino acid sequence.
        /// </summary>
        private static string StripModifications(string sequence)
        {
            if (string.IsNullOrEmpty(sequence)) return "";
            // Remove (UniMod:N) and [UniMod:N] patterns
            var sb = new System.Text.StringBuilder(sequence.Length);
            bool inMod = false;
            char modOpen = ' ';
            foreach (char c in sequence)
            {
                if (!inMod)
                {
                    if (c == '(' || c == '[') { inMod = true; modOpen = c; }
                    else sb.Append(c);
                }
                else
                {
                    char close = modOpen == '(' ? ')' : ']';
                    if (c == close) inMod = false;
                }
            }
            return sb.ToString();
        }
    }
}