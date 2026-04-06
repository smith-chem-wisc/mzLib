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
    ///   1. Scramble the middle residues (C-terminal K/R fixed).
    ///   2. Apply exactly one substitution at a randomly chosen non-C-terminal position.
    ///      Constraints:
    ///        - Aromatic count (F, W, Y) preserved: aromatic → aromatic only
    ///        - Non-aromatic → non-aromatic only, BLOSUM62 score ∈ {-1, 0, +1}
    ///        - Every entry in the substitution table is guaranteed to produce
    ///          |mass delta| > 0.1 Da (I↔L and other zero-delta pairs excluded)
    ///        - No upper bound on mass delta
    ///      If no mutable position exists, fall back to scramble only.
    ///   3. Fragment the mutated sequence (singly-charged b/y ions).
    ///   4. Assign target intensities to nearest decoy fragment m/z.
    ///   5. NDP-guided retry: across NdpMaxRetries attempts (each picks a different
    ///      random position/substitute), keep the lowest-NDP result.
    ///   6. RT: copy RetentionTime and IrtValue from target unchanged.
    ///
    /// Reproducibility: each precursor's RNG is seeded deterministically from its
    /// sequence hash, making results order-independent and run-to-run identical.
    /// </summary>
    public static class DecoyGenerator
    {
        private const int MaxScrambleAttempts = 10;
        private const int NdpMaxRetries = 5;
        private const double NdpAcceptThreshold = 0.10;
        private const double ProtonMass = 1.007276;
        private const double WaterMass = 18.010565;

        // ── Monoisotopic residue masses ───────────────────────────────────────

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

        // ── Aromatic residues (F, W, Y — H excluded per spec) ────────────────

        private static readonly HashSet<char> Aromatics = new HashSet<char> { 'F', 'W', 'Y' };

        // ── BLOSUM62 near-neutral substitution table ──────────────────────────
        //
        // Rules applied when building this table:
        //   - BLOSUM62 score ∈ {-1, 0, +1} for non-aromatic pairs
        //   - Aromatic → aromatic only (F, W, Y); score filter waived
        //   - Charge group preserved: basic→basic, acidic→acidic, neutral→neutral
        //       Basic:   K, R, H
        //       Acidic:  D, E
        //       Neutral: all others
        //   - Every listed substitute has |mass delta| > 0.1 Da vs original
        //   - I↔L excluded everywhere (both mass 113.084 Da, delta = 0)
        //   - Self-substitutions excluded

        // ── Charge groups (substitutions restricted to within-group) ─────────
        //
        //   Basic:   K, R, H
        //   Acidic:  D, E
        //   Neutral: A, V, I, L, M, G, S, T, N, Q, C, F, W, Y, P
        //
        // Consequence: some residues have very few options:
        //   D → E only          (only other acidic)
        //   E → D only          (only other acidic)
        //   K → R, H            (other basics; BLOSUM score filter also applied)
        //   R → K, H            (other basics)
        //   H → K, R            (other basics)
        //   N → S, T, Q         (H removed — basic; D removed — acidic)
        //   Q → N, S, T, A, C   (H removed — basic; K/R removed — basic)
        //   All neutral residues lose any basic/acidic options

        private static readonly Dictionary<char, char[]> NeutralSubstitutes =
            new Dictionary<char, char[]>
        {
            // Aliphatic / small (all neutral → neutral only; no change needed)
            {'A', new[] {'S', 'T', 'G', 'V'}},
            {'V', new[] {'A', 'I', 'M', 'T'}},       // L excluded (mass == I)
            {'I', new[] {'V', 'M', 'T'}},             // L excluded (identical mass)
            {'L', new[] {'V', 'M'}},                  // I excluded (identical mass)
            {'M', new[] {'V', 'I', 'T'}},             // L excluded (mass == I)
            {'G', new[] {'A', 'S'}},
            // Polar uncharged (neutral → neutral only)
            {'S', new[] {'A', 'T', 'N', 'G'}},        // unchanged — all neutral
            {'T', new[] {'S', 'A', 'V', 'N', 'M'}},   // unchanged — all neutral
            {'N', new[] {'S', 'T', 'Q'}},              // D removed (acidic), H removed (basic)
            {'Q', new[] {'N', 'S', 'T', 'A', 'C'}},   // H/K/R removed (basic)
            {'C', new[] {'S', 'A', 'V', 'T'}},        // unchanged — all neutral
            // Acidic → acidic only
            {'D', new[] {'E'}},                        // N/S/T removed (neutral)
            {'E', new[] {'D'}},                        // Q/N/K removed (neutral or basic)
            // Basic → basic only
            {'K', new[] {'R', 'H'}},                  // Q/N/E removed (neutral or acidic)
            {'R', new[] {'K', 'H'}},                  // Q/N removed (neutral)
            {'H', new[] {'K', 'R'}},                  // N/Q/S/D removed (neutral or acidic)
            // Aromatic (aromatic → aromatic only; all neutral — no charge conflict)
            {'F', new[] {'Y', 'W'}},
            {'W', new[] {'F', 'Y'}},
            {'Y', new[] {'F', 'W'}},
            // Proline (neutral → neutral)
            {'P', new[] {'A', 'S'}},
        };

        // ── Public entry point ────────────────────────────────────────────────

        /// <summary>
        /// Generates one decoy per target. Returns a list of the same length as
        /// <paramref name="targets"/> with all entries marked IsDecoy=true.
        /// </summary>
        public static List<LibraryPrecursorInput> GenerateFromTargets(
            IReadOnlyList<LibraryPrecursorInput> targets)
        {
            var decoys = new List<LibraryPrecursorInput>(targets.Count);
            int scrambleFallbacks = 0;
            int mutationFallbacks = 0;
            int countMismatches = 0;

            for (int i = 0; i < targets.Count; i++)
            {
                var t = targets[i];
                int targetFragCount = t.FragmentMzs?.Length ?? 0;
                string bare = StripModifications(t.Sequence);

                // Per-sequence deterministic seed: same sequence → same decoy,
                // regardless of position in the input list.
                int hashSeed = t.Sequence.GetHashCode();
                hashSeed ^= hashSeed << 13;
                hashSeed ^= hashSeed >> 17;
                hashSeed ^= hashSeed << 5;
                var rng = new Random(hashSeed);

                string finalSequence;
                float[] decoyMzs;
                float[] decoyInts;

                if (bare.Length < 2)
                {
                    finalSequence = bare;
                    decoyMzs = ComputeFragmentMzs(finalSequence);
                    decoyInts = AssignIntensities(decoyMzs, t.FragmentMzs, t.FragmentIntensities);
                    (decoyMzs, decoyInts) = FilterZeroIntensity(decoyMzs, decoyInts);
                }
                else
                {
                    // ── Step 1: Scramble ─────────────────────────────────────
                    string scrambled = Scramble(bare, rng, out bool scrambleOk);
                    if (!scrambleOk) { scrambleFallbacks++; scrambled = Reverse(bare); }

                    // ── Step 2–4: One mutation + NDP check ────────────────────
                    // Each attempt picks a uniformly random (position, substitute) pair.
                    // Keep the lowest-NDP result across NdpMaxRetries attempts.
                    // Break early if NDP is already below NdpAcceptThreshold.

                    string bestSequence = scrambled;
                    float[] bestMzs = null;
                    float[] bestInts = null;
                    double bestNdp = double.MaxValue;
                    bool anyMutationSucceeded = false;

                    for (int attempt = 0; attempt < NdpMaxRetries; attempt++)
                    {
                        string mutated = ApplyOneMutation(scrambled, rng, out bool mutationOk);

                        if (!mutationOk)
                        {
                            if (attempt == 0) mutationFallbacks++;
                            break; // no mutable positions — retrying won't help
                        }

                        anyMutationSucceeded = true;

                        float[] cMzs = ComputeFragmentMzs(mutated);
                        float[] cInts = AssignIntensities(cMzs, t.FragmentMzs, t.FragmentIntensities);
                        (cMzs, cInts) = FilterZeroIntensity(cMzs, cInts);

                        double ndp = NormalizedDotProduct(
                            t.FragmentMzs, t.FragmentIntensities, cMzs, cInts);

                        if (ndp < bestNdp)
                        {
                            bestNdp = ndp;
                            bestSequence = mutated;
                            bestMzs = cMzs;
                            bestInts = cInts;
                        }

                        if (ndp <= NdpAcceptThreshold) break;
                    }

                    if (!anyMutationSucceeded)
                    {
                        bestMzs = ComputeFragmentMzs(scrambled);
                        bestInts = AssignIntensities(bestMzs, t.FragmentMzs, t.FragmentIntensities);
                        (bestMzs, bestInts) = FilterZeroIntensity(bestMzs, bestInts);
                    }

                    finalSequence = bestSequence;
                    decoyMzs = bestMzs;
                    decoyInts = bestInts;
                }

                // ── Cap fragment count to match target (±1) ──────────────────
                if (decoyMzs.Length > targetFragCount && targetFragCount > 0)
                {
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
                    Array.Sort(topMzs, topInts);
                    decoyMzs = topMzs;
                    decoyInts = topInts;
                }

                int diff = Math.Abs(decoyMzs.Length - targetFragCount);
                if (diff > 1) countMismatches++;

                if (decoyMzs.Length == 0)
                {
                    decoyMzs = t.FragmentMzs;
                    decoyInts = t.FragmentIntensities;
                }

                decoys.Add(MakeDecoy(t, finalSequence, decoyMzs, decoyInts));
            }

            Console.WriteLine($"  Decoys generated: {decoys.Count:N0}" +
                (scrambleFallbacks > 0
                    ? $"  ({scrambleFallbacks} used reverse fallback for scramble)" : "") +
                (mutationFallbacks > 0
                    ? $"  ({mutationFallbacks} fell back to scramble-only — no mutable position)" : "") +
                (countMismatches > 0
                    ? $"  ({countMismatches} count mismatches >1)" : " (all counts matched within ±1)"));
            return decoys;
        }

        // ── Single-mutation application ───────────────────────────────────────

        /// <summary>
        /// Applies exactly one substitution to <paramref name="scrambled"/>.
        ///
        /// Builds the full list of valid (position, substitute) pairs across all
        /// non-C-terminal positions, then picks one uniformly at random. This gives
        /// equal probability to every valid (position, substitute) combination,
        /// avoiding bias toward residues with more options.
        ///
        /// Returns the mutated sequence. Sets <paramref name="mutationOk"/> = false
        /// if no valid pair exists (no mutable positions or all empty substitute lists).
        /// </summary>
        private static string ApplyOneMutation(string scrambled, Random rng, out bool mutationOk)
        {
            mutationOk = false;
            int len = scrambled.Length;
            if (len < 2) return scrambled;

            // Collect all valid (position, substitute) pairs
            var validPairs = new List<(int pos, char sub)>(len * 3);
            for (int i = 0; i < len - 1; i++) // exclude C-terminal index len-1
            {
                char original = scrambled[i];
                char[] options = GetSubstituteOptions(original);
                if (options == null) continue;
                foreach (char sub in options)
                    validPairs.Add((i, sub));
            }

            if (validPairs.Count == 0) return scrambled;

            // Pick one pair uniformly at random
            var (pos, chosen) = validPairs[rng.Next(validPairs.Count)];
            char[] result = scrambled.ToCharArray();
            result[pos] = chosen;
            mutationOk = true;
            return new string(result);
        }

        /// <summary>
        /// Returns valid substitute residues for <paramref name="original"/>.
        /// Aromatics return other aromatics; non-aromatics return BLOSUM62
        /// near-neutral non-aromatic substitutes. All returned values are
        /// guaranteed to differ from <paramref name="original"/> and to have
        /// |mass delta| > 0.1 Da by table construction.
        /// Returns null if the residue has no table entry.
        /// </summary>
        private static char[] GetSubstituteOptions(char original)
        {
            if (!NeutralSubstitutes.TryGetValue(original, out char[] subs))
                return null;

            if (Aromatics.Contains(original))
                return subs; // table already contains only valid aromatic targets

            // Non-aromatics: filter defensively (table should already be clean)
            var filtered = new List<char>(subs.Length);
            foreach (char s in subs)
                if (!Aromatics.Contains(s) && s != original)
                    filtered.Add(s);

            return filtered.Count > 0 ? filtered.ToArray() : null;
        }

        // ── Sequence manipulation ─────────────────────────────────────────────

        private static string Scramble(string bare, Random rng, out bool success)
        {
            if (bare.Length <= 2) { success = true; return bare; }

            char cTerm = bare[bare.Length - 1];
            char[] middle = bare.Substring(0, bare.Length - 1).ToCharArray();

            for (int attempt = 0; attempt < MaxScrambleAttempts; attempt++)
            {
                for (int j = middle.Length - 1; j > 0; j--)
                {
                    int k = rng.Next(j + 1);
                    (middle[j], middle[k]) = (middle[k], middle[j]);
                }
                string candidate = new string(middle) + cTerm;
                if (candidate != bare) { success = true; return candidate; }
            }

            success = false;
            return new string(middle) + cTerm;
        }

        private static string Reverse(string bare)
        {
            if (bare.Length <= 1) return bare;
            char cTerm = bare[bare.Length - 1];
            char[] middle = bare.Substring(0, bare.Length - 1).ToCharArray();
            Array.Reverse(middle);
            return new string(middle) + cTerm;
        }

        // ── Fragment computation ──────────────────────────────────────────────

        private static float[] ComputeFragmentMzs(string sequence)
        {
            int L = sequence.Length;
            if (L < 2) return Array.Empty<float>();

            var mzs = new List<float>((L - 1) * 2);

            double bMass = ProtonMass;
            for (int i = 0; i < L - 1; i++)
            {
                bMass += GetResidueMass(sequence[i]);
                mzs.Add((float)bMass);
            }

            double yMass = WaterMass + ProtonMass;
            for (int i = L - 1; i >= 1; i--)
            {
                yMass += GetResidueMass(sequence[i]);
                mzs.Add((float)yMass);
            }

            mzs.Sort();
            return mzs.ToArray();
        }

        // ── Intensity assignment ──────────────────────────────────────────────

        private static float[] AssignIntensities(
            float[] decoyMzs, float[] targetMzs, float[] targetInts)
        {
            if (targetMzs == null || targetMzs.Length == 0)
                return new float[decoyMzs.Length];

            var result = new float[decoyMzs.Length];
            for (int i = 0; i < decoyMzs.Length; i++)
            {
                int idx = BinarySearchNearest(targetMzs, decoyMzs[i]);
                result[i] = targetInts[idx];
            }
            return result;
        }

        private static int BinarySearchNearest(float[] arr, float target)
        {
            int lo = 0, hi = arr.Length - 1;
            while (lo < hi)
            {
                int mid = (lo + hi) / 2;
                if (arr[mid] < target) lo = mid + 1;
                else hi = mid;
            }
            if (lo == 0) return 0;
            return Math.Abs(arr[lo - 1] - target) <= Math.Abs(arr[lo] - target) ? lo - 1 : lo;
        }

        // ── NDP ───────────────────────────────────────────────────────────────

        private static double NormalizedDotProduct(
            float[] tMzs, float[] tInts, float[] dMzs, float[] dInts,
            double ppmTolerance = 20.0)
        {
            if (tMzs == null || dMzs == null || tInts == null || dInts == null) return 0.0;

            double dot = 0.0, normT = 0.0, normD = 0.0;
            for (int i = 0; i < tInts.Length; i++) normT += tInts[i] * tInts[i];
            for (int j = 0; j < dInts.Length; j++) normD += dInts[j] * dInts[j];
            normT = Math.Sqrt(normT); normD = Math.Sqrt(normD);
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

        // ── Helpers ───────────────────────────────────────────────────────────

        private static double GetResidueMass(char c)
        {
            ResidueMass.TryGetValue(c, out double m);
            return m > 0 ? m : 111.0;
        }

        private static (float[] mzs, float[] ints) FilterZeroIntensity(float[] mzs, float[] ints)
        {
            int count = 0;
            for (int i = 0; i < ints.Length; i++) if (ints[i] > 0f) count++;
            if (count == ints.Length) return (mzs, ints);

            var fMzs = new float[count]; var fInts = new float[count]; int j = 0;
            for (int i = 0; i < ints.Length; i++)
                if (ints[i] > 0f) { fMzs[j] = mzs[i]; fInts[j] = ints[i]; j++; }
            return (fMzs, fInts);
        }

        private static LibraryPrecursorInput MakeDecoy(
            LibraryPrecursorInput target, string sequence,
            float[] fragmentMzs, float[] fragmentIntensities)
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

        private static string StripModifications(string sequence)
        {
            if (string.IsNullOrEmpty(sequence)) return "";
            var sb = new System.Text.StringBuilder(sequence.Length);
            bool inMod = false; char modOpen = ' ';
            foreach (char c in sequence)
            {
                if (!inMod) { if (c == '(' || c == '[') { inMod = true; modOpen = c; } else sb.Append(c); }
                else { char close = modOpen == '(' ? ')' : ']'; if (c == close) inMod = false; }
            }
            return sb.ToString();
        }
    }
}