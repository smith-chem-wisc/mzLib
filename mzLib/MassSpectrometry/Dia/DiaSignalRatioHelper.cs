// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Phase 13, Prompt 5: Signal ratio deviation features.
    /// 
    /// Compares the observed fragment intensity pattern at the apex scan against
    /// the library spectrum pattern using log2 ratio deviations. Each fragment's
    /// observed fractional intensity is compared to its library fractional intensity.
    /// 
    /// True peptide signals match the library pattern closely (low deviation).
    /// Interfered or random matches show high deviations.
    /// </summary>
    internal static class DiaSignalRatioHelper
    {
        /// <summary>
        /// Computes signal ratio deviation features and populates the result.
        /// 
        /// Algorithm:
        ///   1. Normalize both observed and library intensities to fractional intensities
        ///   2. For each fragment with nonzero signal in both, compute log2(observed_frac / library_frac)
        ///   3. Compute mean-absolute, max-absolute, and std of these signed deviations
        /// 
        /// Requires ≥3 fragments with valid intensities in both observed and library.
        /// </summary>
        /// <param name="apexIntensities">Per-fragment intensities at the apex scan.</param>
        /// <param name="libraryIntensities">Library fragment intensities (parallel array).</param>
        /// <param name="fragmentCount">Number of fragment ions.</param>
        /// <param name="result">DiaSearchResult to populate.</param>
        public static void ComputeSignalRatioFeatures(
            float[] apexIntensities,
            float[] libraryIntensities,
            int fragmentCount,
            DiaSearchResult result)
        {
            if (fragmentCount < 3)
                return;

            // Compute total intensities for normalization
            float obsTotal = 0f;
            float libTotal = 0f;
            for (int f = 0; f < fragmentCount; f++)
            {
                obsTotal += apexIntensities[f];
                libTotal += libraryIntensities[f];
            }

            if (obsTotal <= 0f || libTotal <= 0f)
                return;

            // Compute log2 ratio deviations for fragments with signal in both
            Span<float> deviations = stackalloc float[fragmentCount];
            int nValid = 0;

            for (int f = 0; f < fragmentCount; f++)
            {
                float obs = apexIntensities[f];
                float lib = libraryIntensities[f];

                if (obs <= 0f || lib <= 0f)
                    continue;

                float obsFrac = obs / obsTotal;
                float libFrac = lib / libTotal;

                // log2(observed_fraction / library_fraction)
                float logRatio = MathF.Log2(obsFrac / libFrac);
                deviations[nValid++] = logRatio;
            }

            if (nValid < 3)
                return;

            // Mean absolute deviation
            float sumAbs = 0f;
            float maxAbs = 0f;
            float sumSigned = 0f;
            for (int i = 0; i < nValid; i++)
            {
                float absVal = MathF.Abs(deviations[i]);
                sumAbs += absVal;
                if (absVal > maxAbs) maxAbs = absVal;
                sumSigned += deviations[i];
            }
            result.MeanSignalRatioDev = sumAbs / nValid;
            result.MaxSignalRatioDev = maxAbs;

            // Std of signed deviations (population)
            float mean = sumSigned / nValid;
            float sumSqDev = 0f;
            for (int i = 0; i < nValid; i++)
            {
                float dev = deviations[i] - mean;
                sumSqDev += dev * dev;
            }
            result.StdSignalRatioDev = MathF.Sqrt(sumSqDev / nValid);
        }
    }
}