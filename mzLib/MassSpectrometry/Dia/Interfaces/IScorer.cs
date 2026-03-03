// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Abstraction for scoring similarity between observed and expected fragment intensities.
    /// 
    /// Decoupled from extraction so that:
    ///   - Different scoring methods can be plugged in (dot product, spectral angle, etc.)
    ///   - CPU and GPU implementations can coexist behind the same interface
    ///   - Scoring logic can be tested independently of data loading and extraction
    /// </summary>
    public interface IScorer
    {
        /// <summary>
        /// Computes a similarity score between observed and expected intensity vectors.
        /// Both spans must have the same length.
        /// 
        /// Returns a value in [0, 1] where 1 is perfect match (for most implementations).
        /// </summary>
        /// <param name="observed">Fragment intensities extracted from DIA data.</param>
        /// <param name="expected">Predicted or library fragment intensities.</param>
        float Score(ReadOnlySpan<float> observed, ReadOnlySpan<float> expected);
    }
}
