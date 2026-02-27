// Copyright 2026 mzLib Contributors
// Licensed under the GNU Lesser General Public License v3.0

using System;

namespace MassSpectrometry.Dia
{
    /// <summary>
    /// Describes a single fragment ion to extract from DIA data.
    /// 
    /// This is a value type (struct) used in batch extraction. Each query specifies
    /// a target m/z, a PPM tolerance for matching, a retention time window to restrict
    /// the search, and the isolation window ID that contains this fragment.
    /// 
    /// Keeping this as a struct ensures that batches of queries (FragmentQuery[]) are
    /// stored contiguously in memory with no per-element heap overhead.
    /// </summary>
    public readonly struct FragmentQuery
    {
        /// <summary>Target m/z of the fragment ion to extract.</summary>
        public readonly float TargetMz;

        /// <summary>PPM tolerance for m/z matching.</summary>
        public readonly float TolerancePpm;

        /// <summary>Minimum retention time (minutes) for extraction window.</summary>
        public readonly float RtMin;

        /// <summary>Maximum retention time (minutes) for extraction window.</summary>
        public readonly float RtMax;

        /// <summary>
        /// The DIA isolation window ID that this fragment should be extracted from.
        /// This is the integer window ID assigned by DiaScanIndexBuilder.
        /// </summary>
        public readonly int WindowId;

        /// <summary>
        /// Caller-defined identifier to track which precursor/fragment this query belongs to.
        /// Not used by the extraction engine itself — purely for result correlation.
        /// </summary>
        public readonly int QueryId;

        public FragmentQuery(float targetMz, float tolerancePpm, float rtMin, float rtMax, int windowId, int queryId = 0)
        {
            TargetMz = targetMz;
            TolerancePpm = tolerancePpm;
            RtMin = rtMin;
            RtMax = rtMax;
            WindowId = windowId;
            QueryId = queryId;
        }
    }

    /// <summary>
    /// Result of extracting a single fragment ion across scans in an RT window.
    /// 
    /// Contains the extracted ion chromatogram (XIC) data as offsets into a shared
    /// result buffer. This avoids per-result heap allocation — the actual RT and 
    /// intensity values live in contiguous arrays managed by the extraction engine.
    /// </summary>
    public readonly struct FragmentResult
    {
        /// <summary>The QueryId from the corresponding FragmentQuery, for result correlation.</summary>
        public readonly int QueryId;

        /// <summary>Number of XIC data points extracted for this fragment.</summary>
        public readonly int DataPointCount;

        /// <summary>
        /// Offset into the shared RT result buffer where this fragment's XIC RT values begin.
        /// </summary>
        public readonly int RtBufferOffset;

        /// <summary>
        /// Offset into the shared intensity result buffer where this fragment's XIC intensity values begin.
        /// </summary>
        public readonly int IntensityBufferOffset;

        /// <summary>
        /// Sum of all extracted intensities across the RT window.
        /// Provided as a convenience for quick scoring without re-iterating the XIC.
        /// </summary>
        public readonly float TotalIntensity;

        public FragmentResult(int queryId, int dataPointCount, int rtBufferOffset, int intensityBufferOffset, float totalIntensity)
        {
            QueryId = queryId;
            DataPointCount = dataPointCount;
            RtBufferOffset = rtBufferOffset;
            IntensityBufferOffset = intensityBufferOffset;
            TotalIntensity = totalIntensity;
        }
    }
}
