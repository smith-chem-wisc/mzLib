using System;
using System.Collections.Generic;

namespace MassSpectrometry.Deconvolution.Consensus
{
    /// <summary>
    /// Greedy mass-trace builder, charge-locked.
    ///
    /// For each MS1 scan in input order, every envelope is offered to the
    /// open traces. An envelope joins a trace if (a) charges match, (b) the
    /// trace has been touched within MaxGap+1 scans of the current one, and
    /// (c) the envelope's mass is within tolerance of the trace's anchor.
    /// Best match (closest mass) wins ties. No drift: the anchor is the
    /// first envelope's mass, fixed for the life of the trace.
    ///
    /// Inputs: per-scan envelope lists in scan order, aligned to a parallel
    /// list of <see cref="MsDataScan"/>. Output: every trace (open or closed
    /// at end-of-input) in a flat list. Length-1 singletons are included --
    /// downstream code applies length filters per its own taste.
    /// </summary>
    public static class MassTraceBuilder
    {
        /// <param name="envelopeScorer">
        /// Optional per-envelope quality scorer, called once for each envelope with the scan it
        /// came from. Supplied here rather than computed later because this is the only point in
        /// the pipeline where an envelope and its spectrum are both in hand; recovering the pair
        /// afterwards would mean re-reading the file. When null, scores are recorded as NaN and
        /// every downstream aggregate reports NaN, which callers must treat as "not scored"
        /// rather than as a score of zero.
        /// </param>
        public static List<MassTrace> BuildTraces(
            IReadOnlyList<MsDataScan> ms1Scans,
            IReadOnlyList<IReadOnlyList<IsotopicEnvelope>> perScanEnvelopes,
            double toleranceDa,
            int maxGap,
            Func<IsotopicEnvelope, MsDataScan, double> envelopeScorer = null)
        {
            if (perScanEnvelopes.Count != ms1Scans.Count)
                throw new System.ArgumentException(
                    $"perScanEnvelopes ({perScanEnvelopes.Count}) must align 1:1 with ms1Scans ({ms1Scans.Count}).",
                    nameof(perScanEnvelopes));

            var open = new List<MassTrace>();
            var closed = new List<MassTrace>();
            int nextId = 1;

            for (int scanIdx = 0; scanIdx < ms1Scans.Count; scanIdx++)
            {
                // Retire open traces whose last-touched scan is too old.
                for (int i = open.Count - 1; i >= 0; i--)
                {
                    int gap = scanIdx - open[i].LastScanIndex - 1;
                    if (gap > maxGap)
                    {
                        closed.Add(open[i]);
                        open.RemoveAt(i);
                    }
                }

                foreach (var env in perScanEnvelopes[scanIdx])
                {
                    MassTrace best = null!;
                    double bestDelta = double.MaxValue;
                    foreach (var t in open)
                    {
                        if (t.Charge != env.Charge) continue;
                        // Don't add a second envelope from the same scan to the same trace.
                        if (t.LastScanIndex == scanIdx) continue;
                        double d = System.Math.Abs(env.MonoisotopicMass - t.AnchorMass);
                        if (d <= toleranceDa && d < bestDelta)
                        {
                            best = t;
                            bestDelta = d;
                        }
                    }

                    var entry = (scanIdx,
                                 ms1Scans[scanIdx].OneBasedScanNumber,
                                 ms1Scans[scanIdx].RetentionTime,
                                 env.MonoisotopicMass,
                                 env.TotalIntensity,
                                 envelopeScorer is null
                                     ? double.NaN
                                     : envelopeScorer(env, ms1Scans[scanIdx]));

                    if (best != null)
                    {
                        best.Envelopes.Add(entry);
                    }
                    else
                    {
                        var nt = new MassTrace
                        {
                            Id = nextId++,
                            Charge = env.Charge,
                            AnchorMass = env.MonoisotopicMass,
                        };
                        nt.Envelopes.Add(entry);
                        open.Add(nt);
                    }
                }
            }

            closed.AddRange(open);
            return closed;
        }
    }
}
