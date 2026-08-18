using System.Collections.Generic;
using System.IO;
using System.Linq;
using MassSpectrometry;
using MassSpectrometry.Deconvolution.Consensus;

namespace Test.MassSpectrometryTests.Deconvolution
{
    /// <summary>
    /// Computes the generic per-envelope deconvolution score
    /// (<see cref="DeconvolutionScorer.ScoreEnvelope(IsotopicEnvelope, Chemistry.AverageResidue)"/>:
    /// averagine cosine, ppm error, peak completeness, intensity-ratio consistency — logistic-
    /// combined) for every Classic envelope, then aggregates it per consensus feature and writes
    /// a sidecar TSV aligned to the <c>_ms1.feature</c> row order (FeatureIndex == the file's
    /// sequential Id). This gives a uniform, search-independent feature-quality value for a
    /// computed cutoff. Scored at generation time because the consensus pipeline keeps only
    /// mass+intensity per envelope and drops the peaks the scorer needs.
    /// </summary>
    internal static class ConsensusFeatureScoring
    {
        internal static void WriteScores(
            string scoresPath,
            IReadOnlyList<MassFeature> features,
            IReadOnlyList<List<IsotopicEnvelope>> perScanEnvelopes)
        {
            var model = new Averagine();

            // (scanIndex, mass-key) -> generic envelope score, computed once per envelope.
            var scoreByEnv = new Dictionary<(int, long), double>();
            for (int si = 0; si < perScanEnvelopes.Count; si++)
                foreach (var env in perScanEnvelopes[si])
                {
                    var key = (si, (long)System.Math.Round(env.MonoisotopicMass * 1000.0));
                    if (!scoreByEnv.ContainsKey(key))
                        scoreByEnv[key] = DeconvolutionScorer.ScoreEnvelope(env, model);
                }

            using var w = new StreamWriter(scoresPath);
            w.WriteLine("FeatureIndex\tMass\tRTStart\tRTEnd\tChargeMin\tChargeMax\tChargeCount\tSummedIntensity\tNEnvelopes\tMeanScore\tIntWeightedScore\tMaxScore");
            for (int fi = 0; fi < features.Count; fi++)
            {
                var f = features[fi];
                double sum = 0, wsum = 0, wtot = 0, maxs = double.NegativeInfinity;
                int n = 0;
                foreach (var t in f.Traces)
                    foreach (var e in t.Envelopes)
                        if (scoreByEnv.TryGetValue((e.ScanIndex, (long)System.Math.Round(e.OriginalMass * 1000.0)), out double s))
                        {
                            sum += s; wsum += s * e.Intensity; wtot += e.Intensity;
                            if (s > maxs) maxs = s; n++;
                        }
                double mean = n > 0 ? sum / n : 0;
                double iw = wtot > 0 ? wsum / wtot : mean;
                w.WriteLine($"{fi}\t{f.ConsensusMass:F5}\t{f.RTStart:F4}\t{f.RTEnd:F4}\t{f.Charges.Min()}\t{f.Charges.Max()}\t{f.ChargeCount}\t{f.SummedIntensity:E4}\t{n}\t{mean:F5}\t{iw:F5}\t{(n > 0 ? maxs : 0):F5}");
            }
        }
    }
}
