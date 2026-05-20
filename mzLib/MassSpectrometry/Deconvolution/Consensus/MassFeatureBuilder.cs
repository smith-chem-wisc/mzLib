using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry.Deconvolution.Consensus
{
    /// <summary>
    /// Stitches per-charge <see cref="CorrectedTrace"/> entries back into
    /// cross-charge <see cref="MassFeature"/> objects.
    ///
    /// Algorithm:
    ///   1. Sort the corrected traces by ConsensusMass.
    ///   2. Sweep with a ppm-wide window; for each pair of traces inside
    ///      the window, if (a) their charges differ and (b) their RT
    ///      ranges overlap, union-find them.
    ///   3. Each connected component becomes a <see cref="MassFeature"/>.
    ///
    /// Same-charge co-mass traces are intentionally NOT merged here: if
    /// two traces share a charge and mass, the trace builder would already
    /// have merged them when their scan ranges overlapped. If they didn't
    /// overlap (co-eluting twins or a Phase-2 split), they remain separate
    /// features.
    /// </summary>
    public static class MassFeatureBuilder
    {
        public static List<MassFeature> BuildFeatures(
            IReadOnlyList<CorrectedTrace> corrected,
            double massPpm)
        {
            // Sort traces by consensus mass to enable a single sliding-window pass.
            var sorted = corrected.OrderBy(t => t.ConsensusMass).ToArray();
            int n = sorted.Length;

            // Union-find over trace indices.
            int[] parent = new int[n];
            for (int i = 0; i < n; i++) parent[i] = i;
            int Find(int x)
            {
                while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
                return x;
            }
            void Union(int x, int y)
            {
                int rx = Find(x), ry = Find(y);
                if (rx != ry) parent[rx] = ry;
            }

            for (int i = 0; i < n; i++)
            {
                double anchorMass = sorted[i].ConsensusMass;
                double window = anchorMass * massPpm * 1e-6;
                for (int j = i + 1; j < n; j++)
                {
                    double dm = sorted[j].ConsensusMass - anchorMass;
                    if (dm > window) break;
                    if (sorted[i].Charge == sorted[j].Charge) continue;
                    if (!RtOverlap(sorted[i], sorted[j])) continue;
                    Union(i, j);
                }
            }

            // Collect connected components.
            var byRoot = new Dictionary<int, List<CorrectedTrace>>();
            for (int i = 0; i < n; i++)
            {
                int r = Find(i);
                if (!byRoot.TryGetValue(r, out var bucket))
                {
                    bucket = new List<CorrectedTrace>();
                    byRoot[r] = bucket;
                }
                bucket.Add(sorted[i]);
            }

            // Finalise each component, then assign IDs in a deterministic order. Dictionary
            // enumeration is not order-stable, so without this sort the feature IDs -- and the
            // row order of the written _ms1.feature file -- would vary run-to-run on identical
            // input. Order by (consensus mass, then min charge, then RT start).
            var features = new List<MassFeature>(byRoot.Count);
            foreach (var (_, traces) in byRoot)
            {
                var f = new MassFeature { Traces = traces };
                f.Finalise();
                features.Add(f);
            }
            features.Sort((a, b) =>
            {
                int c = a.ConsensusMass.CompareTo(b.ConsensusMass);
                if (c != 0) return c;
                c = a.Charges.Min().CompareTo(b.Charges.Min());
                if (c != 0) return c;
                return a.RTStart.CompareTo(b.RTStart);
            });
            for (int i = 0; i < features.Count; i++)
                features[i].Id = i + 1;
            return features;
        }

        private static bool RtOverlap(CorrectedTrace a, CorrectedTrace b)
            => a.FirstRT <= b.LastRT && b.FirstRT <= a.LastRT;
    }
}
