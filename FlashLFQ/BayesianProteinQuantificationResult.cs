using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    public class BayesianProteinQuantificationResult
    {
        public readonly string condition1;
        public readonly string condition2;
        public readonly double[] mus;
        public readonly double[] sds;
        public readonly double[] nus;
        public readonly double cutoff;
        public readonly List<double> foldChangeMeasurements;

        public BayesianProteinQuantificationResult(string condition1, string condition2, double[] mus, double[] sds, double[] nus, double cutoff, List<double> fcs)
        {
            this.condition1 = condition1;
            this.condition2 = condition2;
            this.mus = mus;
            this.sds = sds;
            this.nus = nus;
            this.cutoff = cutoff;
            this.foldChangeMeasurements = fcs;
        }

        public double ComputePosteriorErrorProbability()
        {
            return Math.Min(1, ComputeBayesFactor());
        }

        public double ComputeBayesFactor()
        {
            double bayesFactor;

            int numIncreasing = mus.Count(p => p > cutoff);
            int numDecreasing = mus.Count(p => p < -cutoff);

            if (numIncreasing > numDecreasing)
            {
                bayesFactor = (double)mus.Count(p => p < cutoff) / numIncreasing;
            }
            else
            {
                bayesFactor = (double)mus.Count(p => p > -cutoff) / numDecreasing;
            }

            return bayesFactor;
        }
    }
}
