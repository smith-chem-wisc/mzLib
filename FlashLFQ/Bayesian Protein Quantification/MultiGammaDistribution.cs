using MathNet.Numerics.Distributions;
using System;
using System.Collections.Generic;

namespace FlashLFQ
{
    public class MultiGammaDistribution : IContinuousDistribution
    {
        private List<Gamma> Gammas;
        private List<int> Weights;

        public MultiGammaDistribution(List<Gamma> gammas, List<int> weights)
        {
            Gammas = gammas;
            Weights = weights;
        }

        public double Mode => throw new NotImplementedException();

        public double Minimum => throw new NotImplementedException();

        public double Maximum => throw new NotImplementedException();

        public double Mean => throw new NotImplementedException();

        public double Variance => throw new NotImplementedException();

        public double StdDev => throw new NotImplementedException();

        public double Entropy => throw new NotImplementedException();

        public double Skewness => throw new NotImplementedException();

        public double Median => throw new NotImplementedException();

        public Random RandomSource { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

        public double CumulativeDistribution(double x)
        {
            throw new NotImplementedException();
        }

        public double DensityLn(double x)
        {
            throw new NotImplementedException();
        }

        public double Density(double x)
        {
            double pdf = 1;
            int totalWeight = 0;

            for (int i = 0; i < Gammas.Count; i++)
            {
                var gamma = Gammas[i];
                int weight = Weights[i];
                totalWeight += weight;
                
                pdf *= Math.Pow(gamma.Density(x), weight);
            }

            return Math.Pow(pdf, 1.0 / totalWeight);
        }

        public double Sample()
        {
            throw new NotImplementedException();
        }

        public void Samples(double[] values)
        {
            throw new NotImplementedException();
        }

        public IEnumerable<double> Samples()
        {
            throw new NotImplementedException();
        }
    }
}
