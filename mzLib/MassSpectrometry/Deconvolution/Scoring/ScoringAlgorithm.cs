using System;
using System.Collections.Generic;
using System.Dynamic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


namespace MassSpectrometry.Deconvolution.Scoring
{
    public abstract class ScoringAlgorithm
    {
        public IScoreArgs ScoreArguments { get; }

        public ScoringAlgorithm(IScoreArgs scoreArguments)
        {
            ScoreArguments = scoreArguments;
        }
        public abstract double Score();

        #region normalization
        protected double[] NormalizeSquareRootSpectrumSum(double[] spectrum)
        {
            double sqrtSum = spectrum.Select(y => Math.Sqrt(y)).Sum();

            for (int i = 0; i < spectrum.Length; i++)
            {
                spectrum[i] = Math.Sqrt(spectrum[i]) / sqrtSum;
            }
            return spectrum;
        }

        protected double[] NormalizeMostAbundantPeak(double[] spectrum)
        {
            double max = spectrum.Max();

            for (int i = 0; i < spectrum.Length; i++)
            {
                spectrum[i] = spectrum[i] / max;
            }
            return spectrum;
        }

        protected double[] NormalizeSpectrumSum(double[] spectrum)
        {
            double sum = spectrum.Sum();

            for (int i = 0; i < spectrum.Length; i++)
            {
                spectrum[i] = spectrum[i] / sum;
            }
            return spectrum;
        }
        #endregion
    }

    public interface IScoreArgs
    {
    }

    public class MinimalSpectraArgs : IScoreArgs
    {
        public MinimalSpectrum experimentalSpectrum { get; set; }
        public MinimalSpectrum theoreticalSpectgurm { get; set; }
    }

    public class IsotopicEnvelopeArgs : IScoreArgs
    {
        public IsotopicEnvelope ExperimentalEnvelope { get; set; }
        public IsotopicEnvelope TheoreticalEnvelope { get; set; }
    }
}
