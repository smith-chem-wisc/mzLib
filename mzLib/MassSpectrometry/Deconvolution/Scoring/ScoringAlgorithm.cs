using System;
using System.Collections.Generic;
using System.Dynamic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;


namespace MassSpectrometry.Deconvolution.Scoring
{
    public abstract class ScoringAlgorithm
    {
        public PpmTolerance PpmTolerance { get; }

        public ScoringAlgorithm(PpmTolerance tolerance)
        {
            PpmTolerance = tolerance;
        }
        public abstract double GetScore(IScoreArgs args);
    }

    public interface IScoreArgs
    {
    }

    public class MinimalSpectraArgs : IScoreArgs
    {
        public MinimalSpectrum ExperimentalSpectrum { get; set; }
        public MinimalSpectrum TheoreticalSpectrum { get; set; }

        public MinimalSpectraArgs(MinimalSpectrum experimentalSpectrum, MinimalSpectrum theoreticalSpectrum)
        {
            ExperimentalSpectrum = experimentalSpectrum;
            TheoreticalSpectrum = theoreticalSpectrum;
        }
    }

    public class IsotopicEnvelopeArgs : IScoreArgs
    {
        public IsotopicEnvelope ExperimentalEnvelope { get; set; }
        public IsotopicEnvelope TheoreticalEnvelope { get; set; }
    }
}
