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
