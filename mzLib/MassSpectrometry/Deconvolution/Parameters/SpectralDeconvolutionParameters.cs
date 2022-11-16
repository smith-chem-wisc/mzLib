using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry.Proteomics;
using MassSpectrometry.Proteomics.ProteolyticDigestion;
using MathNet.Numerics.Optimization;
using MzLibUtil;

namespace MassSpectrometry.Deconvolution.Parameters
{

    public class SpectralDeconvolutionParameters : DeconvolutionParameters
    {
        public List<Protein> Proteins { get;  }
        public List<Modification> FixedModifications { get; }
        public List<Modification> VariableModifications { get; }
        public DigestionParams DigestionParams { get; }
        public List<SilacLabel> SilacLabels { get; }
        public DoubleRange ScanRange { get; }
        public bool FindTopDownTruncationProducts { get; }
        public int BinsPerDalton { get; }
        public double FineResolutionForIsotopicDistribution { get; }
        public double MinProbabilityForIsotopicDistribution { get; }
        public double AmbiguityThresholdForIsotopicDistribution { get; }
        public int MaxConsecutiveMissedIsotopicPeaks { get;  }
        private bool FindNonDatabasePeaks { get; } // This should be linked to a method that generates Averagine envelopes


        public SpectralDeconvolutionParameters(int minAssumedChargeState, int maxAssumedChargeState,
            double deconvolutionTolerancePpm, List<Protein> proteins, List<Modification> fixedModifications,
            List<Modification> variableModifications, DigestionParams digestionParams,
            List<SilacLabel> silacLabels, bool findTopDownTruncationProducts, double scanMinimumMz, double scanMaximumMz,
            int binsPerDalton = 10, double fineResolutionForIsotopicDistribution = 0.125, double minProbabilityForIsotopicDistribution = 1e-8,
            double ambiguityThresholdForIsotopicDistribution = 0.9, int maxConsecutiveMissedIsotopicPeaks = 1,
            bool findNonDatabasePeaks = false) :
            base(minAssumedChargeState, maxAssumedChargeState, deconvolutionTolerancePpm)
        {
            Proteins = proteins;
            FixedModifications = fixedModifications;
            VariableModifications = variableModifications;
            DigestionParams = digestionParams;
            SilacLabels = silacLabels;
            FindTopDownTruncationProducts = findTopDownTruncationProducts;
            ScanRange = new DoubleRange(scanMinimumMz, scanMaximumMz);
            BinsPerDalton = binsPerDalton;
            FineResolutionForIsotopicDistribution = fineResolutionForIsotopicDistribution;
            MinProbabilityForIsotopicDistribution = minProbabilityForIsotopicDistribution;
            AmbiguityThresholdForIsotopicDistribution = ambiguityThresholdForIsotopicDistribution;
            FindNonDatabasePeaks = findNonDatabasePeaks;
        }
    }
}
