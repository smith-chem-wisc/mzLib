using MassSpectrometry;
using MzLibUtil;
using Proteomics.ProteolyticDigestion;
using Proteomics.RetentionTimePrediction;
using System.Collections.Generic;

namespace FlashLFQ
{
    public class FlashLfqParameters
    {
        public FlashLfqParameters()
        {
            Normalize = false;
            PpmTolerance = 10.0;
            IsotopePpmTolerance = 5.0;
            Integrate = false;
            NumIsotopesRequired = 2;
            IdSpecificChargeState = false;
            UseSharedPeptidesForProteinQuant = false;
            QuantifyAmbiguousPeptides = false;
            Silent = false;
            MaxThreads = -1;

            // IsoTracker settings
            IsoTracker = false;

            // MBR settings
            MatchBetweenRuns = false;
            MbrPpmTolerance = 10.0;
            MaxMbrRtWindow = 1.0;
            RequireMsmsIdInCondition = false;
            MbrQValueThreshold = 0.05;
            DonorCriterion = DonorCriterion.Score;
            DonorQValueThreshold = 0.01;

            // settings for the Bayesian protein quantification engine
            BayesianProteinQuant = false;
            ProteinQuantBaseCondition = null;
            ProteinQuantFoldChangeCutoff = 0.1;
            McmcSteps = 3000;
            McmcBurninSteps = 1000;
            PairedSamples = false;
            RandomSeed = null;
        }

        public bool Normalize { get; set; }
        public double PpmTolerance { get; set; }
        public double IsotopePpmTolerance { get; set; }
        public bool Integrate { get; set; }
        public int NumIsotopesRequired { get; set; }
        public bool IdSpecificChargeState { get; set; }
        public bool UseSharedPeptidesForProteinQuant { get; set; }
        public bool QuantifyAmbiguousPeptides { get; set; }
        public bool Silent { get; set; }
        public int MaxThreads { get; set; }

        //IsoTracker settings
        public bool IsoTracker { get; set; } //Searching parameter for the FlashLFQ engine

        // MBR settings
        public bool MatchBetweenRuns { get; set; }
        public double MbrPpmTolerance { get; set; }
        public double MaxMbrRtWindow { get; set; }
        public bool RequireMsmsIdInCondition { get; set; }
        public double MbrQValueThreshold { get; set; }
        /// <summary>
        /// Specifies how the donor peak for MBR is selected. 
        /// 'Score' selects the donor peak associated with the highest scoring PSM
        /// 'Intensity' selects the donor peak with the max intensity
        /// 'Neighbors' selects the donor peak with the most neighboring peaks
        /// </summary>
        public DonorCriterion DonorCriterion { get; set; }
        public double DonorQValueThreshold { get; set; }

        
        // settings for the Bayesian protein quantification engine
        public bool BayesianProteinQuant { get; set; }
        public string ProteinQuantBaseCondition { get; set; }
        public double ProteinQuantFoldChangeCutoff { get; set; }
        public int McmcSteps { get; set; }
        public int McmcBurninSteps { get; set; }
        public bool PairedSamples { get; set; }
        public int? RandomSeed { get; set; }
    }

    public enum DonorCriterion
    {
        Score,
        Intensity,
        Neighbors
    }
}
