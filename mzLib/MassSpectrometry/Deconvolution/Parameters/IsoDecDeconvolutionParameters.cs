using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class IsoDecDeconvolutionParameters : DeconvolutionParameters
    {
        public override DeconvolutionType DeconvolutionType { get; protected set; } = DeconvolutionType.IsoDecDeconvolution;
        public int PhaseRes { get; protected set; } = 8;
        public int Verbose { get; protected set; } = 0;
        public int PeakWindow { get; protected set; } = 80;
        public float PeakThreshold { get; protected set; } = (float)0.0001;
        public int MinPeaks { get; protected set; } = 3;
        public float Css_Threshold { get; set; } = (float)0.7;
        public float MatchTolerance { get; set; } = (float)5;
        public int MaxShift { get; protected set; } = 3;
        public float[] MzWindow { get; set; } = new float[] { (float)-1.05, (float)2.05 };
        public float[] PlusOneIntWindow { get; protected set; } = new float[] { (float)0.1, (float)0.6 };
        public int KnockdownRounds { get; set; } = 5;
        public float MinScoreDiff { get; set; } = (float)0.1;
        public float MinAreaCovered { get; set; } = (float)0.15;
        public int IsoLength { get; protected set; } = 64;
        public double MassDiffC { get; protected set; } = 1.0033;
        public float AdductMass { get; set; } = (float)1.00727276467;
        public int MinusOneAreasZero { get; set; } = 1;
        public float IsotopeThreshold { get; set; } = (float)0.01;
        public float DataThreshold { get; set; } = (float)0.05;
        public float ZScoreThreshold { get; protected set; } = (float)0.95;
        public bool ReportMulitpleMonoisos { get; set; } = true;

        public IsoDecDeconvolutionParameters(
            int phaseres = 8,
            bool reportmultiplemonoisos = true,
            float css_threshold = (float)0.7, 
            float match_tolerance = (float)5, 
            int maxshift = 3, 
            float[] mzwindow = null, 
            int knockdown_rounds = 5, 
            float min_area_covered = (float)0.20, 
            float relativedatathreshold = (float)0.05)
            : base(1,50,Polarity.Positive)
        {
            this.PhaseRes = phaseres;
            this.ReportMulitpleMonoisos = reportmultiplemonoisos;
            this.Css_Threshold = css_threshold;
            this.MatchTolerance = match_tolerance;
            this.MaxShift = maxshift;
            if (mzwindow != null)
                this.MzWindow = mzwindow;
            else this.MzWindow = new float[] { (float)-1.05, (float)2.05 };
            this.KnockdownRounds = knockdown_rounds;
            this.MinAreaCovered = min_area_covered;
        }
    }
}
