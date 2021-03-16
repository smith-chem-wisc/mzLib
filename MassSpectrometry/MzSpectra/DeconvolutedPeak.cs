using System;
using System.Collections.Generic;
using System.Text;

namespace MassSpectrometry
{
    public class DeconvolutedPeak
    {
        public double ExperimentalMz;
        public double TheoreticalMz;
        public int Charge;
        public double ExperimentalIntensity;
        public double TheoreticalIntensity;
        public double TheoreticalNormalizedAbundance;
        public int IsotopeNumber;

        public DeconvolutedPeak(double experimentalMz, double theorMz, int z, double expIntensity, double theorIntensity, int isotopeNumber, double theoreticalNormalizedAbundance)
        {
            this.ExperimentalMz = experimentalMz;
            this.TheoreticalMz = theorMz;
            this.Charge = z;
            this.ExperimentalIntensity = expIntensity;
            this.TheoreticalIntensity = theorIntensity;
            this.IsotopeNumber = isotopeNumber;
            this.TheoreticalNormalizedAbundance = theoreticalNormalizedAbundance;
        }

        public override string ToString()
        {
            return TheoreticalIntensity.ToString("F1") + " : " + ExperimentalIntensity.ToString("F1");
        }
    }
}
