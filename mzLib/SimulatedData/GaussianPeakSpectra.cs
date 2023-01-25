using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;

namespace SimulatedData
{
    public class GaussianPeakSpectra : SimulatedData
    {
        public double Mean { get; }
        public double Stddev { get; }
        public double IntensityMultiple { get; }
        public GaussianPeakSpectra(double mean, double stddev, double intensityMultiple, 
            int length, double startValue, double spacing) 
            : base(length, startValue, spacing) 
        {
            Mean = mean;
            Stddev = stddev;
            IntensityMultiple = intensityMultiple;
            ApplyElementwise(GaussianFunc, Yarray);
        }
        protected GaussianPeakSpectra()
        {

        }
        protected double GaussianFunc(double d) => IntensityMultiple * Normal.PDF(Mean, Stddev, d);
    }
}
