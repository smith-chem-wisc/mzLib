using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SimulatedData
{
    public class LowFrequencyNoiseParameters
    {
        public int PeakNumberLimitLow { get; }
        public int PeakNumberLimitHigh { get; }
        public double PeakLocationLimitLow { get; }
        public double PeakLocationLimitHigh { get; }
        public double PeakWidthLimitLow { get; }
        public double PeakWidthLimitHigh { get; }
        public double PeakIntensityLimitLow { get; }
        public double PeakIntensityLimitHigh { get; }
        public (double Min, double Max)? ExcludedZone { get; }
        public PeakShapeOptions PeakShapeOptions { get; }

        public LowFrequencyNoiseParameters(int peakNumberLimitLow, int peakNumberLimitHigh,
            double peakLocationLimitLow, double peakLocationLimitHigh,
            double peakWidthLimitLow, double peakWidthLimitHigh,
            double peakIntensityLimitLow, double peakIntensityLimitHigh,
            (double Min, double Max)? excludedZone = null,
            PeakShapeOptions peakShapeOptions = PeakShapeOptions.Gaussian)
        {
            PeakNumberLimitLow = peakNumberLimitLow;
            PeakNumberLimitHigh = peakNumberLimitHigh;
            PeakLocationLimitLow = peakLocationLimitLow;
            PeakLocationLimitHigh = peakLocationLimitHigh;
            PeakWidthLimitLow = peakWidthLimitLow;
            PeakWidthLimitHigh = peakWidthLimitHigh;
            PeakIntensityLimitLow = peakIntensityLimitLow;
            PeakIntensityLimitHigh = peakIntensityLimitHigh;
            ExcludedZone = excludedZone;
            PeakShapeOptions = peakShapeOptions;
        }
    }
}
