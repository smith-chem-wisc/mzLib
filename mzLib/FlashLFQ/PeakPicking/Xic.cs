using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using MathNet.Numerics.Interpolation;
using MzLibUtil;

namespace FlashLFQ.PeakPicking
{
    public class Xic
    {
        public List<IndexedMassSpectralPeak> ImsPeaks;
        public readonly double PeakFindingMass;
        public readonly SpectraFileInfo SpectraFile;
        public CubicSpline Spline;
        public List<Extremum> Extrema { get; private set; }
        public Dictionary<int, DoubleRange> PeakRegions { get; internal set; }
        /// <summary>
        /// If true, this XIC is used as the retention time reference.
        /// All other XICs  are aligned such that they are maximally similar to
        /// this XIC.
        /// </summary>
        public readonly bool ReferenceXic;
        /// <summary>
        /// The number of minutes added to the retention time of each ImsPeak to
        /// align this Xic with the reference XIC
        /// </summary>
        public readonly double RetentionTimeAdjustment;

        public Xic(List<IndexedMassSpectralPeak> imsPeaks, double peakFindingMass, SpectraFileInfo spectraFile, 
            double retentionTimeAdjustment, bool referenceXic = false, CubicSpline spline = null)
        {
            ImsPeaks = imsPeaks;
            PeakFindingMass = peakFindingMass;
            SpectraFile = spectraFile;
            RetentionTimeAdjustment = retentionTimeAdjustment;
            ReferenceXic = referenceXic;
            Spline = spline ?? CubicSpline.InterpolateAkimaSorted(
                imsPeaks.Select(p => p.RetentionTime+retentionTimeAdjustment).ToArray(),
                imsPeaks.Select(p => p.Intensity).ToArray());
            FindExtrema();
        }

        public void FindExtrema()
        {
            if (Extrema.IsNotNullOrEmpty())
            {
                return;
            }

            Extrema = new List<Extremum>();
            double[] firstDerivativeZeros = Spline.StationaryPoints();
            for (int i = 0; i < firstDerivativeZeros.Length; i++)
            {
                ExtremumType type = Spline.Differentiate2(firstDerivativeZeros[i]) > 0 
                    ? ExtremumType.Minimum
                    : ExtremumType.Maximum;
                Extrema.Add(new Extremum(
                    retentionTime: firstDerivativeZeros[i],
                    intensity: (int)Spline.Interpolate(firstDerivativeZeros[i]),
                    type));
            }
        }

    }
}
