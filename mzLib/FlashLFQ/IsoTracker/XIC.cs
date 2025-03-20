using MathNet.Numerics.Interpolation;
using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.IntegralTransforms;
using System.Numerics;

namespace FlashLFQ.IsoTracker
{
    /// <summary>
    /// An XIC should contain a list of every indexed Mass Spectral Peak (imsPeak) that corresponds to the MS1 signal
    /// produced by a given m/z, in a given time range, in one LC-MS Run
    /// </summary>
    public class XIC
    {
        /// <summary>
        /// A list of smoothed intensity values.
        /// </summary>
        public List<double> SmoothedIntensity { get; set; }
        /// <summary>
        /// A list of smoothed retention time values.
        /// </summary>
        public List<double> SmoothedRetentionTime { get; set; }
        /// <summary>
        /// A list of imsPeak objects.
        /// </summary>
        public List<IIndexedMzPeak> Ms1Peaks { get; init; }
        /// <summary>
        /// The peak-finding mz used to detect the imsPeaks. Every peak should have this m/z+- tolerance.
        /// </summary>
        public double PeakFindingMz { get; init; }
        /// <summary>
        /// The spectra file the XIC came from.
        /// </summary>
        public SpectraFileInfo SpectraFile { get; init; }
        /// <summary>
        /// For a collection of XICs, one should be the reference for time alignment. We picked the XIC with the most MS2ID.
        /// </summary>
        public bool Reference { get; set; }
        /// <summary>
        /// A spline will be used for data interpolation.
        /// </summary>
        public LinearSpline LinearSpline { get; private set; }
        /// <summary>
        /// A cubic spline for smoothed data interpolation.
        /// </summary>
        public CubicSpline SmoothedCubicSpline { get; private set; }
        /// <summary>
        /// The time shift from the reference XIC.
        /// </summary>
        public double RtShift { get; private set; }
        /// <summary>
        /// A list of extrema points.
        /// </summary>
        public List<Extremum> Extrema { get; set; }
        public List<Identification> Ids;

        public XIC(List<IIndexedMzPeak> peaks, double peakFindingMass, SpectraFileInfo spectraFile, bool Isreference = false, List<Identification> ids = null, int smoothDegree = 5)
        {

            Ms1Peaks = peaks;
            PeakFindingMz = peakFindingMass;
            SpectraFile = spectraFile;
            Reference = Isreference;
            Ids = ids;
            BuildLinearSpline();
            if (peaks.Count >= 5)
            {
                BuildSmoothedCubicSpline(smoothDegree);
            }
            else
            {
                this.SmoothedCubicSpline = null;
            }

        }

        /// <summary>
        /// Pad the XIC with 5 peaks before the first peak and 5 peaks after the last peak. The intensity of the padded peaks is 0.
        /// </summary>
        /// <returns></returns>
        private List<IIndexedMzPeak> PadPeaks()

        {
            var paddedPeaks = new List<IIndexedMzPeak>();
            var firstPeak = Ms1Peaks[0];
            var lastPeak = Ms1Peaks[Ms1Peaks.Count - 1];
            double gap = (lastPeak.RetentionTime - firstPeak.RetentionTime) / (Ms1Peaks.Count - 1);

            // because we hope to have an odd number of peaks, we have to add the odd number padded peaks



            for (int i = 5; i > 0; i--) //add 4 peaks before the first peak
            {
                paddedPeaks.Add(new IndexedMassSpectralPeak(0, 0, 0, firstPeak.RetentionTime - gap * i));    // not sure about the m/z and index    
            }

            for (int i = 0; i < Ms1Peaks.Count; i++)
            {
                paddedPeaks.Add(Ms1Peaks[i]);
            }

            for (int i = 1; i < 6; i++) //add 5 peaks after the last peak
            {
                paddedPeaks.Add(new IndexedMassSpectralPeak(0, 0, 0, lastPeak.RetentionTime + gap * i));
            }

            return paddedPeaks;
        }

        internal void BuildLinearSpline()
        {
            double[] x = PadPeaks().Select(p => p.RetentionTime).ToArray();
            double[] y = PadPeaks().Select(p => p.Intensity).ToArray();
            this.LinearSpline = LinearSpline.InterpolateSorted(x, y);  // I am not sure what to put in the last parameter
        }

        /// <summary>
        /// calculate the retention time shift among the reference. Then store the value in the RtShift property.
        /// </summary>
        /// <param name="xicToAlign"> The reference XIC</param>
        /// <param name="resolution"> The number of the timePoint for X-correlation </param>
        /// <returns> The retention to shift to align to the reference </returns>
        public double AlignXICs(XIC referenceXIC, int resolution = 1000)
        {
            var referSpline = referenceXIC.LinearSpline;
            var toAlignSpline = this.LinearSpline;

            double timegap = (Ms1Peaks.Last().RetentionTime - Ms1Peaks[0].RetentionTime) / (Ms1Peaks.Count - 1);
            double initialTime = Ms1Peaks[0].RetentionTime - 5.0 * timegap; //after the padding, the first peak move ahead 5 timegap
            double FinalTime = Ms1Peaks.Last().RetentionTime + 5.0 * timegap; //after the padding, the last peak move back 5 timegap
            double time = initialTime;

            // Create two intensity arrays of the reference and target XICs by interpolating the time point. 
            // The number of timePoints depend on the resolution.
            // If the resolution is 1000, the timePoint is 1/1000 = 0.001s.
            Complex[] reference = new Complex[(int)((FinalTime - initialTime) * resolution + 2)];
            Complex[] toAlign = new Complex[(int)((FinalTime - initialTime) * resolution + 2)];
            int index = 0;
            while (time < FinalTime)
            {
                reference[index] = new Complex(referSpline.Interpolate(time), 0);
                toAlign[index] = new Complex(toAlignSpline.Interpolate(time), 0);
                time += (1.0 / resolution);
                index++;
            }

            Fourier.Forward(reference);
            Fourier.Forward(toAlign);
            Complex[] product = new Complex[reference.Length];
            for (int i = 0; i < reference.Length; i++)
            {
                product[i] = reference[i] * Complex.Conjugate(toAlign[i]); //element-wise multiplication
            }
            Fourier.Inverse(product);

            for (int i = 0; i < product.Length / 2; i++) //swap the first half and the second half of the array, ex the timeLine(0,2pi) -> (-pi,pi)
            {
                Complex temp = product[i];
                product[i] = product[product.Length / 2 + i];
                product[product.Length / 2 + i] = temp;
            }

            double maxMagnitude = product.Max(p => p.Magnitude);
            int indexForTheMaxValue = Array.FindIndex(product, p => p.Magnitude == maxMagnitude);
            double rtShift = -(product.Length / 2 - indexForTheMaxValue) * (1.0 / resolution);
            RtShift = rtShift;
            return rtShift;
        }

        /// <summary>
        /// Try to smooth the XIC by averaging the intensity of the points (weight averaging then sum averaging).
        /// Using the smoothedXIC to generate the cubic spline date for Extremum finding.
        /// </summary>
        /// <param name="pointsToAverage"> Should be odds number. The number of points to average for smoothing the XIC </param>
        /// <exception cref="ArgumentException"></exception>
        public void BuildSmoothedCubicSpline(int smoothDegree)
        {
            double[] retentionTime = Ms1Peaks.Select(p => p.RetentionTime).ToArray();
            double[] intensity = Ms1Peaks.Select(p => p.Intensity).ToArray();

            SmoothedIntensity = IntensitySmoothing_weight(intensity, smoothDegree).ToList();
            SmoothedIntensity = IntensitySmoothing_normal(SmoothedIntensity.ToArray(), smoothDegree).ToList();


            SmoothedRetentionTime = retentionTime.ToList();
            this.SmoothedCubicSpline = CubicSpline.InterpolateAkimaSorted(retentionTime, SmoothedIntensity.ToArray());
        }

        /// <summary>
        /// Smooth the intensity of the XIC by averaging the intensity of the points (weight averaging then sum averaging).
        /// </summary>
        /// <param name="intensity"></param>
        /// <param name="pointsToAverage"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public double[] IntensitySmoothing_weight(double[] intensity, int pointsToAverage)
        {
            if (pointsToAverage <= 0)
            {
                throw new ArgumentException("pointsToAverage must be greater than 0");
            }          
            double[] smoothedIntensity = new double[intensity.Length];

            for (int i = 0; i < intensity.Length; i++) //smooth the intensity
            {
                List<double> weightList = new List<double>();

                if (i < pointsToAverage / 2)
                {
                    smoothedIntensity[i] = intensity[i];
                }

                else if (i >= intensity.Length - pointsToAverage / 2)
                {
                    smoothedIntensity[i] = intensity[i];
                }

                else
                {
                    double sum;
                    double weightSum = 0;

                    for (int j = 0; j < pointsToAverage; j++)
                    {
                        weightList.Add(intensity[i + j - pointsToAverage / 2]);
                    }

                    sum = weightList.Sum();

                    for (int j = 0; j < pointsToAverage; j++)
                    {
                        double weightRatio = weightList[j] / sum;
                        weightSum += weightList[j] * weightRatio;
                    }

                    smoothedIntensity[i] = weightSum;
                }

            }
            return smoothedIntensity;
        }

        /// <summary>
        /// Smooth the intensity of the XIC by averaging the intensity of the points (normal averaging).
        /// </summary>
        /// <param name="intensity"></param>
        /// <param name="pointsToAverage"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public double[] IntensitySmoothing_normal(double[] intensity, int pointsToAverage) 
        {
            if (pointsToAverage <= 0)
            {
                throw new ArgumentException("pointsToAverage must be greater than 0");
            }
            double[] smoothedIntensity = new double[intensity.Length];

            for (int i = 0; i < intensity.Length; i++) //smooth the intensity
            {
                List<double> weightList = new List<double>();

                if (i < pointsToAverage / 2)
                {
                    smoothedIntensity[i] = intensity[i];
                }

                else if (i >= intensity.Length - pointsToAverage / 2)
                {
                    smoothedIntensity[i] = intensity[i];
                }

                else
                {
                    double sum = 0;

                    for (int j = 0; j < pointsToAverage; j++)
                    {
                        sum += intensity[i + j - pointsToAverage / 2];
                    }
                    smoothedIntensity[i] = sum / pointsToAverage;
                }

            }
            return smoothedIntensity;
        }

        /// <summary>
        /// Use the second derivative to find the local maximun and minimum points. Notably, at here, the time of the extrema is the time projected to the reference XIC
        /// For example, if the time of the extrema is 10.5, the time of the extrema in the reference XIC is 10.5 + RtShift
        /// The reason to do this is to make the shift is easier to find the shared extrema
        /// </summary>
        public void FindExtrema()
        {
            Extrema = new List<Extremum>();
            if (SmoothedCubicSpline == null)
            {
                return;
            }

            double[] extrePoints = SmoothedCubicSpline.StationaryPoints().Order().ToArray(); //extrePoint is the retentionTime for the point's first derivative is zero
            foreach (var point in extrePoints) 
            {
                if (SmoothedCubicSpline.Differentiate2(point) < 0) //Local Maximun point
                {
                    var intensity_atPoint = LinearSpline.Interpolate(point);
                    Extrema.Add(new Extremum(intensity_atPoint, point + RtShift, ExtremumType.Maximum)); // Store the            
                }
                if (SmoothedCubicSpline.Differentiate2(point) > 0) //Local Minmun point
                {
                    var intensity_atPoint = LinearSpline.Interpolate(point);
                    Extrema.Add(new Extremum(intensity_atPoint, point + RtShift, ExtremumType.Minimum));                   
                }
            }
            Extrema.Sort();
        }
    }
}
