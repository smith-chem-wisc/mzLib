using MathNet.Numerics.Interpolation;
using MathNet.Numerics.Providers.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.IntegralTransforms;
using System.Runtime.CompilerServices;
using System.Numerics;
using Omics.Fragmentation;
using Readers.Generated;
using System.Threading;
using Easy.Common.Extensions;
using System.Timers;
using CsvHelper.Configuration.Attributes;
using Chemistry;

namespace FlashLFQ.Alex_project
{
    public class XIC
    {
        public List<IndexedMassSpectralPeak> Ms1Peaks { get; init; }
        public  double PeakFindingMz { get; init; }
        public  SpectraFileInfo SpectraFile { get; init; }
        public  bool Reference { get; init; }
        public MathNet.Numerics.Interpolation.LinearSpline LinearSpline { get; private set;}
        public CubicSpline SmoothedCubicSpline { get; private set; }
        public double RtShift { get; private set; }
        public List<Extremum> Extrema { get; set; }

        public XIC(List<IndexedMassSpectralPeak> peaks, double peakFindingMass, SpectraFileInfo spectraFile, bool Isreference = false)
        {
            Ms1Peaks = peaks;
            PeakFindingMz = peakFindingMass;
            SpectraFile = spectraFile;
            Reference = Isreference;
            BulidLinearSpline();
        }

        private List<IndexedMassSpectralPeak> PadPeaks() 

        {
            var paddedPeaks = new List<IndexedMassSpectralPeak>();
            var firstPeak = Ms1Peaks[0];
            var lastPeak = Ms1Peaks[Ms1Peaks.Count - 1];
            double gap = (lastPeak.RetentionTime - firstPeak.RetentionTime) / (Ms1Peaks.Count-1);

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

        internal void BulidLinearSpline()
        {
            double[] x = PadPeaks().Select(p => p.RetentionTime).ToArray();
            double[] y = PadPeaks().Select(p => p.Intensity).ToArray();
            this.LinearSpline = LinearSpline.InterpolateSorted(x, y);  // I am not sure what to put in the last parameter
        }


        /// <summary>
        /// calculate the retention time shift between two XICs. Then store the value in the RtShift property
        /// </summary>
        /// <param name="xicToAlign"> The reference XIC</param>
        /// <param name="resolution"> The number of the timepoint</param>
        /// <returns> The retention to shift to align to the reference </returns>
        public double AlignXICs(XIC referenceXIC, int resolution = 100)
        {
            var referSpline = referenceXIC.LinearSpline;
            var toAlignSpline = this.LinearSpline;

            double timegap = (this.Ms1Peaks.Last().RetentionTime - this.Ms1Peaks[0].RetentionTime) / (this.Ms1Peaks.Count - 1);
            double initialTime = this.Ms1Peaks[0].RetentionTime-5.0*timegap; //after the padding, the first peak move ahead 5 timegap
            double FinalTime = this.Ms1Peaks.Last().RetentionTime+5.0*timegap; //after the padding, the last peak move back 5 timegap
            double time = initialTime;

            // create two arrays to store the interpolated values of the two XICs
            Complex[] reference = new Complex[(int)((FinalTime - initialTime)* resolution + 2)];
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

            for (int i = 0; i < product.Length / 2 ; i++) //swap the first half and the second half of the array, ex the timeLine(0,2pi) -> (-pi,pi)
            {
               Complex temp = product[i];
               product[i]= product[product.Length / 2+i]; 
               product[product.Length / 2+i] = temp;               
            }

            double maxMagnitude = product.Max(p => p.Magnitude);
            int indexForTheMaxValue = Array.FindIndex(product, p => p.Magnitude == maxMagnitude);
            double rtShift = -(product.Length/2 - indexForTheMaxValue) * (1.0 / resolution);
            RtShift = rtShift;
            return rtShift;

            // Example
        }

        /// <summary>
        /// Try to smooth the XIC by averaging the intensity of the points
        /// </summary>
        /// <param name="pointsToAverage"> should be odds number. The number of points to average for smoothing the XIC </param>
        /// <exception cref="ArgumentException"></exception>
        public void BuildSmoothedCubicSpline(int pointsToAverage = 3)
        {
            if (pointsToAverage <= 0)
            {
                throw new ArgumentException("pointsToAverage must be greater than 0");
            }

            double[] intensity = Ms1Peaks.Select(p => p.Intensity).ToArray();
            double[] retentionTime = Ms1Peaks.Select(p => p.RetentionTime + RtShift).ToArray();
            double[] smoothedIntensity = new double[intensity.Length];

            for (int i = 0; i < intensity.Length; i++) //smooth the intensity
            {
                if (i < pointsToAverage / 2)
                {
                    smoothedIntensity[i] = 0;
                }

                else if (i >= intensity.Length - pointsToAverage / 2)
                {
                    smoothedIntensity[i] = 0;
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

            this.SmoothedCubicSpline = CubicSpline.InterpolateAkimaSorted(retentionTime, smoothedIntensity);
        }


        /// <summary>
        /// Use the second derivative to find the local maxmun and minimum points
        /// </summary>
        public void FindExtrema()
        {
            Extrema = new List<Extremum>();
            double[] extrePoints = SmoothedCubicSpline.StationaryPoints(); //extrePoint is the retentionTime for the point's first derivative is zero
            foreach (var point in extrePoints) 
            {
                if (SmoothedCubicSpline.Differentiate2(point) < 0) //Local Maximun point
                {
                    var intensity_atPoint = LinearSpline.Interpolate(point - RtShift);
                    Extrema.Add(new Extremum(intensity_atPoint, point, ExtremumType.Maximum));                
                }
                if (SmoothedCubicSpline.Differentiate2(point) > 0) //Local Minmun point
                {
                    var intensity_atPoint = LinearSpline.Interpolate(point - RtShift);
                    Extrema.Add(new Extremum(intensity_atPoint, point, ExtremumType.Minimum));                   
                }
            }

            Extrema.Sort();

        }


    }
}
