using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using MzLibUtil;

namespace MassSpectrometry
{
    public class IsoDecAlgorithm(DeconvolutionParameters deconParameters)
        : DeconvolutionAlgorithm(deconParameters)
    {
        private static string _phaseModelPath;
        static IsoDecAlgorithm()
        {
            _phaseModelPath = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "Deconvolution", "Algorithms", "IsoDecResources",  "phase_model.bin");
        }

        [StructLayout(LayoutKind.Sequential, Pack =1)]
        public struct MatchedPeak
        {
            public float mz;
            public int z;
            public float monoiso;
            public float peakmass;
            public float avgmass;
            public float area;
            public float peakint;
            [MarshalAs(UnmanagedType.ByValArray, SizeConst = 64)]
            public float[] matchedindsiso;
            [MarshalAs(UnmanagedType.ByValArray, SizeConst = 64)]
            public float[] matchedindsexp;
            [MarshalAs(UnmanagedType.ByValArray, SizeConst = 64)]
            public float[] isomz;
            [MarshalAs(UnmanagedType.ByValArray, SizeConst = 64)]
            public float[] isodist;
            [MarshalAs(UnmanagedType.ByValArray, SizeConst = 64)]
            public float[] isomass;
            [MarshalAs(UnmanagedType.ByValArray, SizeConst = 16)]
            public float[] monoisos;
            int startindex;
            int endindex;
            public float startmz;
            public float endmz;
            public float score;
            public int realisolength;
        }

        public struct IsoSettings
        {
            public int phaseres; // Precision of encoding matrix
            public int verbose; // Verbose output
            public int peakwindow; // Peak Detection Window
            public float peakthresh; // Peak Detection Threshold
            public int minpeaks; // Minimum Peaks for an allowed peak
            public float css_thresh; // Minimum cosine similarity score for isotope distribution
            public float matchtol; // Match Tolerance for peak detection in ppm
            public int maxshift; // Maximum shift allowed for isotope distribution
            [MarshalAs(UnmanagedType.ByValArray, SizeConst =2)]
            public float[] mzwindow; // MZ Window for isotope distribution
            [MarshalAs(UnmanagedType.ByValArray, SizeConst = 2)]
            public float[] plusoneintwindow; // Plus One Intensity range. Will be used for charge state 1
            public int knockdown_rounds; // Number of knockdown rounds
            public float min_score_diff; // Minimum score difference for isotope distribution to allow missed monoisotopic peaks
            public float minareacovered; // Minimum area covered by isotope distribution. Use in or with css_thresh
            public int isolength; // Isotope Distribution Length
            public double mass_diff_c; // Mass difference between isotopes
            public float adductmass; // Adduct Mass
            public int minusoneaszero; // Use set the -1 isotope as 0 to help force better alignments
            public float isotopethreshold; // Threshold for isotope distribution. Will remove relative intensities below this.
            public float datathreshold; // Threshold for data. Will remove relative intensities below this relative to max intensity in each cluster
            public float zscore_threshold; //Ratio above which a secondary charge state prediction will be returned.
        }



        [DllImport("Deconvolution/Algorithms/IsoDecResources/isodeclib.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int process_spectrum(double[] mz, float[] intensity, int len, string modelpath, IntPtr matchedpeaks, IsoSettings settings);

        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            var firstIndex = spectrum.GetClosestPeakIndex(range.Minimum);
            var lastIndex = spectrum.GetClosestPeakIndex(range.Maximum);

            var mzs = spectrum.XArray[firstIndex..lastIndex]
                .Select(p => p)
                .ToArray();
            var intensities = spectrum.YArray[firstIndex..lastIndex]
                .Select(p => (float)p)
                .ToArray();

            IntPtr matchedPeaksPtr = Marshal.AllocHGlobal(Marshal.SizeOf(typeof(MatchedPeak)) * intensities.Length);
            IsoSettings settings = DeconParametersToIsoSettings(DeconvolutionParameters as IsoDecDeconvolutionParameters);
            int result = process_spectrum(mzs, intensities, intensities.Length, _phaseModelPath , matchedPeaksPtr, settings);
            if(result > 0)
            {
                MatchedPeak[] matchedpeaks = new MatchedPeak[result];
                for(int i = 0;i<result;i++)
                {
                    matchedpeaks[i] = Marshal.PtrToStructure<MatchedPeak>(matchedPeaksPtr + i * Marshal.SizeOf(typeof(MatchedPeak)));
                }
                return ConvertToIsotopicEnvelopes(DeconvolutionParameters as IsoDecDeconvolutionParameters, matchedpeaks, spectrum);
            }

            else return new List<IsotopicEnvelope>();
        }

        private List<IsotopicEnvelope> ConvertToIsotopicEnvelopes(IsoDecDeconvolutionParameters parameters, MatchedPeak[] matchedpeaks, MzSpectrum spectrum)
        {
            List<IsotopicEnvelope> result = new List<IsotopicEnvelope>();
            int currentId = 0;
            foreach(MatchedPeak peak in matchedpeaks)
            {
                List<(double,double)> peaks = new List<(double,double)> ();
                for (int i = 0; i < peak.realisolength; i++)
                {

                    List<int> indicesWithinTolerance = spectrum.GetPeakIndicesWithinTolerance(peak.isomz[i], new PpmTolerance(5));
                    double maxIntensity = 0;
                    int maxIndex = -1;
                    foreach (int index in indicesWithinTolerance)
                    {
                        if (spectrum.YArray[index] > maxIntensity) { maxIntensity = spectrum.YArray[index]; maxIndex = index; }
                    }
                    if (maxIndex >= 0)
                    {
                        peaks.Add((spectrum.XArray[maxIndex], spectrum.YArray[maxIndex]));
                    }
                    else
                    {
                        peaks.Add((peak.isomz[i], 0));
                    }

                }

                if(parameters.ReportMulitpleMonoisos)
                {
                    foreach (float monoiso in peak.monoisos)
                    {
                        if (monoiso > 0) { result.Add(new IsotopicEnvelope(currentId, peaks, (double)monoiso, peak.z, peak.peakint, peak.score)); }
                        else break;

                    }
                }
                else { result.Add(new IsotopicEnvelope(currentId, peaks, (double)peak.monoiso, peak.z, peak.peakint, peak.score)); }
                currentId++;
            }
            return result;
        }

        public IsoSettings DeconParametersToIsoSettings(IsoDecDeconvolutionParameters parameters)
        {
            IsoSettings result = new IsoSettings();
            result.phaseres = parameters.PhaseRes;
            result.verbose = parameters.Verbose;
            result.peakwindow = parameters.PeakWindow;
            result.peakthresh = parameters.PeakThreshold;
            result.minpeaks = parameters.MinPeaks;
            result.css_thresh = (float)0.7;
            result.matchtol = 5;
            result.maxshift = 3;
            result.mzwindow = [(float)-1.05, (float)2.05];
            result.plusoneintwindow = [(float)0.1, (float)0.6];
            result.knockdown_rounds = 5;
            result.min_score_diff = (float)0.1;
            result.minareacovered = (float)0.15;
            result.isolength = 64;
            result.mass_diff_c = 1.0033;
            //If polarity is positive, adduct is a proton, if negative, it's the loss of a proton.
            if(parameters.Polarity == Polarity.Positive) { result.adductmass = (float)1.007276467; }
            else { result.adductmass = (float)-1.007276467; }
            result.minusoneaszero = 1;
            result.isotopethreshold = (float)0.01;
            result.datathreshold = (float)0.05;
            result.zscore_threshold = (float)0.95;
            return result;
        }
    }
}
