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
        private static readonly string _phaseModelPath;
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
            public float score;
            public int realisolength;
        }

        [DllImport("Deconvolution/Algorithms/IsoDecResources/isodeclib.dll", EntryPoint = "process_spectrum", CallingConvention = CallingConvention.Cdecl)]
        protected static extern int process_spectrum(double[] cmz, float[] cintensity, int c, string fname, IntPtr matchedpeaks, IsoDecDeconvolutionParameters.IsoSettings settings);

        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            var deconParams = DeconvolutionParameters as IsoDecDeconvolutionParameters ?? throw new MzLibException("Deconvolution params and algorithm do not match");

            var firstIndex = spectrum.GetClosestPeakIndex(range.Minimum);
            var lastIndex = spectrum.GetClosestPeakIndex(range.Maximum);

            var mzs = spectrum.XArray[firstIndex..lastIndex]
                .Select(p => p)
                .ToArray();
            var intensities = spectrum.YArray[firstIndex..lastIndex]
                .Select(p => (float)p)
                .ToArray();

            IntPtr matchedPeaksPtr = IntPtr.Zero;
            try
            {
                matchedPeaksPtr = Marshal.AllocHGlobal(Marshal.SizeOf(typeof(MatchedPeak)) * intensities.Length);
                IsoDecDeconvolutionParameters.IsoSettings settings = deconParams.ToIsoSettings();
                int result = process_spectrum(mzs, intensities, intensities.Length, _phaseModelPath, matchedPeaksPtr, settings);
                    if (result <= 0)
                    return Enumerable.Empty<IsotopicEnvelope>();

                // Handle results
                MatchedPeak[] matchedpeaks = new MatchedPeak[result];
                for (int i = 0; i < result; i++)
                {
                    matchedpeaks[i] = Marshal.PtrToStructure<MatchedPeak>(matchedPeaksPtr + i * Marshal.SizeOf(typeof(MatchedPeak)));
                }
                var envelopes = ConvertToIsotopicEnvelopes(deconParams, matchedpeaks, spectrum).ToList();
                return envelopes;
            }
            finally
            {
                if (matchedPeaksPtr != IntPtr.Zero)
                {
                    Marshal.FreeHGlobal(matchedPeaksPtr);
                }
            }
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
                int charge = peak.z;
                if(parameters.Polarity == Polarity.Negative) { charge = -peak.z; }
                if(parameters.ReportMulitpleMonoisos)
                {
                    foreach (float monoiso in peak.monoisos)
                    {
                        if (monoiso > 0) { result.Add(new IsotopicEnvelope(currentId, peaks, (double)monoiso, charge, peak.peakint, peak.score)); }
                    }
                }
                else { result.Add(new IsotopicEnvelope(currentId, peaks, (double)peak.monoiso, charge, peak.peakint, peak.score)); }
                currentId++;
            }
            return result;
        }
    }
}
