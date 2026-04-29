using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using MzLibUtil;

namespace MassSpectrometry
{
    /// <summary>
    /// Performs deconvolution on a single spectrum or region of spectrum using the Isodec algorithm
    /// <remarks>
    /// Isodec only needs to region of interest and does not use surrounding charge states as references.
    /// Isodec can report multiple monoisotopic masses for a single peak if enabled by ReportMultipleMonoisos parameter
    ///     In this case, the resulting isotopic envelopes will have the same precursor ID.
    /// </remarks>
    /// </summary>
    internal class IsoDecAlgorithm : DeconvolutionAlgorithm
    {
        internal IsoDecAlgorithm(DeconvolutionParameters deconParameters) : base(deconParameters)
        {

        }

        /// <summary>
        /// Struct passed by pointer in memory to the Isodec.dll
        /// </summary>
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

        /// <summary>
        /// Calls the Isodec.dll to perform deconvolution on the given spectrum
        /// The Isodec.dll requires three other dll's as dependencies: isogenmass.dll, libmmd.dll, scml_dispmd.dll
        /// </summary>
        /// <param name="cmz"></param>
        /// <param name="cintensity"></param>
        /// <param name="c"></param>
        /// <param name="fname"></param>
        /// <param name="matchedpeaks"></param>
        /// <param name="settings"></param>
        /// <returns></returns>

        [DllImport("isodeclib.dll", EntryPoint = "process_spectrum", CallingConvention = CallingConvention.Cdecl)]
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

            var mpArray = new byte[intensities.Length * Marshal.SizeOf(typeof(MatchedPeak))];
            GCHandle handle = GCHandle.Alloc(mpArray, GCHandleType.Pinned);
            try
            {
                IntPtr matchedPeaksPtr = (IntPtr)handle.AddrOfPinnedObject();
                IsoDecDeconvolutionParameters.IsoSettings settings = deconParams.ToIsoSettings();
                int result = process_spectrum(mzs, intensities, intensities.Length, null, matchedPeaksPtr, settings);
                if (result <= 0)
                    return Enumerable.Empty<IsotopicEnvelope>();

                // Handle results
                MatchedPeak[] matchedpeaks = new MatchedPeak[result];
                for (int i = 0; i < result; i++)
                {
                    matchedpeaks[i] = Marshal.PtrToStructure<MatchedPeak>(matchedPeaksPtr + i * Marshal.SizeOf(typeof(MatchedPeak)));
                }

                return ConvertToIsotopicEnvelopes(deconParams, matchedpeaks, spectrum);
            }
            finally
            {
                handle.Free();
            }
        }

        /// <summary>
        /// Converts the isodec output (MatchedPeak) to IsotopicEnvelope for return
        /// </summary>
        /// <param name="parameters"></param>
        /// <param name="matchedpeaks"></param>
        /// <param name="spectrum"></param>
        /// <returns></returns>
        private List<IsotopicEnvelope> ConvertToIsotopicEnvelopes(IsoDecDeconvolutionParameters parameters, MatchedPeak[] matchedpeaks, MzSpectrum spectrum)
        {
            List<IsotopicEnvelope> result = new List<IsotopicEnvelope>();
            int currentId = 0;
            var tolerance = new PpmTolerance(5);
            foreach(MatchedPeak peak in matchedpeaks)
            {
                List<(double,double)> peaks = new List<(double,double)> ();
                for (int i = 0; i < peak.realisolength; i++)
                {

                    List<int> indicesWithinTolerance = spectrum.GetPeakIndicesWithinTolerance(peak.isomz[i], tolerance);
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
