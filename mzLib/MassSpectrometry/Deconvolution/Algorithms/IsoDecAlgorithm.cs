using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using MathNet.Numerics.Statistics;

namespace MassSpectrometry
{
    public class IsoDecAlgorithm(DeconvolutionParameters deconParameters)
        : DeconvolutionAlgorithm(deconParameters)
    {
        public static string _isoDecDllPath;
        private static string _libmmDllPath;
        private static string _phaseModelPath;
        private static string _scmlDispMdDllPath;
        static IsoDecAlgorithm()
        {
            
            _isoDecDllPath = Path.Combine(new[]
            {
                AppDomain.CurrentDomain.BaseDirectory, "MassSpectrometry", "Deconvolution", "Resources", "isodeclib.dll"
            });

            _libmmDllPath = Path.Combine(new[]
            {
                AppDomain.CurrentDomain.BaseDirectory, "MassSpectrometry", "Deconvolution", "Resources", "libmm.dll"
            });

            _phaseModelPath = Path.Combine(new[]
            {
                AppDomain.CurrentDomain.BaseDirectory, "MassSpectrometry", "Deconvolution", "Resources",
                "phase_model.bin"
            });

            _scmlDispMdDllPath = Path.Combine(new[]
            {
                AppDomain.CurrentDomain.BaseDirectory, "MassSpectrometry", "Deconvolution", "Resources",
                "scml_disp_md.dll"
            });
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
            public int[] matchedindsiso;
            [MarshalAs(UnmanagedType.ByValArray, SizeConst = 64)]
            public int[] matchedindsexp;
            [MarshalAs(UnmanagedType.ByValArray, SizeConst = 128)]
            public float[] isomz;
            [MarshalAs(UnmanagedType.ByValArray, SizeConst = 128)]
            public float[] isodist;
            [MarshalAs(UnmanagedType.ByValArray, SizeConst = 128)]
            public float[] isomass;
            public int startindex;
            public int endindex;
        }


        
        [DllImport(@"C:\Python\UniDec3\unidec\IsoDec\src\isodec\x64\Release\isodeclib.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern int process_spectrum(float[] mz, float[] intensity, int len, string modelpath, IntPtr matchedpeaks);


        public override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            var firstIndex = spectrum.GetClosestPeakIndex(range.Minimum);
            var lastIndex = spectrum.GetClosestPeakIndex(range.Maximum);

            var mzs = spectrum.XArray[firstIndex..lastIndex]
                .Select(p => (float)p)
                .ToArray();
            var intensities = spectrum.YArray[firstIndex..lastIndex]
                .Select(p => (float)p)
                .ToArray();

            IntPtr matchedPeaksPtr = Marshal.AllocHGlobal(Marshal.SizeOf(typeof(MatchedPeak)) * intensities.Length);
            int result = process_spectrum(mzs, intensities, intensities.Length, @"C:\Python\UniDec3\unidec\IsoDec\phase_model.bin", matchedPeaksPtr);
            if(result > 0)
            {
                MatchedPeak[] matchedpeaks = new MatchedPeak[result];
                for(int i = 0;i<result;i++)
                {
                    matchedpeaks[i] = Marshal.PtrToStructure<MatchedPeak>(matchedPeaksPtr + i * Marshal.SizeOf(typeof(MatchedPeak)));
                }
                return ConvertToIsotopicEnvelopes(matchedpeaks, mzs, intensities);
            }

            else return new List<IsotopicEnvelope>();
        }

        private List<IsotopicEnvelope> ConvertToIsotopicEnvelopes(MatchedPeak[] matchedpeaks, float[] mzs, float[] intensities)
        {
            List<IsotopicEnvelope> result = new List<IsotopicEnvelope>();
            foreach(MatchedPeak peak in matchedpeaks)
            {
                List<(double,double)> peaks = new List<(double,double)> ();
                List<double> listofratios = new List<double>();
                for(int i = 0;i<peak.matchedindsexp.Length;i++)
                {
                    if (i > 0 && peak.matchedindsexp[i] == 0)
                        break;
                    else
                    {
                        listofratios.Add(peak.isodist[i] / intensities[peak.startindex + i]);
                        peaks.Add((mzs[peak.startindex + i], intensities[peak.startindex + i]));
                    }
                }
                result.Add(new IsotopicEnvelope(peaks, peak.monoiso, peak.z, peak.peakint, Statistics.StandardDeviation(listofratios), 0));
            }
            return result;
        }
    }
}
