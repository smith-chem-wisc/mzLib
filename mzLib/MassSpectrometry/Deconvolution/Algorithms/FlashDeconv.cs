using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static MassSpectrometry.IsoDecAlgorithm;

namespace MassSpectrometry
{
    internal class FlashDeconv : DeconvolutionAlgorithm
    {
        internal FlashDeconv(DeconvolutionParameters deconParameters) : base(deconParameters)
        {

        }
        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            var t = LogTransformSpectrum(spectrum);

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
            foreach (MatchedPeak peak in matchedpeaks)
            {
                List<(double, double)> peaks = new List<(double, double)>();
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
                if (parameters.Polarity == Polarity.Negative) { charge = -peak.z; }
                if (parameters.ReportMulitpleMonoisos)
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
    
        private MzSpectrum LogTransformSpectrum(MzSpectrum spectrum)
        {
            return new MzSpectrum(spectrum.XArray.Select(x => Math.Log(x)).ToArray(), spectrum.YArray, true);
        }
    }
}
