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
            var logTransformedSpectrum = LogTransformSpectrum(spectrum);
            var acceptibleLogMzDifferences = AllAcceptibleLogMzDifferences();

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
        private List<KeyValuePair<int, double>> AllAcceptibleLogMzDifferences(int lowValue = 1, int highValue = 60)
        {
            return Enumerable.Range(lowValue, highValue)
                    .Select(i => new KeyValuePair<int, double>(i, Math.Log(i)))
                    .ToList();
        }
        private List<(int mzDiffInt,double mzDiff)> Acceptible(int lowValue = 1, int highValue = 60)
        {
            return Enumerable.Range(lowValue, highValue)
            .Select(i => (i, Math.Log(i)))
            .ToList();
        }
        private List<(int, int)> FindIndexedLogMassDifferences(List<(int mzDiffInt, double mzDiff)> acceptibleLogMzDifferences, MzSpectrum spectrum, int selectedLogMzIndex)
        {
            double tolerance = 0.0001;
            List<(int, int)> indexedLogMassDifferences = new ();
            if(selectedLogMzIndex < spectrum.XArray.Length)
            {
                int logMzDiffIndex = 1;
                while (logMzDiffIndex < acceptibleLogMzDifferences.Select(i=>i.mzDiffInt).Max())
                {
                    var i = FindIndexWithinTolerance(spectrum.XArray.ToList(), spectrum.XArray[selectedLogMzIndex] + acceptibleLogMzDifferences[logMzDiffIndex].mzDiff, tolerance);
                    if (i >= 0)
                    {
                        indexedLogMassDifferences.Add((selectedLogMzIndex, i));
                    }
                    logMzDiffIndex++;
                }

            }
            return indexedLogMassDifferences;


        }
        // Finds the index of the first value in a sorted list of doubles matching the target within a given tolerance
        int FindIndexWithinTolerance(List<double> sortedList, double target, double tolerance)
        {
            for (int i = 0; i < sortedList.Count; i++)
            {
                if (Math.Abs(sortedList[i] - target) <= tolerance)
                    return i;
            }
            return -1; // Not found
        }
        // Finds the longest sequence of consecutive integers in a sorted list
        List<int> FindLongestConsecutiveSequence(List<int> sortedList)
        {
            if (sortedList == null || sortedList.Count == 0)
                return new List<int>();

            int maxStart = 0, maxLength = 1;
            int currentStart = 0, currentLength = 1;

            for (int i = 1; i < sortedList.Count; i++)
            {
                if (sortedList[i] == sortedList[i - 1] + 1)
                {
                    currentLength++;
                }
                else
                {
                    if (currentLength > maxLength)
                    {
                        maxLength = currentLength;
                        maxStart = currentStart;
                    }
                    currentStart = i;
                    currentLength = 1;
                }
            }

            // Check at the end
            if (currentLength > maxLength)
            {
                maxLength = currentLength;
                maxStart = currentStart;
            }

            return sortedList.GetRange(maxStart, maxLength);
        }
    }
}
