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
            var myList = new List<List<(int,int)>>();
            for (int i = 0; i < logTransformedSpectrum.XArray.Length; i++)
            {
                var l = FindIndexedLogMassDifferences(acceptibleLogMzDifferences, logTransformedSpectrum, i);
                if (l.Count > 1)
                {
                    myList.Add(l);
                }
            }
            
            return new List<IsotopicEnvelope>();

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
        private List<KeyValuePair<int,double>> AllAcceptibleLogMzDifferences(int lowValue = 2, int highValue = 60)
        {
            return Enumerable.Range(lowValue, highValue)
            .Select(i => new KeyValuePair<int, double>(i, Math.Log(i)))
            .ToList();
        }
        private List<(int, int)> FindIndexedLogMassDifferences(List<KeyValuePair<int, double>> acceptibleLogMzDifferences, MzSpectrum logTransformedSpectrum, int selectedLogMzIndex)
        {
            double tolerance = 0.01;
            //List<(int, int)> indexedLogMassDifferences = new() {(selectedLogMzIndex, 0) };
            List<(int, int)> indexedLogMassDifferences = new();
            foreach (var pair in acceptibleLogMzDifferences)
            {
                double targetLogMz = logTransformedSpectrum.XArray[selectedLogMzIndex] + pair.Value;
                if (targetLogMz > (logTransformedSpectrum.XArray.Max() + tolerance))
                    break;
                else
                {
                    var indexOfTargetMzInTheLogMzArray = FindIndexWithinTolerance(logTransformedSpectrum.XArray, targetLogMz, tolerance);
                    if (indexOfTargetMzInTheLogMzArray >= 0)
                    {
                        indexedLogMassDifferences.Add((indexOfTargetMzInTheLogMzArray, pair.Key));
                    }
                    else
                    {
                        // If no match is found, we can break or continue based on the logic needed
                        // For now, we will just continue to the next pair
                        continue;
                    }
                }
            }
            return indexedLogMassDifferences;
        }
        // Finds the index of the first value in a sorted list of doubles matching the target within a given tolerance
        private int FindIndexWithinTolerance(double[] sortedList, double target, double tolerance)
        {
            for (int i = 0; i < sortedList.Length; i++)
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
