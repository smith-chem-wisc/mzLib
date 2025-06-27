using MzLibUtil;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
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
            var arrayWideDifferences = StandardArrayWideLogMzDifferences(acceptibleLogMzDifferences, logTransformedSpectrum);
            var matches = CountMatchesWithinTolerance(arrayWideDifferences, logTransformedSpectrum.XArray);
            var matchesIntensities = ConvertMatchesToIntensities(matches, logTransformedSpectrum.YArray);
            List<string> result = new List<string>();
            result.Add(string.Join("\t", logTransformedSpectrum.XArray));
            foreach(var match in matches)
            {
                result.Add($"{match.Key}\t{string.Join("\t", match.Value)}");
            }
            foreach (var mi in matchesIntensities)
            {
                result.Add($"{mi.Key}\t{string.Join("\t", mi.Value)}");
            }
            File.WriteAllLines(@"C:\Users\trish\Downloads\result.tsv", result);

            var myList = new List<List<(int,int)>>();
 
            return new List<IsotopicEnvelope>();
        }

        private List<KeyValuePair<int, double[]>> ConvertMatchesToIntensities(List<KeyValuePair<int, int[]>> matches, double[] yArray)
        {
            List < KeyValuePair<int, double[]> > result = new List<KeyValuePair<int, double[]>>();
            foreach (var match in matches)
            {
                double[] intensities = match.Value.Zip(yArray, (i, d) => i * d).ToArray();
                result.Add(new KeyValuePair<int, double[]>(match.Key, intensities));
            }
            return result;
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
            .Select(i => new KeyValuePair<int, double>(i, Math.Log(i) - Math.Log(i-1)))
            .ToList();
        }
        private List<KeyValuePair<int, double[]>> StandardArrayWideLogMzDifferences(List<KeyValuePair<int, double>> acceptibleLogMzDifferences, MzSpectrum logTransformedSpectrum)
        {
            List<KeyValuePair<int, double[]>> runningSum = new List<KeyValuePair<int, double[]>>();
            double sum = 0;
            foreach (var pair in acceptibleLogMzDifferences)
            {
                runningSum.Add(new KeyValuePair<int, double[]>(pair.Key, logTransformedSpectrum.XArray.Select(x => x - pair.Value).ToArray()));
            }
            return runningSum;
        }

        List<KeyValuePair<int, int[]>> CountMatchesWithinTolerance(
            List<KeyValuePair<int, double[]>> arrayWideDifferences,
            double[] logTransformedSpectrumXarray,
            double tolerance = 10.0)
        {
            var result = new List<KeyValuePair<int, int[]>>();

            foreach (var pair in arrayWideDifferences)
            {
                var counts = new int[pair.Value.Length];
                int refIdx = 0;

                for (int i = 0; i < pair.Value.Length; i++)
                {
                    double target = pair.Value[i];
                    int count = 0;

                    // Advance refIdx to the first possible match
                    while (refIdx < logTransformedSpectrumXarray.Length && logTransformedSpectrumXarray[refIdx] < target - LogMzDependentTolerance(logTransformedSpectrumXarray[refIdx]))
                        refIdx++;

                    int tempIdx = refIdx;
                    // Count all matches within tolerance
                    while (tempIdx < logTransformedSpectrumXarray.Length && logTransformedSpectrumXarray[tempIdx] <= target + LogMzDependentTolerance(logTransformedSpectrumXarray[refIdx]))
                    {
                        if (Math.Abs(logTransformedSpectrumXarray[tempIdx] - target) <= LogMzDependentTolerance(logTransformedSpectrumXarray[refIdx]))
                            count++;
                        tempIdx++;
                    }
                    counts[i] = count;
                }

                result.Add(new KeyValuePair<int, int[]>(pair.Key, counts));
            }

            return result;
        }
        private double LogMzDependentTolerance(double logMz, double tolerance = 10.0)
        {
            var m = Math.Exp(logMz);
            var mPlus = m + m * tolerance / 1000000.0;
            var lmPlus = Math.Log(mPlus);
            var newT = lmPlus - logMz;
            return newT;
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
