using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

namespace TopDownEngine.Indexing;

public class ThickIndexView
{
    private const int FineBinsPerDalton = 100;
    private const int ThickBinsPerDalton = 10;
    private const int FineBinsPerThickBin = FineBinsPerDalton / ThickBinsPerDalton;
    private const double HalfThickBinWidth = 0.05;

    private static readonly FieldInfo IndexedPeaksField =
        typeof(IndexingEngine<IndexedMassSpectralPeak>).GetField("IndexedPeaks", BindingFlags.Instance | BindingFlags.NonPublic)
        ?? throw new MzLibException("Unable to locate indexed peak storage on PeakIndexingEngine.");

    private static readonly MethodInfo BinarySearchMethod =
        typeof(IndexingEngine<IndexedMassSpectralPeak>).GetMethod("BinarySearchForIndexedPeak", BindingFlags.Static | BindingFlags.NonPublic)
        ?? throw new MzLibException("Unable to locate BinarySearchForIndexedPeak helper.");

    private readonly List<IndexedMassSpectralPeak>[] _thickIndexedPeaks;

    public ThickIndexView(PeakIndexingEngine fineIndex)
    {
        FineIndex = fineIndex ?? throw new ArgumentNullException(nameof(fineIndex));
        var fineIndexedPeaks = IndexedPeaksField.GetValue(FineIndex) as List<IndexedMassSpectralPeak>[]
            ?? throw new MzLibException("Error: Attempt to construct ThickIndexView before peak indexing was performed.");
        _thickIndexedPeaks = CoarsenBins(fineIndexedPeaks);
    }

    public PeakIndexingEngine FineIndex { get; }

    public IIndexedPeak GetIndexedPeak(double mz, int zeroBasedScanIndex)
    {
        var bins = GetBinsInRange(mz);
        if (bins.Count == 0)
        {
            return null;
        }

        var peakIndicesInBins = bins.Select(b => BinarySearchForIndexedPeak(b, zeroBasedScanIndex)).ToList();
        return GetBestPeakFromBins(bins, mz, zeroBasedScanIndex, peakIndicesInBins);
    }

    public List<IIndexedPeak> GetXic(double mz, double retentionTime, int missedScansAllowed,
        double maxPeakHalfWidth = double.MaxValue, Dictionary<IIndexedPeak, ExtractedIonChromatogram> matchedPeaks = null)
    {
        if (FineIndex.ScanInfoArray == null)
        {
            throw new MzLibException("Error: Attempt to retrieve XIC before peak indexing was performed.");
        }

        int scanIndex = -1;
        foreach (ScanInfo scan in FineIndex.ScanInfoArray)
        {
            if (scan.RetentionTime < retentionTime)
            {
                scanIndex = scan.ZeroBasedScanIndex;
            }
            else
            {
                break;
            }
        }

        return GetXicByScanIndex(mz, scanIndex, missedScansAllowed, maxPeakHalfWidth, matchedPeaks);
    }

    public List<IIndexedPeak> GetXicByScanIndex(double mz, int zeroBasedStartIndex, int missedScansAllowed,
        double maxPeakHalfWidth = double.MaxValue, Dictionary<IIndexedPeak, ExtractedIonChromatogram> matchedPeaks = null)
    {
        if (FineIndex.ScanInfoArray == null)
        {
            throw new MzLibException("Error: Attempt to retrieve XIC before peak indexing was performed.");
        }

        var xic = new List<IIndexedPeak>();
        var allBins = GetBinsInRange(mz);
        if (allBins.Count == 0)
        {
            return xic;
        }

        int[] peakPointerArray = allBins.Select(b => BinarySearchForIndexedPeak(b, zeroBasedStartIndex)).ToArray();
        var initialPeak = GetBestPeakFromBins(allBins, mz, zeroBasedStartIndex, peakPointerArray);
        if (initialPeak != null)
        {
            xic.Add(initialPeak);
        }

        foreach (int direction in new[] { -1, 1 })
        {
            int missedPeaks = 0;
            int currentZeroBasedScanIndex = zeroBasedStartIndex;
            var pointerArrayCopy = new int[peakPointerArray.Length];
            Array.Copy(peakPointerArray, pointerArrayCopy, peakPointerArray.Length);

            while (missedPeaks <= missedScansAllowed)
            {
                currentZeroBasedScanIndex += direction;

                if (currentZeroBasedScanIndex < 0
                    || currentZeroBasedScanIndex > FineIndex.ScanInfoArray.Length - 1
                    || (initialPeak != null
                        && Math.Abs(FineIndex.ScanInfoArray[currentZeroBasedScanIndex].RetentionTime - initialPeak.RetentionTime) > maxPeakHalfWidth))
                {
                    break;
                }

                for (int i = 0; i < pointerArrayCopy.Length; i++)
                {
                    switch (direction)
                    {
                        case -1:
                            do
                            {
                                pointerArrayCopy[i]--;
                            } while (pointerArrayCopy[i] >= 0
                                && allBins[i][pointerArrayCopy[i]].ZeroBasedScanIndex > currentZeroBasedScanIndex - 1);
                            pointerArrayCopy[i]++;
                            break;
                        case 1:
                            while (pointerArrayCopy[i] < allBins[i].Count - 1
                                && allBins[i][pointerArrayCopy[i]].ZeroBasedScanIndex < currentZeroBasedScanIndex)
                            {
                                pointerArrayCopy[i]++;
                            }
                            break;
                    }
                }

                var nextPeak = GetBestPeakFromBins(allBins, mz, currentZeroBasedScanIndex, pointerArrayCopy);
                if (nextPeak == null || (matchedPeaks != null && matchedPeaks.ContainsKey(nextPeak)))
                {
                    missedPeaks++;
                }
                else
                {
                    if (initialPeak == null)
                    {
                        initialPeak = nextPeak;
                    }

                    xic.Add(nextPeak);
                    missedPeaks = 0;
                }
            }
        }

        xic.Sort((x, y) => x.RetentionTime.CompareTo(y.RetentionTime));
        return xic;
    }

    public List<List<IndexedMassSpectralPeak>> GetBinsInRange(double mz)
    {
        int ceilingMz = (int)Math.Ceiling((mz + HalfThickBinWidth) * ThickBinsPerDalton);
        int floorMz = (int)Math.Floor((mz - HalfThickBinWidth) * ThickBinsPerDalton);
        var bins = new List<List<IndexedMassSpectralPeak>>();

        for (int index = floorMz; index <= ceilingMz; index++)
        {
            if (index < 0 || index >= _thickIndexedPeaks.Length)
            {
                continue;
            }

            var bin = _thickIndexedPeaks[index];
            if (bin != null)
            {
                bins.Add(bin);
            }
        }

        return bins;
    }

    public Dictionary<int, float> GetBinMaximumIntensities(float minimumIntensity = 0)
    {
        var binMaximumIntensities = new Dictionary<int, float>();
        for (int i = 0; i < _thickIndexedPeaks.Length; i++)
        {
            var bin = _thickIndexedPeaks[i];
            if (bin == null || bin.Count == 0)
            {
                continue;
            }

            float maxIntensity = 0;
            foreach (var peak in bin)
            {
                if (peak.Intensity > maxIntensity)
                {
                    maxIntensity = peak.Intensity;
                }
            }

            if (maxIntensity >= minimumIntensity)
            {
                binMaximumIntensities[i] = maxIntensity;
            }
        }

        return binMaximumIntensities;
    }

    private static List<IndexedMassSpectralPeak>[] CoarsenBins(List<IndexedMassSpectralPeak>[] fineIndexedPeaks)
    {
        var coarsenedBins = new List<IndexedMassSpectralPeak>[(fineIndexedPeaks.Length + FineBinsPerThickBin - 1) / FineBinsPerThickBin];

        for (int fineBinIndex = 0; fineBinIndex < fineIndexedPeaks.Length; fineBinIndex++)
        {
            var fineBin = fineIndexedPeaks[fineBinIndex];
            if (fineBin == null)
            {
                continue;
            }

            int coarseBinIndex = fineBinIndex / FineBinsPerThickBin;
            coarsenedBins[coarseBinIndex] ??= new List<IndexedMassSpectralPeak>();
            coarsenedBins[coarseBinIndex].AddRange(fineBin);
        }

        foreach (var bin in coarsenedBins.Where(p => p != null))
        {
            bin.Sort((left, right) =>
            {
                int scanComparison = left.ZeroBasedScanIndex.CompareTo(right.ZeroBasedScanIndex);
                return scanComparison != 0 ? scanComparison : left.M.CompareTo(right.M);
            });
        }

        return coarsenedBins;
    }

    private static IIndexedPeak GetBestPeakFromBins(List<List<IndexedMassSpectralPeak>> allBins, double mz,
        int zeroBasedScanIndex, IList<int> peakIndicesInBins)
    {
        IIndexedPeak bestPeak = null;
        for (int i = 0; i < allBins.Count; i++)
        {
            var candidate = GetPeakFromBin(allBins[i], mz, zeroBasedScanIndex, peakIndicesInBins[i]);
            if (candidate == null)
            {
                continue;
            }

            if (bestPeak == null || Math.Abs(candidate.M - mz) < Math.Abs(bestPeak.M - mz))
            {
                bestPeak = candidate;
            }
        }

        return bestPeak;
    }

    private static IIndexedPeak GetPeakFromBin(List<IndexedMassSpectralPeak> bin, double mz, int zeroBasedScanIndex, int peakIndexInBin)
    {
        if (peakIndexInBin < 0 || peakIndexInBin >= bin.Count)
        {
            return null;
        }

        IIndexedPeak bestPeak = null;
        for (int i = peakIndexInBin; i < bin.Count; i++)
        {
            var peak = bin[i];
            if (peak.ZeroBasedScanIndex > zeroBasedScanIndex)
            {
                break;
            }

            if (peak.ZeroBasedScanIndex == zeroBasedScanIndex
                && Math.Abs(peak.M - mz) <= HalfThickBinWidth
                && (bestPeak == null || Math.Abs(peak.M - mz) < Math.Abs(bestPeak.M - mz)))
            {
                bestPeak = peak;
            }
        }

        return bestPeak;
    }

    private static int BinarySearchForIndexedPeak(List<IndexedMassSpectralPeak> indexedPeaks, int zeroBasedScanIndex)
    {
        return (int)BinarySearchMethod.Invoke(null, new object[] { indexedPeaks, zeroBasedScanIndex })!;
    }
}
