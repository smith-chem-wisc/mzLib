using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;

namespace TopDownEngine.Features;

public sealed class CrossRunFeatureMatcher
{
    private static readonly MethodInfo GetBinsInRangeMethod =
        typeof(IndexingEngine<IndexedMassSpectralPeak>).GetMethod("GetBinsInRange", BindingFlags.Instance | BindingFlags.NonPublic)
        ?? throw new MzLibException("Unable to locate GetBinsInRange helper on PeakIndexingEngine.");

    public IReadOnlyList<FeatureGroup> MatchDonorBoxes(
        SpectraFileInfo donorFile,
        IReadOnlyList<FeatureBox> donorBoxes,
        IReadOnlyDictionary<SpectraFileInfo, PeakIndexingEngine> indexesByFile,
        IReadOnlyDictionary<SpectraFileInfo, Func<double, double>> rtWarpsToReference,
        IReadOnlyList<double> nullShifts,
        double pValueThreshold = 0.05)
    {
        if (donorFile == null)
        {
            throw new ArgumentNullException(nameof(donorFile));
        }

        if (donorBoxes == null)
        {
            throw new ArgumentNullException(nameof(donorBoxes));
        }

        if (indexesByFile == null)
        {
            throw new ArgumentNullException(nameof(indexesByFile));
        }

        if (rtWarpsToReference == null)
        {
            throw new ArgumentNullException(nameof(rtWarpsToReference));
        }

        if (nullShifts == null)
        {
            throw new ArgumentNullException(nameof(nullShifts));
        }

        if (!indexesByFile.ContainsKey(donorFile))
        {
            throw new ArgumentException("Missing donor file index.", nameof(indexesByFile));
        }

        if (!rtWarpsToReference.ContainsKey(donorFile))
        {
            throw new ArgumentException("Missing donor file RT warp.", nameof(rtWarpsToReference));
        }

        if (double.IsNaN(pValueThreshold) || double.IsInfinity(pValueThreshold) || pValueThreshold <= 0 || pValueThreshold > 1)
        {
            throw new ArgumentOutOfRangeException(nameof(pValueThreshold), "pValueThreshold must be in (0, 1].");
        }

        List<FeatureGroup> groups = new(donorBoxes.Count);
        foreach (FeatureBox donorBox in donorBoxes)
        {
            Dictionary<SpectraFileInfo, (int KTrue, double PValue, double Intensity)> matches = new();
            matches[donorFile] = (donorBox.PeakCount, 0.0, donorBox.TotalIntensity);

            foreach ((SpectraFileInfo acceptorFile, PeakIndexingEngine acceptorIndex) in indexesByFile)
            {
                if (ReferenceEquals(acceptorFile, donorFile))
                {
                    continue;
                }

                if (!rtWarpsToReference.ContainsKey(acceptorFile))
                {
                    continue;
                }

                Func<double, double> donorToAcceptorRtWarp = BuildDonorToAcceptorWarp(
                    donorFile,
                    acceptorFile,
                    indexesByFile,
                    rtWarpsToReference);

                BinomialMatchScore score = BinomialScorer.ScoreMatch(
                    donorBox,
                    acceptorIndex,
                    donorToAcceptorRtWarp,
                    nullShifts);

                if (score.PValue > pValueThreshold)
                {
                    continue;
                }

                DoubleRange warpedRtRange = WarpRtRange(donorBox.RtRange, donorToAcceptorRtWarp);
                DoubleRange mzRange = new(donorBox.MzRange.Minimum, donorBox.MzRange.Maximum);
                double matchedIntensity = SumIntensityInBox(acceptorIndex, warpedRtRange, mzRange);

                matches[acceptorFile] = (score.KTrue, score.PValue, matchedIntensity);
            }

            groups.Add(new FeatureGroup(donorBox, matches, matches.Count));
        }

        return groups;
    }

    private static Func<double, double> BuildDonorToAcceptorWarp(
        SpectraFileInfo donorFile,
        SpectraFileInfo acceptorFile,
        IReadOnlyDictionary<SpectraFileInfo, PeakIndexingEngine> indexesByFile,
        IReadOnlyDictionary<SpectraFileInfo, Func<double, double>> rtWarpsToReference)
    {
        if (ReferenceEquals(donorFile, acceptorFile))
        {
            return static rt => rt;
        }

        Func<double, double> donorToReference = rtWarpsToReference[donorFile];
        Func<double, double> acceptorToReference = rtWarpsToReference[acceptorFile];
        DoubleRange acceptorRtDomain = GetRtDomain(indexesByFile[acceptorFile]);

        return donorRt =>
        {
            double referenceRt = donorToReference(donorRt);
            return InvertMonotoneWarp(acceptorToReference, acceptorRtDomain, referenceRt);
        };
    }

    private static DoubleRange GetRtDomain(PeakIndexingEngine index)
    {
        if (index.ScanInfoArray == null || index.ScanInfoArray.Length == 0)
        {
            return new DoubleRange(0.0, 0.0);
        }

        double minRt = index.ScanInfoArray.Min(s => s.RetentionTime);
        double maxRt = index.ScanInfoArray.Max(s => s.RetentionTime);
        return new DoubleRange(minRt, maxRt);
    }

    private static double InvertMonotoneWarp(Func<double, double> warp, DoubleRange domain, double targetValue)
    {
        double low = domain.Minimum;
        double high = domain.Maximum;
        if (high <= low)
        {
            return low;
        }

        double lowValue = warp(low);
        double highValue = warp(high);
        bool increasing = highValue >= lowValue;

        if (increasing)
        {
            if (targetValue <= lowValue)
            {
                return low;
            }

            if (targetValue >= highValue)
            {
                return high;
            }

            for (int i = 0; i < 48; i++)
            {
                double mid = 0.5 * (low + high);
                if (warp(mid) < targetValue)
                {
                    low = mid;
                }
                else
                {
                    high = mid;
                }
            }
        }
        else
        {
            if (targetValue >= lowValue)
            {
                return low;
            }

            if (targetValue <= highValue)
            {
                return high;
            }

            for (int i = 0; i < 48; i++)
            {
                double mid = 0.5 * (low + high);
                if (warp(mid) > targetValue)
                {
                    low = mid;
                }
                else
                {
                    high = mid;
                }
            }
        }

        return 0.5 * (low + high);
    }

    private static DoubleRange WarpRtRange(DoubleRange donorRtRange, Func<double, double> rtWarp)
    {
        double warpedMin = rtWarp(donorRtRange.Minimum);
        double warpedMax = rtWarp(donorRtRange.Maximum);
        return new DoubleRange(Math.Min(warpedMin, warpedMax), Math.Max(warpedMin, warpedMax));
    }

    private static double SumIntensityInBox(PeakIndexingEngine index, DoubleRange rtRange, DoubleRange mzRange)
    {
        if (index.ScanInfoArray == null || rtRange.Maximum <= rtRange.Minimum || mzRange.Maximum <= mzRange.Minimum)
        {
            return 0.0;
        }

        HashSet<int> scanIndicesInWindow = new();
        for (int i = 0; i < index.ScanInfoArray.Length; i++)
        {
            ScanInfo scanInfo = index.ScanInfoArray[i];
            if (scanInfo.RetentionTime >= rtRange.Minimum && scanInfo.RetentionTime < rtRange.Maximum)
            {
                scanIndicesInWindow.Add(scanInfo.ZeroBasedScanIndex);
            }
        }

        if (scanIndicesInWindow.Count == 0)
        {
            return 0.0;
        }

        double centerMz = (mzRange.Minimum + mzRange.Maximum) / 2.0;
        double halfWidth = (mzRange.Maximum - mzRange.Minimum) / 2.0;
        AbsoluteTolerance mzTolerance = new(halfWidth);

        var bins = (List<List<IndexedMassSpectralPeak>>)GetBinsInRangeMethod.Invoke(index, new object[] { centerMz, mzTolerance });
        double total = 0.0;
        foreach (List<IndexedMassSpectralPeak> bin in bins)
        {
            foreach (IndexedMassSpectralPeak peak in bin)
            {
                if (scanIndicesInWindow.Contains(peak.ZeroBasedScanIndex)
                    && peak.M >= mzRange.Minimum
                    && peak.M < mzRange.Maximum)
                {
                    total += peak.Intensity;
                }
            }
        }

        return total;
    }
}
