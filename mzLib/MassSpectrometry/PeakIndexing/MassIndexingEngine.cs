using MathNet.Numerics.RootFinding;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class MassIndexingEngine: IndexingEngine <IndexedMass>
    {
        protected override int BinsPerDalton => 1;
        public int MaxMass { get; set; } = 30000;
        public ScanInfo[] ScanInfoArray { get; private set; }

        public MassIndexingEngine()
        {
        }

        public bool IndexPeaks(MsDataScan[] scanArray, DeconvolutionParameters deconParameters, MzRange mzRange = null, double minMass = 0, int minCharge = 1)
        {
            if (scanArray.IsNullOrEmpty() || scanArray.All(p => p == null))
                return false;

            IndexedPeaks = new List<IndexedMass> [MaxMass];
            ScanInfoArray = new ScanInfo[scanArray.Length];
            for (int scanIndex = 0; scanIndex < scanArray.Length; scanIndex++)
            {
                ScanInfoArray[scanIndex] = new ScanInfo(scanArray[scanIndex].OneBasedScanNumber, scanIndex, scanArray[scanIndex].RetentionTime);
                var envelopes = Deconvoluter.Deconvolute(scanArray[scanIndex].MassSpectrum, deconParameters, mzRange);
                foreach (var envelope in envelopes)
                {
                    if (envelope.MonoisotopicMass < minMass || envelope.Charge < minCharge)
                        continue;
                    int roundedMass = (int)Math.Round(envelope.MonoisotopicMass * BinsPerDalton, 0);
                    IndexedPeaks[roundedMass] ??= new List<IndexedMass>();
                    IndexedPeaks[roundedMass].Add(new IndexedMass(envelope, scanArray[scanIndex].RetentionTime, scanIndex, scanArray[scanIndex].MsnOrder));
                }
            }
            if (IndexedPeaks == null || IndexedPeaks.Length == 0)
                return false;
            else
                return true;
        }

        public IndexedMass GetIndexedPeak(double mass, int charge, int zeroBasedScanIndex, PpmTolerance ppmTolerance)
        {
            var bins = GetBinsInRange(mass, ppmTolerance);
            if (bins.Count == 0) return null;
            List<int> peakIndicesInBins = bins.Select(b => BinarySearchForIndexedPeak(b, zeroBasedScanIndex)).ToList();
            return GetBestPeakFromBins(bins, mass, charge, zeroBasedScanIndex, peakIndicesInBins, ppmTolerance);
        }

        public static IndexedMass GetBestPeakFromBins(List<List<IndexedMass>> allBins, double mass, int charge, int zeroBasedScanIndex, IList<int> peakIndicesInBins, PpmTolerance ppmTolerance)
        {
            IndexedMass bestMass = null;
            for (int i = 0; i < allBins.Count; i++)
            {
                var tempPeak = GetPeakFromBin(allBins[i], mass, charge, zeroBasedScanIndex, peakIndicesInBins[i], ppmTolerance);
                if (tempPeak.IsDefaultOrNull()) continue;
                // Check if the mass is within the tolerance, if the charge state is the same, and if it is closer to the target M than the current peak
                if (bestMass.IsDefaultOrNull() || Math.Abs(tempPeak.M - mass) < Math.Abs(bestMass.M - mass))
                {
                    bestMass = tempPeak;
                }
            }
            return bestMass;
        }

        public static IndexedMass GetPeakFromBin(List<IndexedMass> bin, double mass, int charge, int zeroBasedScanIndex, int peakIndexInBin, PpmTolerance ppmTolerance)
        {
            IndexedMass bestPeak = null;
            if (peakIndexInBin < 0 || peakIndexInBin >= bin.Count) return bestPeak;
            for (int i = peakIndexInBin; i < bin.Count; i++)
            {
                var peak = bin[i];

                if (peak.ZeroBasedScanIndex > zeroBasedScanIndex)
                {
                    break;
                }

                if (ppmTolerance.Within(peak.M, mass)
                    && peak.ZeroBasedScanIndex == zeroBasedScanIndex
                    && peak.Charge == charge
                    && (bestPeak.IsDefaultOrNull() || Math.Abs(peak.M - mass) < Math.Abs(bestPeak.M - mass)))
                {
                    bestPeak = peak;
                }
            }
            return bestPeak;
        }

        public List<IndexedMass> GetXic(double mass, int charge, int zeroBasedStartIndex, PpmTolerance ppmTolerance, int missedScansAllowed, double maxPeakHalfWidth = double.MaxValue)
        {
            if (IndexedPeaks == null || ScanInfoArray == null) throw new MzLibException("Error: Attempt to retrieve XIC before peak indexing was performed");

            List<IndexedMass> xic = new List<IndexedMass>();
            var initialMass = GetIndexedPeak(mass, charge, zeroBasedStartIndex, ppmTolerance);

            if (initialMass.IsNotDefaultOrNull())
                xic.Add(initialMass);

            foreach (int direction in new List<int> { -1, 1 })
            {
                int missedPeaks = 0; 
                int currentZeroBasedScanIndex = zeroBasedStartIndex;

                while (missedPeaks <= missedScansAllowed)
                {
                    // increment the scan index we're searching for
                    currentZeroBasedScanIndex += direction;
                    if (currentZeroBasedScanIndex < 0 || currentZeroBasedScanIndex > ScanInfoArray.Length - 1 || Math.Abs(xic.Last().RetentionTime - initialMass.RetentionTime) > maxPeakHalfWidth)
                        break;

                    // Search for the next peak
                    var nextPeak = GetIndexedPeak(mass, charge, currentZeroBasedScanIndex, ppmTolerance);

                    // Add the peak to the XIC or increment the missed peaks
                    if (nextPeak == null)
                        missedPeaks++;
                    else
                    {
                        xic.Add(nextPeak);
                        missedPeaks = 0;
                    }
                }
            }
            // Sort the XIC in place
            xic.Sort((x, y) => x.RetentionTime.CompareTo(y.RetentionTime));
            return xic;
        }
    }
}
