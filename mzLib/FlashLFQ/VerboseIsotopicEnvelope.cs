using Chemistry;
using MathNet.Numerics;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    public class VerboseIsotopicEnvelope : IsotopicEnvelope
    {
        public Dictionary<int, IndexedMassSpectralPeak> PeakDictionary { get; }

        public VerboseIsotopicEnvelope(
          IndexedMassSpectralPeak mostAbundantPeak,
          List<IndexedMassSpectralPeak> allPeaks,
          int chargeState,
          double monoisotopicMass,
          int isotopePpmTolerance = 5,
          double intensity = -1.0)
          : base(mostAbundantPeak, chargeState, intensity > -1.0 ? intensity : allPeaks.Sum(p => p.Intensity))
        {
            this.PeakDictionary = WritePeakDictionary(allPeaks.OrderBy(peak => peak.Mz).ToList(), monoisotopicMass.ToMz(chargeState), isotopePpmTolerance);
            this.RetentionTime = mostAbundantPeak.RetentionTime.Round(4);
        }

        public double RetentionTime { get; }

        public override string ToString() => "+" + this.ChargeState.ToString() + "|" + this.Intensity.ToString("F0") + "|" + this.IndexedPeak.RetentionTime.ToString("F3") + "|" + this.IndexedPeak.ZeroBasedMs1ScanIndex.ToString();

        public static Dictionary<int, IndexedMassSpectralPeak> WritePeakDictionary(
          List<IndexedMassSpectralPeak> peaks,
          double monoisotopicMz,
          int isotopePpmTolerance)
        {
            Dictionary<int, IndexedMassSpectralPeak> dictionary = new Dictionary<int, IndexedMassSpectralPeak>();
            PpmTolerance ppmTolerance = new PpmTolerance((double)isotopePpmTolerance);
            int index = 0;
            int num = 0;
            while (dictionary.Count < peaks.Count)
            {
                if (ppmTolerance.Within(peaks[index].Mz, monoisotopicMz + 1.0033548381 * (double)num))
                {
                    dictionary.Add(num++, peaks[index++]);
                }
                else
                {
                    if (ppmTolerance.GetMinimumValue(monoisotopicMz + 1.0033548381 * (double)num) > peaks[index].Mz)
                        return WritePeakDictionary(peaks, monoisotopicMz, isotopePpmTolerance + 5);
                    ++num;
                }
            }
            return dictionary;
        }

        public static string GetIsotopePeakName((int isotopeNumber, int chargeState) key) => "i" + key.isotopeNumber.ToString() + "z" + key.chargeState.ToString();
    }
}