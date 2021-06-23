using BayesianEstimation;
using Chemistry;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    public class ImputationEngine
    {
        private Dictionary<(SpectraFileInfo, int rt, DoubleRange mzWindow), List<(double hdiStart, double hdiEnd)>> FileAndScanMzRangeToNoise;

        public ImputationEngine()
        {
            FileAndScanMzRangeToNoise = new Dictionary<(SpectraFileInfo, int rt, DoubleRange mzWindow), List<(double hdiStart, double hdiEnd)>>();
        }

        public void CalculateIntensityBaselines(SpectraFileInfo file, MsDataScan[] fileData)
        {
            foreach (MsDataScan scan in fileData.Where(p => p != null))
            {
                List<(double mz, double intensity)> peaks = new List<(double mz, double intensity)>();

                for (double minMz = 0; minMz < scan.ScanWindowRange.Maximum; minMz += 100)
                {
                    peaks.Clear();
                    double maxMz = minMz + 100;

                    int index = scan.MassSpectrum.GetClosestPeakIndex(minMz);

                    for (int i = index; i < scan.MassSpectrum.XArray.Length; i++)
                    {
                        double mz = scan.MassSpectrum.XArray[i];
                        double intensity = scan.MassSpectrum.YArray[i];

                        if (mz > maxMz)
                        {
                            break;
                        }
                        if (mz < minMz || intensity <= 0)
                        {
                            continue;
                        }

                        peaks.Add((mz, intensity));
                    }

                    if (peaks.Count < 20)
                    {
                        continue;
                    }

                    var peakIntensities = peaks.Select(p => p.intensity).ToArray();
                    var hdi = Util.GetHighestDensityInterval(peakIntensities, 0.1);

                    var mzRange = new DoubleRange((int)minMz, (int)maxMz);
                    int roundedRt = (int)scan.RetentionTime;

                    FileAndScanMzRangeToNoise.TryAdd((file, roundedRt, mzRange), new List<(double hdiStart, double hdiEnd)>());

                    FileAndScanMzRangeToNoise[(file, roundedRt, mzRange)].Add(hdi);
                }
            }
        }

        public ChromatographicPeak ImputePeak(SpectraFileInfo file, Identification id, double rt, List<(double mass, double abundance)> theoreticalIsotopeAbundances, int seed)
        {
            var rand = new Random(seed);
            double mz = id.PeakfindingMass.ToMz(id.PrecursorChargeState);

            double mzRangeStart = (int)(mz - (mz % 100));
            double mzRangeEnd = mzRangeStart + 100;

            if (!FileAndScanMzRangeToNoise.TryGetValue((file, (int)rt, new DoubleRange(mzRangeStart, mzRangeEnd)), out var hdiList))
            {
                return null;
            }

            var hdi = hdiList[rand.Next(0, hdiList.Count)];

            double imputedPeakIntensity = rand.Next((int)hdi.hdiStart, (int)hdi.hdiEnd);
            double imputedEnvelopeIntensity = 0;

            // impute isotope peak intensities for the envelope
            for (int i = 0; i < theoreticalIsotopeAbundances.Count; i++)
            {
                imputedEnvelopeIntensity += theoreticalIsotopeAbundances[i].abundance * imputedPeakIntensity;
            }

            var peak = new ChromatographicPeak(id, PeakType.Imputed, file);

            peak.IsotopicEnvelopes.Add(new IsotopicEnvelope(new IndexedMassSpectralPeak(mz, imputedPeakIntensity, 0, rt), id.PrecursorChargeState, imputedEnvelopeIntensity));
            peak.CalculateIntensityForThisFeature(false);

            return peak;
        }
    }
}
