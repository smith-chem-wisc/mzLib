using IO.Thermo;
using MassSpectrometry;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace Deconvolution
{
    internal class Program
    {

        #region Private Methods

        private static void Main(string[] args)
        {
            int? minScan = null;
            int? maxScan = null;
            double deconvolutionTolerancePpm = 10;
            //string filename = @"C:\Users\stepa\Desktop\DeconvolutionStuff\03-01-17_B2A_targeted_td_yeast_fract6_intact.raw";
            //string filename = @"C:\Users\stepa\Desktop\neucode_for_deconv\10-23-15_A_fract5_rep2.raw";
            string filename = @"C:\Users\stepa\Desktop\neucode_for_deconv\10-28-15_S_fract5_rep1.raw";
            int maxAssumedChargeState = 30;
            double intensityRatioLimit = 2;
            Func<IMzPeak, bool> peakFilter = b => (b as ThermoMzPeak).SignalToNoise > 2;

            Loaders.LoadElements("elements2.dat");

            var f = ThermoDynamicData.InitiateDynamicConnection(filename);

            minScan = minScan ?? 1;
            maxScan = maxScan ?? f.NumSpectra;

            var allAggregateGroups = new List<IsotopicEnvelope>[maxScan.Value - minScan.Value + 1];
            Parallel.ForEach(Partitioner.Create(minScan.Value, maxScan.Value + 1), fff =>
            {
                for (int scanIndex = fff.Item1; scanIndex < fff.Item2; scanIndex++)
                {
                    var theScan = f.GetOneBasedScan(scanIndex);
                    Console.WriteLine("Deconvoluting scan " + theScan.OneBasedScanNumber);
                    allAggregateGroups[scanIndex - minScan.Value] = theScan.MassSpectrum.Deconvolute(theScan.ScanWindowRange, maxAssumedChargeState, deconvolutionTolerancePpm, intensityRatioLimit, peakFilter).ToList();
                }
            });

            double aggregationTolerancePpm = 10;

            IEnumerable<AggregatedIsotopeEnvelopes> nice = AggregateIsotopicEnvelopes(allAggregateGroups, aggregationTolerancePpm).ToList();

            Console.WriteLine(string.Join(Environment.NewLine, nice.OrderBy(b => -b.NumPeaks).Take(10)));

            using (System.IO.StreamWriter file =
            new System.IO.StreamWriter(@"out.txt"))
            {
                file.WriteLine(string.Join(Environment.NewLine, nice.OrderBy(b => -b.NumPeaks).Select(b => b.OneLineString())));
            }
        }

        private static IEnumerable<AggregatedIsotopeEnvelopes> AggregateIsotopicEnvelopes(List<IsotopicEnvelope>[] allEnvelopes, double aggregationTolerancePpm)
        {
            List<AggregatedIsotopeEnvelopes> currentListOfGroups = new List<AggregatedIsotopeEnvelopes>();
            for (int zeroBasedScanIndex = 0; zeroBasedScanIndex < allEnvelopes.Length; zeroBasedScanIndex++)
            {
                Console.WriteLine("Aggregating zeroBasedScanIndex: " + zeroBasedScanIndex);
                foreach (var isotopicEnvelope in allEnvelopes[zeroBasedScanIndex])
                {
                    AggregatedIsotopeEnvelopes matchingProteoform = GetMatchingGroup(isotopicEnvelope.monoisotopicMass, currentListOfGroups, aggregationTolerancePpm);

                    if (matchingProteoform == null)
                    {
                        var newGroupScans = new AggregatedIsotopeEnvelopes();
                        newGroupScans.AddEnvelope(isotopicEnvelope, zeroBasedScanIndex);
                        currentListOfGroups.Add(newGroupScans);
                    }
                    else
                    {
                        matchingProteoform.AddEnvelope(isotopicEnvelope, zeroBasedScanIndex);
                    }
                }
                foreach (var ok in currentListOfGroups.Where(b => b.maxScanIndex < zeroBasedScanIndex))
                    yield return ok;
                currentListOfGroups.RemoveAll(b => b.maxScanIndex < zeroBasedScanIndex);
            }
            foreach (var ok in currentListOfGroups)
                yield return ok;
        }

        private static AggregatedIsotopeEnvelopes GetMatchingGroup(double mass, List<AggregatedIsotopeEnvelopes> currentListOfGroups, double aggregationTolerancePpm)
        {
            foreach (var possibleGroup in currentListOfGroups)
            {
                if (Math.Abs(mass - possibleGroup.mass) / possibleGroup.mass * 1e6 <= aggregationTolerancePpm ||
                    Math.Abs(mass + 1.002868314 - possibleGroup.mass) / possibleGroup.mass * 1e6 <= aggregationTolerancePpm ||
                    Math.Abs(mass + 2.005408917 - possibleGroup.mass) / possibleGroup.mass * 1e6 <= aggregationTolerancePpm ||
                    Math.Abs(mass + 3.007841294 - possibleGroup.mass) / possibleGroup.mass * 1e6 <= aggregationTolerancePpm ||
                    Math.Abs(mass - 1.002868314 - possibleGroup.mass) / possibleGroup.mass * 1e6 <= aggregationTolerancePpm ||
                    Math.Abs(mass - 2.005408917 - possibleGroup.mass) / possibleGroup.mass * 1e6 <= aggregationTolerancePpm ||
                    Math.Abs(mass - 3.007841294 - possibleGroup.mass) / possibleGroup.mass * 1e6 <= aggregationTolerancePpm)
                    return possibleGroup;
            }
            return null;
        }

        #endregion Private Methods

    }
}