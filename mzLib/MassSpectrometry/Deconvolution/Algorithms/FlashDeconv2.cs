using Chemistry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry.Deconvolution.Algorithms
{
    internal class FlashDeconv2 : DeconvolutionAlgorithm
    {
        internal FlashDeconv2(DeconvolutionParameters deconParameters) : base(deconParameters)
        {

        }
        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            var logTransformedSpectrum = LogTransformSpectrum(spectrum);
            var logX = logTransformedSpectrum.XArray;
            var logY = logTransformedSpectrum.YArray;
            var acceptibleLogMzDifferencesBetweenNearbyValues = AllAcceptibleLogMzDifferencesForAdjacentValues();
            var groups = FindMatchingGroups(logX, logY, acceptibleLogMzDifferencesBetweenNearbyValues);
            var massIntensityGroups = CreateNeutralMassIntensityGroups(groups);
            return new List<IsotopicEnvelope>();
        }
        private MzSpectrum LogTransformSpectrum(MzSpectrum spectrum, double intensityThresholdForFilter = 0.01)
        {
            var filtered = spectrum.XArray
                .Zip(spectrum.YArray, (x, y) => new { x, y })
                .Where(pair => pair.y > intensityThresholdForFilter)
                .ToArray();

            double[] xArray = filtered.Select(pair => Math.Log(pair.x)).ToArray();
            double[] yArray = filtered.Select(pair => pair.y).ToArray();

            return new MzSpectrum(xArray, yArray, true);
        }
        private List<double> AllAcceptibleLogMzDifferencesForAdjacentValues(int lowValue = 1, int highValue = 60)
        {
            return Enumerable.Range(lowValue, highValue)
            .Select(i => Math.Log(i+1)-Math.Log(i))
            .ToList();
        }
        // Finds all groups in logTransformedXArray where the differences between consecutive elements
        // (from high to low) match (within a tolerance) a subsequence of acceptibleLogMzDifferences.
        // Returns a list of (X[], Y[]) for each group found.
        private static List<(double[] X, double[] Y, int[] ChargeState)> FindMatchingGroups(
            double[] logTransformedXArray,
            double[] yArray,
            List<double> acceptibleLogMzDifferences)
        {
            var results = new List<(double[], double[], int[])>();
            int n = logTransformedXArray.Length;

            // Sort by descending m/z (i.e., descending log(m/z))
            var sorted = logTransformedXArray
                .Select((val, idx) => (val, idx))
                .OrderByDescending(x => x.val)
                .ToList();

            var sortedLogX = sorted.Select(x => x.val).ToArray();
            var sortedY = sorted.Select(x => yArray[x.idx]).ToArray();
            var sortedIndices = sorted.Select(x => x.idx).ToArray();

            for (int start = 0; start < n - 1; start++)
            {
                var groupIndices = new List<int> { start };
                List<int> groupCharges = new ();
                int prevIdx = start;
                double longRangeExpectedDiff = 0;
                double firstValue = sortedLogX[start];
                for (int p = 0; p < acceptibleLogMzDifferences.Count; p++)
                {
                    double expectedDiff = acceptibleLogMzDifferences[p];
                    double prevValue = sortedLogX[prevIdx];
                    double tolerance = LogMzDependentTolerance(prevValue);
                    bool found = false;
                    for (int nextIdx = prevIdx + 1; nextIdx < n; nextIdx++)
                    {
                        double diff = prevValue - sortedLogX[nextIdx]; // Correct direction: high to low
                        if (Math.Abs(diff - expectedDiff) <= tolerance)
                        {
                            longRangeExpectedDiff += diff; // Accumulate the difference
                            double longRangeDiff = firstValue - sortedLogX[nextIdx]; // Calculate the long-range expected difference
                            if (Math.Abs(longRangeDiff - longRangeExpectedDiff) <= tolerance)
                            {
                                if(groupCharges.Count == 0) 
                                {
                                    groupCharges.Add(p + 1); // Charge states are 1-indexed
                                    groupCharges.Add(p + 2); // Charge states are 1-indexed
                                }
                                else
                                {
                                    groupCharges.Add(p + 2); // Charge states are 1-indexed
                                }
                                groupIndices.Add(nextIdx);
                                prevIdx = nextIdx;
                                found = true;
                                break;
                            }
                        }
                    }
                    if (!found)
                        break; // Stop if the next expected difference is not found
                }

                if (groupIndices.Count > 1)
                {
                    var groupX = groupIndices.Select(i => sortedLogX[i]).ToArray();
                    var groupY = groupIndices.Select(i => sortedY[i]).ToArray();
                    results.Add((groupX, groupY, groupCharges.ToArray()));
                }
            }
            return results;
        }
        // Removes any group from 'groups' where the double[] X is a subset of any other double[] X
        private static List<(double[] X, double[] Y)> RemoveSubsetGroups(List<(double[] X, double[] Y)> groups)
        {
            // Sort groups by length of X (ascending)
            var sortedGroups = groups
                .Select((g, idx) => (g, idx))
                .OrderBy(x => x.g.X.Length)
                .ToList();

            var toRemove = new HashSet<int>();
            int groupCount = sortedGroups.Count;

            for (int i = 0; i < groupCount; i++)
            {
                var (group, idx) = sortedGroups[i];
                var xArray = group.X;

                for (int j = i + 1; j < groupCount; j++)
                {
                    var (otherGroup, otherIdx) = sortedGroups[j];
                    var otherX = otherGroup.X;

                    // Check if every value in xArray has a match in otherX within tolerance
                    bool isSubset = xArray.All(x =>
                        otherX.Any(ox => Math.Abs(ox - x) <= LogMzDependentTolerance(x))
                    );

                    if (isSubset && xArray.Length < otherX.Length)
                    {
                        toRemove.Add(idx);
                        break;
                    }
                }
            }

            // Return groups not marked for removal, preserving original order
            return groups
                .Where((g, idx) => !toRemove.Contains(idx))
                .ToList();
        }
        // Creates new groups from filteredGroups where each value in double[] X is transformed by the inverse natural log (Math.Exp)
        private static List<(double[] X, double[] Y)> TransformGroupsToExpX(List<(double[] X, double[] Y)> filteredGroups)
        {
            return filteredGroups
                .Select(g => (X: g.X.Select(Math.Exp).ToArray(), Y: g.Y))
                .ToList();
        }

        private static List<(double[] neutralMass, double[] intensity)> CreateNeutralMassIntensityGroups(
            List<(double[] X, double[] Y, int[] ChargeState)> groups)
        {
            var result = new List<(double[] neutralMass, double[] intensity)>();

            foreach (var group in groups)
            {
                int len = group.X.Length;
                var neutralMass = new double[len];
                for (int i = 0; i < len; i++)
                {
                    neutralMass[i] = NeutralMassFromLogMz(group.X[i], group.ChargeState[i]);
                }
                result.Add((neutralMass, group.Y));
            }

            return result;
        }
        private static double LogMzDependentTolerance(double logMz, double tolerance = 250.0)
        {
            var m = Math.Exp(logMz);
            var mPlus = m + m * tolerance / 1000000.0;
            var lmPlus = Math.Log(mPlus);
            var newT = lmPlus - logMz;
            return newT;
            //return 0.0001;
        }

        private static double NeutralMassFromLogMz(double logmz, int chargeState)
        {
            return Math.Exp(logmz).ToMass(chargeState);
        }
    }
}
