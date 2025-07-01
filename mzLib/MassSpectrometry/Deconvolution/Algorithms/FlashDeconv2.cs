using Chemistry;
using MassSpectrometry.Deconvolution.Parameters;
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
            FilterMassIntensityGroupsByPpmTolerance(
                massIntensityGroups,
                out var likelyCorrect,
                out var likelyIncorrect,
                correctPpmTolerance: 25,
                incorrectPpmTolerance: 250,
                correctFraction: 0.7);
            var neutralMassIntensityGroups = GetMostCommonNeutralMassAndSummedIntensity(
                likelyCorrect,
                ppmTolerance: 25);

            var neutralMassSpectrum = new NeutralMassSpectrum(
                neutralMassIntensityGroups.Select(g => g.mostCommonNeutralMass.ToMz(1)).ToArray(),
                    neutralMassIntensityGroups.Select(g => g.summedIntensity).ToArray(),
                    Enumerable.Repeat(1, neutralMassIntensityGroups.Count).ToArray(),
                    shouldCopy: true);
            FlashDeconvDeconvolutionParameters deconvolutionParameters = (FlashDeconvDeconvolutionParameters)DeconvolutionParameters;
            var k = Deconvoluter.Deconvolute(neutralMassSpectrum, deconvolutionParameters, range);
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
        
        public static void FilterMassIntensityGroupsByPpmTolerance(
            IEnumerable<(double[] neutralMass, double[] intensity)> massIntensityGroups,
            out List<(double[] neutralMass, double[] intensity)> likelyCorrect,
            out List<(double[] neutralMass, double[] intensity)> likelyIncorrect,
            double correctPpmTolerance = 25,
            double incorrectPpmTolerance = 250,
            double correctFraction = 0.7)
        {
            likelyCorrect = new List<(double[] neutralMass, double[] intensity)>();
            likelyIncorrect = new List<(double[] neutralMass, double[] intensity)>();

            foreach (var group in massIntensityGroups)
            {
                var masses = group.neutralMass;
                if (masses.Length < 2)
                {
                    // Not enough data to judge, treat as correct by default
                    likelyCorrect.Add(group);
                    continue;
                }

                int closeCount = 0;
                int farCount = 0;
                int totalPairs = 0;

                // Compare all pairs
                for (int i = 0; i < masses.Length; i++)
                {
                    for (int j = i + 1; j < masses.Length; j++)
                    {
                        double ppm = Math.Abs(masses[i] - masses[j]) / ((masses[i] + masses[j]) / 2.0) * 1e6;
                        if (ppm <= correctPpmTolerance)
                            closeCount++;
                        if (ppm > incorrectPpmTolerance)
                            farCount++;
                        totalPairs++;
                    }
                }

                // If most pairs are close, it's likely correct
                if (totalPairs == 0 || (closeCount >= correctFraction * totalPairs && farCount < (1 - correctFraction) * totalPairs))
                    likelyCorrect.Add(group);
                else if (farCount > (1 - correctFraction) * totalPairs)
                    likelyIncorrect.Add(group);
                else
                    likelyIncorrect.Add(group); // ambiguous, treat as incorrect
            }
        }
        public static List<(double mostCommonNeutralMass, double summedIntensity)> GetMostCommonNeutralMassAndSummedIntensity(
            IEnumerable<(double[] neutralMass, double[] intensity)> likelyCorrectGroups,
            double ppmTolerance = 25)
        {
            var results = new List<(double mostCommonNeutralMass, double summedIntensity)>();

            foreach (var group in likelyCorrectGroups)
            {
                var masses = group.neutralMass;
                var intensities = group.intensity;

                // Cluster masses within ppmTolerance
                var clusters = new List<List<int>>(); // Each cluster is a list of indices

                for (int i = 0; i < masses.Length; i++)
                {
                    bool added = false;
                    for (int c = 0; c < clusters.Count; c++)
                    {
                        // Compare to first mass in cluster
                        double refMass = masses[clusters[c][0]];
                        double ppm = Math.Abs(masses[i] - refMass) / refMass * 1e6;
                        if (ppm <= ppmTolerance)
                        {
                            clusters[c].Add(i);
                            added = true;
                            break;
                        }
                    }
                    if (!added)
                    {
                        clusters.Add(new List<int> { i });
                    }
                }

                // Find the largest cluster (the mode)
                var modeCluster = clusters.OrderByDescending(cl => cl.Count).First();
                // Use the average mass of the mode cluster as the representative value
                double mostCommonNeutralMass = modeCluster.Select(idx => masses[idx]).Average();
                double summedIntensity = modeCluster.Select(idx => intensities[idx]).Sum();

                results.Add((mostCommonNeutralMass, summedIntensity));
            }

            return results;
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

        //internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        //{
        //    throw new NotImplementedException();
        //}
    }
}
