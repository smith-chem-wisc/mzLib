using Chemistry;
using MzLibUtil;
using System;
using System.Collections.Generic;
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
            var acceptibleLogMzDifferences = AllAcceptibleLogMzDifferences();
            var groups = FindMatchingGroups(logX, logY, acceptibleLogMzDifferences);
            var filteredGroups = RemoveSubsetGroups(groups);
            var transformedFilteredGroups = TransformGroupsToExpX(filteredGroups);
            var mzValuesWithIntensitiesAndChargeStates = AssignChargeStatesToGroups(transformedFilteredGroups, acceptibleLogMzDifferences);
            var massIntensityGroups = CreateMassIntensityGroups(mzValuesWithIntensitiesAndChargeStates);
            return new List<IsotopicEnvelope>();
        }
        private MzSpectrum LogTransformSpectrum(MzSpectrum spectrum)
        {
            return new MzSpectrum(spectrum.XArray.Select(x => Math.Log(x)).ToArray(), spectrum.YArray, true);
        }
        private List<double> AllAcceptibleLogMzDifferences(int lowValue = 1, int highValue = 60)
        {
            return Enumerable.Range(lowValue, highValue)
            .Select(i => Math.Log(i+1)-Math.Log(i))
            .ToList();
        }
        // Finds all groups in logTransformedXArray where the differences between consecutive elements
        // (from high to low) match (within a tolerance) a subsequence of acceptibleLogMzDifferences.
        // Returns a list of (X[], Y[]) for each group found.
        private static List<(double[] X, double[] Y)> FindMatchingGroups(
            double[] logTransformedXArray,
            double[] yArray,
            List<double> acceptibleLogMzDifferences)
        {
            var results = new List<(double[], double[])>();
            int n = logTransformedXArray.Length;
            var reversedPattern = acceptibleLogMzDifferences.AsEnumerable().Reverse().ToList();

            for (int start = 0; start < n - 1; start++)
            {
                var groupIndices = new List<int> { start };
                int prevIdx = start;
                int pIdx = 0;

                while (pIdx < reversedPattern.Count)
                {
                    bool found = false;
                    double expectedDiff = reversedPattern[pIdx];
                    double prevValue = logTransformedXArray[prevIdx];
                    double tolerance = LogMzDependentTolerance(prevValue);

                    // Search for the next index that matches the expected difference within tolerance
                    for (int nextIdx = prevIdx + 1; nextIdx < n; nextIdx++)
                    {
                        double diff = logTransformedXArray[nextIdx] - prevValue;
                        if (Math.Abs(diff - expectedDiff) <= tolerance)
                        {
                            groupIndices.Add(nextIdx);
                            prevIdx = nextIdx;
                            found = true;
                            break;
                        }
                    }

                    pIdx++;
                    if (!found)
                    {
                        // Allow missing pattern value, just move to next pattern
                        continue;
                    }
                }

                if (groupIndices.Count > 1)
                {
                    var groupX = groupIndices.Select(i => logTransformedXArray[i]).ToArray();
                    var groupY = groupIndices.Select(i => yArray[i]).ToArray();
                    results.Add((groupX, groupY));
                }
            }
            return results;
        }
        // Removes any group from 'groups' where the double[] X is a subset of any other double[] X
        private static List<(double[] X, double[] Y)> RemoveSubsetGroups(List<(double[] X, double[] Y)> groups)
        {
            var result = new List<(double[] X, double[] Y)>();
            var groupCount = groups.Count;

            for (int i = 0; i < groupCount; i++)
            {
                var xSet = new HashSet<double>(groups[i].X);
                bool isSubset = false;

                for (int j = 0; j < groupCount; j++)
                {
                    if (i == j) continue;
                    var otherSet = new HashSet<double>(groups[j].X);
                    if (xSet.IsSubsetOf(otherSet) && xSet.Count < otherSet.Count)
                    {
                        isSubset = true;
                        break;
                    }
                }

                if (!isSubset)
                    result.Add(groups[i]);
            }

            return result;
        }
        // Creates new groups from filteredGroups where each value in double[] X is transformed by the inverse natural log (Math.Exp)
        private static List<(double[] X, double[] Y)> TransformGroupsToExpX(List<(double[] X, double[] Y)> filteredGroups)
        {
            return filteredGroups
                .Select(g => (X: g.X.Select(Math.Exp).ToArray(), Y: g.Y))
                .ToList();
        }
        // Assigns charge states to each value in double[] X for each group, based on the pattern of log differences
        private static List<(double[] X, double[] Y, int[] chargeStates)> AssignChargeStatesToGroups(
            List<(double[] X, double[] Y)> transformedFilteredGroups,
            List<double> acceptibleLogMzDifferences)
        {
            // Map: log(n) => n for n in [2, 60]
            var logToInt = Enumerable.Range(2, 59)
                .ToDictionary(n => Math.Log(n), n => n);

            var result = new List<(double[] X, double[] Y, int[] chargeStates)>();

            foreach (var group in transformedFilteredGroups)
            {
                int len = group.X.Length;
                var chargeStates = new int[len];

                // Start with charge state 1 at the highest X (first, since X is sorted high-to-low in pattern)
                chargeStates[0] = 1;

                for (int i = 1; i < len; i++)
                {
                    double diff = Math.Log(group.X[i - 1]) - Math.Log(group.X[i]);
                    int chargePrev = chargeStates[i - 1];
                    int matchedCharge = -1;

                    // Find the integer n such that diff ≈ log(n) within tolerance
                    foreach (var kvp in logToInt)
                    {
                        double tolerance = LogMzDependentTolerance(Math.Log(group.X[i - 1]));
                        if (Math.Abs(diff - kvp.Key) <= tolerance)
                        {
                            matchedCharge = chargePrev + 1;
                            break;
                        }
                    }

                    // If not found, keep previous charge or set to -1 (unknown)
                    chargeStates[i] = matchedCharge > 0 ? matchedCharge : -1;
                }

                result.Add((group.X, group.Y, chargeStates));
            }

            return result;
        }
        private static double LogMzDependentTolerance(double logMz, double tolerance = 10.0)
        {
            var m = Math.Exp(logMz);
            var mPlus = m + m * tolerance / 1000000.0;
            var lmPlus = Math.Log(mPlus);
            var newT = lmPlus - logMz;
            return newT;
        }
        // Creates new groups where mass is computed using NeutralMassFromMz for each X and chargeState
        private List<(double[] mass, double[] intensity)> CreateMassIntensityGroups(
            List<(double[] X, double[] Y, int[] chargeStates)> mzValuesWithIntensitiesAndChargeStates)
        {
            var result = new List<(double[] mass, double[] intensity)>();

            foreach (var group in mzValuesWithIntensitiesAndChargeStates)
            {
                int len = group.X.Length;
                var mass = new double[len];
                for (int i = 0; i < len; i++)
                {
                    mass[i] = NeutralMassFromMz(group.X[i], group.chargeStates[i]);
                }
                result.Add((mass, group.Y));
            }

            return result;
        }
        private static double NeutralMassFromMz(double mz, int chargeState)
        {
            return mz.ToMass(chargeState);
        }
    }
}
