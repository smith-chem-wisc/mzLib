using Chemistry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;


namespace MassSpectrometry.Deconvolution.Algorithms
{
    internal class FlashDeconv2 : DeconvolutionAlgorithm
    {
        private static readonly ClassicDeconvolutionParameters ToEnvelopeDeconParams = new ClassicDeconvolutionParameters(1, 1, 10, 3);
        private static readonly ClassicDeconvolutionAlgorithm ToEnvelopeDeconvoluter = new(ToEnvelopeDeconParams);
        internal FlashDeconv2(DeconvolutionParameters deconParameters) : base(deconParameters)
        {

        }
        // Deconvolute method: Performs deconvolution of the input MzSpectrum within the specified MzRange.
        // The method extracts isotopic envelopes by transforming, grouping, filtering, and summarizing spectral features.
        internal override IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum, MzRange range)
        {
            // 1. Log-transform the spectrum to linearize charge state spacing and filter low-intensity peaks.
            var logTransformedSpectrum = LogTransformSpectrum(spectrum);

            // 3. Find groups of peaks that match expected charge state patterns in log(m/z) space.
            var matchingGroups = FindMatchingGroups(
                logTransformedSpectrum.XArray,
                logTransformedSpectrum.YArray);

            // 4. Remove groups that are subsets of larger groups to avoid redundancy.
            var filteredGroups = RemoveSubsetGroups(
                matchingGroups.Select(g => (g.X, g.Y, g.ChargeState)).ToList());

            // 5. Transform the filtered groups' X values back from log(m/z) to m/z space.
            var expTransformedGroups = TransformGroupsToExpX(filteredGroups);

            // 6. Create neutral mass/intensity groups from the charge state groups.
            var neutralMassIntensityGroups = CreateNeutralMassIntensityGroups(expTransformedGroups);

            // 7. Filter the neutral mass/intensity groups into likely correct and incorrect groups based on ppm tolerance.
            FilterMassIntensityGroupsByPpmTolerance(
                neutralMassIntensityGroups,
                out var likelyCorrectGroups,
                out var likelyIncorrectGroups);

            // 8. For each likely correct group, determine the most common neutral mass (mode) and sum the corresponding intensities.
            var summarizedEnvelopes = GetMostCommonNeutralMassAndSummedIntensity(likelyCorrectGroups);

            var nms = new MzSpectrum(
                summarizedEnvelopes.Select(x => x.mostCommonNeutralMass.ToMz(1)).ToArray(),
                summarizedEnvelopes.Select(x => x.summedIntensity).ToArray(),
                true);

            var envelopes = ToEnvelopeDeconvoluter.Deconvolute(nms, nms.Range);
            return envelopes;
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
        private static List<double> AllAcceptibleLogMzDifferencesForAdjacentValues(int lowValue = 1, int highValue = 60)
        {
            var a = Enumerable.Range(lowValue, highValue - lowValue + 1).ToArray();
            var b = a.Select(i => (double)i).ToArray();
            var k = b.Select(i => Math.Log(i)).ToArray();

            List<double> acceptibleLogMzDifferences = new List<double>();
            acceptibleLogMzDifferences.Add(k[0]);
            for (int i = 1; i < k.Length; i++)
            {
                acceptibleLogMzDifferences.Add(k[i] - k[i - 1] - Math.Log(1.00007898996332));
            }
            return acceptibleLogMzDifferences;
        }

        // Finds all groups in logTransformedXArray where the differences between consecutive elements
        // (from high to low) match (within a tolerance) a subsequence of acceptibleLogMzDifferences.
        // Returns a list of (X[], Y[], ChargeState[]) for each group found.
        private static List<(double[] X, double[] Y, int[] ChargeState)> FindMatchingGroups(
            double[] logTransformedXArray,
            double[] yArray)
        {
            var results = new List<(double[], double[], int[])>();
            int logMzArrayLength = logTransformedXArray.Length;

            // Sort peaks by descending m/z (i.e., descending log(m/z))
            var sorted = logTransformedXArray
                .Select((val, idx) => (val, idx))
                .OrderByDescending(x => x.val)
                .ToArray();

            var allExperimentalSortedLogMz = sorted.Select(x => x.val).ToArray();
            var sortedY = sorted.Select(x => yArray[x.idx]).ToArray();

            // Iterate through each peak as a potential group start
            for (int indexOfFirstPeakInPotentialSeries = 0; indexOfFirstPeakInPotentialSeries < logMzArrayLength - 1; indexOfFirstPeakInPotentialSeries++)
            {
                int minChargeState = 1;
                int maxChargeState = 60;
                for (int lowChargeState = minChargeState; lowChargeState < maxChargeState; lowChargeState++)
                {
                    var allPotentialTargetMzsAndCharges = GetTargetLogMzsAndCharges(allExperimentalSortedLogMz[indexOfFirstPeakInPotentialSeries], lowChargeState, maxChargeState);
                    var indiciesOfFoundPotentialTargetMzs = FindMatchingSearchForIndices(allExperimentalSortedLogMz, allPotentialTargetMzsAndCharges.targetLogMzs, LogMzDependentTolerance(allExperimentalSortedLogMz[indexOfFirstPeakInPotentialSeries]));
                    if(indiciesOfFoundPotentialTargetMzs.Any() && indiciesOfFoundPotentialTargetMzs.Length > 1)
                    {
                        var longestSeriesOfConsecutiveChargesWithMinimumGap = LongestSubarrayWithMaxGap(SelectByIndices(allPotentialTargetMzsAndCharges.charges,indiciesOfFoundPotentialTargetMzs), 1);
                        if (longestSeriesOfConsecutiveChargesWithMinimumGap.Length > 1)
                        {
                            var indiciesOfLongestSeriesOfConsecutiveChargesInTargetArray = FindMatchingSearchForIndices(longestSeriesOfConsecutiveChargesWithMinimumGap, allPotentialTargetMzsAndCharges.charges, 0.1);
                            var groupX = FindMatchingValuesFromSearchIn(allExperimentalSortedLogMz,SelectByIndices(allPotentialTargetMzsAndCharges.targetLogMzs,indiciesOfLongestSeriesOfConsecutiveChargesInTargetArray), 0.01);
                            int[] indices = groupX
                                .Select(sf => Array.FindIndex(allExperimentalSortedLogMz, si => Math.Abs(si - sf) <= 0.0001))
                                .Where(idx => idx != -1)
                                .ToArray();
                            var z = RemoveDuplicates(indices);
                            var groupY = SelectByIndices(sortedY,z);
                            var charges = longestSeriesOfConsecutiveChargesWithMinimumGap;
                            results.Add((groupX, groupY, charges));

                        }
                    }
                }
            }
            return results;
        }
        public static int[] RemoveDuplicates(int[] input)
        {
            return input.Distinct().ToArray();
        }
        public static T[] SelectByIndices<T>(T[] source, int[] indices)
        {
            T[] result = new T[indices.Length];
            for (int i = 0; i < indices.Length; i++)
            {
                result[i] = source[indices[i]];
            }
            return result;
        }
        public static int[] LongestSubarrayWithMaxGap(int[] input, int maxDifference)
        {
            if (input == null || input.Length == 0)
                return new int[0];

            List<int> current = new List<int> { input[0] };
            List<int> best = new List<int>(current);

            for (int i = 1; i < input.Length; i++)
            {
                if (Math.Abs(input[i] - input[i - 1]) <= maxDifference)
                {
                    current.Add(input[i]);
                }
                else
                {
                    if (current.Count > best.Count)
                        best = new List<int>(current);
                    current = new List<int> { input[i] };
                }
            }

            if (current.Count > best.Count)
                best = current;

            return best.ToArray();
        }
        public static double[] FindMatchingValuesFromSearchIn(double[] searchIn, double[] searchFor, double tolerance)
        {
            var matchedValues = new List<double>();

            for (int i = 0; i < searchFor.Length; i++)
            {
                double searchForValue = searchFor[i];
                double closestValue = double.NaN;
                double smallestDiff = double.MaxValue;

                for (int j = 0; j < searchIn.Length; j++)
                {
                    double diff = Math.Abs(searchIn[j] - searchForValue);
                    if (diff <= tolerance && diff < smallestDiff)
                    {
                        smallestDiff = diff;
                        closestValue = searchIn[j];
                    }
                }

                if (!double.IsNaN(closestValue))
                {
                    matchedValues.Add(closestValue);
                }
            }

            return matchedValues.ToArray();
        }
        // Finds indices in searchIn that match values in searchFor within a specified tolerance.
        // Returns an array of indices corresponding to searchFor values that were found in searchIn.
        public static int[] FindMatchingSearchForIndices<T>(T[] searchIn, T[] searchFor, double tolerance) where T : IConvertible
        {
            var matchedIndices = new List<int>();

            for (int i = 0; i < searchFor.Length; i++)
            {
                int closestIndex = -1;
                double smallestDiff = double.MaxValue;
                double searchForValue = Convert.ToDouble(searchFor[i]);
                for (int j = 0; j < searchIn.Length; j++)
                {
                    double searchInValue = Convert.ToDouble(searchIn[j]);
                    double diff = Math.Abs(searchInValue - searchForValue);
                    if (diff <= tolerance && diff < smallestDiff)
                    {
                        smallestDiff = diff;
                        closestIndex = i; // store index from searchFor
                    }
                }
                if (closestIndex != -1)
                {
                    matchedIndices.Add(closestIndex);
                }
            }

            return matchedIndices.ToArray();
        }
        public static (double[] targetLogMzs, int[] charges) GetTargetLogMzsAndCharges(double firstMz, int lowestCharge = 1, int highestCharge = 60)
        {
            double[] mzs = new double[highestCharge - lowestCharge + 1];
            int[] charges = Enumerable.Range(lowestCharge, highestCharge - lowestCharge + 1).ToArray();
            var logDiffs = AllAcceptibleLogMzDifferencesForAdjacentValues(lowestCharge, highestCharge);

            double previousMz = firstMz;
            mzs[0] = firstMz;
            for (int i = 1; i < mzs.Length; i++)
            {
                mzs[i] = previousMz - logDiffs[i];
                previousMz -= logDiffs[i];
            }
            return (mzs, charges);
        }
        
        // Removes any group from 'groups' where the double[] X is a subset of any other double[] X
        private static List<(double[] X, double[] Y, int[] charges)> RemoveSubsetGroups(List<(double[] X, double[] Y, int[] charges)> groups)
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
        private static List<(double[] X, double[] Y, int[] chargeStates)> TransformGroupsToExpX(List<(double[] X, double[] Y, int[] charges)> filteredGroups)
        {
            return filteredGroups
                .Select(g => (X: g.X.Select(Math.Exp).ToArray(), Y: g.Y, g.charges))
                .ToList();
        }

        // Filters mass/intensity groups into likely correct and likely incorrect groups based on neutral mass agreement.
        // A group is likely correct if most neutralMass values are within a small ppm tolerance of each other.
        // A group is likely incorrect if most neutralMass values differ by more than a larger ppm tolerance.
        public static void FilterMassIntensityGroupsByPpmTolerance(
            IEnumerable<(double[] neutralMass, double[] intensity)> massIntensityGroups,
            out List<(double[] neutralMass, double[] intensity)> likelyCorrect,
            out List<(double[] neutralMass, double[] intensity)> likelyIncorrect,
            double correctPpmTolerance = 80,
            double incorrectPpmTolerance = 85,
            double correctFraction = 0.7)
        {
            // Initialize output lists
            likelyCorrect = new List<(double[] neutralMass, double[] intensity)>();
            likelyIncorrect = new List<(double[] neutralMass, double[] intensity)>();

            // Process each group
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

                // Compare all pairs of neutralMass values in the group
                for (int i = 0; i < masses.Length; i++)
                {
                    for (int j = i + 1; j < masses.Length; j++)
                    {
                        // Calculate the ppm difference between the two masses
                        double ppm = Math.Abs(masses[i] - masses[j]) / ((masses[i] + masses[j]) / 2.0) * 1e6;
                        if (ppm <= correctPpmTolerance)
                            closeCount++;
                        if (ppm > incorrectPpmTolerance)
                            farCount++;
                        totalPairs++;
                    }
                }

                // If most pairs are close, consider the group likely correct
                if (totalPairs == 0 || (closeCount >= correctFraction * totalPairs && farCount < (1 - correctFraction) * totalPairs))
                    likelyCorrect.Add(group);
                // If most pairs are far apart, consider the group likely incorrect
                else if (farCount > (1 - correctFraction) * totalPairs)
                    likelyIncorrect.Add(group);
                // Otherwise, treat as incorrect (ambiguous case)
                else
                    likelyIncorrect.Add(group);
            }
        }
        // For each group, finds the most common neutral mass (mode) within a specified ppm tolerance
        // and sums the intensities of the peaks that belong to this mode cluster.
        // Returns a list of (mostCommonNeutralMass, summedIntensity) for all groups.
        public static List<(double mostCommonNeutralMass, double summedIntensity)> GetMostCommonNeutralMassAndSummedIntensity(
            IEnumerable<(double[] neutralMass, double[] intensity)> likelyCorrectGroups,
            double ppmTolerance = 2)
        {
            var results = new List<(double mostCommonNeutralMass, double summedIntensity)>();

            // Process each group individually
            foreach (var group in likelyCorrectGroups)
            {
                var masses = group.neutralMass;
                var intensities = group.intensity;

                // Cluster neutral masses: each cluster contains indices of masses within ppmTolerance
                var clusters = new List<List<int>>();

                for (int i = 0; i < masses.Length; i++)
                {
                    bool added = false;
                    // Try to add the mass to an existing cluster
                    for (int c = 0; c < clusters.Count; c++)
                    {
                        double refMass = masses[clusters[c][0]];
                        double ppm = Math.Abs(masses[i] - refMass) / refMass * 1e6;
                        if (ppm <= ppmTolerance)
                        {
                            clusters[c].Add(i);
                            added = true;
                            break;
                        }
                    }
                    // If not close to any cluster, start a new cluster
                    if (!added)
                    {
                        clusters.Add(new List<int> { i });
                    }
                }

                // Find the largest cluster (the mode)
                var modeCluster = clusters.OrderByDescending(cl => cl.Count).First();
                // Calculate the average mass of the mode cluster as the representative value
                double mostCommonNeutralMass = modeCluster.Select(idx => masses[idx]).Average();
                // Sum the intensities of the mode cluster
                double summedIntensity = modeCluster.Select(idx => intensities[idx]).Sum();

                // Add the result for this group
                results.Add((mostCommonNeutralMass, summedIntensity));
            }

            return results.OrderBy(x => x.mostCommonNeutralMass).ToList();
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
                    neutralMass[i] = group.X[i].ToMass(group.ChargeState[i]);
                }
                result.Add((neutralMass, group.Y));
            }

            return result;
        }
        private static double LogMzDependentTolerance(double logMz, double tolerance = 5000)
        {
            var m = Math.Exp(logMz);
            var mPlus = m + m * tolerance / 1000000.0;
            var lmPlus = Math.Log(mPlus);
            var newT = lmPlus - logMz;
            //return newT;
            return 0.0001;
        }
    }
}
