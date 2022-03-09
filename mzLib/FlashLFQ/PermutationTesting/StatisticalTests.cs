using Accord.Statistics.Testing;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    /// <summary>
    /// This class contains the statsitical tests used in the program
    /// </summary>
    public class StatisticalTests
    {
        /// <summary>
        /// This method is used to generate all possible combinations of choosing 2 indicies from
        /// the given indicies
        /// </summary>
        public List<List<int>> GenerateAllCombinationsOfTwoIndices(List<double> indices)
        {
            List<List<int>> allTwoIndicesCombination = new List<List<int>>();
            for (int i = 0; i < indices.Count; i++)
            {
                for (int j = i + 1; j < indices.Count; j++)
                {
                    allTwoIndicesCombination.Add(new List<int> { i, j });
                }
            }
            return allTwoIndicesCombination;
        }

        /// <summary>
        /// Calculates the standard deviation of intensity values of the protein
        /// </summary>
        public double CalculateProteinIntensityValuesStandardDeviation(List<double> intensityValues, double intensityValuesMean)
        {
            double intensityValuesStdDev = 0;

            foreach (double intensityVal in intensityValues)
            {
                intensityValuesStdDev += (intensityVal - intensityValuesMean) * (intensityVal - intensityValuesMean);
            }

            intensityValuesStdDev /= intensityValues.Count;
            intensityValuesStdDev = Math.Sqrt(intensityValuesStdDev);

            return intensityValuesStdDev;
        }


        /// <summary>
        /// Calculates the N Value for each protein based on its intensity values in two conditions.
        /// the permutationTesting boolean is used to determine if we are caluclating permuted or observed N Values.
        /// 
        /// The N value is a metric which combines the pValue and the magnitude change (through the LogFoldChange) of
        /// the protein amongst the two conditions in order to enable significance classifcation of the protein based
        /// on the same principle of Volcano Plots (Volcano plots allow to identify signifcantly different proteins in large data 
        /// sets comprising of replicate data by plotting the signifance statistic on the y axis, which is the pValue here,
        /// and the magnitude chanbge on the x axis, which is the logFoldChange here)
        /// </summary>
        public void GetNValueUsingTTest(List<double> proteinFirstConditionIntensityValues, List<double> proteinSecondConditionIntensityValues,
            List<double> nValues, List<double> pValues, List<double> logFoldChange, double sOValue)
        {
            double mean1 = proteinFirstConditionIntensityValues.Average();
            double mean2 = proteinSecondConditionIntensityValues.Average();
            double stdDev1 = CalculateProteinIntensityValuesStandardDeviation(proteinFirstConditionIntensityValues, mean1);
            double stdDev2 = CalculateProteinIntensityValuesStandardDeviation(proteinSecondConditionIntensityValues, mean2);
            double variance1 = stdDev1 * stdDev1;
            double variance2 = stdDev2 * stdDev2;

            // the f-test is used to determine whether the variance of the two intensityvalue popluations are equal
            // the significant variable gets whether the null hypothesis can be rejected
            bool significant = new FTest(variance1, variance2, proteinFirstConditionIntensityValues.Count - 1, proteinSecondConditionIntensityValues.Count - 1).Significant;

            // Create two tailed t test to get p values
            TwoSampleTTest ttest = new TwoSampleTTest(mean1, variance1, proteinFirstConditionIntensityValues.Count,
                mean2, variance2, proteinSecondConditionIntensityValues.Count, !significant);

            double pValue = ttest.PValue;
            double logpValue = -Math.Log10(pValue);
            double logfoldChange = mean2 - mean1;

            // compute N-Value for protein
            double nValue = (logpValue * (logfoldChange * logfoldChange - sOValue * sOValue)) / ((logfoldChange) * (logfoldChange));

                nValues.Add(nValue);
                logFoldChange.Add(logfoldChange);
                pValues.Add(pValue);          
        }

        /// <summary>
        /// Generates permuted N Value's for each protein based on its intensity values in two conditions using Permutation Testing
        /// </summary>
        public void GetNValueUsingPermutationtests(List<double> proteinFirstConditionIntensityValues, List<double> proteinSecondConditionIntensityValues,
             List<double> permutedNValues, double sOValue)
        {
            List<List<int>> allTwoIndiciesCombinationsFromFirstCondition = GenerateAllCombinationsOfTwoIndices(proteinFirstConditionIntensityValues);
            List<List<int>> allTwoIndiciesCombinationsFromSecondCondition = GenerateAllCombinationsOfTwoIndices(proteinSecondConditionIntensityValues);

            int count = 0;
            foreach (var twoIndiciesCombinationEntryFromFirstCondition in allTwoIndiciesCombinationsFromFirstCondition)
            {
                foreach (var twoIndiciesCombinationEntryFromSecondCondition in allTwoIndiciesCombinationsFromSecondCondition)
                {
                    // these are the new arrays which will be made after swapping intensity values between the two conditions
                    List<double> swappedFirstConditionIntensityValues = new List<double>();
                    List<double> swappedSecondConditionIntensityValues = new List<double>();

                    int[] indiciesToSwapFromFirstCondition = new int[2];
                    int[] indiciesToSwapFromSecondCondition = new int[2];
                    int removeIndiciesFirstConditionTracker = 0;
                    int removeIndiciesSecondConditionTracker = 0;

                    // store the indices, corresponding to intensity values, to be swapped from first condition
                    foreach (var index in twoIndiciesCombinationEntryFromFirstCondition)
                    {
                        indiciesToSwapFromFirstCondition[removeIndiciesFirstConditionTracker] = index;
                        removeIndiciesFirstConditionTracker++;
                    }

                    // store the indices, corresponding to intensity values, to be swapped from second condition
                    foreach (var index in twoIndiciesCombinationEntryFromSecondCondition)
                    {
                        indiciesToSwapFromSecondCondition[removeIndiciesSecondConditionTracker] = index;
                        removeIndiciesSecondConditionTracker++;
                    }

                    // add the intensity values to be swapped from first condition into the second condition
                    for (int j = 0; j < indiciesToSwapFromFirstCondition.Count(); j++)
                    {
                        swappedSecondConditionIntensityValues.Add(proteinFirstConditionIntensityValues[indiciesToSwapFromFirstCondition[j]]);
                    }

                    // add the intensity values to be swapped from second condition into the first condition
                    for (int j = 0; j < indiciesToSwapFromSecondCondition.Count(); j++)
                    {
                        swappedFirstConditionIntensityValues.Add(proteinSecondConditionIntensityValues[indiciesToSwapFromSecondCondition[j]]);
                    }

                    // now we add the remaining intensity values from the first condition to the swappedFirstCondition Array
                    for (int j = 0; j < proteinFirstConditionIntensityValues.Count(); j++)
                    {
                        if (indiciesToSwapFromFirstCondition.Contains(j)) continue;
                        swappedFirstConditionIntensityValues.Add(proteinFirstConditionIntensityValues[j]);
                    }

                    // now we add the remaining intensity values from the second condition to the swappedSecondCondition Array
                    for (int j = 0; j < proteinSecondConditionIntensityValues.Count(); j++)
                    {
                        if (indiciesToSwapFromSecondCondition.Contains(j)) continue;
                        swappedSecondConditionIntensityValues.Add(proteinSecondConditionIntensityValues[j]);
                    }

                    // at this stage we have the newly made swapped arrays with mixture of groups.
                    // need to proceed with T tests for these groups to generate permuted p values.
                    GetNValueUsingTTest(swappedFirstConditionIntensityValues, swappedSecondConditionIntensityValues, permutedNValues,
                        null, null, sOValue, true);
                }
                count++;
                if (count == 2) break;
            }
        }

        /// <summary>
        /// Used to determine the N value threhold for FDR purposes. It does so combining the permuted N values 
        /// and observed N values in descending order and determines when (permuted N values)/(observed N Values) > Target FDR.
        /// </summary>
        public int GetNumberOfSignificantChangers(List<double> observedNValues, List<double> permutedNValues, double FDR)
        {
            //sort n-values
            List<double> sortedObservedNValues = observedNValues.OrderByDescending(x => x).ToList();
            permutedNValues = permutedNValues.OrderByDescending(x => x).ToList();

            int numObserved = 0;
            int numPermuted = 0;
            int observedIndex = 0;
            int permutedIndex = 0;

            while (observedIndex < sortedObservedNValues.Count || permutedIndex < permutedNValues.Count)
            {
                if (observedIndex < sortedObservedNValues.Count && permutedIndex < permutedNValues.Count)
                {
                    if (sortedObservedNValues[observedIndex] > permutedNValues[permutedIndex])
                    {
                        numObserved++;
                        observedIndex++;
                    }
                    else if (sortedObservedNValues[observedIndex] < permutedNValues[permutedIndex])
                    {
                        numPermuted++;
                        permutedIndex++;
                        if (numObserved == 0 || ((numPermuted*1d / numObserved) > FDR))
                        {
                            return numObserved;
                        }
                    }
                    else
                    {
                        numObserved++;
                        numPermuted++;
                        if ((numPermuted*1d / numObserved) > FDR)
                        {
                            return numObserved;
                        }
                        observedIndex++;
                        permutedIndex++;
                    }
                }
                else if (permutedIndex < permutedNValues.Count)
                {
                    numPermuted++;
                    permutedIndex++;
                    if (numObserved == 0 || (numPermuted*1d / numObserved) > FDR)
                    {
                        return numObserved;
                    }
                }
                else //if (observedIndex < sortedObservedNValues.Count)
                {
                    return sortedObservedNValues.Count;
                }
            }

            return numObserved;
        }
    }
}