using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    public class PermutationTestingEngine
    {
        private FlashLfqResults QuantResults;
        private const double MAX_Q_VALUE = 0.05;
        //TODO add functionality for user-specified parameters

        public PermutationTestingEngine(FlashLfqResults quantResults)
        {
            bool conditionsAreDefined = quantResults.SpectraFiles.All(p => !string.IsNullOrWhiteSpace(p.Condition));
            if (!conditionsAreDefined)
            {
                throw new MzLibException("Conditions must be defined to run the Bayesian protein quantification");
            }
            QuantResults = quantResults;
        }

        public void Run()
        {
            //get all conditions
            List<string> conditions = QuantResults.SpectraFiles.Select(x => x.Condition).Distinct().ToList();

            //get proteins
            ProteinGroup[] proteins = QuantResults.ProteinGroups.Values.ToArray();

            //iterate through all possible comparisons (i.e. if 3 conditions, do AB, AC, BC)
            for (int conditionOneIndex = 0; conditionOneIndex < conditions.Count - 1; conditionOneIndex++)
            {
                for (int conditionTwoIndex = conditionOneIndex + 1; conditionTwoIndex < conditions.Count; conditionTwoIndex++)
                {
                    //get two conditions to be compared
                    string firstCondition = conditions[conditionOneIndex];
                    string secondCondition = conditions[conditionTwoIndex];

                    //get all files for these conditions
                    List<SpectraFileInfo> conditionOneFiles = QuantResults.SpectraFiles.Where(x => x.Condition.Equals(firstCondition)).ToList();
                    List<SpectraFileInfo> conditionTwoFiles = QuantResults.SpectraFiles.Where(x => x.Condition.Equals(secondCondition)).ToList();
                    List<SpectraFileInfo> allFilesInComparison = conditionOneFiles.Concat(conditionTwoFiles).ToList();
                    int maxNumValidValues = allFilesInComparison.Count;

                    //determine number of valid values for each protein in the comparison
                    int[] numValidValuesPerProtein = new int[proteins.Length];
                    for(int i=0; i<proteins.Length; i++)
                    {
                        ProteinGroup protein = proteins[i];
                        numValidValuesPerProtein[i] = allFilesInComparison.Count(file => protein.GetIntensity(file) > 0);                  
                    }
                    Array.Sort(numValidValuesPerProtein, proteins);

                    //starting values for iteration
                    double sOValue = 0.1;
                    double meanFraction = 0.1;
                    int minNumValidValues = 3;

                    //best set of parameters
                    int maxSignificantCount = 0;
                    int bestMinNumValidValues = 0;
                    double bestMeanFraction = 0;
                    double bestS0Value = 0;

                    //results of best set
                    List<ProteinGroup> proteinsFromBestSetOfParams = new List<ProteinGroup>();
                    List<double> bestNValues = new List<double>();
                    List<double> bestpValues = new List<double>();
                    List<double> bestFoldChanges = new List<double>();


                    ///BEGIN ITERATING OVER ALL SETS OF PARAMETERS
                    int validValueIndex = 0; //used to record the index of proteins with min valid values
                    for (; minNumValidValues <= maxNumValidValues; minNumValidValues++)
                    {
                        //remove proteins which don't meet the required number of valid values
                        for(; validValueIndex < proteins.Length; validValueIndex++)
                        {
                            if(numValidValuesPerProtein[validValueIndex] >= minNumValidValues)
                            {
                                break;
                            }
                        }
                        List<ProteinGroup> proteinsForThisComparison = proteins.SubArray(validValueIndex, proteins.Length - validValueIndex).ToList();

                        //we need to impute missing data. This requires finding the distribution of protein intensities for each file.
                        foreach(ProteinGroup protein in proteinsForThisComparison)
                        {
                            //get the sum of intensity values
                        }

                        for (; meanFraction < 1; meanFraction += 0.3) //FIXME: maybe useless, validate
                        {                                
                            // imputes missing intensity values for each protein
                            ImputationProcess imputationProcess = new ImputationProcess();
                            List<ProteinGroup> imputedProteins = imputationProcess.RunImputationProcess(proteinsForThisComparison, allFilesInComparison, meanFraction);

                            for (; sOValue < 1; sOValue += 0.1)
                            {
                                // Declaring variables which will be generated after T-Tests and Permutation Tests
                                List<double> observedNValues = new List<double>(); // will store observed N values
                                List<double> observedPValues = new List<double>(); // will store observed P values
                                List<double> observedLogFoldChange = new List<double>(); // will store observed Log Fold Change values
                                List<double> permutedNValues = new List<double>(); // will store permuted N values
                                StatisticalTests statisticalTests = new StatisticalTests();

                                // Compute observed and permuted N Values for each protein using T Tests and Permutation Testing
                                for (int i = 0; i < imputedProteins.Count; i++)
                                {
                                    ProteinGroup protein = imputedProteins[i];
                                    List<double> conditionOneIntensities = new List<double>();
                                    List<double> conditionTwoIntensities = new List<double>();

                                    // get the protein's intensity values corresponding to the chosen pair of conditions
                                    foreach(SpectraFileInfo file in conditionOneFiles)
                                    {
                                        conditionOneIntensities.Add(protein.GetIntensity(file));
                                    }    
                                    foreach(SpectraFileInfo file in conditionTwoFiles)
                                    {
                                        conditionTwoIntensities.Add(protein.GetIntensity(file));
                                    }

                                    // Compute observed N Values with the chosen pair of conditions using T-Tests and
                                    // store in observedNValues array
                                    statisticalTests.GetNValueUsingTTest(conditionOneIntensities, conditionTwoIntensities, observedNValues, observedPValues, observedLogFoldChange, sOValue);
                                    

                                    // Compute permuted N Values with the chosen pair of conditions using T-Tests and 
                                    // store in permutedNValues array
                                    statisticalTests.GetNValueUsingPermutationtests(conditionOneIntensities, conditionTwoIntensities, permutedNValues, sOValue);
                                }

                                // Usually more permuted values than observed, so resize the permuted N values list to be the same size as the observed list, emulating a target-decoy comparison
                                ResizePermutedArray(permutedNValues, permutedNValues.Count - observedNValues.Count);

                                // get number of significant proteins
                                int numSigChangers = statisticalTests.GetNumberOfSignificantChangers(observedNValues, permutedNValues, MAX_Q_VALUE);
                                if(numSigChangers>maxSignificantCount)
                                {
                                    maxSignificantCount = numSigChangers;
                                    bestMinNumValidValues = minNumValidValues;
                                    bestMeanFraction = meanFraction;
                                    bestS0Value = sOValue;
                                    proteinsFromBestSetOfParams = imputedProteins;
                                    bestNValues = observedNValues;
                                    bestpValues = observedPValues;
                                    bestFoldChanges = observedLogFoldChange;
                                }
                            }
                        }
                    }

                    //assign the significance results from the best set of parameters
                    //TODO
                }
            }
        }

        /// <summary>
        /// Used to resize the permuted N values containing array so that its size is equal 
        /// to the original N values array size. Sample evenly from the distribution of permuted N-values.
        /// </summary>
        public void ResizePermutedArray(List<double> permutedNValues, int sizeDifference)
        {
            if (sizeDifference != 0)
            {
                //need to evenly remove some values
                permutedNValues.Sort();
                double trackElements = (Convert.ToDouble(permutedNValues.Count) - 1) / sizeDifference;
                double loopMaintainer = trackElements;

                while (trackElements < permutedNValues.Count)
                {
                    int roundedNumber = (int)Math.Round(trackElements);
                    permutedNValues[roundedNumber] = -1;
                    trackElements += loopMaintainer;
                }

                permutedNValues.RemoveAll(value => value == -1);
            }
        }
    }
}
