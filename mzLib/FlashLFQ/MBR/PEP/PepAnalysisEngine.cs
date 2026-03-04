using Microsoft.ML;
using Microsoft.ML.Data;
using Microsoft.ML.Trainers.FastTree;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Omics;
using System.Collections;
using System.Security.Policy;
using System.Text.RegularExpressions;
using System.Reflection;

namespace FlashLFQ.PEP
{
    public class PepAnalysisEngine
    {
        public double PipScoreCutoff;

        private static int _randomSeed = 42;

        /// <summary>
        /// This method contains the hyper-parameters that will be used when training the machine learning model
        /// </summary>
        /// <returns> Options object to be passed in to the FastTree constructor </returns>
        public FastTreeBinaryTrainer.Options BGDTreeOptions =>
            new FastTreeBinaryTrainer.Options
            {
                NumberOfThreads = 1,
                NumberOfTrees = 100,
                MinimumExampleCountPerLeaf = 10,
                NumberOfLeaves = 20,
                LearningRate = 0.2,
                LabelColumnName = "Label",
                FeatureColumnName = "Features",
                Seed = _randomSeed,
                FeatureSelectionSeed = _randomSeed,
                RandomStart = false,
                UnbalancedSets = true
            };
        
        public List<MbrChromatographicPeak> Peaks { get;  }
        public string OutputFolder { get; set; }
        public int MaxThreads { get; set; }
        public double PepTrainingFraction { get; set; }

        public PepAnalysisEngine(List<MbrChromatographicPeak> peaks, string outputFolder, int maxThreads, double pepTrainingFraction = 0.25)
        {
            Peaks = peaks;
            OutputFolder = outputFolder;
            MaxThreads = maxThreads;
            PepTrainingFraction = pepTrainingFraction;
        }

        public string ComputePEPValuesForAllPeaks()
        {
            string[] trainingVariables = ChromatographicPeakData.trainingInfos["standard"];

            #region Construct Donor Groups
            // this is target peak not target peptide
            List<DonorGroup> donors= new();
            foreach(var donorGroup in Peaks
                .Where(peak => peak.DetectionType == DetectionType.MBR)
                .OrderByDescending(peak => peak.MbrScore)
                .GroupBy(peak => peak.Identifications.First())) //Group by donor peptide
            {
                var donorId = donorGroup.Key;
                var targetAcceptors = donorGroup.Where(peak => !peak.RandomRt).ToList();
                var decoyAcceptors = donorGroup.Where(peak => peak.RandomRt).ToList();
                donors.Add(new DonorGroup(donorId, targetAcceptors, decoyAcceptors));
            }

            // Fix the order
            donors = OrderDonorGroups(donors);

            var peakScores = donors.SelectMany(donor => donor.Select(p => p.MbrScore)).OrderByDescending(score => score).ToList();
            PipScoreCutoff = peakScores[(int)Math.Floor(peakScores.Count * PepTrainingFraction)]; //Select the top N percent of all peaks, only use those as positive examples

            MLContext mlContext = new MLContext(_randomSeed);
            //the number of groups used for cross-validation is hard-coded at three. Do not change this number without changing other areas of effected code.
            const int numGroups = 3;

            List<int>[] donorGroupIndices = GetDonorGroupIndices(donors, numGroups, PipScoreCutoff);

            #endregion

            #region Create Groups and Model
            IEnumerable<ChromatographicPeakData>[] ChromatographicPeakDataGroups = new IEnumerable<ChromatographicPeakData>[numGroups];
            for (int i = 0; i < numGroups; i++)
            {
                ChromatographicPeakDataGroups[i] = CreateChromatographicPeakData(donors, donorGroupIndices[i], MaxThreads);

                if (!ChromatographicPeakDataGroups[i].Any(p => p.Label == true) 
                    || !ChromatographicPeakDataGroups[i].Any(p => p.Label == false))
                {
                    return "Posterior error probability analysis failed. This can occur for small data sets when some sample groups are missing positive or negative training examples.";
                }
            }

            TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[] trainedModels 
                = new TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[numGroups];

            var trainer = mlContext.BinaryClassification.Trainers.FastTree(BGDTreeOptions);
            var pipeline = mlContext.Transforms.Concatenate("Features", trainingVariables)
                .Append(trainer);

            List<CalibratedBinaryClassificationMetrics> allMetrics = new List<CalibratedBinaryClassificationMetrics>();

            #endregion

            #region Training and Cross Validation First iteration

            for (int groupIndexNumber = 0; groupIndexNumber < numGroups; groupIndexNumber++)
            {

                List<int> allGroupIndexes = Enumerable.Range(0, numGroups).ToList();
                allGroupIndexes.RemoveAt(groupIndexNumber);

                //concat doesn't work in a loop, therefore I had to hard code the concat to group 3 out of 4 lists. if the const int numGroups value is changed, then the concat has to be changed accordingly.
                IDataView dataView = mlContext.Data.LoadFromEnumerable(
                    ChromatographicPeakDataGroups[allGroupIndexes[0]]
                    .Concat(ChromatographicPeakDataGroups[allGroupIndexes[1]]));

                trainedModels[groupIndexNumber] = pipeline.Fit(dataView);
                var myPredictions = trainedModels[groupIndexNumber].Transform(mlContext.Data.LoadFromEnumerable(ChromatographicPeakDataGroups[groupIndexNumber]));
                CalibratedBinaryClassificationMetrics metrics = mlContext.BinaryClassification.Evaluate(data: myPredictions, labelColumnName: "Label", scoreColumnName: "Score");

                //Parallel operation of the following code requires the method to be stored and then read, once for each thread
                //if not output directory is specified, the model cannot be stored, and we must force single-threaded operation
                if (OutputFolder != null)
                {
                    mlContext.Model.Save(trainedModels[groupIndexNumber], dataView.Schema, Path.Combine(OutputFolder, "model.zip"));
                }

                Compute_PEP_For_All_Peaks(donors, donorGroupIndices[groupIndexNumber], mlContext, trainedModels[groupIndexNumber], OutputFolder, MaxThreads);

                allMetrics.Add(metrics);
            }

            #endregion
            #region Iterative Training

            for(int trainingIteration = 0; trainingIteration < 9; trainingIteration++)
            {
                ChromatographicPeakDataGroups = new IEnumerable<ChromatographicPeakData>[numGroups];
                for (int i = 0; i < numGroups; i++)
                {
                    ChromatographicPeakDataGroups[i] = CreateChromatographicPeakDataIteration(donors, donorGroupIndices[i], MaxThreads);

                    if (!ChromatographicPeakDataGroups[i].Any(p => p.Label == true)
                        || !ChromatographicPeakDataGroups[i].Any(p => p.Label == false))
                    {
                        return "Posterior error probability analysis failed. This can occur for small data sets when some sample groups are missing positive or negative training examples.";
                    }
                }

                for (int groupIndexNumber = 0; groupIndexNumber < numGroups; groupIndexNumber++)
                {
                    List<int> allGroupIndexes = Enumerable.Range(0, numGroups).ToList();
                    allGroupIndexes.RemoveAt(groupIndexNumber);

                    IDataView dataView = mlContext.Data.LoadFromEnumerable(
                        ChromatographicPeakDataGroups[allGroupIndexes[0]]
                        .Concat(ChromatographicPeakDataGroups[allGroupIndexes[1]]));

                    trainedModels[groupIndexNumber] = pipeline.Fit(dataView);
                    var myPredictions = trainedModels[groupIndexNumber].Transform(mlContext.Data.LoadFromEnumerable(ChromatographicPeakDataGroups[groupIndexNumber]));
                    CalibratedBinaryClassificationMetrics metrics = mlContext.BinaryClassification.Evaluate(data: myPredictions, labelColumnName: "Label", scoreColumnName: "Score");

                    //Parallel operation of the following code requires the method to be stored and then read, once for each thread
                    //if not output directory is specified, the model cannot be stored, and we must force single-threaded operation
                    if (OutputFolder != null)
                    {
                        mlContext.Model.Save(trainedModels[groupIndexNumber], dataView.Schema, Path.Combine(OutputFolder, "model.zip"));
                    }

                    Compute_PEP_For_All_Peaks(donors, donorGroupIndices[groupIndexNumber], mlContext, trainedModels[groupIndexNumber], OutputFolder, MaxThreads);

                    allMetrics.Add(metrics);
                }
            }
            #endregion

            return AggregateMetricsForOutput(allMetrics);
        }

        public static List<DonorGroup> OrderDonorGroups(List<DonorGroup> donors)
        {
            return donors.OrderByDescending(donor => donor.TargetAcceptors.Count)
                .ThenByDescending(donor => donor.DecoyAcceptors.Count)
                .ThenByDescending(donor => donor.BestTargetMbrScore)
                .ToList();
        }

        //we add the indexes of the targets and decoys to the groups separately in the hope that we'll get at least one target and one decoy in each group.
        //then training can possibly be more successful.
        public static List<int>[] GetDonorGroupIndices(List<DonorGroup> donors, int numGroups, double scoreCutoff)
        {
            List<int>[] groupsOfIndices = new List<int>[numGroups];
            for (int i = 0; i < numGroups; i++)
            {
                groupsOfIndices[i] = new List<int>();
            }

            int myIndex = 0;

            while (myIndex < donors.Count)
            {
                int subIndex = 0;
                while (subIndex < numGroups && myIndex < donors.Count)
                {
                    groupsOfIndices[subIndex].Add(myIndex);

                    subIndex++;
                    myIndex++;
                }
            }

            EqualizeDonorGroupIndices(donors, groupsOfIndices, scoreCutoff, numGroups);
            
            return groupsOfIndices;
        }

        /// <summary>
        /// Equalizes partitions used for cross-validation. The goal is to have the same number of targets and decoys in each partition
        /// </summary>
        /// <param name="donors"> List of all DonorGroups to be classified</param>
        /// <param name="groupsOfIndices">An array of lists. Each list contains the indices of donor groups for a given partition </param>
        /// <param name="scoreCutoff"> The MBR Score cutoff that determines which MBR target peaks will be used as positive training examples</param>
        /// /// <param name="numGroups">Number of groups used for cross-validation, default = 3 </param>
        public static void EqualizeDonorGroupIndices(List<DonorGroup> donors, List<int>[] groupsOfIndices, double scoreCutoff, int numGroups = 3)
        {
            HashSet<int> swappedDonors = new HashSet<int>(); // Keep track of everything we've swapped so we don't swap it again
            // Outer loop iterates over the groups of indices (partitions) three times
            // after each inner loop iterations, the number of ttargtes and decoys in each adjacent group is equal, but commonly group 1 and 3 will have a different number
            // of targets and decoys. Looping three times should resolve this
            for (int i = 0; i < numGroups*3 - 1; i++)
            {
                int groupA = i % numGroups;
                int groupB = (i + 1) % numGroups;
                int targetsA = 0;
                int targetsB = 0;
                int decoysA = 0;
                int decoysB = 0;
                foreach (int index in groupsOfIndices[groupA])
                {
                    targetsA += donors[index].TargetAcceptors.Count(peak => peak.MbrScore >= scoreCutoff);
                    decoysA += donors[index].DecoyAcceptors.Count;
                }
                foreach (int index in groupsOfIndices[groupB])
                {
                    targetsB += donors[index].TargetAcceptors.Count(peak => peak.MbrScore >= scoreCutoff);
                    decoysB += donors[index].DecoyAcceptors.Count;
                }

                bool stuck = false;
                int outerIterations = 0;
                int minIndex = groupsOfIndices[groupA].Min();

                // Calculate the difference in targets and decoys between the two groups
                int targetSurplus = targetsA - targetsB;
                int decoySurplus = decoysA - decoysB;

                while ((Math.Abs(targetSurplus) > 1 | Math.Abs(decoySurplus) > 1) && !stuck && outerIterations < 3)
                {
                    bool swapped = false;
                    outerIterations++;

                    int innerIterations = 0;
                    // start from the bottom of group 1, trying to swap peaks.
                    // If group 1 has more targets than group 2, we want to swap groups to equalize the number of targets in each group
                    while (Math.Abs(targetSurplus) > 1 & !stuck & innerIterations < 3)
                    {
                        innerIterations++;
                        swapped = false;
                        // Traverse the list of donor indices in descending order, looking for a good candidate to swap
                        foreach (int donorIndexA in groupsOfIndices[groupA].Where(idx => !swappedDonors.Contains(idx)).OrderByDescending(idx => idx))
                        {
                            int donorIndexATargetCount = donors[donorIndexA].TargetAcceptors.Count(peak => peak.MbrScore > scoreCutoff);
                            switch (targetSurplus > 0)
                            {
                                case true: // i.e., too many targets
                                    if (donorIndexATargetCount < 1) continue; // No targets to swap
                                    foreach (int donorIndexB in groupsOfIndices[groupB].Where(idx => !swappedDonors.Contains(idx)).OrderByDescending(idx => idx))
                                    {
                                        if (donors[donorIndexB].TargetAcceptors.Count(peak => peak.MbrScore > scoreCutoff) < donorIndexATargetCount)
                                        {
                                            GroupSwap(donors, groupsOfIndices, donorIndexA, donorIndexB, groupA, groupB,
                                                scoreCutoff, swappedDonors, ref targetSurplus, ref decoySurplus);
                                            swapped = true;
                                            break;
                                        }
                                    }
                                    break;
                                case false: // i.e., too few targets
                                    foreach (int donorIndexB in groupsOfIndices[groupB].Where(idx => !swappedDonors.Contains(idx)).OrderByDescending(idx => idx))
                                    {
                                        if (donors[donorIndexB].TargetAcceptors.Count(peak => peak.MbrScore > scoreCutoff) > donorIndexATargetCount)
                                        {
                                            GroupSwap(donors, groupsOfIndices, donorIndexA, donorIndexB, groupA, groupB,
                                                scoreCutoff, swappedDonors, ref targetSurplus, ref decoySurplus);
                                            swapped = true;
                                            break;
                                        }
                                    }
                                    break;
                            }

                            // If we reach the index of the list of donorGroups, set stuck to true so that the outer loop will break
                            if (donorIndexA == minIndex)
                            {
                                stuck = true;
                                break;
                            }
                            if (swapped)
                                break;

                        }
                    }

                    innerIterations = 0;
                    // Now we'll do the decoys
                    while (Math.Abs(decoySurplus) > 1 & !stuck & innerIterations < 3)
                    {
                        innerIterations++;
                        swapped = false;
                        foreach (int donorIndexA in groupsOfIndices[groupA].Where(idx => !swappedDonors.Contains(idx)).OrderByDescending(idx => idx))
                        {
                            int donorIndexADecoyCount = donors[donorIndexA].DecoyAcceptors.Count();
                            switch (decoySurplus > 0)
                            {
                                case true: // i.e., too many decoys
                                    if (donorIndexADecoyCount < 1) continue; // No decoys to swap
                                    foreach (int donorIndexB in groupsOfIndices[groupB].Where(idx => !swappedDonors.Contains(idx)).OrderByDescending(idx => idx))
                                    {
                                        if (donors[donorIndexB].DecoyAcceptors.Count() < donorIndexADecoyCount)
                                        {
                                            GroupSwap(donors, groupsOfIndices, donorIndexA, donorIndexB, groupA, groupB,
                                                scoreCutoff, swappedDonors, ref targetSurplus, ref decoySurplus);
                                            swapped = true;
                                            break;
                                        }
                                    }
                                    break;
                                case false: // i.e., too few decoys
                                    foreach (int donorIndexB in groupsOfIndices[groupB].Where(idx => !swappedDonors.Contains(idx)).OrderByDescending(idx => idx))
                                    {
                                        if (donors[donorIndexB].DecoyAcceptors.Count() > donorIndexADecoyCount)
                                        {
                                            GroupSwap(donors, groupsOfIndices, donorIndexA, donorIndexB, groupA, groupB,
                                                scoreCutoff, swappedDonors, ref targetSurplus, ref decoySurplus);
                                            swapped = true;
                                            break;
                                        }
                                    }
                                    break;
                            }

                            // If we reach the index of the list of donorGroups, set stuck to true so that the outer loop will break
                            if (donorIndexA == minIndex)
                            {
                                stuck = true;
                                break;
                            }
                            if (swapped)
                                break;
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Takes in a list of donor groups and a list of indices for each group, and swaps two groups of indices
        /// Updates the targetSurplus and decoySurplus variables
        /// Updates the swappedDonors hash set to keep track of which donors have been swapped
        /// This is done to equalize the number of targets and decoys in each paritition for cross validation
        /// </summary>
        public static void GroupSwap(
            List<DonorGroup> donors, 
            List<int>[] groupsOfIndices, 
            int donorIndexA, 
            int donorIndexB,
            int groupsOfIndicesIndexA,
            int groupsOfIndicesIndexB,
            double scoreCutoff,
            HashSet<int> swappedDonors,
            ref int targetSurplus,
            ref int decoySurplus)
        {
            // Multiply by two because the surplus is the difference between the two groups
            // So removing one peak from one group and adding it to the other group is a difference of two
            targetSurplus += 2 * (
                donors[donorIndexB].TargetAcceptors.Count(peak => peak.MbrScore >= scoreCutoff) -
                donors[donorIndexA].TargetAcceptors.Count(peak => peak.MbrScore >= scoreCutoff));
            decoySurplus += 2 * (
                donors[donorIndexB].DecoyAcceptors.Count -
                donors[donorIndexA].DecoyAcceptors.Count);

            groupsOfIndices[groupsOfIndicesIndexA].Add(donorIndexB);
            groupsOfIndices[groupsOfIndicesIndexA].Remove(donorIndexA);
            groupsOfIndices[groupsOfIndicesIndexB].Add(donorIndexA);
            groupsOfIndices[groupsOfIndicesIndexB].Remove(donorIndexB);
        }

        /// <summary>
        /// Creates chromatographic peak data that will be used to train the machine learning model
        /// Classifies peaks as positive or negative training examples
        /// Positive training examples are peaks with MBR scores above the 25th percentile,
        /// Negative training examples are peaks with random retention times
        /// </summary>
        /// <param name="donors">The list of donor groups.</param>
        /// <param name="donorIndices">The list of donor indices.</param>
        /// <param name="maxThreads">The maximum number of threads.</param>
        /// <returns>The enumerable of chromatographic peak data.</returns>
        public IEnumerable<ChromatographicPeakData> CreateChromatographicPeakData(List<DonorGroup> donors, List<int> donorIndices, int maxThreads)
        {
            object ChromatographicPeakDataListLock = new object();
            List<ChromatographicPeakData> ChromatographicPeakDataList = new List<ChromatographicPeakData>();
            int[] threads = Enumerable.Range(0, maxThreads).ToArray();

            List<double> pipScores = new();
            foreach(int i in donorIndices)
            {
                pipScores.AddRange(donors[i].Select(peak => peak.MbrScore));
            }
            pipScores.Sort((a, b) => b.CompareTo(a)); // This is a descending sort
            double groupSpecificPipScoreCutoff = pipScores[(int)Math.Floor(pipScores.Count * 0.25)];

            Parallel.ForEach(Partitioner.Create(0, donorIndices.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    List<ChromatographicPeakData> localChromatographicPeakDataList = new List<ChromatographicPeakData>();
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        var donor = donors[donorIndices[i]];
                        foreach (var peak in donor)
                        {
                            ChromatographicPeakData newChromatographicPeakData = new ChromatographicPeakData();
                            if (peak.RandomRt)
                            {
                                newChromatographicPeakData = CreateOneChromatographicPeakDataEntry(peak, label: false);
                                localChromatographicPeakDataList.Add(newChromatographicPeakData);
                            }
                            else if (!peak.RandomRt & peak.MbrScore >= groupSpecificPipScoreCutoff)
                            {
                                newChromatographicPeakData = CreateOneChromatographicPeakDataEntry(peak, label: true);
                                localChromatographicPeakDataList.Add(newChromatographicPeakData);
                            }
                        }
                    }
                    lock (ChromatographicPeakDataListLock)
                    {
                        ChromatographicPeakDataList.AddRange(localChromatographicPeakDataList);
                    }
                });

            ChromatographicPeakData[] pda = ChromatographicPeakDataList.ToArray();

            return pda.AsEnumerable();
        }

        /// <summary>
        /// Creates chromatographic peak data, but uses PEP values instead of MBR scores to select the positive training examples
        /// </summary>
        /// <param name="donors">The list of donor groups.</param>
        /// <param name="donorIndices">The list of donor indices.</param>
        /// <param name="maxThreads">The maximum number of threads.</param>
        /// <returns>The enumerable of chromatographic peak data.</returns>
        public IEnumerable<ChromatographicPeakData> CreateChromatographicPeakDataIteration(List<DonorGroup> donors, List<int> donorIndices, int maxThreads)
        {
            object ChromatographicPeakDataListLock = new object();
            List<ChromatographicPeakData> ChromatographicPeakDataList = new List<ChromatographicPeakData>();
            int[] threads = Enumerable.Range(0, maxThreads).ToArray();

            List<double> peps = new();
            foreach (int i in donorIndices)
            {
                peps.AddRange(donors[i].Select(peak => peak.MbrPep ?? 1));
            }
            peps.Sort();
            double groupSpecificPepCutoff = peps[(int)Math.Floor(peps.Count * 0.25)];

            Parallel.ForEach(Partitioner.Create(0, donorIndices.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    List<ChromatographicPeakData> localChromatographicPeakDataList = new List<ChromatographicPeakData>();
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        var donor = donors[donorIndices[i]];
                        foreach (var peak in donor)
                        {
                            ChromatographicPeakData newChromatographicPeakData = new ChromatographicPeakData();
                            if (peak.RandomRt)
                            {
                                newChromatographicPeakData = CreateOneChromatographicPeakDataEntry(peak, label: false);
                                localChromatographicPeakDataList.Add(newChromatographicPeakData);
                            }
                            else if (!peak.RandomRt & peak.MbrPep <= groupSpecificPepCutoff)
                            {
                                newChromatographicPeakData = CreateOneChromatographicPeakDataEntry(peak, label: true);
                                localChromatographicPeakDataList.Add(newChromatographicPeakData);
                            }
                        }
                    }
                    lock (ChromatographicPeakDataListLock)
                    {
                        ChromatographicPeakDataList.AddRange(localChromatographicPeakDataList);
                    }
                });

            ChromatographicPeakData[] pda = ChromatographicPeakDataList.ToArray();

            return pda.AsEnumerable();
        }

        public static void Compute_PEP_For_All_Peaks(
            List<DonorGroup> donors,
            List<int> donorIndices,
            MLContext mLContext,
            TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>> trainedModel,
            string outputFolder, int maxThreads)
        {
            object lockObject = new object();

            //the trained model is not threadsafe. Therefore, to use the same model for each thread saved the model to disk. Then each thread reads its own copy of the model back from disk.
            //If there is no output folder specified, then this can't happen. We set maxthreads eqaul to one and use the model that gets passed into the method.
            if (String.IsNullOrEmpty(outputFolder))
            {
                maxThreads = 1;
            }

            Parallel.ForEach(Partitioner.Create(0, donorIndices.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {

                    ITransformer threadSpecificTrainedModel;
                    if (maxThreads == 1)
                    {
                        threadSpecificTrainedModel = trainedModel;
                    }
                    else
                    {
                        threadSpecificTrainedModel = mLContext.Model.Load(Path.Combine(outputFolder, "model.zip"), out DataViewSchema savedModelSchema);
                    }

                    // one prediction engine per thread, because the prediction engine is not thread-safe
                    var threadPredictionEngine = mLContext.Model.CreatePredictionEngine<ChromatographicPeakData, TruePositivePrediction>(threadSpecificTrainedModel);

                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        DonorGroup donor = donors[donorIndices[i]];

                        foreach(MbrChromatographicPeak peak in donor)
                        {
                            ChromatographicPeakData pd = CreateOneChromatographicPeakDataEntry(peak, label: !peak.RandomRt);
                            var pepValuePrediction = threadPredictionEngine.Predict(pd);
                            peak.MbrPep = 1 - pepValuePrediction.Probability;
                        }
                    }
                });
        }

        public static string AggregateMetricsForOutput(List<CalibratedBinaryClassificationMetrics> allMetrics)
        {
            List<double> accuracy = allMetrics.Select(m => m.Accuracy).ToList();
            List<double> areaUnderRocCurve = allMetrics.Select(m => m.AreaUnderRocCurve).ToList();
            List<double> areaUnderPrecisionRecallCurve = allMetrics.Select(m => m.AreaUnderPrecisionRecallCurve).ToList();
            List<double> F1Score = allMetrics.Select(m => m.F1Score).ToList();
            List<double> logLoss = allMetrics.Select(m => m.LogLoss).ToList();
            List<double> logLossReduction = allMetrics.Select(m => m.LogLossReduction).ToList();
            List<double> positivePrecision = allMetrics.Select(m => m.PositivePrecision).ToList();
            List<double> positiveRecall = allMetrics.Select(m => m.PositiveRecall).ToList();
            List<double> negativePrecision = allMetrics.Select(m => m.NegativePrecision).ToList();
            List<double> negativeRecall = allMetrics.Select(m => m.NegativeRecall).ToList();

            // log-loss can stochastically take on a value of infinity.
            // correspondingly, log-loss reduction can be negative infinity.
            // when this happens for one or more of the metrics, it can lead to uninformative numbers.
            // so, unless they are all infinite, we remove them from the average. If they are all infinite, we report that.

            logLoss.RemoveAll(x => x == Double.PositiveInfinity);
            logLossReduction.RemoveAll(x => x == Double.NegativeInfinity);

            double logLossAverage = Double.PositiveInfinity;
            double logLossReductionAverage = Double.NegativeInfinity;

            if ((logLoss != null) && (logLoss.Any()))
            {
                logLossAverage = logLoss.Average();
            }

            if ((logLossReduction != null) && (logLossReduction.Any()))
            {
                logLossReductionAverage = logLossReduction.Average();
            }

            StringBuilder s = new StringBuilder();
            s.AppendLine();
            s.AppendLine("************************************************************");
            s.AppendLine("*       Metrics for Determination of PEP Using Binary Classification      ");
            s.AppendLine("*-----------------------------------------------------------");
            s.AppendLine("*       Accuracy:  " + accuracy.Average().ToString());
            s.AppendLine("*       Area Under Curve:  " + areaUnderRocCurve.Average().ToString());
            s.AppendLine("*       Area under Precision recall Curve:  " + areaUnderPrecisionRecallCurve.Average().ToString());
            s.AppendLine("*       F1Score:  " + F1Score.Average().ToString());
            s.AppendLine("*       LogLoss:  " + logLossAverage.ToString());
            s.AppendLine("*       LogLossReduction:  " + logLossReductionAverage.ToString());
            s.AppendLine("*       PositivePrecision:  " + positivePrecision.Average().ToString());
            s.AppendLine("*       PositiveRecall:  " + positiveRecall.Average().ToString());
            s.AppendLine("*       NegativePrecision:  " + negativePrecision.Average().ToString());
            s.AppendLine("*       NegativeRecall:  " + negativeRecall.Average().ToString());
            s.AppendLine("************************************************************");
            return s.ToString();
        }

        public static ChromatographicPeakData CreateOneChromatographicPeakDataEntry(MbrChromatographicPeak peak,bool label)
        {

            peak.PepPeakData = new ChromatographicPeakData
            {
                PpmErrorScore = (float)peak.PpmScore,
                IntensityScore = (float)peak.IntensityScore,
                RtScore = (float)peak.RtScore,
                ScanCountScore = (float)peak.ScanCountScore,
                IsotopicDistributionScore = (float)peak.IsotopicDistributionScore,

                PpmErrorRaw = (float)Math.Abs(peak.MassError),
                IntensityRaw = (float)Math.Log2(peak.Intensity),
                RtPredictionErrorRaw = (float)Math.Abs(peak.RtPredictionError),
                ScanCountRaw = (float)peak.IsotopicEnvelopes.Count,
                IsotopicPearsonCorrelation = (float)(peak.IsotopicPearsonCorrelation),

                Label = label,

            };

            return peak.PepPeakData;
        }
    }
}