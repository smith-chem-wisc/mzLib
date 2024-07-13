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

namespace FlashLFQ.PEP
{
    public static class PEP_Analysis_Cross_Validation
    {
        public static double PipScoreCutoff;

        public static string ComputePEPValuesForAllPeaks(ref List<ChromatographicPeak> peaks, string outputFolder, int maxThreads, double pepTrainingFraction = 0.25)
        {
            string[] trainingVariables = ChromatographicPeakData.trainingInfos["standard"];

            // Combine the groups and ensure that the order is always stable.
            peaks = peaks.OrderBy(p => p.GetHashCode()).ToList();

            var peakScores = peaks.Select(peak => peak.MbrScore).OrderByDescending(score => score).ToList();
            PipScoreCutoff = peakScores[(int)Math.Floor(peakScores.Count * pepTrainingFraction)]; //Select the top N percent of all peaks, only use those as positive examples

            MLContext mlContext = new MLContext();
            //the number of groups used for cross-validation is hard-coded at four. Do not change this number without changing other areas of effected code.
            const int numGroups = 3;

            List<int>[] peakGroupIndices = Get_Peak_Group_Indices(peaks, numGroups);
            IEnumerable<ChromatographicPeakData>[] ChromatographicPeakDataGroups = new IEnumerable<ChromatographicPeakData>[numGroups];
            for (int i = 0; i < numGroups; i++)
            {
                ChromatographicPeakDataGroups[i] = CreateChromatographicPeakData(peaks, peakGroupIndices[i], maxThreads);
            }

            TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[] trainedModels 
                = new TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[numGroups];

            var trainer = mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features", numberOfTrees: 100);
            var pipeline = mlContext.Transforms.Concatenate("Features", trainingVariables)
                .Append(trainer);

            List<CalibratedBinaryClassificationMetrics> allMetrics = new List<CalibratedBinaryClassificationMetrics>();

            bool allSetsContainPositiveAndNegativeTrainingExamples = true;
            int groupNumber = 0;
            while (allSetsContainPositiveAndNegativeTrainingExamples == true && groupNumber < numGroups)
            {
                if (ChromatographicPeakDataGroups[groupNumber].Where(p => p.Label == true).Count() == 0 
                    || ChromatographicPeakDataGroups[groupNumber].Where(p => p.Label == false).Count() == 0)
                {
                    allSetsContainPositiveAndNegativeTrainingExamples = false;
                }
                groupNumber++;
            }

            if (allSetsContainPositiveAndNegativeTrainingExamples)
            {
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
                    if (outputFolder != null)
                    {
                        mlContext.Model.Save(trainedModels[groupIndexNumber], dataView.Schema, Path.Combine(outputFolder, "model.zip"));
                    }

                    Compute_Peak_PEP(peaks, peakGroupIndices[groupIndexNumber], mlContext, trainedModels[groupIndexNumber], outputFolder, maxThreads);

                    allMetrics.Add(metrics);
                }

                return AggregateMetricsForOutput(allMetrics);
            }
            else
            {
                return "Posterior error probability analysis failed. This can occur for small data sets when some sample groups are missing positive or negative training examples.";
            }
        }

        public static string ComputePEPValuesForAllPeaksIterations(ref List<ChromatographicPeak> peaks, string outputFolder, int maxThreads, double pepTrainingFraction = 0.25)
        {
            string[] trainingVariables = ChromatographicPeakData.trainingInfos["standard"];

            // Combine the groups and ensure that the order is always stable.
            peaks = peaks.OrderBy(p => p.GetHashCode()).ToList();

            var peakScores = peaks.Where(peak => peak.PipPep != null).Select(peak => (double)peak.PipPep).OrderBy(pep => pep).ToList();
            PipScoreCutoff = peakScores[(int)Math.Floor(peakScores.Count * pepTrainingFraction)]; //Select the top N percent of all peaks, only use those as positive examples

            MLContext mlContext = new MLContext();
            //the number of groups used for cross-validation is hard-coded at four. Do not change this number without changing other areas of effected code.
            const int numGroups = 3;

            List<int>[] peakGroupIndices = Get_Peak_Group_Indices_Iteration(peaks, numGroups);
            IEnumerable<ChromatographicPeakData>[] ChromatographicPeakDataGroups = new IEnumerable<ChromatographicPeakData>[numGroups];
            for (int i = 0; i < numGroups; i++)
            {
                ChromatographicPeakDataGroups[i] = CreateChromatographicPeakDataIteration(peaks, peakGroupIndices[i], maxThreads);
            }

            TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[] trainedModels
                = new TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[numGroups];

            var trainer = mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features", numberOfTrees: 100);
            var pipeline = mlContext.Transforms.Concatenate("Features", trainingVariables)
                .Append(trainer);

            List<CalibratedBinaryClassificationMetrics> allMetrics = new List<CalibratedBinaryClassificationMetrics>();

            bool allSetsContainPositiveAndNegativeTrainingExamples = true;
            int groupNumber = 0;
            while (allSetsContainPositiveAndNegativeTrainingExamples == true && groupNumber < numGroups)
            {
                if (ChromatographicPeakDataGroups[groupNumber].Where(p => p.Label == true).Count() == 0
                    || ChromatographicPeakDataGroups[groupNumber].Where(p => p.Label == false).Count() == 0)
                {
                    allSetsContainPositiveAndNegativeTrainingExamples = false;
                }
                groupNumber++;
            }

            if (allSetsContainPositiveAndNegativeTrainingExamples)
            {
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
                    if (outputFolder != null)
                    {
                        mlContext.Model.Save(trainedModels[groupIndexNumber], dataView.Schema, Path.Combine(outputFolder, "model.zip"));
                    }

                    Compute_Peak_PEP(peaks, peakGroupIndices[groupIndexNumber], mlContext, trainedModels[groupIndexNumber], outputFolder, maxThreads);

                    allMetrics.Add(metrics);
                }

                return AggregateMetricsForOutput(allMetrics);
            }
            else
            {
                return "Posterior error probability analysis failed. This can occur for small data sets when some sample groups are missing positive or negative training examples.";
            }
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

        public static void Compute_Peak_PEP(
            List<ChromatographicPeak> peaks,
            List<int> peakIndicies,
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

            Parallel.ForEach(Partitioner.Create(0, peakIndicies.Count),
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
                        ChromatographicPeak peak = peaks[peakIndicies[i]];

                        if (peak != null)
                        {
                            ChromatographicPeakData pd = CreateOneChromatographicPeakDataEntry(peak, label: !peak.RandomRt) ;
                            var pepValuePrediction = threadPredictionEngine.Predict(pd);
                            peak.PipPep = 1 - pepValuePrediction.Probability;

                            //A score is available using the variable pepvaluePrediction.Score
                        }
                    }
                });
        }

        //we add the indexes of the targets and decoys to the groups separately in the hope that we'll get at least one target and one decoy in each group.
        //then training can possibly be more successful.
        public static List<int>[] Get_Peak_Group_Indices(List<ChromatographicPeak> peaks, int numGroups)
        {
            //Dictionary<string, ChromatographicPeak> bestTargetPeaks = peaks
            //    .Where(p => !p.RandomRt)
            //    .GroupBy(p => p.Identifications.First().ModifiedSequence)
            //    .ToDictionary(g => g.Key, g => g.MaxBy(p => p.MbrScore));

            //Dictionary<string, ChromatographicPeak> bestDecoyPeaks = peaks
            //    .Where(p => p.RandomRt)
            //    .GroupBy(p => p.Identifications.First().ModifiedSequence)
            //    .ToDictionary(g => g.Key, g => g.MaxBy(p => p.MbrScore));

            List<int>[] groupsOfIndicies = new List<int>[numGroups];
            for (int i = 0; i < numGroups; i++)
            {
                groupsOfIndicies[i] = new List<int>();
            }

            // Targets and Decoys are labeled and used for training, remaining peaks are not
            List<int> targetPeakIndices = new List<int>();
            List<int> decoyPeakIndices = new List<int>();
            List<int> remainingPeakIndices = new List<int>();

            for (int i = 0; i < peaks.Count; i++)
            {
                if (peaks[i].RandomRt)
                    //&& bestDecoyPeaks.TryGetValue(peaks[i].Identifications.First().ModifiedSequence, out ChromatographicPeak bestDecoyPeak)
                    //&& peaks[i] == bestDecoyPeak)
                {
                    decoyPeakIndices.Add(i);
                }
                else if (!peaks[i].RandomRt
                    && peaks[i].MbrScore >= PipScoreCutoff )
                    //&& bestTargetPeaks.TryGetValue(peaks[i].Identifications.First().ModifiedSequence, out ChromatographicPeak bestTargetPeak)
                    //&& peaks[i] == bestTargetPeak)
                {
                    targetPeakIndices.Add(i);
                }
                else
                {
                    remainingPeakIndices.Add(i);
                }
            }

            int myIndex = 0;

            while (myIndex < decoyPeakIndices.Count)
            {
                int subIndex = 0;
                while (subIndex < numGroups && myIndex < decoyPeakIndices.Count)
                {
                    groupsOfIndicies[subIndex].Add(decoyPeakIndices[myIndex]);
                    subIndex++;
                    myIndex++;
                }
            }

            myIndex = 0;

            while (myIndex < targetPeakIndices.Count)
            {
                int subIndex = 0;
                while (subIndex < numGroups && myIndex < targetPeakIndices.Count)
                {
                    groupsOfIndicies[subIndex].Add(targetPeakIndices[myIndex]);
                    subIndex++;
                    myIndex++;
                }
            }

            myIndex = 0;

            while (myIndex < remainingPeakIndices.Count)
            {
                int subIndex = 0;
                while (subIndex < numGroups && myIndex < remainingPeakIndices.Count)
                {
                    groupsOfIndicies[subIndex].Add(remainingPeakIndices[myIndex]);
                    subIndex++;
                    myIndex++;
                }
            }

            return groupsOfIndicies;
        }

        //we add the indexes of the targets and decoys to the groups separately in the hope that we'll get at least one target and one decoy in each group.
        //then training can possibly be more successful.
        public static List<int>[] Get_Peak_Group_Indices_Iteration(List<ChromatographicPeak> peaks, int numGroups)
        {
            //Dictionary<string, ChromatographicPeak> bestTargetPeaks = peaks
            //    .Where(p => !p.RandomRt)
            //    .GroupBy(p => p.Identifications.First().ModifiedSequence)
            //    .ToDictionary(g => g.Key, g => g.MinBy(p => p.PipPep));

            //Dictionary<string, ChromatographicPeak> bestDecoyPeaks = peaks
            //    .Where(p => p.RandomRt)
            //    .GroupBy(p => p.Identifications.First().ModifiedSequence)
            //    .ToDictionary(g => g.Key, g => g.MinBy(p => p.PipPep));

            List<int>[] groupsOfIndicies = new List<int>[numGroups];
            for (int i = 0; i < numGroups; i++)
            {
                groupsOfIndicies[i] = new List<int>();
            }

            // Targets and Decoys are labeled and used for training, remaining peaks are not
            List<int> targetPeakIndices = new List<int>();
            List<int> decoyPeakIndices = new List<int>();
            List<int> remainingPeakIndices = new List<int>();

            for (int i = 0; i < peaks.Count; i++)
            {
                if (peaks[i].RandomRt)
                    //&& bestDecoyPeaks.TryGetValue(peaks[i].Identifications.First().ModifiedSequence, out ChromatographicPeak bestDecoyPeak)
                    //&& peaks[i] == bestDecoyPeak)
                {
                    decoyPeakIndices.Add(i);
                }
                else if (!peaks[i].RandomRt
                    && peaks[i].PipPep <= PipScoreCutoff)
                    //&& bestTargetPeaks.TryGetValue(peaks[i].Identifications.First().ModifiedSequence, out ChromatographicPeak bestTargetPeak)
                    //&& peaks[i] == bestTargetPeak)
                {
                    targetPeakIndices.Add(i);
                }
                else
                {
                    remainingPeakIndices.Add(i);
                }
            }

            int myIndex = 0;

            while (myIndex < decoyPeakIndices.Count)
            {
                int subIndex = 0;
                while (subIndex < numGroups && myIndex < decoyPeakIndices.Count)
                {
                    groupsOfIndicies[subIndex].Add(decoyPeakIndices[myIndex]);
                    subIndex++;
                    myIndex++;
                }
            }

            myIndex = 0;

            while (myIndex < targetPeakIndices.Count)
            {
                int subIndex = 0;
                while (subIndex < numGroups && myIndex < targetPeakIndices.Count)
                {
                    groupsOfIndicies[subIndex].Add(targetPeakIndices[myIndex]);
                    subIndex++;
                    myIndex++;
                }
            }

            myIndex = 0;

            while (myIndex < remainingPeakIndices.Count)
            {
                int subIndex = 0;
                while (subIndex < numGroups && myIndex < remainingPeakIndices.Count)
                {
                    groupsOfIndicies[subIndex].Add(remainingPeakIndices[myIndex]);
                    subIndex++;
                    myIndex++;
                }
            }

            return groupsOfIndicies;
        }


        public static IEnumerable<ChromatographicPeakData> CreateChromatographicPeakData(List<ChromatographicPeak> peaks, List<int> peakIndicies, int maxThreads)
        {
            object ChromatographicPeakDataListLock = new object();
            List<ChromatographicPeakData> ChromatographicPeakDataList = new List<ChromatographicPeakData>();
            int[] threads = Enumerable.Range(0, maxThreads).ToArray();

            Parallel.ForEach(Partitioner.Create(0, peakIndicies.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    List<ChromatographicPeakData> localChromatographicPeakDataList = new List<ChromatographicPeakData>();
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        var peak = peaks[peakIndicies[i]];
                        ChromatographicPeakData newChromatographicPeakData = new ChromatographicPeakData();

                        bool label;

                        if (peak.RandomRt)
                        {
                            label = false;
                            newChromatographicPeakData = CreateOneChromatographicPeakDataEntry(peak, label);
                            localChromatographicPeakDataList.Add(newChromatographicPeakData);
                        }
                        else if (!peak.RandomRt && peak.MbrScore >= PipScoreCutoff)
                        {
                            label = true;
                            newChromatographicPeakData = CreateOneChromatographicPeakDataEntry(peak, label);
                            localChromatographicPeakDataList.Add(newChromatographicPeakData);
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

        public static IEnumerable<ChromatographicPeakData> CreateChromatographicPeakDataIteration(List<ChromatographicPeak> peaks, List<int> peakIndicies, int maxThreads)
        {
            object ChromatographicPeakDataListLock = new object();
            List<ChromatographicPeakData> ChromatographicPeakDataList = new List<ChromatographicPeakData>();
            int[] threads = Enumerable.Range(0, maxThreads).ToArray();

            Parallel.ForEach(Partitioner.Create(0, peakIndicies.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    List<ChromatographicPeakData> localChromatographicPeakDataList = new List<ChromatographicPeakData>();
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        var peak = peaks[peakIndicies[i]];
                        ChromatographicPeakData newChromatographicPeakData = new ChromatographicPeakData();

                        bool label;

                        if (peak.RandomRt)
                        {
                            label = false;
                            newChromatographicPeakData = CreateOneChromatographicPeakDataEntry(peak, label);
                            localChromatographicPeakDataList.Add(newChromatographicPeakData);
                        }
                        else if (!peak.RandomRt && peak.PipPep <= PipScoreCutoff)
                        {
                            label = true;
                            newChromatographicPeakData = CreateOneChromatographicPeakDataEntry(peak, label);
                            localChromatographicPeakDataList.Add(newChromatographicPeakData);
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


        public static ChromatographicPeakData CreateOneChromatographicPeakDataEntry(ChromatographicPeak peak,bool label)
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