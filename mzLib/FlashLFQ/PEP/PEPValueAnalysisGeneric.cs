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
using Microsoft.FSharp.Collections;

namespace FlashLFQ.PEP
{
    public class DonorGroup : IEnumerable<ChromatographicPeak>
    {
        public Identification DonorId { get; }
        public List<ChromatographicPeak> Acceptors { get; }
        public bool IsPeakDecoy { get; }

        public DonorGroup(Identification donorId, List<ChromatographicPeak> acceptors, bool isPeakDecoy)
        {
            DonorId = donorId;
            Acceptors = acceptors;
            IsPeakDecoy = isPeakDecoy;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(DonorId.GetHashCode(), IsPeakDecoy.GetHashCode());
        }

        public ChromatographicPeak BestPeakByScore => Acceptors.MaxBy(acceptor => acceptor.MbrScore);

        // The null coalescing is kinda bad here, really it should throw an exception if the value is null
        public ChromatographicPeak BestPeakByPep => Acceptors.MinBy(acceptor => acceptor.PipPep ?? 1);

        public double BestMbrScore => Acceptors.Max(acceptor => acceptor.MbrScore);
        public double BestPep => Acceptors.Min(acceptor => acceptor.PipPep ?? 1);

        public IEnumerator<ChromatographicPeak> GetEnumerator()
        {
            return Acceptors.GetEnumerator();
        }

        System.Collections.IEnumerator System.Collections.IEnumerable.GetEnumerator()
        {
            return Acceptors.GetEnumerator();
        }
    }


    public static class PEP_Analysis_Cross_Validation
    {
        public static double PipScoreCutoff;

        public static string ComputePEPValuesForAllPeaks(ref List<ChromatographicPeak> peaks, string outputFolder, int maxThreads, double pepTrainingFraction = 0.25)
        {
            string[] trainingVariables = ChromatographicPeakData.trainingInfos["standard"];

            #region Construct Donor Groups
            // this is target peak not target peptide
            List<DonorGroup> donors= new();
            foreach(var donorGroup in peaks.Where(peak => peak.IsMbrPeak & !peak.RandomRt).GroupBy(peak => peak.Identifications.First())) //Group by donor peptide
            {
                var donorId = donorGroup.Key;
                var acceptors = donorGroup.ToList();
                donors.Add(new DonorGroup(donorId, acceptors, isPeakDecoy: false));
            }

            foreach (var donorGroup in peaks.Where(peak => peak.IsMbrPeak & peak.RandomRt).GroupBy(peak => peak.Identifications.First())) //Group by donor peptide
            {
                var donorId = donorGroup.Key;
                var acceptors = donorGroup.ToList();
                donors.Add(new DonorGroup(donorId, acceptors, isPeakDecoy: true));
            }

            // Fix the order
            donors = donors.OrderByDescending(donor => donor.GetHashCode()).ToList();

            var peakScores = donors.Select(donor => donor.BestMbrScore).OrderByDescending(score => score).ToList();
            PipScoreCutoff = peakScores[(int)Math.Floor(peakScores.Count * pepTrainingFraction)]; //Select the top N percent of all peaks, only use those as positive examples

            MLContext mlContext = new MLContext();
            //the number of groups used for cross-validation is hard-coded at four. Do not change this number without changing other areas of effected code.
            const int numGroups = 3;

            List<int>[] donorGroupIndices = Get_Donor_Group_Indices(donors, numGroups);

            #endregion

            #region Create Groups and Model
            IEnumerable<ChromatographicPeakData>[] ChromatographicPeakDataGroups = new IEnumerable<ChromatographicPeakData>[numGroups];
            for (int i = 0; i < numGroups; i++)
            {
                ChromatographicPeakDataGroups[i] = CreateChromatographicPeakData(donors, donorGroupIndices[i], maxThreads);

                if (!ChromatographicPeakDataGroups[i].Any(p => p.Label == true) 
                    || !ChromatographicPeakDataGroups[i].Any(p => p.Label == false))
                {
                    return "Posterior error probability analysis failed. This can occur for small data sets when some sample groups are missing positive or negative training examples.";
                }
            }

            TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[] trainedModels 
                = new TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[numGroups];

            var trainer = mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features", numberOfTrees: 100);
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
                if (outputFolder != null)
                {
                    mlContext.Model.Save(trainedModels[groupIndexNumber], dataView.Schema, Path.Combine(outputFolder, "model.zip"));
                }

                Compute_PEP_For_All_Peaks(donors, donorGroupIndices[groupIndexNumber], mlContext, trainedModels[groupIndexNumber], outputFolder, maxThreads);

                allMetrics.Add(metrics);
            }

            #endregion
            #region Iterative Training

            for(int trainingIteration = 0; trainingIteration < 4; trainingIteration++)
            {
                ChromatographicPeakDataGroups = new IEnumerable<ChromatographicPeakData>[numGroups];
                for (int i = 0; i < numGroups; i++)
                {
                    ChromatographicPeakDataGroups[i] = CreateChromatographicPeakDataIteration(donors, donorGroupIndices[i], maxThreads);

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
                    if (outputFolder != null)
                    {
                        mlContext.Model.Save(trainedModels[groupIndexNumber], dataView.Schema, Path.Combine(outputFolder, "model.zip"));
                    }

                    Compute_PEP_For_All_Peaks(donors, donorGroupIndices[groupIndexNumber], mlContext, trainedModels[groupIndexNumber], outputFolder, maxThreads);

                    allMetrics.Add(metrics);
                }
            }
            #endregion

            return AggregateMetricsForOutput(allMetrics);
        }

        //we add the indexes of the targets and decoys to the groups separately in the hope that we'll get at least one target and one decoy in each group.
        //then training can possibly be more successful.
        public static List<int>[] Get_Donor_Group_Indices(List<DonorGroup> donors, int numGroups)
        {
            List<int>[] groupsOfIndicies = new List<int>[numGroups];
            for (int i = 0; i < numGroups; i++)
            {
                groupsOfIndicies[i] = new List<int>();
            }

            // Targets and Decoys are labeled and used for training, remaining peaks are not
            List<int> targetIndices = new List<int>();
            List<int> decoyIndices = new List<int>();
            List<int> remainingIndices = new List<int>();

            for (int i = 0; i < donors.Count; i++)
            {
                if (donors[i].IsPeakDecoy)
                {
                    decoyIndices.Add(i);
                }
                else if (!donors[i].IsPeakDecoy
                    && donors[i].BestMbrScore >= PipScoreCutoff)
                {
                    targetIndices.Add(i);
                }
                else
                {
                    remainingIndices.Add(i);
                }
            }

            int myIndex = 0;

            while (myIndex < decoyIndices.Count)
            {
                int subIndex = 0;
                while (subIndex < numGroups && myIndex < decoyIndices.Count)
                {
                    groupsOfIndicies[subIndex].Add(decoyIndices[myIndex]);
                    subIndex++;
                    myIndex++;
                }
            }

            myIndex = 0;

            while (myIndex < targetIndices.Count)
            {
                int subIndex = 0;
                while (subIndex < numGroups && myIndex < targetIndices.Count)
                {
                    groupsOfIndicies[subIndex].Add(targetIndices[myIndex]);
                    subIndex++;
                    myIndex++;
                }
            }

            myIndex = 0;

            while (myIndex < remainingIndices.Count)
            {
                int subIndex = 0;
                while (subIndex < numGroups && myIndex < remainingIndices.Count)
                {
                    groupsOfIndicies[subIndex].Add(remainingIndices[myIndex]);
                    subIndex++;
                    myIndex++;
                }
            }

            return groupsOfIndicies;
        }

        public static IEnumerable<ChromatographicPeakData> CreateChromatographicPeakData(List<DonorGroup> donors, List<int> donorIndices, int maxThreads)
        {
            object ChromatographicPeakDataListLock = new object();
            List<ChromatographicPeakData> ChromatographicPeakDataList = new List<ChromatographicPeakData>();
            int[] threads = Enumerable.Range(0, maxThreads).ToArray();

            Parallel.ForEach(Partitioner.Create(0, donorIndices.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    List<ChromatographicPeakData> localChromatographicPeakDataList = new List<ChromatographicPeakData>();
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        var donor = donors[donorIndices[i]];
                        ChromatographicPeakData newChromatographicPeakData = new ChromatographicPeakData();

                        bool label;

                        if (donor.IsPeakDecoy)
                        {
                            label = false;
                            newChromatographicPeakData = CreateOneChromatographicPeakDataEntry(donor.BestPeakByScore, label);
                            localChromatographicPeakDataList.Add(newChromatographicPeakData);
                        }
                        else if (!donor.IsPeakDecoy && donor.BestMbrScore >= PipScoreCutoff)
                        {
                            label = true;
                            newChromatographicPeakData = CreateOneChromatographicPeakDataEntry(donor.BestPeakByScore, label);
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

                        foreach(ChromatographicPeak peak in donor.Acceptors)
                        {
                            ChromatographicPeakData pd = CreateOneChromatographicPeakDataEntry(peak, label: !peak.RandomRt);
                            var pepValuePrediction = threadPredictionEngine.Predict(pd);
                            peak.PipPep = 1 - pepValuePrediction.Probability;
                        }
                    }
                });
        }

        public static IEnumerable<ChromatographicPeakData> CreateChromatographicPeakDataIteration(List<DonorGroup> donors, List<int> donorIndices, int maxThreads)
        {
            object ChromatographicPeakDataListLock = new object();
            List<ChromatographicPeakData> ChromatographicPeakDataList = new List<ChromatographicPeakData>();
            int[] threads = Enumerable.Range(0, maxThreads).ToArray();

            List<double> peps = donors.Select(donors => donors.BestPep).OrderBy(pep => pep).ToList();
            double groupSpecificPepCutoff = peps[(int)Math.Floor(peps.Count * 0.25)];

            Parallel.ForEach(Partitioner.Create(0, donorIndices.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    List<ChromatographicPeakData> localChromatographicPeakDataList = new List<ChromatographicPeakData>();
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        var donor = donors[donorIndices[i]];
                        ChromatographicPeakData newChromatographicPeakData = new ChromatographicPeakData();

                        bool label;
                        if (donor.IsPeakDecoy)
                        {
                            label = false;
                            newChromatographicPeakData = CreateOneChromatographicPeakDataEntry(donor.BestPeakByPep, label);
                            localChromatographicPeakDataList.Add(newChromatographicPeakData);
                        }
                        else if (!donor.IsPeakDecoy && donor.BestPep <= groupSpecificPepCutoff)
                        {
                            label = true;
                            newChromatographicPeakData = CreateOneChromatographicPeakDataEntry(donor.BestPeakByPep, label);
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