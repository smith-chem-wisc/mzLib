using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    internal class StatisticalAnalysisEngine
    {
        private FlashLFQResults Results;
        private readonly double PValueForCallingSignificance;
        private readonly double MinimumFoldChange;

        public StatisticalAnalysisEngine(FlashLFQResults results, double pValueForCallingSignificance, double minimumFoldChange)
        {
            Results = results;
            PValueForCallingSignificance = pValueForCallingSignificance;
            MinimumFoldChange = minimumFoldChange;
        }

        public void PerformStatisticalAnalysis()
        {
            // calculate intensities for proteins/peptides
            Results.CalculatePeptideResults(true);
            Results.CalculatePeptideResults(false);
            Results.CalculateProteinResults();

            List<string> output = new List<string>();
            Dictionary<string, List<double[]>> conditionToSamples = new Dictionary<string, List<double[]>>();
            List<Peptide> peptides = Results.peptideModifiedSequences.Values.ToList();
            int samples = 0;

            // sum intensities by biorep (bioreps are the samples for t-testing)
            var conditions = Results.spectraFiles.GroupBy(p => p.Condition);

            foreach (var condition in conditions)
            {
                List<double[]> samplesForThisCondition = new List<double[]>();
                var bioreps = condition.GroupBy(p => p.BiologicalReplicate);

                foreach (var biorep in bioreps)
                {
                    samples++;
                    double[] biorepIntensities = new double[peptides.Count];

                    var fractions = biorep.GroupBy(p => p.Fraction);

                    foreach (var fraction in fractions)
                    {
                        // the technical reps are averaged, so they are interchangeable; just use the first one
                        SpectraFileInfo file = fraction.First();

                        for (int p = 0; p < peptides.Count; p++)
                        {
                            biorepIntensities[p] += peptides[p].GetIntensity(file);
                        }
                    }

                    samplesForThisCondition.Add(biorepIntensities);
                }

                conditionToSamples.Add(condition.Key, samplesForThisCondition);
            }

            // log-transform intensities
            foreach (var condition in conditionToSamples)
            {
                foreach (var biorep in condition.Value)
                {
                    for (int p = 0; p < biorep.Length; p++)
                    {
                        if (biorep[p] > 0)
                        {
                            biorep[p] = Math.Log(biorep[p], 2);
                        }
                    }
                }
            }

            // perform t tests between conditions
            // we're assuming two samples and that the first condition is the one to t-test against
            var cond = conditionToSamples.ToList();

            for (int c = 1; c < cond.Count; c++)
            {
                // calculate p-value for each peptide
                for (int p = 0; p < peptides.Count; p++)
                {
                    // get the intensity of the peptide in the bioreps for the two conditions
                    var cond1Samples = cond[0].Value.Select(biorep => biorep[p]).ToList();
                    var cond2Samples = cond[c].Value.Select(biorep => biorep[p]).ToList();

                    int degreesOfFreedom = cond1Samples.Count + cond2Samples.Count - 2;

                    // calculate test statistic
                    double avgCondition1 = cond1Samples.Average();
                    double avgCondition2 = cond2Samples.Average();

                    double varianceCond1 = Statistics.Variance(cond1Samples);
                    double varianceCond2 = Statistics.Variance(cond2Samples);

                    double pooledVariance = ((cond1Samples.Count - 1) * varianceCond1 + (cond2Samples.Count - 1) * varianceCond2) / (degreesOfFreedom);

                    double testStatistic = (avgCondition1 - avgCondition2) / Math.Sqrt(pooledVariance / cond1Samples.Count + pooledVariance / cond2Samples.Count);

                    // calculate p-value
                    double pValue = double.NaN;
                    if (degreesOfFreedom > 0 && (!cond1Samples.Any(v => v == 0) && !cond2Samples.Any(v => v == 0)))
                    {
                        pValue = 2.0 * (1.0 - StudentT.CDF(0, 1, degreesOfFreedom, Math.Abs(testStatistic)));
                    }

                    // TODO: this cuts off at log fold change, not fold change
                    bool isSignificantChanger = (pValue <= PValueForCallingSignificance) && (Math.Abs(avgCondition1 - avgCondition2) >= MinimumFoldChange);

                    output.Add(peptides[p].Sequence + "\t" +
                        string.Join("\t", cond1Samples.Select(v => Math.Pow(2, v))) + "\t" +
                        string.Join("\t", cond2Samples.Select(v => Math.Pow(2, v))) + "\t" +
                        string.Join("\t", cond1Samples.Select(v => v)) + "\t" +
                        string.Join("\t", cond2Samples.Select(v => v)) + "\t" +
                        "\t" + (avgCondition2 - avgCondition1) + "\t" +
                        pValue + "\t" +
                        isSignificantChanger);
                }

                // TODO: apply Benjamini-Hochberg correction
            }

            // DEBUG write output
            //System.IO.File.WriteAllLines(@"C:\Data\InSilicoFractionatedYeast\StatisticalTesting.tsv", output);
        }
    }
}