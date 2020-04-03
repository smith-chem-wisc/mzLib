using BayesianEstimation;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    public class PairedProteinQuantResult : ProteinQuantificationEngineResult
    {
        public readonly List<(Peptide peptide, List<double> foldChanges)> PeptideFoldChangeMeasurements;

        public PairedProteinQuantResult(ProteinGroup protein, List<Peptide> peptides, string controlCondition, string treatmentCondition,
            bool useSharedPeptides, FlashLfqResults results, Dictionary<(Peptide, string, int), (double, DetectionType)> peptideToSampleQuantity)
            : base(protein, peptides, controlCondition, treatmentCondition)
        {
            throw new NotImplementedException();
            //PeptideFoldChangeMeasurements = GetPeptideFoldChanges(useSharedPeptides, results, peptideToSampleQuantity);
        }

        /// <summary>
        /// Estimates the fold-change between conditions for a protein, given its constituent peptides' fold change measurements.
        /// </summary>
        public void EstimateProteinFoldChange(int? randomSeed, int nBurnin, int n, double nullHypothesisCutoff = double.NaN)
        {
            throw new NotImplementedException();
            //    IsStatisticallyValid = this.PeptideFoldChangeMeasurements.Count > 1;

            //    bool skepticalPrior = !double.IsNaN(nullHypothesisCutoff);

            //    var res = ProteinQuantificationEngine.FitProteinQuantModel(
            //        measurements: PeptideFoldChangeMeasurements.SelectMany(p => p.foldChanges.Select(v => new Datum(v, 1))).ToList(), //TODO: weight
            //        skepticalPrior: skepticalPrior,
            //        paired: true,
            //        randomSeed: randomSeed,
            //        burnin: nBurnin,
            //        n: n,
            //        nullHypothesisInterval: nullHypothesisCutoff,
            //        sdPrior: null); //TODO

            //    if (!skepticalPrior)
            //    {
            //        this.FoldChangePointEstimate = res.mus.Median();
            //        this.StandardDeviationPointEstimate = res.sds.Median();
            //        this.NuPointEstimate = res.sds.Median();
            //    }
            //    else
            //    {
            //        this.NullHypothesisInterval = nullHypothesisCutoff;
            //        this.CalculatePosteriorErrorProbability(res.mus);
            //    }
        }

        //public override string ToString()
        //{
        //    int nPeptides = PeptideFoldChangeMeasurements.Count;
        //    int nMeasurements = PeptideFoldChangeMeasurements.SelectMany(p => p.foldChanges).Count();
        //    var measurementsString = string.Join(",", PeptideFoldChangeMeasurements.Select(p => p.peptide.Sequence + ":" + string.Join(";", p.foldChanges.Select(v => v.ToString("F3")))));

        //    if (measurementsString.Length > 32000)
        //    {
        //        measurementsString = "Too long for Excel";
        //    }

        //    return
        //        Protein.ProteinGroupName + "\t" +
        //        Protein.GeneName + "\t" +
        //        Protein.Organism + "\t" +
        //        ControlCondition + "\t" +
        //        TreatmentCondition + "\t" +
        //        NullHypothesisInterval.Value + "\t" +
        //        FoldChangePointEstimate + "\t" +
        //        StandardDeviationPointEstimate + "\t" +
        //        ControlConditionIntensity + "\t" +
        //        TreatmentConditionIntensity + "\t" +
        //        nPeptides + "\t" +
        //        nMeasurements + "\t" +
        //        measurementsString + "\t" +
        //        PosteriorErrorProbability + "\t" +
        //        FalseDiscoveryRate + "\t\t";
        //}

        public static new string TabSeparatedHeader()
        {
            throw new NotImplementedException();
            //    return
            //        "Protein Group" + "\t" +
            //        "Gene" + "\t" +
            //        "Organism" + "\t" +
            //        "Control Condition" + "\t" +
            //        "Treatment Condition" + "\t" +
            //        "Log2 Fold-Change Cutoff" + "\t" +
            //        "Protein Log2 Fold-Change" + "\t" +
            //        "Standard Deviation of Peptide Log2 Fold-Changes" + "\t" +
            //        "Protein Intensity in Control Condition" + "\t" +
            //        "Protein Intensity in Treatment Condition" + "\t" +
            //        "Number of Peptides" + "\t" +
            //        "Number of Fold-Change Measurements" + "\t" +
            //        "List of Fold-Change Measurements Grouped by Peptide" + "\t" +
            //        "Posterior Error Probability" + "\t" +
            //        "False Discovery Rate" + "\t\t";
        }

        /// <summary>
        /// Computes a list of fold-change measurements between the constituent peptides of this protein between the control and treatment condition.
        /// </summary>
        //private List<(Peptide peptide, List<double> foldChanges)> GetPeptideFoldChanges(bool useSharedPeptides, FlashLfqResults flashLfqResults,
        //    Dictionary<(Peptide, string, int), (double, DetectionType)> PeptideToSampleQuantity)
        //{
        //    List<(Peptide, List<double>)> allPeptideFoldChanges = new List<(Peptide, List<double>)>();

        //    int numSamples = flashLfqResults.SpectraFiles.Where(p => p.Condition == ControlCondition).Max(p => p.BiologicalReplicate) + 1;

        //    foreach (Peptide peptide in Peptides)
        //    {
        //        if (!peptide.UseForProteinQuant || (!useSharedPeptides && peptide.ProteinGroups.Count > 1))
        //        {
        //            continue;
        //        }

        //        List<double> peptideFoldChanges = new List<double>();

        //        for (int sample = 0; sample < numSamples; sample++)
        //        {
        //            if (PeptideToSampleQuantity.TryGetValue((peptide, ControlCondition, sample), out var controlIntensity)
        //                && PeptideToSampleQuantity.TryGetValue((peptide, TreatmentCondition, sample), out var treatmentIntensity))
        //            {
        //                double? foldChange = GetLogFoldChange(controlIntensity.Item1, treatmentIntensity.Item1);

        //                if (foldChange != null)
        //                {
        //                    peptideFoldChanges.Add(foldChange.Value);
        //                    NMeasurements++;
        //                }
        //            }
        //        }

        //        allPeptideFoldChanges.Add((peptide, peptideFoldChanges));
        //    }

        //    return allPeptideFoldChanges;
        //}

        /// <summary>
        /// Computes the log-fold change between two intensities. If there is an error (e.g., one of the intensities is zero), 
        /// null is returned.
        /// </summary>
        //private double? GetLogFoldChange(double intensity1, double intensity2)
        //{
        //    double logFoldChange = Math.Log(intensity2, 2) - Math.Log(intensity1, 2);

        //    if (!double.IsNaN(logFoldChange) && !double.IsInfinity(logFoldChange))
        //    {
        //        return logFoldChange;
        //    }

        //    return null;
        //}
    }
}
