using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Accord.Math.Decompositions;
using Accord.Math;
using MathNet.Numerics.Statistics;
using System.Collections.Concurrent;

namespace FlashLFQ
{
    /// <summary>
    /// This is the "advanced" protein quantification engine used by FlashLFQ. It weights peptides by how well they co-vary with other peptides assigned to the same protein.
    /// The protein intensity is simply a weighted average of the peptide intensities. For a peptide to be used in protein quant, it must have a weight of at least 0.5.
    /// This uses a C# translation of a portion of Diffacto, which in turn is a python translation of an R program called FARMS.
    /// 
    /// Diffacto: 
    ///    Zhang B, Pirmoradian M, Zubarev R, and Käll L (2017). Covariation
    ///    of Peptide Abundances Accurately Reflects Protein Concentration
    ///    Differences. Molecular and Cellular Proteomics, mcp-O117,
    ///    http://www.mcponline.org/content/16/5/936.full.
    /// FARMS: 
    ///    Hochreiter S, Clevert D, and Obermayer K (2006). A new summarization
    ///    method for affymetrix probe level data. Bioinformatics, 22(8),
    ///    http://bioinformatics.oxfordjournals.org/cgi/content/abstract/22/8/943.
    /// </summary>
    public class ProteinQuantificationEngine
    {
        private const double MIN_WEIGHT = 0.5;
        private const double MU = 0.1;
        private const double MIN_NOISE = 1e-4;
        private const double ALPHA = 0.1;
        private const int MAX_ITER = 1000;
        private readonly FlashLfqResults results;
        private readonly int maxThreads;

        /// <summary>
        /// Constructs the protein quantification engine
        /// </summary>
        public ProteinQuantificationEngine(FlashLfqResults results, int maxThreads)
        {
            this.maxThreads = maxThreads;
            this.results = results;
        }

        /// <summary>
        /// Runs the protein quantification engine
        /// </summary>
        public void Run()
        {
            // link proteins to peptides
            Dictionary<ProteinGroup, List<Peptide>> proteinsToPeptides = new Dictionary<ProteinGroup, List<Peptide>>();
            foreach (var peptide in results.PeptideModifiedSequences)
            {
                foreach (var protein in peptide.Value.proteinGroups)
                {
                    if (proteinsToPeptides.TryGetValue(protein, out var peptides))
                    {
                        peptides.Add(peptide.Value);
                    }
                    else
                    {
                        proteinsToPeptides.Add(protein, new List<Peptide> { peptide.Value });
                    }
                }
            }

            var proteinList = proteinsToPeptides.ToList();

            Parallel.ForEach(Partitioner.Create(0, proteinList.Count), new ParallelOptions { MaxDegreeOfParallelism = maxThreads }, partitionRange =>
            {
                for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                {
                    // peptides must have at least one non-zero intensity measurement to be used for this protein quant analysis
                    List<Peptide> peptides = new List<Peptide>();

                    for (int p = 0; p < proteinList[i].Value.Count; p++)
                    {
                        foreach (var file in results.SpectraFiles.Where(f => f.TechnicalReplicate == 0))
                        {
                            if (proteinList[i].Value[p].GetIntensity(file) > 0)
                            {
                                peptides.Add(proteinList[i].Value[p]);
                                break;
                            }
                        }
                    }

                    // if this protein has no valid peptides (i.e., all missing values) its intensity is 0
                    if (!peptides.Any())
                    {
                        for (int s = 0; s < results.SpectraFiles.Count; s++)
                        {
                            proteinList[i].Key.SetIntensity(results.SpectraFiles[s], 0);
                        }

                        lock (results.ProteinGroups)
                        {
                            results.ProteinGroups.Add(proteinList[i].Key.ProteinGroupName, proteinList[i].Key);
                        }

                        continue;
                    }

                    // put intensities from the samples into a 2d (peptides x samples) array
                    double[][] intensityArray = GetIntensityArray(peptides);

                    // subtract reference sample and log-transform
                    for (int p = 0; p < intensityArray.Length; p++)
                    {
                        double avg = Math.Log(intensityArray[p].Average(), 2);

                        intensityArray[p] = intensityArray[p].Select(v => Math.Log(v, 2) - avg).ToArray();
                    }

                    // remove NaN/infinity
                    for (int p = 0; p < peptides.Count; p++)
                    {
                        for (int s = 0; s < intensityArray[p].Length; s++)
                        {
                            if (double.IsNaN(intensityArray[p][s]) || double.IsInfinity(intensityArray[p][s]))
                            {
                                intensityArray[p][s] = 0;
                            }
                        }
                    }

                    double[] weights;

                    if (peptides.Count > 1)
                    {
                        // calculate peptide weights with FARMS
                        var weightsAndNoise = FastFarms(readouts: intensityArray, mu: MU, weight: ALPHA, max_iter: MAX_ITER, force_iter: false, min_noise: MIN_NOISE);
                        weights = weightsAndNoise.Item1;
                    }
                    else
                    {
                        weights = new double[] { 1.0 };
                    }

                    double[] proteinIntensitiesPerFile = new double[results.SpectraFiles.Count];

                    // calculate weighted-average peptide intensity
                    for (int p = 0; p < peptides.Count; p++)
                    {
                        if (weights[p] >= MIN_WEIGHT && !double.IsNaN(weights[p]))
                        {
                            for (int s = 0; s < proteinIntensitiesPerFile.Length; s++)
                            {
                                proteinIntensitiesPerFile[s] += (weights[p] * peptides[p].GetIntensity(results.SpectraFiles[s]));
                            }
                        }
                    }

                    // store results
                    for (int s = 0; s < results.SpectraFiles.Count; s++)
                    {
                        proteinList[i].Key.SetIntensity(results.SpectraFiles[s], proteinIntensitiesPerFile[s]);
                    }

                    lock (results.ProteinGroups)
                    {
                        results.ProteinGroups.Add(proteinList[i].Key.ProteinGroupName, proteinList[i].Key);
                    }
                }
            });
        }

        /// <summary>
        /// Bayesian Factor Analysis for Proteomics Summarization
        ///    A C# translation of function "generateExprVal.method.farms" from
        ///    Bioconductor FARMS.
        ///    [http://www.bioconductor.org/packages/release/bioc/html/farms.html]
        ///    [http://www.bioinf.jku.at/publications/papers/farms/supplementary.ps]
        ///
        /// Reference:
        ///    Hochreiter S, Clevert D, and Obermayer K (2006). A new summarization
        ///    method for affymetrix probe level data. Bioinformatics, 22(8),
        ///    http://bioinformatics.oxfordjournals.org/cgi/content/abstract/22/8/943.
        ///
        /// Inputs:
        ///    probes: peptide abundance array(N peptides, M samples) in log scale.
        ///    weight: Hyperparameter(backscale factor) value in the range of[0, 1]
        ///             which determines the influence of the prior.
        ///    mu:     Hyperparameter value which allows to quantify different aspects
        ///             of potential prior knowledge.A value near zero assumes that
        ///             most genes do not contain a signal, and introduces a bias for
        ///            loading matrix elements near zero.
        /// </summary>
        public static Tuple<double[], double> FastFarms(double[][] readouts, double mu = 0, double weight = 0.5, int max_iter = 1000, bool force_iter = false, double min_noise = 1e-4, double fill_nan = 0.0)
        {
            // subtract average intensity
            for (int p = 0; p < readouts.Length; p++)
            {
                double iAvg = readouts[p].Average();
                readouts[p] = readouts[p].Select(i => i - iAvg).ToArray();
            }

            // calculate std dev
            double[] xsd = readouts.Select(p => Statistics.PopulationStandardDeviation(p)).ToArray();

            for (int p = 0; p < xsd.Length; p++)
            {
                if (xsd[p] < min_noise)
                {
                    xsd[p] = 1.0;
                }
            }

            // divide by std dev
            for (int p = 0; p < readouts.Length; p++)
            {
                readouts[p] = readouts[p].Select(i => i / xsd[p]).ToArray();
            }

            // handle infinity/NaN
            for (int p = 0; p < readouts.Length; p++)
            {
                for (int f = 0; f < readouts[p].Length; f++)
                {
                    if (double.IsNaN(readouts[p][f]) || double.IsInfinity(readouts[p][f]))
                    {
                        readouts[p][f] = fill_nan;
                    }
                }
            }

            // calculate covariance matrix
            double[,] C = GetCovarianceMatrix(readouts);

            // positive definite
            for (int p = 0; p < C.GetLength(0); p++)
            {
                for (int f = 0; f < C.GetLength(1); f++)
                {
                    if (C[p, f] < 0)
                        C[p, f] = 0;
                }
            }

            // singular value decomposition
            var svd = new SingularValueDecomposition(C);
            double[,] u = svd.LeftSingularVectors;
            double[] s = svd.Diagonal;
            double[,] v = svd.RightSingularVectors;
            v = Matrix.Transpose(v);

            // min noise
            for (int i = 0; i < s.Length; i++)
            {
                if (s[i] < min_noise)
                {
                    s[i] = min_noise;
                }
            }

            var sDiag = Matrix.Diagonal(s);
            C = Matrix.Dot(Matrix.Dot(u, sDiag), v);

            // initiation
            double[] λ = Matrix.Diagonal(C).Select(i => Math.Sqrt(i * 0.75)).ToArray();
            double[] ψ = new double[λ.Length];

            for (int i = 0; i < λ.Length; i++)
            {
                ψ[i] = Matrix.Diagonal(C)[i] - Math.Pow(λ[i], 2);
            }

            double[] old_psi = new double[ψ.Length];
            ψ.CopyTo(old_psi);

            double alpha = weight * readouts.Length;
            double E = 1.0;

            for (int i = 0; i < max_iter; i++)
            {
                // E step
                double[] φ = Elementwise.Divide(λ, ψ);
                double a = 1 + Matrix.Dot(λ, φ);
                double[] η = φ.Select(j => j / a).ToArray();
                double[] ζ = Matrix.Dot(C, η);
                E = 1 - Matrix.Dot(η, λ) + Matrix.Dot(η, ζ);

                // M step
                λ = Elementwise.Divide(ζ, Elementwise.Add(Elementwise.Multiply(ψ, alpha), E));

                var temp1 = Elementwise.Multiply(ζ, λ);
                var temp2 = Elementwise.Multiply(ψ, alpha);
                var temp3 = Elementwise.Multiply(temp2, λ);
                var temp4 = Elementwise.Subtract(mu, λ);
                var temp5 = Elementwise.Multiply(temp3, temp4);

                ψ = Elementwise.Add(Elementwise.Subtract(Matrix.Diagonal(C), temp1), temp5);

                for (int j = 0; j < ψ.Length; j++)
                {
                    if (ψ[j] < Math.Pow(min_noise, 2))
                    {
                        ψ[j] = Math.Pow(min_noise, 2);
                    }
                }

                if (!force_iter)
                {
                    if (Elementwise.Subtract(ψ, old_psi).Max(p => Math.Abs(p)) / old_psi.Max() < min_noise / 10)
                    {
                        break;
                    }
                }

                ψ.CopyTo(old_psi);
            }

            double[] loading = Elementwise.Multiply(λ, Math.Sqrt(E));
            var k = Elementwise.Divide(loading, ψ);

            double max = loading.Max();
            double[] weights = loading.Select(p => p / max).ToArray();
            double noise = 1.0 / (1.0 + Matrix.Dot(loading, k));

            // return weights, noise
            return new Tuple<double[], double>(weights, noise);
        }

        /// <summary>
        /// Builds a 2d matrix of covariance data 
        /// (how each peptide co-varies with the other peptides)
        /// </summary>
        private static double[,] GetCovarianceMatrix(double[][] data)
        {
            int n = data.Length;
            double[][] cov = new double[n][];

            for (int i = 0; i < n; i++)
            {
                cov[i] = new double[n];

                for (int j = 0; j < n; j++)
                {
                    cov[i][j] = Statistics.PopulationCovariance(data[i], data[j]);
                }
            }

            double[,] ret = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    ret[i, j] = cov[i][j];
                }
            }

            return ret;
        }

        /// <summary>
        /// Converts the list of peptides w/ intensities to a (peptide x intensity)
        /// 2d array structure for use in FARMS
        /// </summary>
        private double[][] GetIntensityArray(List<Peptide> peptides)
        {
            // only use the first tech rep for calculating peptide weights... could change this later
            var spectraFiles = results.SpectraFiles.Where(p => p.TechnicalReplicate == 0).ToList();

            double[][] intensityArray = new double[peptides.Count][];

            if (spectraFiles.Max(p => p.Fraction == 0))
            {
                // non-fractionated data
                for (int p = 0; p < peptides.Count; p++)
                {
                    intensityArray[p] = new double[spectraFiles.Count];

                    for (int s = 0; s < spectraFiles.Count; s++)
                    {
                        intensityArray[p][s] = peptides[p].GetIntensity(spectraFiles[s]);
                    }
                }
            }
            else
            {
                // fractionated data; need to sum by biorep before entering it into the array
                var cond = spectraFiles.Select(v => v.Condition).Distinct().ToList();
                Dictionary<Tuple<string, int>, int> conditionAndBiorepToSampleNumber = new Dictionary<Tuple<string, int>, int>();

                int sampleNumber = 0;
                for (int c = 0; c < cond.Count; c++)
                {
                    var bioreps = spectraFiles.Where(v => v.Condition == cond[c]).Select(v => v.BiologicalReplicate).Distinct();

                    foreach (var biorep in bioreps)
                    {
                        conditionAndBiorepToSampleNumber.Add(new Tuple<string, int>(cond[c], biorep), sampleNumber);
                        sampleNumber++;
                    }
                }

                for (int p = 0; p < peptides.Count; p++)
                {
                    intensityArray[p] = new double[sampleNumber + 1];

                    for (int s = 0; s < spectraFiles.Count; s++)
                    {
                        SpectraFileInfo spectraFile = spectraFiles[s];
                        sampleNumber = conditionAndBiorepToSampleNumber[new Tuple<string, int>(spectraFile.Condition, spectraFile.BiologicalReplicate)];

                        intensityArray[p][sampleNumber] += peptides[p].GetIntensity(spectraFiles[s]);
                    }
                }
            }

            return intensityArray;
        }
    }
}
