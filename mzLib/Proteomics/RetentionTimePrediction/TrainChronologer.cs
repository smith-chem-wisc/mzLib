using Proteomics.PSM;
using System;
using System.Collections.Generic;
using System.Linq;
using TorchSharp;
using TorchSharp.Modules;

namespace Proteomics.RetentionTimePrediction
{
    public static class TrainChronologer
    {
        public static torch.Tensor Tensorize(PsmFromTsv psm)
        {
            if (psm.BaseSeq.Length <= 50)
            {
                var dictionary = new Dictionary<(char, string), int>()
                {
                    { ('A', ""), 1  }, //'Alanine
                    { ('C', ""), 2  }, //'Cysteine
                    { ('D', ""), 3  }, //'Aspartate
                    { ('E', ""), 4  }, //'Glutamate
                    { ('F', ""), 5  }, //'Phenylalaline
                    { ('G', ""), 6  }, //'Glycine
                    { ('H', ""), 7  }, //'Histidine
                    { ('I', ""), 8  }, //'Isoleucine
                    { ('K', ""), 9  }, //'Lysine
                    { ('L', ""), 10 }, //'Leucine
                    { ('M', ""), 11 }, //'Methionine
                    { ('N', ""), 12 }, //'Asparagine
                    { ('P', ""), 13 }, //'Proline
                    { ('Q', ""), 14 }, //'Glutamine
                    { ('R', ""), 15 }, //'Argenine
                    { ('S', ""), 16 }, //'Serine
                    { ('T', ""), 17 }, //'Threonine
                    { ('V', ""), 18 }, //'Valine
                    { ('W', ""), 19 }, //'Tryptophane
                    { ('Y', ""), 20 }, //'Tyrosine
                    //{ ('C', "Carbamidomethyl on C"), 21 }, //'Carbamidomethyl
                    //{ ('M', "Oxidation on M"), 22 }, //'Oxidized
                    ////_residueWithModToTensorInt.Add(('C',null),23);//'S - carbamidomethylcysteine
                    //{ ('E', "Glu to PyroGlu"), 24 }, //'Pyroglutamate
                    //{ ('S', "Phosphorylation on S"), 25 }, //'Phosphoserine
                    //{ ('T', "Phosphorylation on T"), 26 }, //'Phosphothreonine
                    //{ ('Y', "Phosphorylation on Y"), 27 }, //'Phosphotyrosine
                    //{ ('K', "Accetylation on K"), 28 }, //'Acetylated
                    //{ ('K', "Succinylation on K"), 29 }, //'Succinylated
                    //{ ('K', "Ubiquitination on K"), 30 }, //'Ubiquitinated
                    //{ ('K', "Methylation on K"), 31 }, //'Monomethyl
                    //{ ('K', "Dimethylation on K"), 32 }, //'Dimethyl
                    //{ ('K', "Trimethylation on K"), 33 }, //'Trimethyl
                    //{ ('R', "Methylation on R"), 34 }, //'Monomethyl
                    //{ ('R', "Dimethylation on R"), 35 }, //'Dimethyl
                };

                var tensor = torch.zeros(1, 52, torch.ScalarType.Int64);

                tensor[0][0] = 38; //C-terminus
                tensor[0][51] = 44; //N-terminus

                for (int i = 1; i < psm.BaseSeq.Length - 1; i++) //base sequence for the moment
                {
                    tensor[0][i] = dictionary[(psm.BaseSeq[i], "")];
                }

                return tensor;
            }

            return torch.ones(1, 52, torch.ScalarType.Int64);
        }

        //todo: implement this method training_fuctions.py
        public static (List<(torch.Tensor, double)>, List<(torch.Tensor, double)>) RetentionTimeToTensorDatabase(
            List<PsmFromTsv> dataFiles, int seed, double validationFraction)
        {
            var trainTestDb = new Dictionary<string, List<(torch.Tensor, double)>>()
            {
                {"train", new List<(torch.Tensor, double)>()},
                {"test", new List<(torch.Tensor, double)>()}
            };

            var allData = new List<(torch.Tensor, double)>();

            var sources = new HashSet<string>();

            foreach (var dataFile in dataFiles)
            {
                if (dataFile.DecoyContamTarget.Equals("T"))
                {
                    var db =
                        (dataFile.FileNameWithoutExtension, dataFile.RetentionTime, dataFile.BaseSeq); //base seq for the moment
                    var tensor = Tensorize(dataFile);
                    
                    if(tensor.Equals(torch.ones(1,52,torch.ScalarType.Int64)))
                        continue;

                    if (tensor[0][0].item<Int64>().Equals((Int64)38))
                    {
                        allData.Add((tensor, db.RetentionTime.Value)); //todo: add encoded sequence tensor
                        sources.Add(db.FileNameWithoutExtension);
                    }

                    var a = 0;
                }
            }

            var trainingSet = allData.Take((int)(allData.Count * (1 - validationFraction))).ToList();
            var testSet = allData.Skip((int)(allData.Count * (1 - validationFraction))).ToList();

            trainTestDb["train"] = trainingSet;
            trainTestDb["test"] = testSet;

            return (trainTestDb["train"], trainTestDb["test"]);
        }

        /// <summary>
        /// LogLLoss module for Chronologer.
        ///
        /// https://github.com/searlelab/chronologer/blob/main/src/chronologer/training_functions.py
        /// </summary>
        internal class LogLLoss : torch.nn.Module<torch.Tensor, torch.Tensor, torch.Tensor>
        {
            public LogLLoss(int numberOfSources,
                ReturnDistribution family = ReturnDistribution.Laplace,
                double fdr = 0.01) : base(nameof(LogLLoss))
            {
                _numberOfSources = numberOfSources;
                _family = family;
                _fdr = fdr;

                _sourceScale = torch.nn.Linear(_numberOfSources, 1, false);
                torch.nn.init.constant_(_sourceScale.weight, 10.0);


                RegisterComponents();

            }
            public override torch.Tensor forward(torch.Tensor prediction, torch.Tensor target)
            {
                var scale = _sourceScale.forward(target).clamp(1e-7);
                var distribution = new Laplace(center: target, scale: scale);
                var logL_loss = -1 * distribution.LogL(prediction);
                return logL_loss.mean();//todo: not finished, its a place holder
            }

            private int _numberOfSources { get; set; }
            private double _fdr { get; set; }
            private ReturnDistribution _family { get; set; }
            private Linear _sourceScale { get; set; }
        }

        public enum ReturnDistribution
        {
            Gaussian,
            Laplace,
            Gumbel
        }

        internal class Gaussian
        {
            public Gaussian(torch.Tensor data = null, torch.Tensor center = null,
                torch.Tensor scale = null)
            {
                if (data.NumberOfElements == 0)
                {
                    if (center.GetType() == typeof(torch.Tensor))
                    {
                        _mu = center;
                    }
                    else
                    {
                        _mu = torch.full(new long[] { 1 }, center.ToScalar());
                    }

                    if (scale.GetType() == typeof(torch.Tensor))
                    {
                        _sigma = scale;
                    }
                    else
                    {
                        _sigma = torch.full(new long[] { 1 }, scale.ToScalar());
                    }
                }
                else
                {
                    _mu = torch.mean(data);
                    _sigma = torch.std(data);
                }
            }

            public torch.Tensor CDF(torch.Tensor x)
            {
                return 0.5 * (1 + torch.erf((x - _mu) / (_sigma * Math.Sqrt(2))));
            }

            public torch.Tensor PPF(torch.Tensor q)
            {
                return _mu + _sigma * Math.Sqrt(2) * torch.erfinv(2 * q - 1);
            }

            public torch.Tensor LogL(torch.Tensor x)
            {
                return -1 * (torch.log(_sigma * Math.Sqrt(2 * Math.PI)) + 0.5 * torch.pow((x - _mu), 1) / _sigma);
            }
            private torch.Tensor _mu { get; set; }
            private torch.Tensor _sigma { get; set; }
        }

        internal class Laplace
        {
            public Laplace(torch.Tensor data = null, torch.Tensor center = null,
                torch.Tensor scale = null)
            {
                if (data.NumberOfElements == 0)
                {
                    if (center.GetType() == typeof(torch.Tensor))
                    {
                        _mu = center;
                    }
                    else
                    {
                        _mu = torch.full(new long[] { 1 }, center.ToScalar());
                    }

                    if (scale.GetType() == typeof(torch.Tensor))
                    {
                        _b = scale;
                    }
                    else
                    {
                        _b = torch.full(new long[] { 1 }, scale.ToScalar());
                    }
                }
                else
                {
                    _mu = torch.mean(data);
                    _b = torch.mean(torch.abs(data - _mu));
                }
            }


            public torch.Tensor CDF(torch.Tensor x)
            {
                if (x.less_equal(_mu).Equals(true))
                {
                    return 0.5 * torch.exp((x - _mu) / _b);
                }

                return 1 - 0.5 * torch.exp((-1 * (x - _mu)) / _b);
            }

            public torch.Tensor PPF(torch.Tensor q)
            {
                if (q.less_equal(0.5).Equals(true))
                {
                    return _mu + _b * torch.log(2 * q);
                }
                else
                {
                    return _mu - _b * torch.log(2 - 2 * q);
                }
            }
            public torch.Tensor LogL(torch.Tensor x)
            {
                return -1 * (torch.log(2 * _b) + torch.abs(x - _mu) / _b);
            }
            private torch.Tensor _mu { get; set; }
            private torch.Tensor _b { get; set; }
        }

        internal class Gumbel
        {
            public Gumbel(torch.Tensor data)
            {
                var mean = torch.mean(data);
                var std = torch.std(data);
                _beta = std * Math.Sqrt(6) / Math.PI;
                _mu = mean - 0.57721 * _beta;
            }

            private torch.Tensor CDF(torch.Tensor x)
            {
                return torch.exp(-1 * torch.exp(-1 * (x - _mu) / _beta));
            }

            private torch.Tensor PPF(torch.Tensor q)
            {
                return _mu - _beta * torch.log(-1 * torch.log(q));
            }

            private torch.Tensor _beta { get; set; }
            private torch.Tensor _mu { get; set; }
        }
    }
}
