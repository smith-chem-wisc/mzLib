using System;
using System.Collections.Generic;
using TorchSharp;
using TorchSharp.Modules;

namespace Proteomics.RetentionTimePrediction
{
    public static class TrainChronologer
    {
        //todo: implement this method training_fuctions.py
        public static (Dictionary<string, torch.utils.data.Dataset>, int) RetentionTimeToTensorDatabase(
            int[] dataFiles)
        {
            var dbs = new Dictionary<string, List<object>>()
            {
                {"train", new List<object>()},
                {"test", new List<object>()}
            };

            var sources = new HashSet<string>();

            foreach (var dataFile in dataFiles)
            {

            }
            throw new NotImplementedException();
        }
            
        /// <summary>
        /// LogLLoss module for Chronologer.
        ///
        /// https://github.com/searlelab/chronologer/blob/main/src/chronologer/training_functions.py
        /// </summary>
        internal class LogLLoss : torch.nn.Module<torch.Tensor, torch.Tensor>
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
            public override torch.Tensor forward(torch.Tensor input)
            {
                return torch.nn.functional.log_softmax(input, 1);
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
            

            private torch.Tensor CDF(torch.Tensor x)
            {
                if (x.less_equal(_mu).Equals(true))
                {
                    return 0.5 * torch.exp((x - _mu) / _b);
                }
                
                return 1 - 0.5 * torch.exp((-1 * (x - _mu)) / _b);
            }

            private torch.Tensor PPF(torch.Tensor q)
            {
                if (q.less_equal(0.5).Equals(true))
                {
                    return _mu + _b * torch.log(2 * q);
                }
                else
                {
                    return _mu - _b * torch.log(2-2*q);
                }
            }
            private torch.Tensor LogL(torch.Tensor x)
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
            private torch.Tensor _mu {get; set; }
        }
    }
}
