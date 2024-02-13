using System;
using System.Collections.Generic;
using TorchSharp;

namespace Proteomics.RetentionTimePrediction.Chronologer
{
    public class CustomDataset : torch.utils.data.Dataset
    {
        public List<(torch.Tensor, double)> Data { get; set; }
        public override long Count => Data.Count;
        public CustomDataset(List<(torch.Tensor, double)> data)
        {
            Data = data;
        }
        public override Dictionary<string, torch.Tensor> GetTensor(long index)
        {
            return new Dictionary<string, torch.Tensor>()
            {
                {"Retention Time on File", Data[(int) index].Item2},
                {"Encoded Sequence", Data[(int) index].Item1}
            };
        }
    }
}
