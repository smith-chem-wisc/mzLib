using TorchSharp;

namespace MachineLearning
{
    /// <summary>
    /// Abstract class for one dimensional input size datasets for use in DataLoader.
    /// </summary>
    public class TorchDataset : torch.utils.data.Dataset
    {
        public List<(torch.Tensor, double)> Data { get; set; }

        public override long Count => this.Data.Count;

        public TorchDataset(List<(torch.Tensor, double)> data)
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
