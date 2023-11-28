using TorchSharp;

namespace MachineLearning
{
    /// <summary>
    /// Abstract class for one dimensional input size datasets for use in DataLoader.
    /// </summary>
    public class TorchDataset : torch.utils.data.Dataset
    {
        private string _targetName;
        private string _labelName;
        public List<(torch.Tensor, double)> Data { get; set; }

        public override long Count => this.Data.Count;

        public TorchDataset(List<(torch.Tensor, double)> data, string targetName = "Target", string labelName = "Label")
        {
            Data = data;
            _targetName = targetName;
            _labelName = labelName;
        }

        public override Dictionary<string, torch.Tensor> GetTensor(long index)
        {
            return new Dictionary<string, torch.Tensor>()
                {
                    {_labelName, Data[(int) index].Item2},
                    {_targetName, Data[(int) index].Item1}
                };
        }
    }
}
