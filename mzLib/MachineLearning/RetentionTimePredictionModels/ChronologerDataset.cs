using Proteomics.PSM;
using TorchSharp;

namespace MachineLearning.RetentionTimePredictionModels
{
    public class ChronologerDataset : TorchDataset
    {
        public ChronologerDataset(List<(torch.Tensor, double)> data) : base(data)
        {
        }

        public override Dictionary<string, torch.Tensor> GetTensor(long index)
        {
            return new Dictionary<string, torch.Tensor>()
            {
                {"Retention Time on File", Data[(int) index].Item2},
                {"Encoded Sequence", Data[(int) index].Item1}
            };
        }

        /// <summary>
        /// Converts PsmFromTsv list into tensors and splits the data into training and test sets.
        /// </summary>
        /// <param name="dataFiles"></param>
        /// <param name="validationFraction"></param>
        /// <param name="dictionary"></param>
        /// <returns></returns>
        public static (List<(torch.Tensor, double)>, List<(torch.Tensor, double)>) ChronologerRetentionTimeToTensorDatabaseSplit(
            List<PsmFromTsv> dataFiles, double validationFraction, Dictionary<(char, string), int> dictionary)
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

                    var tensor = Chronologer.Tensorize(dataFile, dictionary);

                    if (tensor.Equals(torch.ones(1, 52, torch.ScalarType.Int64)))
                        continue;

                    if (tensor[0][0].item<Int64>().Equals((Int64)38))
                    {
                        allData.Add((tensor, db.RetentionTime.Value)); //todo: add encoded sequence tensor
                        sources.Add(db.FileNameWithoutExtension);
                    }
                }
            }

            var trainingSet = allData.Take((int)(allData.Count * (1 - validationFraction))).ToList();
            var testSet = allData.Skip((int)(allData.Count * (1 - validationFraction))).ToList();

            trainTestDb["train"] = trainingSet;
            trainTestDb["test"] = testSet;

            return (trainTestDb["train"], trainTestDb["test"]);
        }
    }
}
