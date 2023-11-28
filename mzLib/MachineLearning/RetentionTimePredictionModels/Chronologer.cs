using System.Data;
using Easy.Common.Extensions;
using Proteomics.PSM;
using Proteomics.RetentionTimePrediction;
using System.Diagnostics;
using TorchSharp;
using TorchSharp.Modules;
using static UsefulProteomicsDatabases.ChronologerDictionary;
using System.Collections.Generic;

namespace MachineLearning.RetentionTimePredictionModels
{
    /// <summary>
    /// Chronologer is a deep learning model for highly accurate prediction of peptide C18 retention times (reported in % ACN).
    /// Chronologer was trained on a new large harmonized database of > 2.6 million retention time observations
    /// (2.25M unique peptides) constructed from 11 community datasets
    /// and natively supports prediction of 17 different modification types.
    /// With only a few observations of a new modification type (> 10 peptides),
    /// Chronologer can be easily re-trained to predict up to 10 user supplied modifications.
    ///
    /// Damien Beau Wilburn, Ariana E. Shannon, Vic Spicer, Alicia L. Richards, Darien Yeung, Danielle L. Swaney, Oleg V. Krokhin, Brian C. Searle
    /// bioRxiv 2023.05.30.542978; doi: https://doi.org/10.1101/2023.05.30.542978
    /// 
    /// https://github.com/searlelab/chronologer
    ///
    /// Licensed under Apache License 2.0
    /// 
    /// </summary>
    public class Chronologer : DeepTorch
    {
        public DataLoader? TrainDataLoader { get; set; }
        public DataLoader? TestDataLoader { get; set; }
        public override TorchDataset? TrainingDataset { get; set; }
        public override TorchDataset? TestingDataset { get; set; }
        public Dictionary<(char TargetAA, string ModificationId), int> ModificationDictionary { get; set; }

        public Chronologer(Verbosity trainingVerbosity, Dictionary<(char TargetAA, string ModificationId), int> dictionary, string preTrainedModelPath = null,
            bool evalMode = true, TypeOfDictionary dict = TypeOfDictionary.Chronologer) : base(trainingVerbosity, preTrainedModelPath, evalMode)
        {
           ModificationDictionary = GetChronologerDictionary(dict);
        }

        public override void CreateDataSet(List<object> data, double validationFraction, int batchSize)
        {
            if (data.Select(x => x) is List<PsmFromTsv> psmList)
            {
                var trainTestDb = new Dictionary<string, List<(torch.Tensor, double)>>()
                {
                    { "train", new List<(torch.Tensor, double)>() },
                    { "test", new List<(torch.Tensor, double)>() }
                };

                var allData = new List<(torch.Tensor, double)>();

                var sources = new HashSet<string>();

                foreach (var dataFile in psmList)
                {
                    if (dataFile.DecoyContamTarget.Equals("T"))
                    {
                        var db =
                            (dataFile.FileNameWithoutExtension, dataFile.RetentionTime,
                                dataFile.BaseSeq); //base seq for the moment

                        var tensor = Tensorize(dataFile);

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

                TrainingDataset = new TorchDataset(trainTestDb["train"], "Encoded Sequence", "Retention Time on File");
                TestingDataset = new TorchDataset(trainTestDb["test"], "Encoded Sequence", "Retention Time on File");

                // (_trainingData, _testData) = (trainTestDb["train"], trainTestDb["test"]);

                CreateDataLoader(batchSize);
            }
            else
            {
                throw new System.ArgumentException("Object is not of type List<PsmFromTsv>", nameof(data));
            }
        }

        protected override void CreateDataLoader(int batchSize)
        {
            TrainDataLoader = new DataLoader(TrainingDataset, batchSize, shuffle: true, drop_last: true);
            TestDataLoader = new DataLoader(TestingDataset, batchSize, shuffle: true, drop_last: true);
        }

        public override void Train(string modelSavingPath, List<PsmFromTsv> trainingData,
            Dictionary<(char, string), int> dictionary, DeviceType device, float validationFraction, int batchSize, int epochs)
        {
            if (this.training.Equals(false))
            {
                this.train(true);
            }

            //Change input embedding layer to the size of the dictionary
            this.seq_embed = torch.nn.Embedding(dictionary.Count, 64, 0).to(device);
            RegisterComponents();
            this.to(device); //moves the model to GPU

            var (train, test) =
                ChronologerDataset.ChronologerRetentionTimeToTensorDatabaseSplit(trainingData, validationFraction, dictionary);

            var trainDataSet = new CustomDataset(train);
            var testDataSet = new CustomDataset(test);

            var trainLoader = torch.utils.data.DataLoader(trainDataSet, batchSize, shuffle: true, drop_last: true);
            var testLoader = torch.utils.data.DataLoader(testDataSet, batchSize, shuffle: true, drop_last: true);

            var lossFunction = torch.nn.L1Loss();

            var parameters = new List<Parameter>();
            this.parameters().ForEach(x => parameters.Add(x));
            lossFunction.parameters().ForEach(x => parameters.Add(x));

            var optimizer = torch.optim.Adam(parameters, 0.001);
            var scheduler = torch.optim.lr_scheduler.StepLR(optimizer, 25, 0.1, verbose: true);
            //training loop
            var scores = new List<float>();

            for (int i = 0; i < epochs; i++)
            {
                var testLoss = 0.0;

                // Debug.WriteLine($"Epoch {i + 1} of {epochs}");
                var epochLoss = 0.0;
                var runningLoss = 0.0;
                foreach (var batch in trainLoader)
                {
                    var batchX = batch["Encoded Sequence"];
                    var batchY = batch["Retention Time on File"];

                    batchX = batchX.to(device);
                    batchY = batchY.to(device);
                    for (int j = 0; j < batchX.size(0); j++)
                    {
                        var x = batchX[j];
                        double y = batchY[j].item<double>();
                        //
                        // Debug.WriteLine(x.ToString(TensorStringStyle.Julia));
                        // Debug.WriteLine(y);
                        // Debug.WriteLine(batchPass.ToString(TensorStringStyle.Julia));

                        var output = this.forward(x);
                        var loss = lossFunction.forward(output[0],
                            torch.tensor(new[] { (float)y }).to(device));

                        // Debug.WriteLine(loss.item<float>());
                        optimizer.zero_grad();
                        loss.backward();
                        optimizer.step();

                        //statistics per batch
                        runningLoss += loss.item<float>() * batchSize;
                        scores.Add(loss.item<float>());
                    }
                }
                scheduler.step();
                //statistics per epoch
                epochLoss = runningLoss / train.Count;
                Debug.WriteLine($"Epoch {i + 1} loss: {epochLoss}");
            }

            //saving the model
            this.eval();
            this.train(false);
            this.save(modelSavingPath);
        }

        public override torch.Tensor forward(torch.Tensor input)
        {
            input = seq_embed.forward(input).transpose(1, -1);

            var residual = input.clone();
            input = conv_layer_1.forward(input); //renet_block
            input = norm_layer_1.forward(input);
            input = relu.forward(input);
            input = conv_layer_2.forward(input);
            input = norm_layer_2.forward(input);
            input = relu.forward(input);
            input = term_block.forward(input);
            input = residual + input;
            input = relu.forward(input);

            residual = input.clone();
            input = conv_layer_4.forward(input);//renet_block
            input = norm_layer_4.forward(input);
            input = relu.forward(input);
            input = conv_layer_5.forward(input);
            input = norm_layer_5.forward(input);
            input = relu.forward(input);
            input = term_block.forward(input);
            input = residual + input;
            input = relu.forward(input);

            residual = input.clone();
            input = conv_layer_7.forward(input);//renet_block
            input = norm_layer_7.forward(input);
            input = term_block.forward(input);
            input = relu.forward(input);
            input = conv_layer_8.forward(input);
            input = norm_layer_8.forward(input);
            input = relu.forward(input);
            input = term_block.forward(input);
            input = residual + input;
            input = relu.forward(input);

            input = dropout.forward(input);
            input = flatten.forward(input);
            input = output.forward(input);

            return input;
        }

        protected override void EarlyStop()
        {
            throw new NotImplementedException();
        }

        protected override void RegisterComponents()
        {
            throw new NotImplementedException();
        }


        public override torch.Tensor Tensorize(object toTensoize)
        {
            if (toTensoize is PsmFromTsv psm)
            {
                return Tensorize(psm, ModificationDictionary);
            }
            else
            {
                throw new System.ArgumentException("Object is not of type PsmFromTsv", nameof(toTensoize));
            }
        }

        


        public static torch.Tensor Tensorize(PsmFromTsv psm, Dictionary<(char, string), int> dictionary,
            double qValueFiltering = 0.01, string ambiguityLevel = "1")
        {
            var fullSequence = psm.FullSequence.Split(new[] { '[', ']' })
                 .Where(x => !x.Equals("")).ToArray();

            //section to remoce type of mod from the sequence
            //i.e. Common Fixed:Carbamidomethyl and only stay with the target aa after the mod
            for (int i = 0; i < fullSequence.Length; i++)
            {
                if (fullSequence[i].Contains(" "))
                {
                    var tempSeq = fullSequence[i].Split(':');
                    fullSequence[i] = tempSeq[1];
                }
            }

            if (psm.BaseSeq.Length <= 50 &&
                psm.QValue <= qValueFiltering &&
                psm.AmbiguityLevel == 1.ToString()) //Chronologer only takes sequences of length 50 or less
            {
                var tensor = torch.zeros(1, 52, torch.ScalarType.Int64);

                tensor[0][0] = 38; //C-terminus
                var tensorCounter = 1; //skips the first element which is the C-terminus in the tensor
                char modID = ' '; //takes the target aa from inside the loop to hold it for the next loop
                bool mod = false; //assumes that the next loop is not a mod 

                if (fullSequence[0].Contains(" "))
                {
                    mod = true;
                    modID = fullSequence[0][fullSequence[0].Length - 1];
                }

                foreach (var subString in fullSequence)
                {
                    //if mod, enter
                    if (mod)
                    {
                        var key = (char.ToUpper(modID), subString);
                        if (dictionary.ContainsKey(key))
                        {
                            tensor[0][tensorCounter] = dictionary[key];
                            tensorCounter = tensorCounter + 1;
                        }
                        mod = false; //next iteration is not a mod
                        continue;
                    }

                    if (!subString.Contains(" "))
                    {
                        //without mods
                        for (int i = 0; i < subString.Length - 1; i++)
                        {
                            tensor[0][tensorCounter] = dictionary[(char.ToUpper(subString[i]), "")];
                            tensorCounter = tensorCounter + 1;
                        }
                        //next loop will be a mod
                        mod = true;

                        //save target aa for next loop
                        modID = subString[subString.Length - 1];
                    }
                    else
                    {
                        for (int i = 0; i < subString.Length; i++)
                        {
                            tensor[0][tensorCounter] = dictionary[(char.ToUpper(subString[i]), "")];
                            tensorCounter = tensorCounter + 1;
                        }
                        //next loop will be a mod
                        mod = true;

                    }
                }

                tensor[0][tensorCounter] = 44; //N-terminus

                Debug.WriteLine(tensor.ToString(TensorStringStyle.Julia));
                return tensor;
            }

            //if sequence is longer than 50, return a tensor of ones, quick way to remove it later from the dataset
            return torch.ones(1, 52, torch.ScalarType.Int64);
        }

        //All Modules (shortcut modules are for loading the weights only)
        private Embedding seq_embed = torch.nn.Embedding(55, 64, 0);
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_1 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 1);
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_2 = torch.nn.Conv1d(64, 64, 7, Padding.Same, dilation: 1);
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_3 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 1); //shortcut
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_4 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 2);
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_5 = torch.nn.Conv1d(64, 64, 7, Padding.Same, dilation: 2);
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_6 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 2); //shortcut
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_7 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 3);
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_8 = torch.nn.Conv1d(64, 64, 7, Padding.Same, dilation: 3);
        private torch.nn.Module<torch.Tensor, torch.Tensor> conv_layer_9 = torch.nn.Conv1d(64, 64, 1, Padding.Same, dilation: 3); //shortcut
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_1 = torch.nn.BatchNorm1d(64);
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_2 = torch.nn.BatchNorm1d(64);
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_3 = torch.nn.BatchNorm1d(64); //shortcut
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_4 = torch.nn.BatchNorm1d(64);
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_5 = torch.nn.BatchNorm1d(64);
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_6 = torch.nn.BatchNorm1d(64); //shortcut
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_7 = torch.nn.BatchNorm1d(64);
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_8 = torch.nn.BatchNorm1d(64);
        private torch.nn.Module<torch.Tensor, torch.Tensor> norm_layer_9 = torch.nn.BatchNorm1d(64);
        private torch.nn.Module<torch.Tensor, torch.Tensor> term_block = torch.nn.Identity();
        private torch.nn.Module<torch.Tensor, torch.Tensor> relu = torch.nn.ReLU(true);
        private torch.nn.Module<torch.Tensor, torch.Tensor> dropout = torch.nn.Dropout(0.01);
        private torch.nn.Module<torch.Tensor, torch.Tensor> flatten = torch.nn.Flatten(1);
        private torch.nn.Module<torch.Tensor, torch.Tensor> output = torch.nn.Linear(52 * 64, 1);
    }
}
