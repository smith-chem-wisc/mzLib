using Easy.Common.Extensions;
using Proteomics.PSM;
using System.Data;
using System.Diagnostics;
using MathNet.Numerics.Statistics;
using SkiaSharp;
using TorchSharp;
using TorchSharp.Modules;
using TorchSharp.Utils;
using static UsefulProteomicsDatabases.DictionaryBuilder;

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
        public DataLoader? ValidationDataLoader { get; set; }
        public DataLoader? TestDataLoader { get; set; }
        public override TorchDataset? TrainingDataset { get; set; }
        public override TorchDataset? ValidationDataset { get; set; }
        public override TorchDataset? TestingDataset { get; set; }
        public Dictionary<(char TargetAA, string ModificationId), int> ModificationDictionary { get; set; }

        public Chronologer(Dictionary<(char TargetAA, string ModificationId), int> dictionary,
            string preTrainedModelPath = null, bool evalMode = true, TypeOfDictionary dict = TypeOfDictionary.Chronologer)
            : base(preTrainedModelPath, evalMode)
        {
            ModificationDictionary = GetChronologerDictionary(dict);
        }

        public override void CreateDataSet(List<PsmFromTsv> data, float validationFraction, float testingFraction, int batchSize)
        {
            // if (data.Select(x => x).ToList() is List<PsmFromTsv> psmList)
            // {
                var trainTestDb = new Dictionary<string, List<(torch.Tensor, double)>>()
                {
                    { "train", new List<(torch.Tensor, double)>() },
                    {"validation", new List<(torch.Tensor, double)>()},
                    { "test", new List<(torch.Tensor, double)>() }
                };

                var allData = new List<(torch.Tensor, double)>();

                var sources = new HashSet<string>();

                foreach (var dataFile in data)
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

                allData = allData.Randomize().ToList();

                var trainingSet = allData.Take((int)(allData.Count * (1 - validationFraction))).ToList();

                var validationSet = allData.Skip((int)(allData.Count * (1 - validationFraction)))
                    .Take((int)(allData.Count * (1 - testingFraction))).ToList();

                var testSet = allData.Skip((int)(allData.Count * (1 - validationFraction))).ToList();

                trainTestDb["train"] = trainingSet;
                trainTestDb["validation"] = validationSet;
                trainTestDb["test"] = testSet;

                TrainingDataset = new TorchDataset(trainTestDb["train"], "Encoded Sequence", "Retention Time on File");
                ValidationDataset = new TorchDataset(trainTestDb["validation"], "Encoded Sequence", "Retention Time on File");
                TestingDataset = new TorchDataset(trainTestDb["test"], "Encoded Sequence", "Retention Time on File");

                CreateDataLoader(batchSize);
            }
            // else
            // {
            //     throw new System.ArgumentException("Object is not of type List<PsmFromTsv>", nameof(data));
            // }

        protected override void CreateDataLoader(int batchSize)
        {
            TrainDataLoader = new DataLoader(TrainingDataset, batchSize, shuffle: true, drop_last: true);
            ValidationDataLoader = new DataLoader(ValidationDataset, batchSize, shuffle: true, drop_last: true);
            TestDataLoader = new DataLoader(TestingDataset, batchSize, shuffle: true, drop_last: true);
        }

        public override void Train(string modelSavingPath, List<PsmFromTsv> trainingData,
            Dictionary<(char, string), int> dictionary, DeviceType device, float validationFraction = 0.1f, float testingFraction = 0.1f,
            int batchSize = 64, int epochs = 50, int patience = 5)
        {
            double bestScore = double.MaxValue;
            if (this.training.Equals(false))
            {
                this.train(true);
            }

            //Change input embedding layer to the size of the dictionary
            this.seq_embed = torch.nn.Embedding(55, 64, 0).to(device);
            RegisterComponents();
            this.to(device); //moves the model to GPU


            //training DataLoader
            var trainLoader = TrainDataLoader;
            //Validation DataLoader todo: add validation set
            var validationLoader = ValidationDataLoader;
            //Testing DataLoader
            var testLoader = TestDataLoader;

            var lossFunction = torch.nn.L1Loss();

            var parameters = new List<Parameter>();
            this.parameters().ForEach(x => parameters.Add(x));
            // lossFunction.parameters().ForEach(x => parameters.Add(x)); //todo might not be needed

            //Optimizer
            var optimizer = torch.optim.Adam(parameters, 0.001);

            var scores = new List<float>();
            
            //training loop
            for (int i = 0; i < epochs; i++)
            {
                TrainingMode();

                var testLoss = 0.0;

                // Debug.WriteLine($"Epoch {i + 1} of {epochs}");
                var epochLoss = 0.0;
                var runningLoss = 0.0;

                //Tranining
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

                        var output = this.forward(x);
                        var loss = lossFunction.forward(output[0],
                            torch.tensor(new[] { (float)y }).to(device));

                        optimizer.zero_grad();
                        loss.backward();
                        optimizer.step();

                        //statistics per batch
                        runningLoss += loss.item<float>() * batchSize;
                        scores.Add(loss.item<float>());
                    }
                }

                //Validation
                var validationScore = Validate(validationLoader, lossFunction, device);

                //check for ealy stopping
                if (EarlyStop(validationScore, bestScore, patience, out bestScore, out patience))
                    break;

                //Scheduler step
                _scheduler.step();

                //statistics per epoch
                epochLoss = runningLoss / TrainingDataset.Count;
                Debug.WriteLine($"Epoch {i + 1} loss: {epochLoss}");
            }

            //Evaluate model performance
            var (accuracy, predictions, labels) = Test(TestDataLoader, lossFunction, device);

            //saving the model
            SaveModel(modelSavingPath);

            //Model Stats
            ModelPerformance(modelSavingPath, bestScore, accuracy, predictions, labels);
        }

        protected override double Validate(DataLoader? validationDataLoader, torch.nn.Module<torch.Tensor, torch.Tensor, torch.Tensor> criterion,
            DeviceType device)
        {
            EvaluationMode();
            double totalLoss = 0.0;
            using (torch.no_grad())
            {
                foreach (var batch in validationDataLoader)
                {
                    var batchX = batch["Encoded Sequence"];
                    var batchY = batch["Retention Time on File"];

                    batchX = batchX.to(device);
                    batchY = batchY.to(device);

                    for(int i = 0; i < batchX.size(0); i++)
                    {
                        var output = this.Predict(batchX[i]);
                        var loss = criterion.forward(output[0], batchY[i]);

                        totalLoss += loss.item<double>() * batchX.size(0);
                    }
                }
            }
            var averageLoss = totalLoss / validationDataLoader.Count;

            return averageLoss;
        }

        protected override (double, List<double>, List<double>) Test(DataLoader? testingDataLoader,
            torch.nn.Module<torch.Tensor, torch.Tensor, torch.Tensor> criterion,
            DeviceType device)
        {
            EvaluationMode();
            double totalLoss = 0.0;
            List<double> predictions = new List<double>();
            List<double> labels = new List<double>();
            using (torch.no_grad())
            {
                foreach (var batch in testingDataLoader)
                {
                    var batchX = batch["Encoded Sequence"];
                    var batchY = batch["Retention Time on File"];

                    batchX = batchX.to(device);
                    batchY = batchY.to(device);

                    for (int i = 0; i < batchX.size(0); i++)
                    {
                        var output = this.Predict(batchX[i]);
                        var loss = criterion.forward(output[0], batchY[i]);

                        predictions.Add(output.item<float>());
                        labels.Add(batchY.item<double>());
                        totalLoss += loss.item<double>() * batchX.size(0);

                        predictions.Add(output.item<float>());
                        labels.Add(batchY.item<double>());
                        totalLoss += loss.item<double>() * batchX.size(0);
                    }
                }
            }
            var averageLoss = totalLoss / testingDataLoader.Count;

            return (averageLoss, predictions, labels);
        }

        public override torch.Tensor Tensorize(object toTensoize)
        {
            if (toTensoize is PsmFromTsv psm)
            {
                return Tensorize(psm, ModificationDictionary);
            }
            else
            {
                throw new System.ArgumentException("Object is not supported ", nameof(toTensoize));
            }
        }

        /// <summary>
        /// Tensorizes PsmFromTsv objects for use in the Chronologer model.
        /// </summary>
        /// <param name="psm"></param>
        /// <param name="dictionary"></param>
        /// <param name="qValueFiltering"></param>
        /// <param name="ambiguityLevel"></param>
        /// <returns></returns>
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

                //N terminal mods set mod to true

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
                }

                tensor[0][tensorCounter] = 44; //N-terminus

                Debug.WriteLine(tensor.ToString(TensorStringStyle.Julia));
                return tensor;
            }

            //if sequence is longer than 50, return a tensor of ones, quick way to remove it later from the dataset
            return torch.ones(1, 52, torch.ScalarType.Int64);
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
