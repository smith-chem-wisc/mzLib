using Easy.Common.Extensions;
using Proteomics.PSM;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using TorchSharp;
using TorchSharp.Modules;

namespace Proteomics.RetentionTimePrediction
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
    public class Chronologer : torch.nn.Module<torch.Tensor, torch.Tensor>
    {
        public Chronologer() : this(Path.Combine(AppDomain.CurrentDomain.BaseDirectory,
            "RetentionTimePrediction",
            "Chronologer_20220601193755_TorchSharp.dat"))
        {

        }

        /// <summary>
        /// Initializes a new instance of the Chronologer model class with pre-trained weights from the paper
        /// Deep learning from harmonized peptide libraries enables retention time prediction of diverse post
        /// translational modifications paper (https://github.com/searlelab/chronologer).
        /// Eval mode is set to true and training mode is set to false by default.
        ///
        /// Please use .Predict() for using the model, not .forward(). 
        /// </summary>
        /// <param name="weightsPath"></param>
        /// <param name="evalMode"></param>
        public Chronologer(string weightsPath = null, bool evalMode = true) : base(nameof(Chronologer))
        {
            RegisterComponents();
            if (weightsPath != null)
                LoadWeights(weightsPath);//loads weights from the file

            if (evalMode)
            {
                this.eval();
                this.train(false);
            }
        }

        //Without pre-trained weights
        public Chronologer(bool trained = false) : base(nameof(Chronologer))
        {
            RegisterComponents();
        }

        /// <summary>
        /// Do not use for inferring. Use .Predict() instead.
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public override torch.Tensor forward(torch.Tensor x)
        {
            var input = seq_embed.forward(x).transpose(1, -1);

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

        /// <summary>
        /// Loads pre-trained weights from the file Chronologer_20220601193755_TorchSharp.dat.
        /// </summary>
        /// <param name="weightsPath"></param>
        private void LoadWeights(string weightsPath)
        {
            //load weights from the file
            this.load(weightsPath, true);
        }

        /// <summary>
        /// Predicts the retention time of the input peptide sequence. The input must be a torch.Tensor of shape (1, 52).
        /// </summary>
        /// <param name="input"></param>
        /// <returns></returns>
        public torch.Tensor Predict(torch.Tensor input)
        {
            return this.call(input);
        }

        public void Train(string savingPath, List<PsmFromTsv> trainingData, Dictionary<(char, string),
                int> dictionary, double validationFraction = 0.2,
            int seed = 2447, string device = "cpu", string startModelPath = null,
            double intialBatchScaler = 1.0, int batchSize = 64, int epochs = 100)
        {
            if (this.training.Equals(false))
            {
                this.train(true);
            }

            //loads pre-trained weights if startModelPath is not null
            if (startModelPath != null)
            {
                this.load(startModelPath, true);
            }
            //Change input embedding layer to the size of the dictionary
            this.seq_embed = torch.nn.Embedding(dictionary.Count, 64, 0).to(DeviceType.CUDA);
            RegisterComponents();
            this.to(DeviceType.CUDA); //moves the model to GPU


            var (train, test) =
                TrainChronologer.RetentionTimeToTensorDatabase(trainingData, seed, validationFraction, dictionary);

            var trainDataSet = new CustomDataset(train);
            var testDataSet = new CustomDataset(test);

            var trainLoader = torch.utils.data.DataLoader(trainDataSet, batchSize, shuffle: true, drop_last: true);
            var testLoader = torch.utils.data.DataLoader(testDataSet, batchSize, shuffle: true, drop_last: true);
            var lossFunction = torch.nn.L1Loss();
            // var lossFunction = new TrainChronologer.LogLLoss(52, 
            //     TrainChronologer.ReturnDistribution.Laplace, 0.01);

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

                    batchX = batchX.to(DeviceType.CUDA);
                    batchY = batchY.to(DeviceType.CUDA);
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
                            torch.tensor(new[] { (float)y }).to(DeviceType.CUDA));

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
            this.save(savingPath);
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
