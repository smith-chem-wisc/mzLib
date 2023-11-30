using MathNet.Numerics.Statistics;
using Proteomics.PSM;
using TorchSharp;
using TorchSharp.Modules;

namespace MachineLearning
{
    /// <summary>
    /// Abstract class for one dimensional input size TorchSharp models.
    /// </summary>
    public abstract class DeepTorch : torch.nn.Module<torch.Tensor, torch.Tensor>
    {
        public DeepTorch(string preTrainedModelPath = null, bool evalMode = true,
            DeviceType device = DeviceType.CPU)
            : base(nameof(DeepTorch))
        {
            RegisterComponents();

            if (preTrainedModelPath != null)
                LoadPreTrainedModel(preTrainedModelPath);

            if (evalMode)
                EvaluationMode();
            else
                TrainingMode();

            this.to(device);
        }

        /// <summary>
        /// Defines the computation performed at every call. Do not use for inferring, use .Predict() instead.
        /// </summary>
        /// <param name="input"></param>
        /// <returns></returns>
        public abstract override torch.Tensor forward(torch.Tensor input);

        public virtual torch.Tensor Predict(torch.Tensor input)
        {
            return this.call(input);
        }

        /// <summary>
        /// Defines the behavior of the module during training phase.
        /// </summary>
        public abstract void Train(string modelSavingPath, List<PsmFromTsv> trainingData,
            Dictionary<(char, string), int> dictionary, DeviceType device, float validationFraction,
            float testingFraction, int batchSize, int epochs, int patience);

        /// <summary>
        /// Defines the behavior of the model when loading pre-trained weights.
        /// </summary>
        /// <param name="weightsPath"></param>
        /// <param name="strict"></param>
        public virtual void LoadPreTrainedModel(string weightsPath, bool strict = true)
        {
            this.load(weightsPath, strict);
        }

        /// <summary>
        /// Stipulates an early stop condition to prevent overfitting.
        /// </summary>
        protected virtual bool EarlyStop(double score, double currentBestScore, int currentPatience,
            out double bestScore, out int patience)
        {
            if (score < currentBestScore && currentPatience > 0)
            {
                bestScore = score;
                patience = currentPatience - 1;
                return false;
            }

            bestScore = currentBestScore;
            patience = 0;
            return true;

        }

        /// <summary>
        /// Generates a checkpoint to save the model at indicated step.
        /// </summary>
        protected virtual void Checkpoint(string checkPointPath, int epoch)
        {

            Directory.CreateDirectory(checkPointPath);
            Directory.CreateDirectory(checkPointPath + "/Model" + epoch);
            SaveModel(checkPointPath + "/Model" + epoch + "checkpointModel");
        }

        protected virtual void ModelPerformance(string savingPath, double bestScore, double accuracy,
            List<double> predictions,
            List<double> labels)
        {
            using (var writter = new System.IO.StreamWriter(savingPath + ".txt"))
            {
                writter.WriteLine($"Best score: {bestScore}");
                writter.WriteLine($"Accuracy: {accuracy}");
                writter.WriteLine($"R^2: " + Correlation.Pearson(predictions, labels));
                writter.WriteLine($"Predictions: {predictions}");
                writter.WriteLine($"Labels: {labels}");
            }
        }
        /// <summary>
        /// Learning Rate Decay.
        /// </summary>
        protected virtual torch.optim.lr_scheduler.LRScheduler _scheduler =>
            torch.optim.lr_scheduler.StepLR(new Adam(this.parameters()), 25, 0.1);

        /// <summary>
        /// Sets model into training mode.
        /// </summary>
        public virtual void TrainingMode()
        {
            this.train(true);
        }

        /// <summary>
        /// Sets model into evaluation mode.
        /// </summary>
        public virtual void EvaluationMode()
        {
            this.eval();
            this.train(false);
        }

        /// <summary>
        /// Saves the trained weights as a .dat file (TorchSharp format).
        /// </summary>
        /// <param name="modelSavingPath"></param>
        public virtual void SaveModel(string modelSavingPath)
        {
            EvaluationMode();

            if (modelSavingPath.EndsWith(".dat"))
                this.save(modelSavingPath);
            else
                this.save(modelSavingPath + ".dat");
        }

        protected abstract double Validate(DataLoader? validationDataLoader,
            torch.nn.Module<torch.Tensor, torch.Tensor, torch.Tensor> criterion, DeviceType device);
        protected abstract (double, List<double>, List<double>) Test(DataLoader? testingDataLoader,
                       torch.nn.Module<torch.Tensor, torch.Tensor, torch.Tensor> criterion, DeviceType device);
        public abstract TorchDataset? TrainingDataset { get; set; }
        public abstract TorchDataset? TestingDataset { get; set; }
        public abstract TorchDataset? ValidationDataset { get; set; }
        public abstract torch.Tensor Tensorize(object toTensoize);

        // public abstract Dictionary<string, torch.Tensor> GetTensor(long index);

        // sets nullable filed Dataset to a new instance of a TorchDataset
        public abstract void CreateDataSet(List<PsmFromTsv> data, float validationFraction, float testingFraction, int batchSize);

        protected abstract void CreateDataLoader(int batchSize);
    }
}
