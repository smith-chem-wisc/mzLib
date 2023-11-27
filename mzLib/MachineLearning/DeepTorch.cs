using Proteomics.PSM;
using TorchSharp;

namespace MachineLearning
{
    /// <summary>
    /// Abstract class for one dimensional input size TorchSharp models.
    /// </summary>
    public abstract class DeepTorch : torch.nn.Module<torch.Tensor, torch.Tensor>
    {
        public Verbosity TrainingVerbosity { get; set; } //todo: Figure out a way to manipulate level of verbosity

        public DeepTorch(Verbosity trainingVerbosity, string preTrainedModelPath = null, bool evalMode = true)
            : base(nameof(DeepTorch))
        {
            RegisterComponents();

            TrainingVerbosity = trainingVerbosity;

            if (preTrainedModelPath != null)
                LoadPreTrainedModel(preTrainedModelPath);

            if (evalMode)
                EvaluationMode();
            else
                TrainingMode();

        }

        /// <summary>
        /// Defines the computation performed at every call.
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
        public abstract void Train();

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
        protected virtual void EarlyStop()
        {

        }

        /// <summary>
        /// Generates a checkpoint to save the model at indicated step.
        /// </summary>
        protected virtual void Checkpoint()
        {

        }

        /// <summary>
        /// Learning Rate Decay setter.
        /// </summary>
        protected virtual void LearningRateScheduler()
        {

        }

        /// <summary>
        /// Sets model into training mode.
        /// </summary>
        public void TrainingMode()
        {
            this.train(true);
        }

        /// <summary>
        /// Sets model into evaluation mode.
        /// </summary>
        public void EvaluationMode()
        {
            this.eval();
            this.train(false);
        }

        /// <summary>
        /// Tensorizes the data for use in the model. Use the keyword new to hide the base class method and
        /// implement your own method.
        /// </summary>
        /// <returns></returns>
        /// <exception cref="NotImplementedException"></exception>
        public static torch.Tensor Tensorize(PsmFromTsv psm, Dictionary<(char, string), int> dictionary,
            TensorEncodingSchema tensorEncoding, double qValueFiltering = 0.01, string ambiguityLevel = "1")
        {
            throw new NotImplementedException();
        }
    }

    public class asdasdas : DeepTorch
    {
        public asdasdas(Verbosity trainingVerbosity, string preTrainedModelPath = null, bool evalMode = true)
            : base(trainingVerbosity, preTrainedModelPath, evalMode) { }

        public override torch.Tensor forward(torch.Tensor input)
        {
            throw new NotImplementedException();
        }

        public override void Train()
        {
            throw new NotImplementedException();
        }
    }
}
