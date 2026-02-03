namespace PredictionClients.Koina.Interfaces
{
    // TODO: This interface needs to be reworked/refined. Maybe turn into a master abstract class instead?
    public interface IKoinaModelIO
    {
        public string ModelName { get; }
        public int MaxBatchSize { get; }
        public Task RunInferenceAsync();
    }
}
