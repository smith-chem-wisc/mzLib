namespace PredictionClients.Koina.Interfaces
{
    public interface IKoinaModelIO
    {
        public string ModelName { get; }
        public int MaxBatchSize { get; }
        public List<Dictionary<string, object>> ToBatchedRequests();
        public Task RunInferenceAsync();
    }
}
