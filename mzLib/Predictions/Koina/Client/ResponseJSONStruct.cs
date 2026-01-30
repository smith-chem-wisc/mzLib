using Newtonsoft.Json;

namespace Predictions.Koina.Client
{
    public class ResponseJSONStruct
    {
        [JsonProperty("id")]
        public string Id { get; set; } = string.Empty;

        [JsonProperty("model_name")]
        public string ModelName { get; set; } = string.Empty;

        [JsonProperty("model_version")]
        public string ModelVersion { get; set; } = string.Empty;

        [JsonProperty("parameters")]
        public Dictionary<string, object> Parameters { get; set; } = new();

        [JsonProperty("outputs")]
        public List<OutputJSONStruct> Outputs { get; set; } = new();
    }
}
