using Newtonsoft.Json;

namespace Predictions.Koina.SupportedModels.Prosit2020IntensityHCD
{
    public class ResponseJSONStruct
    {
        [JsonProperty("id")]
        public string Id { get; set; }

        [JsonProperty("model_name")]
        public string ModelName { get; set; }

        [JsonProperty("model_version")]
        public string ModelVersion { get; set; }

        [JsonProperty("parameters")]
        public Dictionary<string, object> Parameters { get; set; }

        [JsonProperty("outputs")]
        public List<OutputJSONStruct> Outputs { get; set; }
    }
}
