using Newtonsoft.Json;

namespace Predictions.Koina.SupportedModels.Prosit2020IntensityHCD
{
    public class OutputJSONStruct
    {
        [JsonProperty("name")]
        public string Name { get; set; } = string.Empty;

        [JsonProperty("datatype")]
        public string Datatype { get; set; } = string.Empty;

        [JsonProperty("shape")]
        public List<int> Shape { get; set; } = new();

        [JsonProperty("data")]
        public List<object> Data { get; set; } = new();
    }
}
