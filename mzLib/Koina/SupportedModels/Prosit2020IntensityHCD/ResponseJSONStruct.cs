using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Koina.SupportedModels.Prosit2020IntensityHCD
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
        public List<Koina.SupportedModels.Prosit2020IntensityHCD.OutputJSONStruct> Outputs { get; set; }
    }
}
