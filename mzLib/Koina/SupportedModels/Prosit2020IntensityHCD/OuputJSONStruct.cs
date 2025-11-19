using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Koina.SupportedModels.Prosit2020IntensityHCD
{
    public class OutputJSONStruct
    {
        [JsonProperty("name")]
        public string Name { get; set; }

        [JsonProperty("datatype")]
        public string Datatype { get; set; }

        [JsonProperty("shape")]
        public List<int> Shape { get; set; }

        [JsonProperty("data")]
        public List<object> Data { get; set; }
    }
}
