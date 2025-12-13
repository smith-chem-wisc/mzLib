using Newtonsoft.Json;
using System.Text;

namespace Predictions.Koina.Client
{
    public class HTTP
    {
        public static readonly string _ModelsURL = "https://koina.wilhelmlab.org:443/v2/models/";
        public readonly HttpClient Client;

        public HTTP(int timeoutInMinutes = 1)
        {
            Client = new HttpClient { Timeout = TimeSpan.FromMinutes(timeoutInMinutes) };
        }

        public async Task<string> InferenceRequest(string modelName, Dictionary<string, object> request)
        {
            var json = JsonConvert.SerializeObject(request);
            using var content = new StringContent(json, Encoding.UTF8, "application/json");
            using var response = await Client.PostAsync($"{_ModelsURL}{modelName}/infer", content);

            if (!response.IsSuccessStatusCode)
            {
                var errorContent = await response.Content.ReadAsStringAsync();
                throw new HttpRequestException(
                    $"Request failed with status {(int)response.StatusCode} {response.ReasonPhrase}: {errorContent}");
            }

            // Stream instead of buffering
            using var stream = await response.Content.ReadAsStreamAsync();
            using var reader = new StreamReader(stream);
            return await reader.ReadToEndAsync(); // No buffer limit
        }

        public void TestConnection()
        {
            using var head = new HttpRequestMessage(HttpMethod.Head, "https://koina.wilhelmlab.org/");
            using var headResp = Client.Send(head);
            if (!headResp.IsSuccessStatusCode &&
                headResp.StatusCode != System.Net.HttpStatusCode.MethodNotAllowed)
            {
                throw new HttpRequestException($"Koina unreachable: {(int)headResp.StatusCode} {headResp.ReasonPhrase}");
            }
        }
    }
}
