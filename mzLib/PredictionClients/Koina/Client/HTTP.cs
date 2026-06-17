using Newtonsoft.Json;
using System.Text;

namespace PredictionClients.Koina.Client
{
    /// <summary>
    /// Wraps a single shared <see cref="HttpClient"/> to avoid socket exhaustion across many
    /// prediction sessions. The client has no timeout; each request is bounded by a caller-supplied token.
    /// </summary>
    public static class HTTP
    {
        public static readonly string ModelsURL = "https://koina.wilhelmlab.org:443/v2/models/";

        private static readonly HttpClient Client = new() { Timeout = Timeout.InfiniteTimeSpan };

        public static async Task<string> InferenceRequest(string modelName, Dictionary<string, object> request, CancellationToken cancellationToken = default)
        {
            var json = JsonConvert.SerializeObject(request);
            using var content = new StringContent(json, Encoding.UTF8, "application/json");
            using var response = await Client.PostAsync($"{ModelsURL}{modelName}/infer", content, cancellationToken);

            if (!response.IsSuccessStatusCode)
            {
                var errorContent = await response.Content.ReadAsStringAsync(cancellationToken);
                throw new HttpRequestException(
                    $"Request failed with status {(int)response.StatusCode} {response.ReasonPhrase}: {errorContent}");
            }

            // Stream instead of buffering
            using var stream = await response.Content.ReadAsStreamAsync(cancellationToken);
            using var reader = new StreamReader(stream);
            return await reader.ReadToEndAsync(cancellationToken);
        }

        public static async Task TestConnectionAsync(CancellationToken cancellationToken = default)
        {
            using var head = new HttpRequestMessage(HttpMethod.Head, "https://koina.wilhelmlab.org/");
            using var headResp = await Client.SendAsync(head, cancellationToken);
            if (!headResp.IsSuccessStatusCode &&
                headResp.StatusCode != System.Net.HttpStatusCode.MethodNotAllowed)
            {
                throw new HttpRequestException($"Koina unreachable: {(int)headResp.StatusCode} {headResp.ReasonPhrase}");
            }
        }
    }
}
