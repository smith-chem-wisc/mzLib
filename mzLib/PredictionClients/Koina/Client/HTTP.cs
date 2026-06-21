using Newtonsoft.Json;
using System.Text;

namespace PredictionClients.Koina.Client
{
    /// <summary>
    /// Wraps a single shared <see cref="HttpClient"/> to avoid socket exhaustion across many
    /// prediction sessions. Each request must supply a <see cref="CancellationToken"/> that bounds it
    /// to the caller's session deadline (the token is required, not optional). A generous fixed
    /// <see cref="HttpClient.Timeout"/> is kept only as a coarse backstop so that a caller passing
    /// <see cref="CancellationToken.None"/> still cannot hang a request indefinitely.
    /// </summary>
    public static class HTTP
    {
        public static readonly string ModelsURL = "https://koina.wilhelmlab.org:443/v2/models/";

        // Coarse per-request (per-POST) backstop only. The precise, batch-size-scaled bound is the
        // caller-supplied token (see the CancellationTokenSource in each model's AsyncThrottledPredictor).
        // A single batch chunk completes in seconds, so this never fires for a healthy request; it only
        // guarantees a stalled connection cannot block forever when a token never fires.
        private static readonly HttpClient Client = new() { Timeout = TimeSpan.FromMinutes(10) };

        public static async Task<string> InferenceRequest(string modelName, Dictionary<string, object> request, CancellationToken cancellationToken)
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

        public static async Task TestConnectionAsync(CancellationToken cancellationToken)
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
