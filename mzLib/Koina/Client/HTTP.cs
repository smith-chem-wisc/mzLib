using Koina.SupportedModels;
using Newtonsoft.Json;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Net.Http;
using System.Net.Http.Json;
using System.Text;
using System.Threading.Tasks;
using static System.Net.WebRequestMethods;

namespace Koina.Client
{
    public class HTTP
    {
        public static readonly string _ModelsURL = "https://koina.wilhelmlab.org:443/v2/models/";
        public readonly HttpClient Client;

        public HTTP()
        {
            Client = new HttpClient { Timeout = TimeSpan.FromSeconds(15) };
        }

        public async Task<string> InferenceRequest(string modelName, Dictionary<string, object> request)
        {
            try
            {
                TestConnection();

                var json = JsonConvert.SerializeObject(request);
                using var response = await Client.PostAsync($"{_ModelsURL}{modelName}/infer", new StringContent(json, Encoding.UTF8, "application/json"));

                string payload = await response.Content.ReadAsStringAsync();

                if (response.StatusCode == System.Net.HttpStatusCode.NotFound)
                {
                    throw new HttpRequestException($"404 Not Found for {modelName}. Body:\n{payload}");
                }
                response.EnsureSuccessStatusCode();
                return payload;
            }
            catch (Exception ex)
            {
                throw new HttpRequestException("Inference request failed: " + ex.Message);
            }
        }

        public void TestConnection()
        {
            try
            {
                using var head = new HttpRequestMessage(HttpMethod.Head, "https://koina.wilhelmlab.org/");
                using var headResp = Client.Send(head);
                if (!headResp.IsSuccessStatusCode &&
                    headResp.StatusCode != System.Net.HttpStatusCode.MethodNotAllowed)
                {
                    Assert.Inconclusive($"Koina unreachable: {(int)headResp.StatusCode} {headResp.ReasonPhrase}");
                }
            }
            catch (Exception ex)
            {
                Assert.Inconclusive("Network not available: " + ex.Message);
            }
        }
    }
}
