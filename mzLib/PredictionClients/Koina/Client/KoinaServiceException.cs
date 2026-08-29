using Newtonsoft.Json.Linq;

namespace PredictionClients.Koina.Client
{
    /// <summary>
    /// Koina accepted the request and then failed to run the model.
    ///
    /// This is worth its own type because the transport cannot tell you about it. The server answers
    /// <c>400 Bad Request</c> for BOTH "your request is malformed" and "my GPU fell over mid-inference",
    /// putting the difference only in the response body:
    ///
    /// <code>
    /// {"error":"in ensemble 'Altimeter_2024_intensities', PyTorch execute failure: ...
    ///  RuntimeError: CUDA error: CUBLAS_STATUS_EXECUTION_FAILED when calling
    ///  `cublasSgemmStridedBatched(...)`"}
    /// </code>
    ///
    /// A caller that treats those the same either retries a request that will never work, or gives up
    /// on a request that would succeed on the next call. The distinction is drawn once, here, where
    /// the body is in hand -- rather than by every caller pattern-matching an exception message.
    ///
    /// Derives from <see cref="HttpRequestException"/> so that existing <c>catch</c> clauses keep
    /// working unchanged.
    /// </summary>
    public class KoinaServiceException : HttpRequestException
    {
        public KoinaServiceException(string message, string serverError)
            : base(message)
        {
            ServerError = serverError;
        }

        /// <summary>
        /// The server's own description of what went wrong -- the <c>error</c> field of the response
        /// body when it parses as JSON, otherwise the raw body. Never null; empty when the server sent
        /// nothing.
        /// </summary>
        public string ServerError { get; }

        /// <summary>
        /// Markers that identify a fault in Koina's own inference stack rather than a rejection of what
        /// we sent. Matched case-insensitively against <see cref="ServerError"/>.
        /// </summary>
        /// <remarks>
        /// Deliberately narrow, and it must stay that way. Everything here names the SERVER's execution
        /// failing: a Torch interpreter fault, a CUDA or cuBLAS call failing, a model that would not
        /// load, the GPU running out of memory. None of them can be provoked by a malformed request.
        ///
        /// The safe direction is to under-match. An unrecognised error is treated as our problem, so a
        /// genuine regression in how we build a request still fails the caller loudly. Widening this to
        /// something generic -- "error", "failed", "invalid" -- would silence exactly the failures a
        /// test suite exists to catch.
        /// </remarks>
        private static readonly string[] ServiceFaultMarkers =
        [
            "pytorch execute failure",
            "torchscript",
            "cuda error",
            "cublas_",
            "cudnn_",
            "out of memory",
            "failed to load model",
            "model is not ready",
            "inference server is not live"
        ];

        /// <summary>
        /// True when <paramref name="responseBody"/> describes Koina failing to run the model, as
        /// opposed to rejecting the request.
        /// </summary>
        /// <remarks>
        /// Takes the raw body rather than a parsed object so that a body which is not JSON at all -- an
        /// HTML error page from a proxy, say -- is still classified rather than throwing.
        /// </remarks>
        public static bool IsServiceFault(string responseBody) =>
            ServiceFaultMarkers.Any(marker =>
                ExtractServerError(responseBody).Contains(marker, StringComparison.OrdinalIgnoreCase));

        /// <summary>
        /// The <c>error</c> field of a JSON body, or the body itself when it is not JSON or carries no
        /// such field. Returns an empty string for a null or blank body.
        /// </summary>
        public static string ExtractServerError(string responseBody)
        {
            if (string.IsNullOrWhiteSpace(responseBody)) return string.Empty;

            try
            {
                if (JToken.Parse(responseBody) is JObject body
                    && body["error"]?.Type == JTokenType.String)
                {
                    return body["error"]!.Value<string>() ?? string.Empty;
                }
            }
            catch (Newtonsoft.Json.JsonReaderException)
            {
                // Not JSON. The raw body is the best description available, so fall through to it
                // rather than discarding the one piece of evidence the server sent.
            }

            return responseBody;
        }
    }
}
