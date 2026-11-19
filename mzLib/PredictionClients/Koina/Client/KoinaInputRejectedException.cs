using System.Net;

namespace PredictionClients.Koina.Client
{
    /// <summary>
    /// Koina ran the model and the model refused something in the request: its preprocessing or its
    /// own code raised on one of the inputs.
    ///
    /// Koina answers <c>400 Bad Request</c> for the whole request when that happens, so one bad input
    /// costs every other input sent with it. The real case, from Prosit_2025_irt_40PTM given one
    /// peptide with E[UNIMOD:954]:
    ///
    /// <code>
    /// {"error":"in ensemble 'Prosit_2025_irt_40PTM', Failed to process the request(s) for model
    ///  instance 'Prosit_Preprocess_ac_loss_0_1', message: KeyError: \"E_954: no attachment group
    ///  for residue 'E'\"\n\nAt:\n  /models/repo/Prosit_Preprocess_ac_loss/1/zero_shot.py(77): ..."}
    /// </code>
    ///
    /// The other inputs are worth sending again without the bad one, and this type is what lets the
    /// prediction models tell that case apart: they split a request refused this way until they find
    /// the input the model will not take (see KoinaModelBase.RequestPredictionsAsync).
    ///
    /// Two other 400s look similar and are NOT this, and splitting either would only multiply the
    /// requests:
    ///   * a malformed request, such as <c>{"error":"unexpected inference input 'peptide_sequenes'
    ///     for model ..."}</c>. Triton refuses it before any model runs, so the marker below is
    ///     absent. Every request we build would carry the same fault, so it stays a plain
    ///     <see cref="HttpRequestException"/> and fails loudly;
    ///   * a service fault, such as a CUDA error, which is a <see cref="KoinaServiceException"/> and
    ///     is checked first.
    ///
    /// Derives from <see cref="HttpRequestException"/> so that existing <c>catch</c> clauses keep
    /// working unchanged.
    /// </summary>
    public class KoinaInputRejectedException : HttpRequestException
    {
        public KoinaInputRejectedException(string message, string serverError)
            : base(message, null, HttpStatusCode.BadRequest)
        {
            ServerError = serverError;
        }

        /// <summary>
        /// The server's own description of what it refused -- the <c>error</c> field of the response
        /// body when it parses as JSON, otherwise the raw body. Never null.
        /// </summary>
        public string ServerError { get; }

        /// <summary>
        /// Text that Triton's Python backend puts in front of any exception a model's own code raises
        /// while it handles a request. Matched case-insensitively against the server error.
        /// </summary>
        /// <remarks>
        /// Narrow on purpose, like <see cref="KoinaServiceException"/>'s markers. An unrecognised 400
        /// stays a plain failure, so the safe direction is still to under-match: missing a refusal only
        /// means the whole request fails, as it did before this type existed.
        /// </remarks>
        private const string ModelRaisedMarker = "failed to process the request(s) for model instance";

        /// <summary>
        /// True when a response with <paramref name="statusCode"/> and <paramref name="responseBody"/>
        /// means the model refused one of the inputs: a 400, not a service fault, and carrying the
        /// marker of a model's own code raising.
        /// </summary>
        /// <remarks>
        /// Only a 400. A 5xx is an outage whatever its body says, and must not be split.
        /// </remarks>
        public static bool IsInputRejection(int statusCode, string responseBody) =>
            statusCode == (int)HttpStatusCode.BadRequest
            && !KoinaServiceException.IsServiceFault(responseBody)
            && KoinaServiceException.ExtractServerError(responseBody)
                .Contains(ModelRaisedMarker, StringComparison.OrdinalIgnoreCase);
    }
}
