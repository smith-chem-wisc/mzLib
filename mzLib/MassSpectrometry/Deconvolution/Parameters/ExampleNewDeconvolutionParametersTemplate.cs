using System;
using System.Diagnostics.CodeAnalysis;

namespace MassSpectrometry
{
    [ExcludeFromCodeCoverage]
    public class ExampleNewDeconvolutionParametersTemplate : DeconvolutionParameters
    {
        public override DeconvolutionType DeconvolutionType { get; protected set; }
            = DeconvolutionType.ExampleNewDeconvolutionTemplate;

        public ExampleNewDeconvolutionParametersTemplate(int minCharge, int maxCharge,
            Polarity polarity = Polarity.Positive)
            : base(minCharge, maxCharge, polarity)
        {
        }

        /// <summary>
        /// Replace this with a cached clone of these parameters configured for decoy deconvolution.
        /// See <see cref="ClassicDeconvolutionParameters.ToDecoyParameters"/> for a reference implementation.
        /// </summary>
        public override DeconvolutionParameters ToDecoyParameters() =>
            throw new NotImplementedException(
                "Implement ToDecoyParameters() for your new deconvolution method. " +
                "See ClassicDeconvolutionParameters for the recommended pattern.");
    }
}