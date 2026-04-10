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

        // This algorithm does not yet support decoy deconvolution.
        public override DeconvolutionParameters? ToDecoyParameters() => null;
    }
}