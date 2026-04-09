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

        // ToDecoyParameters() is not overridden — the base class default (returns null)
        // signals that this algorithm does not yet support decoy deconvolution.
    }
}