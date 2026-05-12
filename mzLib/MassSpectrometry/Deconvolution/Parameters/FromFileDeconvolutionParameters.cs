using System;
using System.Diagnostics.CodeAnalysis;

namespace MassSpectrometry
{
    [ExcludeFromCodeCoverage]
    public class FromFileDeconvolutionParameters : DeconvolutionParameters
    {
        private string _resultPath;

        public override DeconvolutionType DeconvolutionType { get; protected set; } = DeconvolutionType.FromFile;
        public FromFileDeconvolutionParameters(string resultPath, int minCharge, int maxCharge, Polarity polarity = Polarity.Positive)
            : base(minCharge, maxCharge, polarity)
        {
            _resultPath = resultPath;

            // Load Results

            // Create some sort of internal representation. 
        }

        // This algorithm does not yet support decoy deconvolution.
        public override DeconvolutionParameters? ToDecoyParameters() => null;
    }
}