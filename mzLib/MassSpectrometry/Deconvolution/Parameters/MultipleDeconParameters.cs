using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;

namespace MassSpectrometry;
public class MultipleDeconParameters : DeconvolutionParameters
{
    public override DeconvolutionType DeconvolutionType { get; protected set; } = DeconvolutionType.Multiple;

    public DeconvolutionParameters[] Parameters { get; protected set; }

    public MultipleDeconParameters(IEnumerable<DeconvolutionParameters> deconParameters, int minCharge, int maxCharge, Polarity polarity = Polarity.Positive, AverageResidue averageResidueModel = null, double expectedIsotopeSpacing = Constants.C13MinusC12) : base(minCharge, maxCharge, polarity, averageResidueModel, expectedIsotopeSpacing)
    {
        Parameters = deconParameters.ToArray();
    }

    private volatile DeconvolutionParameters? _decoyParams = null;
    private readonly object _decoyParamsLock = new();

    public override DeconvolutionParameters ToDecoyParameters()
    {
        if (_decoyParams is not null)
            return _decoyParams;

        lock (_decoyParamsLock)
        {
            if (_decoyParams is not null)
                return _decoyParams;

            var decoySubParams = Parameters.Select(p => p.ToDecoyParameters()).ToArray();

            if (decoySubParams.Any(p => p is null))
                return null; // If any of the individual parameter sets don't support decoys, then we can't do decoys for the whole thing

            return _decoyParams = new MultipleDeconParameters(decoySubParams, MinAssumedChargeState, MaxAssumedChargeState, Polarity, new DecoyAveragine(AverageResidueModel, DecoyAveragine.DefaultDecoyIsotopeSpacing, ExpectedIsotopeSpacing), DecoyAveragine.DefaultDecoyIsotopeSpacing);
        }
    }
}
