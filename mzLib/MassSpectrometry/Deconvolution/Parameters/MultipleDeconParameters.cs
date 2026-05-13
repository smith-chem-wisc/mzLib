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

    public override DeconvolutionParameters ToDecoyParameters()
    {
        var decoyParams = Parameters.Select(p => p.ToDecoyParameters()).ToArray();

        if (decoyParams.Any(p => p is null))
            return null; // If any of the individual parameter sets don't support decoys, then we can't do decoys for the whole thing

        return new MultipleDeconParameters(decoyParams, MinAssumedChargeState, MaxAssumedChargeState, Polarity, new DecoyAveragine(AverageResidueModel, DecoyAveragine.DefaultDecoyIsotopeSpacing, ExpectedIsotopeSpacing), DecoyAveragine.DefaultDecoyIsotopeSpacing);
    }
}
