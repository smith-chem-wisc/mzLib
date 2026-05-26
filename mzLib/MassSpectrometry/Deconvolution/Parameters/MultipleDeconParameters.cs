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
        ArgumentNullException.ThrowIfNull(deconParameters);
        Parameters = deconParameters.ToArray();
        if (Parameters.Length == 0)
            throw new ArgumentException(
                "At least one sub-deconvolution parameter set is required.", nameof(deconParameters));
    }

    #region IEquatable<MultipleDeconParameters>

    protected override bool EqualProperties(DeconvolutionParameters other)
    {
        var o = (MultipleDeconParameters)other;
        return Parameters.SequenceEqual(o.Parameters);
    }

    protected override void AddHashCodes(HashCode hash)
    {
        foreach (var p in Parameters)
            hash.Add(p);
    }

    public override MultipleDeconParameters Clone()
    {
        return new MultipleDeconParameters(
            Parameters.Select(p => (DeconvolutionParameters)p.Clone()),
            MinAssumedChargeState, MaxAssumedChargeState,
            Polarity, AverageResidueModel, ExpectedIsotopeSpacing)
        {
            UseGenericScore = UseGenericScore
        };
    }

    #endregion

    private volatile DeconvolutionParameters? _decoyParams = null;
    private volatile bool _decoyComputed;
    private readonly object _decoyParamsLock = new();

    public override DeconvolutionParameters ToDecoyParameters()
    {
        if (_decoyComputed)
            return _decoyParams;

        lock (_decoyParamsLock)
        {
            if (_decoyComputed)
                return _decoyParams;

            var decoySubParams = Parameters.Select(p => p.ToDecoyParameters()).ToArray();

            if (decoySubParams.Any(p => p is null))
            {
                _decoyComputed = true;
                return null;
            }

            _decoyParams = new MultipleDeconParameters(decoySubParams, MinAssumedChargeState, MaxAssumedChargeState, Polarity, new DecoyAveragine(AverageResidueModel, DecoyAveragine.DefaultDecoyIsotopeSpacing, ExpectedIsotopeSpacing), DecoyAveragine.DefaultDecoyIsotopeSpacing);
            _decoyComputed = true;
            return _decoyParams;
        }
    }
}
