using Chemistry;
using MzLibUtil;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    public static class Deconvoluter
    {
        public static IEnumerable<IsotopicEnvelope> Deconvolute(MsDataScan scan,
            DeconvolutionParameters deconvolutionParameters, MzRange rangeToGetPeaksFrom = null)
        {
            return Deconvolute(scan.MassSpectrum, deconvolutionParameters, rangeToGetPeaksFrom);
        }

        public static IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum,
            DeconvolutionParameters deconvolutionParameters, MzRange rangeToGetPeaksFrom = null)
        {
            rangeToGetPeaksFrom ??= spectrum.Range;
            if (spectrum is NeutralMassSpectrum newt)
                return DeconvoluteNeutralMassSpectrum(newt, rangeToGetPeaksFrom);
            DeconvolutionAlgorithm deconAlgorithm = CreateAlgorithm(deconvolutionParameters);
            return deconAlgorithm.Deconvolute(spectrum, rangeToGetPeaksFrom);
        }

        public static (List<IsotopicEnvelope> Targets, List<IsotopicEnvelope> Decoys)
            DeconvoluteWithDecoys(MsDataScan scan, DeconvolutionParameters parameters,
                MzRange rangeToGetPeaksFrom = null)
            => DeconvoluteWithDecoys(scan.MassSpectrum, parameters, rangeToGetPeaksFrom);

        public static (List<IsotopicEnvelope> Targets, List<IsotopicEnvelope> Decoys)
            DeconvoluteWithDecoys(MzSpectrum spectrum, DeconvolutionParameters parameters,
                MzRange rangeToGetPeaksFrom = null)
        {
            var targets = Deconvolute(spectrum, parameters, rangeToGetPeaksFrom).ToList();

            DeconvolutionParameters decoyParams = parameters.ToDecoyParameters();
            var decoys = Deconvolute(spectrum, decoyParams, rangeToGetPeaksFrom).ToList();

            return (targets, decoys);
        }

        private static DeconvolutionAlgorithm CreateAlgorithm(DeconvolutionParameters parameters)
        {
            return parameters.DeconvolutionType switch
            {
                DeconvolutionType.ClassicDeconvolution => new ClassicDeconvolutionAlgorithm(parameters),
                DeconvolutionType.ExampleNewDeconvolutionTemplate => new ExampleNewDeconvolutionAlgorithmTemplate(parameters),
                DeconvolutionType.IsoDecDeconvolution => new IsoDecAlgorithm(parameters),
                _ => throw new MzLibException("DeconvolutionType not yet supported")
            };
        }

        private static IEnumerable<IsotopicEnvelope> DeconvoluteNeutralMassSpectrum(
            NeutralMassSpectrum neutralSpectrum, MzRange range)
        {
            for (int i = 0; i < neutralSpectrum.XArray.Length; i++)
            {
                double neutralMass = neutralSpectrum.XArray[i];
                double intensity = neutralSpectrum.YArray[i];
                int chargeState = neutralSpectrum.Charges[i];
                if (range.Contains(neutralMass.ToMz(chargeState)))
                    yield return new IsotopicEnvelope(neutralMass, intensity, chargeState);
            }
        }
    }
}