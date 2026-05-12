using Chemistry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace MassSpectrometry
{
    /// <summary>
    /// Context class for all deconvolution
    /// </summary>
    public static class Deconvoluter
    {
        /// <summary>
        /// Static deconvolution of an MsDataScan that does not require Deconvoluter construction
        /// </summary>
        /// <param name="scan">scan to deconvolute</param>
        /// <param name="deconvolutionParameters">decon parameters to use, also determines type of deconvolution used</param>
        /// <param name="rangeToGetPeaksFrom">Range of peaks to deconvolute, if null, will deconvolute entire spectra</param>
        /// <returns></returns>
        public static IEnumerable<IsotopicEnvelope> Deconvolute(MsDataScan scan,
            DeconvolutionParameters deconvolutionParameters, MzRange rangeToGetPeaksFrom = null)
        {
            if (deconvolutionParameters.DeconvolutionType == DeconvolutionType.FromFile && !(rangeToGetPeaksFrom is MzRtRange))
            {
                // If the deconvolution type is FromFile, we expect the range to get peaks from to be an MzRtRange. If it is not provided, we will create one using the scan's retention time and the full m/z range of the spectrum
                rangeToGetPeaksFrom = rangeToGetPeaksFrom is null 
                    ? new MzRtRange(scan.MassSpectrum.Range, scan.RetentionTime) 
                    : new MzRtRange(rangeToGetPeaksFrom, scan.RetentionTime);
            }                

            return Deconvolute(scan.MassSpectrum, deconvolutionParameters, rangeToGetPeaksFrom);
        }

        /// <summary>
        /// Static deconvolution of an MzSpectrum that does not require Deconvoluter construction
        /// </summary>
        /// <param name="spectrum">spectrum to deconvolute</param>
        /// <param name="deconvolutionParameters">decon parameters to use, also determines type of deconvolution used</param>
        /// <param name="rangeToGetPeaksFrom">Range of peaks to deconvolute, if null, will deconvolute entire spectra</param>
        /// <returns></returns>
        public static IEnumerable<IsotopicEnvelope> Deconvolute(MzSpectrum spectrum,
            DeconvolutionParameters deconvolutionParameters, MzRange rangeToGetPeaksFrom = null)
        {
            rangeToGetPeaksFrom ??= spectrum.Range;

            if (deconvolutionParameters.DeconvolutionType == DeconvolutionType.FromFile && !(rangeToGetPeaksFrom is MzRtRange))
                throw new ArgumentException("If DeconvolutionType is FromFile, rangeToGetPeaksFrom must be an MzRtRange");

            // Short circuit deconvolution if it is called on a neutral mass spectrum
            if (spectrum is NeutralMassSpectrum newt)
                return DeconvoluteNeutralMassSpectrum(newt, rangeToGetPeaksFrom);

            // set deconvolution algorithm
            DeconvolutionAlgorithm deconAlgorithm = CreateAlgorithm(deconvolutionParameters);

            // Delegate deconvolution to the algorithm
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

            var decoyParams = parameters.ToDecoyParameters();
            if (decoyParams is null)
            {
                throw new InvalidOperationException(
                    $"DeconvoluteWithDecoys requires decoy support, but {parameters.GetType().Name} " +
                    $"does not implement ToDecoyParameters(). Override ToDecoyParameters() in your " +
                    $"parameter class to enable decoy deconvolution.");
            }
            var decoys = Deconvolute(spectrum, decoyParams, rangeToGetPeaksFrom).ToList();

            return (targets, decoys);
        }

        /// <summary>
        /// Factory method to create the correct deconvolution algorithm from the parameters
        /// </summary>
        /// <param name="parameters"></param>
        /// <returns></returns>
        /// <exception cref="MzLibException"></exception>
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

        /// <summary>
        /// Returns all peaks in the neutral mass spectrum as an isotopic envelope with a single peak
        /// </summary>
        /// <param name="neutralSpectrum"></param>
        /// <param name="range"></param>
        /// <returns></returns>
        private static IEnumerable<IsotopicEnvelope> DeconvoluteNeutralMassSpectrum(
            NeutralMassSpectrum neutralSpectrum, MzRange range)
        {
            for (int i = 0; i < neutralSpectrum.XArray.Length; i++)
            {
                double neutralMass = neutralSpectrum.XArray[i];
                double intensity = neutralSpectrum.YArray[i];
                int chargeState = neutralSpectrum.Charges[i];

                if (range.Contains(neutralMass.ToMz(chargeState)))
                {
                    yield return new IsotopicEnvelope(neutralMass, intensity, chargeState);
                }
            }
        }
    }
}
