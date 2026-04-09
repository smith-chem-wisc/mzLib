using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MzLibUtil;

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
            // set any specific deconvolution parameters found only in the MsDataScan
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

            // Short circuit deconvolution if it is called on a neutral mass spectrum
            if (spectrum is NeutralMassSpectrum newt)
                return DeconvoluteNeutralMassSpectrum(newt, rangeToGetPeaksFrom);

            // set deconvolution algorithm 
            DeconvolutionAlgorithm deconAlgorithm = CreateAlgorithm(deconvolutionParameters);

            // Delegate deconvolution to the algorithm
            return deconAlgorithm.Deconvolute(spectrum, rangeToGetPeaksFrom);
        }

        /// <summary>
        /// Deconvolutes a spectrum to produce target envelopes, then runs a second
        /// decoy deconvolution pass with a shifted isotope distance (0.9444 Da instead
        /// of 1.003355 Da) to produce a null distribution for FDR estimation.
        ///
        /// The decoy pass uses the same algorithm, charge range, mass range, and
        /// tolerance as the target pass — only the isotope spacing changes. This ensures
        /// decoy envelopes face the same quality filters as targets, so the decoy score
        /// distribution reflects the algorithm's false-positive rate rather than a
        /// systematically easier or harder problem.
        /// </summary>
        /// <param name="scan">Scan to deconvolute.</param>
        /// <param name="parameters">
        /// Target deconvolution parameters. These are not modified; a clone is made
        /// internally for the decoy pass.
        /// </param>
        /// <param name="rangeToGetPeaksFrom">
        /// m/z range to deconvolute. If null, the full spectrum range is used.
        /// </param>
        /// <returns>
        /// A tuple of two fully-enumerated lists: target envelopes and decoy envelopes.
        /// The decoy envelopes have <see cref="DeconvolutionParameters.IsDecoyRun"/>
        /// set to true on the parameters used to generate them, but the envelopes
        /// themselves are identical in structure to target envelopes.
        /// </returns>
        public static (List<IsotopicEnvelope> Targets, List<IsotopicEnvelope> Decoys)
            DeconvoluteWithDecoys(
                MsDataScan scan,
                DeconvolutionParameters parameters,
                MzRange rangeToGetPeaksFrom = null)
            => DeconvoluteWithDecoys(scan.MassSpectrum, parameters, rangeToGetPeaksFrom);

        /// <summary>
        /// Deconvolutes a spectrum to produce target envelopes, then runs a second
        /// decoy deconvolution pass with a shifted isotope distance (0.9444 Da instead
        /// of 1.003355 Da) to produce a null distribution for FDR estimation.
        /// </summary>
        /// <param name="spectrum">Spectrum to deconvolute.</param>
        /// <param name="parameters">
        /// Target deconvolution parameters. Not modified; a clone is used for the decoy pass.
        /// </param>
        /// <param name="rangeToGetPeaksFrom">
        /// m/z range to deconvolute. If null, the full spectrum range is used.
        /// </param>
        public static (List<IsotopicEnvelope> Targets, List<IsotopicEnvelope> Decoys)
            DeconvoluteWithDecoys(
                MzSpectrum spectrum,
                DeconvolutionParameters parameters,
                MzRange rangeToGetPeaksFrom = null)
        {
            var targets = Deconvolute(spectrum, parameters, rangeToGetPeaksFrom).ToList();

            DeconvolutionParameters decoyParams =
                DeconvolutionDecoyGenerator.MakeDecoyParameters(parameters);
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
                //TODO: add FLASHDeconvolution here when it gets implemented
                //DeconvolutionType.FLASHDeconvolution => new FLASHDeconvolutionAlgorithm(parameters),
                _ => throw new MzLibException("DeconvolutionType not yet supported")
            };
        }

        /// <summary>
        /// Returns all peaks in the neutral mass spectrum as an isotopic envelope with a single peak
        /// </summary>
        /// <param name="neutralSpectrum"></param>
        /// <param name="range"></param>
        /// <returns></returns>
        private static IEnumerable<IsotopicEnvelope> DeconvoluteNeutralMassSpectrum(NeutralMassSpectrum neutralSpectrum, MzRange range)
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
