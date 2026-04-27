using System.Collections.Generic;
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
        /// Deconvolutes a scan and additionally computes the generic deconvolution score for
        /// every yielded envelope. Each envelope's <see cref="IsotopicEnvelope.GenericScore"/>
        /// is set before it is yielded; the algorithm-specific <see cref="IsotopicEnvelope.Score"/>
        /// is unchanged.
        /// </summary>
        /// <remarks>
        /// If you already have an <see cref="IsotopicEnvelope"/> from a prior deconvolution call
        /// and want to score it without re-running the algorithm, use the
        /// <see cref="IsotopicEnvelopeExtensions.GetOrComputeGenericScore(IsotopicEnvelope, DeconvolutionParameters)"/>
        /// extension method.
        /// </remarks>
        public static IEnumerable<IsotopicEnvelope> DeconvoluteWithGenericScoring(MsDataScan scan,
            DeconvolutionParameters deconvolutionParameters, MzRange rangeToGetPeaksFrom = null)
        {
            return DeconvoluteWithGenericScoring(scan.MassSpectrum, deconvolutionParameters, rangeToGetPeaksFrom);
        }

        /// <summary>
        /// Deconvolutes a spectrum and additionally computes the generic deconvolution score for
        /// every yielded envelope. Each envelope's <see cref="IsotopicEnvelope.GenericScore"/>
        /// is set before it is yielded; the algorithm-specific <see cref="IsotopicEnvelope.Score"/>
        /// is unchanged.
        /// </summary>
        /// <remarks>
        /// Generic scoring adds a per-envelope feature-extraction pass (cosine vs Averagine,
        /// ppm error, completeness, ratio consistency). Use this overload when you need a
        /// score directly comparable across deconvolution algorithms; otherwise call
        /// <see cref="Deconvolute(MzSpectrum, DeconvolutionParameters, MzRange)"/> to avoid the cost.
        ///
        /// If you already have an <see cref="IsotopicEnvelope"/> from a prior deconvolution call
        /// and want to score it without re-running the algorithm, use the
        /// <see cref="IsotopicEnvelopeExtensions.GetOrComputeGenericScore(IsotopicEnvelope, DeconvolutionParameters)"/>
        /// extension method.
        /// </remarks>
        public static IEnumerable<IsotopicEnvelope> DeconvoluteWithGenericScoring(MzSpectrum spectrum,
            DeconvolutionParameters deconvolutionParameters, MzRange rangeToGetPeaksFrom = null)
        {
            AverageResidue model = deconvolutionParameters.AverageResidueModel;
            foreach (var envelope in Deconvolute(spectrum, deconvolutionParameters, rangeToGetPeaksFrom))
            {
                if (envelope.GenericScore == null)
                {
                    double score = DeconvolutionScorer.ScoreEnvelope(envelope, model);
                    envelope.SetGenericScore(score);
                }
                yield return envelope;
            }
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
