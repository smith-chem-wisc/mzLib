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

        /// <summary>
        /// Pairs each MS1 deconvolution feature with the MS2 scan(s) that selected it
        /// for fragmentation, returning one (MS2 scan, isotopic envelope) pair per match.
        /// </summary>
        /// <param name="ms1Features">
        /// Per-charge deconvolution results from any producer — mzLib readers for FlashDeconv
        /// or TopFD <c>.ms1.feature</c> output, or mzLib's own whole-file deconvolution.
        /// </param>
        /// <param name="msDataFile">The MS data file containing the MS2 scans to pair against.</param>
        /// <returns>
        /// One <c>(MsDataScan, IsotopicEnvelope)</c> pair per (feature, MS2-scan) match. A
        /// match requires:
        /// <list type="bullet">
        ///   <item><description>MS2 scan retention time inside
        ///   [<see cref="ISingleChargeMs1Feature.RetentionTimeStart"/>,
        ///    <see cref="ISingleChargeMs1Feature.RetentionTimeEnd"/>].</description></item>
        ///   <item><description>Feature m/z inside <see cref="MsDataScan.IsolationRange"/>.</description></item>
        /// </list>
        /// </returns>
        /// <remarks>
        /// Pure join. No cross-charge consensus, no off-by-one correction, no harmonic filtering —
        /// those concerns belong to the deconvolution producer upstream. Pairing is restricted to
        /// <see cref="MsDataScan.MsnOrder"/> == 2: MS3 isolation targets an MS2 fragment, not an
        /// MS1 precursor. Chimeras (multiple features matching one MS2) emit one pair per match;
        /// a single feature spanning multiple MS2 scans likewise emits one per scan.
        ///
        /// The returned envelope is built via <see cref="IsotopicEnvelope(double, double, int)"/>
        /// from the feature's <see cref="ISingleChargeMs1Feature.Mz"/>,
        /// <see cref="ISingleChargeMs1Feature.Charge"/>, and
        /// <see cref="ISingleChargeMs1Feature.Intensity"/>. The envelope's
        /// <c>Peaks</c> list contains a single synthetic entry — <c>NumberOfIsotopes</c> from the
        /// input is not surfaced here because the per-peak m/z and intensity are unknown for
        /// external reader output. Callers needing real peak lists should consume mzLib's own
        /// whole-file deconvolution output instead.
        /// </remarks>
        public static IEnumerable<(MsDataScan Ms2Scan, IsotopicEnvelope PrecursorEnvelope)>
            PairPrecursorsToMs2(IEnumerable<ISingleChargeMs1Feature> ms1Features, MsDataFile msDataFile)
        {
            var ms2Scans = new List<MsDataScan>();
            foreach (var scan in msDataFile.GetMsDataScans())
            {
                if (scan.MsnOrder != 2) continue;
                if (scan.IsolationRange == null) continue;
                ms2Scans.Add(scan);
            }

            foreach (var feat in ms1Features)
            {
                foreach (var scan in ms2Scans)
                {
                    if (scan.RetentionTime < feat.RetentionTimeStart) continue;
                    if (scan.RetentionTime > feat.RetentionTimeEnd) continue;
                    if (feat.Mz < scan.IsolationRange.Minimum) continue;
                    if (feat.Mz > scan.IsolationRange.Maximum) continue;

                    var envelope = new IsotopicEnvelope(feat.Mz.ToMass(feat.Charge), feat.Intensity, feat.Charge);
                    yield return (scan, envelope);
                }
            }
        }
    }
}
