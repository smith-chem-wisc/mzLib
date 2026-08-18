using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FlashLFQ
{
    /// <summary>
    /// An <see cref="IsotopicEnvelope"/> that also keeps the individual isotope peaks it was summed from, indexed by
    /// isotope number. The base class keeps only the summed intensity and the most abundant peak, which is all
    /// quantification needs; this keeps the rest so that the envelope can be written out peak by peak.
    ///
    /// Built only when FlashLFQ is asked for verbose peaks, because retaining every isotope peak of every envelope of
    /// every charge state of every peak costs memory proportional to the raw data rather than to the results.
    /// </summary>
    public class VerboseIsotopicEnvelope : IsotopicEnvelope
    {
        /// <summary>
        /// The isotope peaks this envelope was summed from, keyed by isotope number: 0 is the monoisotopic peak, 1
        /// the first C13 isotopologue, and so on. Isotopes that were looked for and not observed are absent rather
        /// than present holding null, and an isotope number can be negative when the peak finding mass sits above
        /// the monoisotopic mass.
        /// </summary>
        public IReadOnlyDictionary<int, IIndexedPeak> PeakDictionary { get; }

        /// <summary>
        /// The retention time of the most abundant peak, which is the scan this envelope was observed in.
        /// </summary>
        public double RetentionTime => IndexedPeak.RetentionTime;

        public VerboseIsotopicEnvelope(IIndexedPeak mostAbundantPeak, IReadOnlyList<IIndexedPeak> isotopologuePeaks,
            int chargeState, double intensity, double pearsonCorrelation, double monoisotopicMass,
            Tolerance isotopeTolerance)
            : base(mostAbundantPeak, chargeState, intensity, pearsonCorrelation)
        {
            PeakDictionary = BuildPeakDictionary(isotopologuePeaks, monoisotopicMass, chargeState, isotopeTolerance);
        }

        /// <summary>
        /// Assigns each observed isotope peak the isotope number whose expected m/z it sits closest to, keeping only
        /// peaks that fall inside the tolerance the peaks were found with in the first place.
        /// </summary>
        /// <remarks>
        /// The isotope number is computed from the m/z rather than searched for, so this is one pass over the peaks
        /// with no tolerance widening and no way to fail to terminate. Where two peaks claim the same isotope number
        /// - which the spacing makes possible only just inside the tolerance - the more intense one wins, since that
        /// is the one peak finding would have traced.
        /// </remarks>
        internal static Dictionary<int, IIndexedPeak> BuildPeakDictionary(IReadOnlyList<IIndexedPeak> peaks,
            double monoisotopicMass, int chargeState, Tolerance isotopeTolerance)
        {
            var dictionary = new Dictionary<int, IIndexedPeak>();
            if (peaks == null || peaks.Count == 0 || chargeState == 0)
            {
                return dictionary;
            }

            double monoisotopicMz = monoisotopicMass.ToMz(chargeState);
            double isotopeSpacingMz = Constants.C13MinusC12 / chargeState;

            foreach (IIndexedPeak peak in peaks)
            {
                int isotopeNumber = (int)Math.Round((peak.M - monoisotopicMz) / isotopeSpacingMz,
                    MidpointRounding.AwayFromZero);

                if (!isotopeTolerance.Within(peak.M, monoisotopicMz + isotopeNumber * isotopeSpacingMz))
                {
                    continue;
                }

                if (!dictionary.TryGetValue(isotopeNumber, out IIndexedPeak existing)
                    || peak.Intensity > existing.Intensity)
                {
                    dictionary[isotopeNumber] = peak;
                }
            }

            return dictionary;
        }

        /// <summary>
        /// The column label for one isotope peak of one charge state, e.g. i0z2 for the monoisotopic peak at 2+.
        /// </summary>
        public static string GetIsotopePeakName((int IsotopeNumber, int ChargeState) key)
        {
            return "i" + key.IsotopeNumber.ToString(System.Globalization.CultureInfo.InvariantCulture)
                + "z" + key.ChargeState.ToString(System.Globalization.CultureInfo.InvariantCulture);
        }

        public override string ToString()
        {
            return base.ToString() + "|" + PeakDictionary.Count + " isotopes";
        }
    }
}
