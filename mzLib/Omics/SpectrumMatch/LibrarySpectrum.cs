using System.Text;
using Easy.Common.Extensions;
using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using Omics.Fragmentation;

namespace Omics.SpectrumMatch
{
    public class LibrarySpectrum : MzSpectrum
    {
        public string Sequence { get; set; }
        public double? RetentionTime { get; set; }
        public double PrecursorMz { get; set; }
        public int ChargeState { get; set; }
        public List<MatchedFragmentIon> MatchedFragmentIons { get; set; }
        public bool IsDecoy { get; set; }

        public virtual string Name
        {
            get { return Sequence + "/" + ChargeState; }
        }

        public LibrarySpectrum(string sequence, double precursorMz, int chargeState, List<MatchedFragmentIon> peaks, double? rt, bool isDecoy = false) 
            : base(peaks.Select(p => p.Mz).ToArray(), peaks.Select(p => p.Intensity).ToArray(), false)
        {
            Sequence = sequence;
            PrecursorMz = precursorMz;
            MatchedFragmentIons = peaks;
            ChargeState = chargeState;
            IsDecoy = isDecoy;
            RetentionTime = rt;
            Array.Sort(XArray, YArray);
        }

        /// <summary>
        /// This function enables the spectrum angle to be computed between an individual experimental spectrum and the loaded library spectrum within MetaDraw
        /// </summary>
        /// <param name="librarySpectrum"></param>
        /// <returns></returns>
        public string CalculateSpectralAngleOnTheFly(List<MatchedFragmentIon> spectrumMatchFragments)
        {
            if (!spectrumMatchFragments.Any())
            {
                return "N/A";
            }

            if (spectrumMatchFragments.IsNotNullOrEmpty()) { }
            SpectralSimilarity spectraComparison = new SpectralSimilarity(
                spectrumMatchFragments.Select(f => f.Mz).ToArray(),
                spectrumMatchFragments.Select(f => f.Intensity).ToArray(),
                MatchedFragmentIons.Select(f => f.Mz).ToArray(),
                MatchedFragmentIons.Select(f => f.Intensity).ToArray(),
                SpectralSimilarity.SpectrumNormalizationScheme.MostAbundantPeak,
                toleranceInPpm: 20,
                allPeaks: true);
            double? spectralContrastAngle = spectraComparison.SpectralContrastAngle();

            return spectralContrastAngle == null
                ? "N/A"
                : ((double)spectralContrastAngle).ToString("F4");
        }

        public override string ToString()
        {
            StringBuilder spectrum = new StringBuilder();
            spectrum.Append("Name: " + Name);
            spectrum.Append("\nMW: " + PrecursorMz);
            spectrum.Append("\nComment: ");
            spectrum.Append("Parent=" + PrecursorMz);
            spectrum.Append(" RT=" + RetentionTime);
            spectrum.Append("\nNum peaks: " + MatchedFragmentIons.Count);

            double maxIntensity = MatchedFragmentIons.Select(b => b.Intensity).Max();

            foreach (MatchedFragmentIon matchedIon in MatchedFragmentIons)
            {
                double intensityFraction = matchedIon.Intensity / maxIntensity;

                string neutralLoss = null;
                if (matchedIon.NeutralTheoreticalProduct.NeutralLoss != 0)
                {
                    neutralLoss = "-" + matchedIon.NeutralTheoreticalProduct.NeutralLoss;
                }

                spectrum.Append("\n" + matchedIon.Mz + "\t" + intensityFraction + "\t" + "\"" +
                    matchedIon.NeutralTheoreticalProduct.ProductType.ToString() +
                    matchedIon.NeutralTheoreticalProduct.FragmentNumber.ToString() + "^" +
                    matchedIon.Charge + neutralLoss + "/" + 0 + "ppm" + "\"");
            }

            return spectrum.ToString();
        }
    }
}
