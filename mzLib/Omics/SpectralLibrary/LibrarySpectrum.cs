using System.Globalization;
using System.Text;
using Chemistry;
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

        public virtual string Name => Sequence + "/" + ChargeState;

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
                var product = matchedIon.NeutralTheoreticalProduct;

                string neutralLoss = null;
                if (product.NeutralLoss != 0)
                {
                    neutralLoss = "-" + product.NeutralLoss.ToString("F2", CultureInfo.InvariantCulture);
                }

                string ionAnnotation;
                if (product.IsInternalFragment)
                {
                    // Internal fragment format: bIb[31-34]^1/0ppm
                    ionAnnotation = $"{product.ProductType}I{product.SecondaryProductType}[{product.FragmentNumber}-{product.SecondaryFragmentNumber}]^{matchedIon.Charge}{neutralLoss}/0ppm";
                }
                else
                {
                    // Standard ion format: b3^1/0ppm or b3^1-97.98/0ppm
                    ionAnnotation = $"{product.ProductType}{product.FragmentNumber}^{matchedIon.Charge}{neutralLoss}/0ppm";
                }

                spectrum.Append("\n" + matchedIon.Mz.ToString(CultureInfo.InvariantCulture) + "\t" +
                    intensityFraction.ToString(CultureInfo.InvariantCulture) + "\t" + "\"" + ionAnnotation + "\"");
            }

            return spectrum.ToString();
        }

        // For decoy library spectrum generation, we use the predicted m/z value of the decoy sequence and we use the decoy's corresponding target's library spectrum's intensity values as decoy's intensities
        public static List<MatchedFragmentIon> GetDecoyLibrarySpectrumFromTargetByReverse(LibrarySpectrum targetSpectrum, List<Product> decoyPeptideTheorProducts)
        {
            var decoyFragmentIons = new List<MatchedFragmentIon>();
            foreach (var targetIon in targetSpectrum.MatchedFragmentIons)
            {
                foreach (var decoyPeptideTheorIon in decoyPeptideTheorProducts)
                {
                    if (targetIon.NeutralTheoreticalProduct.ProductType == decoyPeptideTheorIon.ProductType && targetIon.NeutralTheoreticalProduct.FragmentNumber == decoyPeptideTheorIon.FragmentNumber)
                    {
                        double decoyFragmentMz = decoyPeptideTheorIon.NeutralMass.ToMz(targetIon.Charge);
                        Product temProduct = decoyPeptideTheorIon;
                        decoyFragmentIons.Add(new MatchedFragmentIon(temProduct, decoyFragmentMz, targetIon.Intensity, targetIon.Charge));
                    }
                }
            }
            return decoyFragmentIons;
        }
    }
}
