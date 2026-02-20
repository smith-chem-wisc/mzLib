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
            spectrum.Append("\nMW: " + PrecursorMz.ToString(CultureInfo.InvariantCulture));
            spectrum.Append("\nComment: ");
            spectrum.Append("Parent=" + PrecursorMz.ToString(CultureInfo.InvariantCulture));
            spectrum.Append(" RT=" + RetentionTime?.ToString(CultureInfo.InvariantCulture));
            spectrum.Append("\nNum peaks: " + MatchedFragmentIons.Count);

            double maxIntensity = MatchedFragmentIons.Select(b => b.Intensity).Max();

            foreach (MatchedFragmentIon matchedIon in MatchedFragmentIons)
            {
                double intensityFraction = matchedIon.Intensity / maxIntensity;
                var product = matchedIon.NeutralTheoreticalProduct;

                string neutralLoss = null;
                if (product.NeutralLoss != 0)
                {
                    neutralLoss = "-" + product.NeutralLoss.ToString(CultureInfo.InvariantCulture);
                }

                string ionAnnotation;
                if (product.IsInternalFragment)
                {
                    // Internal fragment format: bIb[31-34]^1/0ppm
                    ionAnnotation = $"{product.ProductType}I{product.SecondaryProductType}[{product.FragmentNumber}-{product.SecondaryFragmentNumber}]^{matchedIon.Charge}{neutralLoss}/0ppm";
                }
                else
                {
                    // Standard ion format: b3^1/0ppm or b3^1-97.976895573/0ppm
                    ionAnnotation = $"{product.ProductType}{product.FragmentNumber}^{matchedIon.Charge}{neutralLoss}/0ppm";
                }

                spectrum.Append("\n" + matchedIon.Mz.ToString(CultureInfo.InvariantCulture) + "\t" +
                    intensityFraction.ToString(CultureInfo.InvariantCulture) + "\t\"" + ionAnnotation + "\"");
            }

            return spectrum.ToString();
        }

        /// <summary>
        /// Generates a decoy library spectrum from a target spectrum by matching fragment ions.
        /// For each target ion, finds the corresponding decoy theoretical product with matching properties
        /// and creates a new matched fragment ion using the decoy's m/z but preserving the target's intensity.
        /// </summary>
        /// <remarks>
        /// Matching logic:
        /// <list type="bullet">
        ///   <item><description>Standard ions: Match by ProductType and FragmentNumber</description></item>
        ///   <item><description>Internal fragment ions: Additionally require matching SecondaryProductType and SecondaryFragmentNumber</description></item>
        /// </list>
        /// This prevents false matches between standard and internal ions (e.g., b2 vs bIb[2-5]).
        /// </remarks>
        /// <param name="targetSpectrum">The target library spectrum containing matched fragment ions to reverse</param>
        /// <param name="decoyPeptideTheorProducts">Theoretical products from the decoy peptide sequence</param>
        /// <returns>A list of matched fragment ions for the decoy spectrum, preserving target intensities with decoy m/z values</returns>

        public static List<MatchedFragmentIon> GetDecoyLibrarySpectrumFromTargetByReverse(LibrarySpectrum targetSpectrum, List<Product> decoyPeptideTheorProducts)
        {
            var decoyFragmentIons = new List<MatchedFragmentIon>();
            foreach (var targetIon in targetSpectrum.MatchedFragmentIons)
            {
                foreach (var decoyPeptideTheorIon in decoyPeptideTheorProducts)
                {
                    var targetProduct = targetIon.NeutralTheoreticalProduct;

                    // Check primary properties
                    bool primaryMatch = targetProduct.ProductType == decoyPeptideTheorIon.ProductType
                        && targetProduct.FragmentNumber == decoyPeptideTheorIon.FragmentNumber;

                    if (!primaryMatch) continue;

                    // For internal fragments, also check secondary properties
                    if (targetProduct.IsInternalFragment || decoyPeptideTheorIon.IsInternalFragment)
                    {
                        // Both must be internal fragments with matching secondary properties
                        if (targetProduct.IsInternalFragment != decoyPeptideTheorIon.IsInternalFragment)
                            continue;
                        if (targetProduct.SecondaryProductType != decoyPeptideTheorIon.SecondaryProductType)
                            continue;
                        if (targetProduct.SecondaryFragmentNumber != decoyPeptideTheorIon.SecondaryFragmentNumber)
                            continue;
                    }

                    double decoyFragmentMz = decoyPeptideTheorIon.NeutralMass.ToMz(targetIon.Charge);
                    decoyFragmentIons.Add(new MatchedFragmentIon(decoyPeptideTheorIon, decoyFragmentMz, targetIon.Intensity, targetIon.Charge));
                }
            }
            return decoyFragmentIons;
        }
    }
}
