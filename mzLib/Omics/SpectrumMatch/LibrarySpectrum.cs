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

        public string Name
        {
            get { return Sequence + "/" + ChargeState; }
        }

        public LibrarySpectrum(string fullSequence, double precursorMz, int chargeState, List<MatchedFragmentIon> peaks, double? rt, bool isDecoy = false) 
            : base(peaks.Select(p => p.Mz).ToArray(), peaks.Select(p => p.Intensity).ToArray(), false)
        {
            Sequence = fullSequence;
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

        public string ToFraggerLibraryString(string proteinAccession, string geneName)
        {
            StringBuilder spectrum = new StringBuilder();
            //for now we ignore any ions with neutral loss
            foreach (MatchedFragmentIon matchedIon in MatchedFragmentIons.Where(mfi=>mfi.NeutralTheoreticalProduct.NeutralLoss == 0))
            {
                StringBuilder spectrumRow = new StringBuilder();
                spectrumRow.Append(PrecursorMz + "\t"); // PrecursorMz
                spectrumRow.Append(matchedIon.Mz + "\t"); // ProductMz
                spectrumRow.Append(matchedIon.Annotation.Replace('+', '^') + "\t"); // Annotation
                spectrumRow.Append(proteinAccession + "\t"); // ProteinId
                spectrumRow.Append(geneName + "\t"); // GeneName
                spectrumRow.Append(BaseSequenceFromFullSequence() + "\t"); // PeptideSequence
                spectrumRow.Append(Sequence + "\t"); // ModifiedPeptideSequence
                spectrumRow.Append(ChargeState + "\t"); // PrecursorCharge
                spectrumRow.Append(matchedIon.Intensity + "\t"); // LibraryIntensity
                spectrumRow.Append(RetentionTime + "\t"); // NormalizedRetentionTime
                spectrumRow.Append("" + "\t"); // PrecursorIonMobility
                spectrumRow.Append(matchedIon.NeutralTheoreticalProduct.ProductType + "\t"); // FragmentType
                spectrumRow.Append(matchedIon.Charge + "\t"); // FragmentCharge
                spectrumRow.Append(matchedIon.NeutralTheoreticalProduct.FragmentNumber + "\t"); // FragmentSeriesNumber
                spectrumRow.Append("" + "\t"); // FragmentLossType
                spectrumRow.Append(RetentionTime + "\t"); // AverageExperimentalRetentionTime
                spectrumRow.Append("sp|" + proteinAccession + "|" + geneName + "\t"); // AllMappedProteins
                spectrumRow.Append(geneName + "\t"); // AllMappedGenes
                spectrumRow.Append(1); // Proteotypic
                 

                spectrum.Append(spectrumRow + "\n");
            }

            return spectrum.ToString();
        }
        public string FraggerLibraryHeader()
        {
            return "PrecursorMz\tProductMz\tAnnotation\tProteinId\tGeneName\tPeptideSequence\tModifiedPeptideSequence\tPrecursorCharge\tLibraryIntensity\tNormalizedRetentionTime\tPrecursorIonMobility\tFragmentType\tFragmentCharge\tFragmentSeriesNumber\tFragmentLossType\tAverageExperimentalRetentionTime\tAllMappedProteins\tAllMappedGenes\tProteotypic";
        }
        public string BaseSequenceFromFullSequence()
        {
            StringBuilder result = new StringBuilder();
            int bracketLevel = 0;

            foreach (char c in Sequence)
            {
                if (c == '[')
                {
                    bracketLevel++;
                }
                else if (c == ']')
                {
                    if (bracketLevel > 0)
                    {
                        bracketLevel--;
                    }
                }
                else if (bracketLevel == 0)
                {
                    result.Append(c);
                }
            }

            return result.ToString();
        }
        public override string ToString()
        {
            StringBuilder stringBuilder = new StringBuilder();
            stringBuilder.Append("Name: " + Name);
            stringBuilder.Append("\nMW: " + PrecursorMz);
            stringBuilder.Append("\nComment: ");
            stringBuilder.Append("Parent=" + PrecursorMz);
            stringBuilder.Append(" RT=" + RetentionTime);
            stringBuilder.Append("\nNum peaks: " + MatchedFragmentIons.Count);
            double num = MatchedFragmentIons.Select((MatchedFragmentIon b) => b.Intensity).Max();
            foreach (MatchedFragmentIon matchedFragmentIon in MatchedFragmentIons)
            {
                double num2 = matchedFragmentIon.Intensity / num;
                string text = null;
                if (matchedFragmentIon.NeutralTheoreticalProduct.NeutralLoss != 0.0)
                {
                    text = "-" + matchedFragmentIon.NeutralTheoreticalProduct.NeutralLoss;
                }

                stringBuilder.Append("\n" + matchedFragmentIon.Mz + "\t" + num2 + "\t\"" + matchedFragmentIon.NeutralTheoreticalProduct.ProductType.ToString() + matchedFragmentIon.NeutralTheoreticalProduct.FragmentNumber + "^" + matchedFragmentIon.Charge + text + "/" + 0 + "ppm\"");
            }

            return stringBuilder.ToString();
        }
    }
}
