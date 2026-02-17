using Chemistry;
using MassSpectrometry;
using Omics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

namespace Readers.InternalIons
{
    /// <summary>
    /// Extracts internal fragment ion features from PSMs for downstream analysis.
    /// </summary>
    public static class InternalFragmentFeatureExtractor
    {
        private static readonly Regex InternalFragmentRegex = new(@"[yb]I[yb]\[(\d+)-(\d+)\]",
            RegexOptions.IgnoreCase | RegexOptions.Compiled);

        /// <summary>
        /// Extracts InternalFragmentIon objects from a list of PSMs.
        /// </summary>
        /// <param name="psms">List of PSMs to process.</param>
        /// <param name="msDataFile">The raw data file containing the spectra.</param>
        /// <param name="defaultCollisionEnergy">Default collision energy to use if scan metadata is unavailable.</param>
        /// <returns>List of extracted InternalFragmentIon objects.</returns>
        public static List<InternalFragmentIon> ExtractFromPsms(
            List<PsmFromTsv> psms,
            MsDataFile msDataFile,
            double defaultCollisionEnergy = double.NaN)
        {
            if (psms == null || psms.Count == 0)
                return new List<InternalFragmentIon>();

            var results = new List<InternalFragmentIon>();
            var scanLookup = BuildScanLookup(msDataFile);

            foreach (var psm in psms)
            {
                if (psm.MatchedIons == null || psm.MatchedIons.Count == 0)
                    continue;

                var internalIons = psm.MatchedIons
                    .Where(ion => ion.IsInternalFragment)
                    .ToList();

                if (internalIons.Count == 0)
                    continue;

                double basePeakIntensity = psm.MatchedIons.Max(ion => ion.Intensity);
                if (basePeakIntensity <= 0)
                    basePeakIntensity = 1.0;

                MsDataScan? scan = null;
                if (scanLookup.TryGetValue(psm.Ms2ScanNumber, out var foundScan))
                    scan = foundScan;

                // Try to get collision energy from scan metadata, fall back to default
                double collisionEnergy = double.NaN;
                if (scan?.HcdEnergy != null && TryParseCollisionEnergy(scan.HcdEnergy, out double ce))
                {
                    collisionEnergy = ce;
                }
                else
                {
                    collisionEnergy = defaultCollisionEnergy;
                }

                // Extract all internal fragments for this PSM
                var psmInternalFragments = new List<InternalFragmentIon>();

                foreach (var ion in internalIons)
                {
                    var internalFragment = ExtractSingleInternalFragment(
                        psm, ion, basePeakIntensity, collisionEnergy, scan);

                    if (internalFragment != null)
                        psmInternalFragments.Add(internalFragment);
                }

                // Post-process: Mark isobaric ambiguous ions
                MarkIsobaricAmbiguousIons(psmInternalFragments);

                results.AddRange(psmInternalFragments);
            }

            return results;
        }

        /// <summary>
        /// Groups ions by theoretical mass (rounded to 4 decimal places) and marks
        /// any group with more than one ion as isobaric ambiguous.
        /// </summary>
        private static void MarkIsobaricAmbiguousIons(List<InternalFragmentIon> ions)
        {
            if (ions.Count <= 1)
                return;

            var massGroups = ions.GroupBy(ion => Math.Round(ion.TheoreticalMass, 4));

            foreach (var group in massGroups)
            {
                if (group.Count() > 1)
                {
                    foreach (var ion in group)
                    {
                        ion.IsIsobaricAmbiguous = true;
                    }
                }
            }
        }

        private static bool TryParseCollisionEnergy(string hcdEnergy, out double collisionEnergy)
        {
            collisionEnergy = double.NaN;
            if (string.IsNullOrWhiteSpace(hcdEnergy))
                return false;

            string cleaned = hcdEnergy.Replace("@", "").Replace("HCD", "").Replace("hcd", "").Trim();
            return double.TryParse(cleaned, out collisionEnergy);
        }

        private static InternalFragmentIon? ExtractSingleInternalFragment(
            PsmFromTsv psm,
            MatchedFragmentIon ion,
            double basePeakIntensity,
            double collisionEnergy,
            MsDataScan? scan)
        {
            string baseSequence = psm.BaseSeq ?? psm.FullSequence ?? string.Empty;

            var product = ion.NeutralTheoreticalProduct;
            var (internalSequence, startResidue, endResidue) = ParseInternalFragmentFromProduct(
                product, baseSequence, ion.Annotation);

            if (string.IsNullOrEmpty(internalSequence) || startResidue <= 0 || endResidue <= 0)
                return null;

            double observedMass = ion.Mz.ToMass(ion.Charge);
            double theoreticalMass = product.NeutralMass;
            double localRank = CalculateLocalIntensityRank(scan, ion.Mz, ion.Intensity);

            char nTermFlank = startResidue > 1
                ? baseSequence[startResidue - 2]
                : '-';
            char cTermFlank = endResidue < baseSequence.Length
                ? baseSequence[endResidue]
                : '-';

            return new InternalFragmentIon
            {
                PeptideSequence = baseSequence,
                InternalSequence = internalSequence,
                StartResidue = startResidue,
                EndResidue = endResidue,
                TheoreticalMass = theoreticalMass,
                ObservedMass = observedMass,
                NormalizedIntensity = ion.Intensity / basePeakIntensity,
                LocalIntensityRank = localRank,
                PrecursorCharge = psm.PrecursorCharge,
                CollisionEnergy = collisionEnergy,
                PeptidePEP = psm.PEP,
                PeptideScore = psm.Score,
                NTerminalFlankingResidue = nTermFlank,
                CTerminalFlankingResidue = cTermFlank,
                IsDecoy = psm.DecoyContamTarget?.Contains('D') ?? false,
                SourceFile = psm.FileNameWithoutExtension ?? string.Empty,
                ScanNumber = psm.Ms2ScanNumber.ToString(),
                IsIsobaricAmbiguous = false
            };
        }

        private static (string internalSequence, int startResidue, int endResidue) ParseInternalFragmentFromProduct(
            Product product,
            string peptideSequence,
            string annotation)
        {
            if (product != null && product.IsInternalFragment)
            {
                int start = product.FragmentNumber;
                int end = product.SecondaryFragmentNumber;

                if (start > 0 && end > 0 && start <= end && end <= peptideSequence.Length)
                {
                    string subseq = peptideSequence.Substring(start - 1, end - start + 1);
                    return (subseq, start, end);
                }
            }

            var match = InternalFragmentRegex.Match(annotation ?? string.Empty);
            if (match.Success &&
                int.TryParse(match.Groups[1].Value, out int startPos) &&
                int.TryParse(match.Groups[2].Value, out int endPos))
            {
                if (startPos > 0 && endPos > 0 && startPos <= endPos && endPos <= peptideSequence.Length)
                {
                    string subseq = peptideSequence.Substring(startPos - 1, endPos - startPos + 1);
                    return (subseq, startPos, endPos);
                }
            }

            return (string.Empty, 0, 0);
        }

        private static double CalculateLocalIntensityRank(MsDataScan? scan, double targetMz, double targetIntensity)
        {
            if (scan?.MassSpectrum == null)
                return double.NaN;

            double windowSize = 100.0;
            double minMz = targetMz - windowSize;
            double maxMz = targetMz + windowSize;

            var spectrum = scan.MassSpectrum;
            var mzArray = spectrum.XArray;
            var intensityArray = spectrum.YArray;

            if (mzArray == null || intensityArray == null || mzArray.Length == 0)
                return double.NaN;

            int higherCount = 0;

            for (int i = 0; i < mzArray.Length; i++)
            {
                if (mzArray[i] >= minMz && mzArray[i] <= maxMz)
                {
                    if (intensityArray[i] > targetIntensity)
                        higherCount++;
                }
            }

            return higherCount + 1;
        }

        private static Dictionary<int, MsDataScan> BuildScanLookup(MsDataFile msDataFile)
        {
            var lookup = new Dictionary<int, MsDataScan>();

            if (msDataFile?.GetAllScansList() == null)
                return lookup;

            foreach (var scan in msDataFile.GetAllScansList())
            {
                if (scan != null && !lookup.ContainsKey(scan.OneBasedScanNumber))
                    lookup[scan.OneBasedScanNumber] = scan;
            }

            return lookup;
        }
    }
}