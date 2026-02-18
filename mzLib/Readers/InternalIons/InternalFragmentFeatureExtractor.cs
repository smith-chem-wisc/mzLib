using Chemistry;
using MassSpectrometry;
using Omics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

namespace Readers.InternalIons
{
    public static class InternalFragmentFeatureExtractor
    {
        private static readonly Regex InternalFragmentRegex = new(@"[yb]I[yb]\[(\d+)-(\d+)\]",
            RegexOptions.IgnoreCase | RegexOptions.Compiled);

        private static bool _warnedAboutMissingMatchedIons = false;

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

                double collisionEnergy = double.NaN;
                if (scan?.HcdEnergy != null && TryParseCollisionEnergy(scan.HcdEnergy, out double ce))
                {
                    collisionEnergy = ce;
                }
                else
                {
                    collisionEnergy = defaultCollisionEnergy;
                }

                string fullModifiedSequence = psm.FullSequence ?? string.Empty;
                var modificationsByPosition = ParseModificationPositions(fullModifiedSequence);
                string baseSequence = psm.BaseSeq ?? psm.FullSequence ?? string.Empty;
                int peptideLength = baseSequence.Length;

                // Build lookup for b and y ions
                var bIonLookup = BuildTerminalIonLookup(psm.MatchedIons, "b", basePeakIntensity);
                var yIonLookup = BuildTerminalIonLookup(psm.MatchedIons, "y", basePeakIntensity);

                var psmInternalFragments = new List<InternalFragmentIon>();

                foreach (var ion in internalIons)
                {
                    var internalFragment = ExtractSingleInternalFragment(
                        psm, ion, basePeakIntensity, collisionEnergy, scan,
                        fullModifiedSequence, modificationsByPosition,
                        peptideLength, bIonLookup, yIonLookup);

                    if (internalFragment != null)
                        psmInternalFragments.Add(internalFragment);
                }

                MarkIsobaricAmbiguousIons(psmInternalFragments);
                results.AddRange(psmInternalFragments);
            }

            return results;
        }

        /// <summary>
        /// Builds a lookup dictionary for terminal (b or y) ions.
        /// Key: ion number, Value: normalized intensity (max across charge states).
        /// </summary>
        private static Dictionary<int, double> BuildTerminalIonLookup(
            List<MatchedFragmentIon> matchedIons,
            string ionType,
            double basePeakIntensity)
        {
            var lookup = new Dictionary<int, double>();

            foreach (var ion in matchedIons)
            {
                if (ion.IsInternalFragment)
                    continue;

                var product = ion.NeutralTheoreticalProduct;
                if (product == null)
                    continue;

                // Check if this is the target ion type
                string annotation = ion.Annotation ?? string.Empty;
                bool isTargetType = false;

                if (ionType == "b" && (product.ProductType == ProductType.b ||
                    annotation.StartsWith("b", StringComparison.OrdinalIgnoreCase)))
                {
                    isTargetType = true;
                }
                else if (ionType == "y" && (product.ProductType == ProductType.y ||
                    annotation.StartsWith("y", StringComparison.OrdinalIgnoreCase)))
                {
                    isTargetType = true;
                }

                if (!isTargetType)
                    continue;

                int ionNumber = product.FragmentNumber;
                double normalizedIntensity = ion.Intensity / basePeakIntensity;

                // Keep the maximum intensity for this ion number (across charge states)
                if (!lookup.ContainsKey(ionNumber) || lookup[ionNumber] < normalizedIntensity)
                {
                    lookup[ionNumber] = normalizedIntensity;
                }
            }

            return lookup;
        }

        private static Dictionary<int, string> ParseModificationPositions(string fullModifiedSequence)
        {
            var mods = new Dictionary<int, string>();
            if (string.IsNullOrEmpty(fullModifiedSequence))
                return mods;

            int residuePosition = 0;
            int i = 0;

            while (i < fullModifiedSequence.Length)
            {
                char c = fullModifiedSequence[i];

                if (c == '[')
                {
                    int closeBracket = fullModifiedSequence.IndexOf(']', i);
                    if (closeBracket > i)
                    {
                        string modContent = fullModifiedSequence.Substring(i + 1, closeBracket - i - 1);
                        if (residuePosition > 0)
                        {
                            mods[residuePosition] = modContent;
                        }
                        i = closeBracket + 1;
                    }
                    else
                    {
                        i++;
                    }
                }
                else if (char.IsLetter(c) && char.IsUpper(c))
                {
                    residuePosition++;
                    i++;
                }
                else
                {
                    i++;
                }
            }

            return mods;
        }

        private static string GetModificationsInRange(
            Dictionary<int, string> modificationsByPosition,
            int startResidue,
            int endResidue)
        {
            var modsInRange = modificationsByPosition
                .Where(kvp => kvp.Key >= startResidue && kvp.Key <= endResidue)
                .OrderBy(kvp => kvp.Key)
                .Select(kvp => $"{kvp.Value} at position {kvp.Key}")
                .ToList();

            return string.Join("; ", modsInRange);
        }

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
            MsDataScan? scan,
            string fullModifiedSequence,
            Dictionary<int, string> modificationsByPosition,
            int peptideLength,
            Dictionary<int, double> bIonLookup,
            Dictionary<int, double> yIonLookup)
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

            string modsInFragment = GetModificationsInRange(modificationsByPosition, startResidue, endResidue);

            // ========== B/Y ION CORRELATION ==========
            // For internal fragment at positions [startResidue, endResidue]:
            // - The N-terminal cleavage produces b_(startResidue-1)
            // - The C-terminal cleavage produces y_(peptideLength - endResidue)

            int bIonNumber = startResidue - 1;
            int yIonNumber = peptideLength - endResidue;

            double bIonIntensity = 0.0;
            bool hasMatchedB = false;
            if (bIonNumber > 0 && bIonLookup.TryGetValue(bIonNumber, out double bInt))
            {
                bIonIntensity = bInt;
                hasMatchedB = true;
            }

            double yIonIntensity = 0.0;
            bool hasMatchedY = false;
            if (yIonNumber > 0 && yIonLookup.TryGetValue(yIonNumber, out double yInt))
            {
                yIonIntensity = yInt;
                hasMatchedY = true;
            }

            double byProductScore = bIonIntensity * yIonIntensity;

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
                IsIsobaricAmbiguous = false,
                FullModifiedSequence = fullModifiedSequence,
                ModificationsInInternalFragment = modsInFragment,
                BIonIntensityAtNTerm = bIonIntensity,
                YIonIntensityAtCTerm = yIonIntensity,
                HasMatchedBIonAtNTerm = hasMatchedB,
                HasMatchedYIonAtCTerm = hasMatchedY,
                BYProductScore = byProductScore
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