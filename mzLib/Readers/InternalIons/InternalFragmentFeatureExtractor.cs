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

        private static readonly Dictionary<char, double> KyteDoolittle = new()
        {
            {'I',  4.5}, {'V',  4.2}, {'L',  3.8}, {'F',  2.8},
            {'C',  2.5}, {'M',  1.9}, {'A',  1.8}, {'G', -0.4},
            {'W', -0.9}, {'T', -0.7}, {'S', -0.8}, {'Y', -1.3},
            {'P', -1.6}, {'H', -3.2}, {'D', -3.5}, {'N', -3.5},
            {'E', -3.5}, {'Q', -3.5}, {'K', -3.9}, {'R', -4.5}
        };

        public static List<InternalFragmentIon> ExtractFromPsms(
            List<PsmFromTsv> psms, MsDataFile msDataFile, double defaultCollisionEnergy = double.NaN)
        {
            if (psms == null || psms.Count == 0) return new List<InternalFragmentIon>();

            var results = new List<InternalFragmentIon>();
            var scanLookup = BuildScanLookup(msDataFile);

            foreach (var psm in psms)
            {
                if (psm.MatchedIons == null || psm.MatchedIons.Count == 0) continue;

                var internalIons = psm.MatchedIons.Where(ion => ion.IsInternalFragment).ToList();
                if (internalIons.Count == 0) continue;

                double basePeakIntensity = psm.MatchedIons.Max(ion => ion.Intensity);
                if (basePeakIntensity <= 0) basePeakIntensity = 1.0;

                MsDataScan? scan = null;
                if (scanLookup.TryGetValue(psm.Ms2ScanNumber, out var foundScan)) scan = foundScan;

                double tic = GetTotalIonCurrent(scan, psm.MatchedIons);

                double collisionEnergy = double.NaN;
                if (scan?.HcdEnergy != null && TryParseCollisionEnergy(scan.HcdEnergy, out double ce))
                    collisionEnergy = ce;
                else
                    collisionEnergy = defaultCollisionEnergy;

                string fullModifiedSequence = psm.FullSequence ?? string.Empty;
                var modificationsByPosition = ParseModificationPositions(fullModifiedSequence);
                string baseSequence = psm.BaseSeq ?? psm.FullSequence ?? string.Empty;
                int peptideLength = baseSequence.Length;

                var bIonLookup = BuildTerminalIonLookupTic(psm.MatchedIons, "b", tic);
                var yIonLookup = BuildTerminalIonLookupTic(psm.MatchedIons, "y", tic);

                var psmInternalFragments = new List<InternalFragmentIon>();

                foreach (var ion in internalIons)
                {
                    var internalFragment = ExtractSingleInternalFragment(
                        psm, ion, basePeakIntensity, tic, collisionEnergy, scan,
                        fullModifiedSequence, modificationsByPosition,
                        peptideLength, bIonLookup, yIonLookup, baseSequence);

                    if (internalFragment != null) psmInternalFragments.Add(internalFragment);
                }

                MarkIsobaricAmbiguousIons(psmInternalFragments);
                results.AddRange(psmInternalFragments);
            }

            return results;
        }

        private static double GetTotalIonCurrent(MsDataScan? scan, List<MatchedFragmentIon> matchedIons)
        {
            if (scan?.TotalIonCurrent > 0) return scan.TotalIonCurrent;
            if (scan?.MassSpectrum?.YArray != null && scan.MassSpectrum.YArray.Length > 0)
                return scan.MassSpectrum.YArray.Sum();
            double matchedSum = matchedIons.Sum(ion => ion.Intensity);
            return matchedSum > 0 ? matchedSum : 1.0;
        }

        private static Dictionary<int, double> BuildTerminalIonLookupTic(
            List<MatchedFragmentIon> matchedIons, string ionType, double tic)
        {
            var lookup = new Dictionary<int, double>();
            if (tic <= 0) return lookup;

            foreach (var ion in matchedIons)
            {
                if (ion.IsInternalFragment) continue;
                var product = ion.NeutralTheoreticalProduct;
                if (product == null) continue;

                string annotation = ion.Annotation ?? string.Empty;
                bool isTargetType = (ionType == "b" && (product.ProductType == ProductType.b ||
                    annotation.StartsWith("b", StringComparison.OrdinalIgnoreCase))) ||
                    (ionType == "y" && (product.ProductType == ProductType.y ||
                    annotation.StartsWith("y", StringComparison.OrdinalIgnoreCase)));

                if (!isTargetType) continue;

                int ionNumber = product.FragmentNumber;
                double ticNormalized = ion.Intensity / tic;

                if (!lookup.ContainsKey(ionNumber) || lookup[ionNumber] < ticNormalized)
                    lookup[ionNumber] = ticNormalized;
            }
            return lookup;
        }

        private static int CountBasicResidues(string sequence, int startIdx, int endIdx)
        {
            int count = 0;
            for (int i = startIdx; i <= endIdx && i < sequence.Length; i++)
            {
                char c = sequence[i];
                if (c == 'K' || c == 'R' || c == 'H') count++;
            }
            return count;
        }

        private static InternalFragmentIon? ExtractSingleInternalFragment(
            PsmFromTsv psm, MatchedFragmentIon ion, double basePeakIntensity, double tic,
            double collisionEnergy, MsDataScan? scan, string fullModifiedSequence,
            Dictionary<int, string> modificationsByPosition, int peptideLength,
            Dictionary<int, double> bIonLookup, Dictionary<int, double> yIonLookup, string baseSequence)
        {
            var product = ion.NeutralTheoreticalProduct;
            var (internalSequence, startResidue, endResidue) = ParseInternalFragmentFromProduct(
                product, baseSequence, ion.Annotation);

            if (string.IsNullOrEmpty(internalSequence) || startResidue <= 0 || endResidue <= 0)
                return null;

            double observedMass = ion.Mz.ToMass(ion.Charge);
            double theoreticalMass = product.NeutralMass;
            double localRank = CalculateLocalIntensityRank(scan, ion.Mz, ion.Intensity);

            char nTermFlank = startResidue > 1 ? baseSequence[startResidue - 2] : '-';
            char cTermFlank = endResidue < baseSequence.Length ? baseSequence[endResidue] : '-';

            string modsInFragment = GetModificationsInRange(modificationsByPosition, startResidue, endResidue);

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
            int distanceFromCTerm = peptideLength - endResidue;
            double maxTerminalIonIntensity = Math.Max(bIonIntensity, yIonIntensity);
            bool hasBothTerminalIons = hasMatchedB && hasMatchedY;

            int basicInBSpan = startResidue > 1 ? CountBasicResidues(baseSequence, 0, startResidue - 2) : 0;
            int basicInYSpan = endResidue < peptideLength ? CountBasicResidues(baseSequence, endResidue, peptideLength - 1) : 0;

            // New scorer features
            bool isProlineAtInternalNTerm = !string.IsNullOrEmpty(internalSequence) && internalSequence[0] == 'P';
            bool isTerminalRescue = distanceFromCTerm <= 3;
            double nTermFlankHydrophobicity = KyteDoolittle.TryGetValue(nTermFlank, out double hScore) ? hScore : 0.0;

            return new InternalFragmentIon
            {
                PeptideSequence = baseSequence,
                InternalSequence = internalSequence,
                StartResidue = startResidue,
                EndResidue = endResidue,
                TheoreticalMass = theoreticalMass,
                ObservedMass = observedMass,
                NormalizedIntensity = ion.Intensity / basePeakIntensity,
                TicNormalizedIntensity = tic > 0 ? ion.Intensity / tic : 0,
                TotalIonCurrent = tic,
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
                BYProductScore = byProductScore,
                DistanceFromCTerm = distanceFromCTerm,
                MaxTerminalIonIntensity = maxTerminalIonIntensity,
                HasBothTerminalIons = hasBothTerminalIons,
                BasicResiduesInBIonSpan = basicInBSpan,
                BasicResiduesInYIonSpan = basicInYSpan,
                IsProlineAtInternalNTerminus = isProlineAtInternalNTerm,
                IsTerminalRescue = isTerminalRescue,
                NTerminalFlankingHydrophobicity = nTermFlankHydrophobicity
            };
        }

        #region Helper Methods

        private static Dictionary<int, string> ParseModificationPositions(string fullModifiedSequence)
        {
            var mods = new Dictionary<int, string>();
            if (string.IsNullOrEmpty(fullModifiedSequence)) return mods;
            int residuePosition = 0, i = 0;
            while (i < fullModifiedSequence.Length)
            {
                char c = fullModifiedSequence[i];
                if (c == '[')
                {
                    int close = fullModifiedSequence.IndexOf(']', i);
                    if (close > i)
                    {
                        if (residuePosition > 0) mods[residuePosition] = fullModifiedSequence.Substring(i + 1, close - i - 1);
                        i = close + 1;
                    }
                    else i++;
                }
                else if (char.IsLetter(c) && char.IsUpper(c)) { residuePosition++; i++; }
                else i++;
            }
            return mods;
        }

        private static string GetModificationsInRange(Dictionary<int, string> mods, int start, int end) =>
            string.Join("; ", mods.Where(kvp => kvp.Key >= start && kvp.Key <= end)
                .OrderBy(kvp => kvp.Key).Select(kvp => $"{kvp.Value} at position {kvp.Key}"));

        private static void MarkIsobaricAmbiguousIons(List<InternalFragmentIon> ions)
        {
            if (ions.Count <= 1) return;
            foreach (var group in ions.GroupBy(ion => Math.Round(ion.TheoreticalMass, 4)))
                if (group.Count() > 1)
                    foreach (var ion in group) ion.IsIsobaricAmbiguous = true;
        }

        private static bool TryParseCollisionEnergy(string hcdEnergy, out double ce)
        {
            ce = double.NaN;
            if (string.IsNullOrWhiteSpace(hcdEnergy)) return false;
            return double.TryParse(hcdEnergy.Replace("@", "").Replace("HCD", "").Replace("hcd", "").Trim(), out ce);
        }

        private static (string, int, int) ParseInternalFragmentFromProduct(Product product, string peptideSeq, string annotation)
        {
            if (product != null && product.IsInternalFragment)
            {
                int s = product.FragmentNumber, e = product.SecondaryFragmentNumber;
                if (s > 0 && e > 0 && s <= e && e <= peptideSeq.Length)
                    return (peptideSeq.Substring(s - 1, e - s + 1), s, e);
            }
            var match = InternalFragmentRegex.Match(annotation ?? "");
            if (match.Success && int.TryParse(match.Groups[1].Value, out int sp) && int.TryParse(match.Groups[2].Value, out int ep))
                if (sp > 0 && ep > 0 && sp <= ep && ep <= peptideSeq.Length)
                    return (peptideSeq.Substring(sp - 1, ep - sp + 1), sp, ep);
            return ("", 0, 0);
        }

        private static double CalculateLocalIntensityRank(MsDataScan? scan, double mz, double intensity)
        {
            if (scan?.MassSpectrum == null) return double.NaN;
            var spec = scan.MassSpectrum;
            if (spec.XArray == null || spec.YArray == null) return double.NaN;
            int higher = 0;
            for (int i = 0; i < spec.XArray.Length; i++)
                if (spec.XArray[i] >= mz - 100 && spec.XArray[i] <= mz + 100 && spec.YArray[i] > intensity) higher++;
            return higher + 1;
        }

        private static Dictionary<int, MsDataScan> BuildScanLookup(MsDataFile msDataFile)
        {
            var lookup = new Dictionary<int, MsDataScan>();
            if (msDataFile?.GetAllScansList() == null) return lookup;
            foreach (var scan in msDataFile.GetAllScansList())
                if (scan != null && !lookup.ContainsKey(scan.OneBasedScanNumber))
                    lookup[scan.OneBasedScanNumber] = scan;
            return lookup;
        }

        #endregion
    }
}