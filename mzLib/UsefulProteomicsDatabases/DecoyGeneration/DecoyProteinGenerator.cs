using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Omics.BioPolymer;
using Omics.Modifications;
using Omics;

namespace UsefulProteomicsDatabases
{
    public static class DecoyProteinGenerator
    {
        public static List<Protein> GenerateDecoys(List<Protein> proteins, DecoyType decoyType, int maxThreads = -1, string decoyIdentifier = "DECOY")
        {
            return decoyType switch
            {
                DecoyType.None => new List<Protein>(),
                DecoyType.Reverse => GenerateReverseDecoys(proteins, maxThreads, decoyIdentifier),
                DecoyType.Slide => GenerateSlideDecoys(proteins, maxThreads, decoyIdentifier),
                _ => throw new ArgumentException("Decoy type " + decoyType + " is not implemented.")
            };
        }

        private static List<Protein> GenerateReverseDecoys(List<Protein> proteins, int maxThreads = -1, string decoyIdentifier = "DECOY")
        {
            List<Protein> decoyProteins = new();
            Parallel.ForEach(proteins, new ParallelOptions { MaxDegreeOfParallelism = maxThreads }, protein =>
            {
                // Reverse sequence (keep initiator M if present)
                char[] sequenceArray = protein.BaseSequence.ToCharArray();
                bool startsWithM = protein.BaseSequence.StartsWith("M", StringComparison.Ordinal);
                int[] positionMapping = GeneratePositionMapping(protein.BaseSequence, startsWithM);
                if (startsWithM)
                {
                    Array.Reverse(sequenceArray, 1, sequenceArray.Length - 1);
                }
                else
                {
                    Array.Reverse(sequenceArray);
                }
                string reversedSequence = new(sequenceArray);

                // Reverse consensus (non‑variant) sequence
                char[] nonVariantArray = protein.ConsensusVariant.BaseSequence.ToCharArray();
                int[] consensusPositionMapping = GeneratePositionMapping(protein.ConsensusVariant.BaseSequence, startsWithM);
                if (protein.ConsensusVariant.BaseSequence.StartsWith("M", StringComparison.Ordinal))
                {
                    Array.Reverse(nonVariantArray, 1, nonVariantArray.Length - 1);
                }
                else
                {
                    Array.Reverse(nonVariantArray);
                }

                // Reverse mods
                Dictionary<int, List<Modification>> decoyModifications = GetReversedModifications(protein, startsWithM);

                // Reverse proteolysis products
                List<TruncationProduct> decoyPP = new();
                foreach (var pp in protein.TruncationProducts)
                {
                    if (startsWithM)
                    {
                        decoyPP.Add(new TruncationProduct(pp.OneBasedBeginPosition, pp.OneBasedEndPosition, $"{decoyIdentifier} {pp.Type}"));
                    }
                    else
                    {
                        decoyPP.Add(new TruncationProduct(
                            protein.BaseSequence.Length - pp.OneBasedEndPosition + 1,
                            protein.BaseSequence.Length - pp.OneBasedBeginPosition + 1,
                            $"{decoyIdentifier} {pp.Type}"));
                    }
                }

                // Reverse disulfide bonds
                List<DisulfideBond> decoyDisulfides = new();
                foreach (var bond in protein.DisulfideBonds)
                {
                    if (startsWithM)
                    {
                        decoyDisulfides.Add(new DisulfideBond(
                            bond.OneBasedBeginPosition == 1 ? 1 : protein.BaseSequence.Length - bond.OneBasedEndPosition + 2,
                            protein.BaseSequence.Length - bond.OneBasedBeginPosition + 2,
                            $"{decoyIdentifier} {bond.Description}"));
                    }
                    else
                    {
                        decoyDisulfides.Add(new DisulfideBond(
                            protein.BaseSequence.Length - bond.OneBasedEndPosition + 1,
                            protein.BaseSequence.Length - bond.OneBasedBeginPosition + 1,
                            $"{decoyIdentifier} {bond.Description}"));
                    }
                }

                // Reverse splice sites
                List<SpliceSite> decoySpliceSites = new();
                foreach (var spliceSite in protein.SpliceSites)
                {
                    if (startsWithM && spliceSite.OneBasedBeginPosition == 1 && spliceSite.OneBasedEndPosition == 1)
                    {
                        decoySpliceSites.Add(new SpliceSite(1, 1, $"{decoyIdentifier} {spliceSite.Description}"));
                    }
                    else if (startsWithM && spliceSite.OneBasedBeginPosition == 1)
                    {
                        int end = protein.BaseSequence.Length - spliceSite.OneBasedBeginPosition + 1;
                        int begin = end - spliceSite.OneBasedEndPosition + spliceSite.OneBasedBeginPosition;
                        decoySpliceSites.Add(new SpliceSite(begin, end, $"{decoyIdentifier} {spliceSite.Description}"));
                    }
                    else if (startsWithM)
                    {
                        decoySpliceSites.Add(new SpliceSite(
                            protein.BaseSequence.Length - spliceSite.OneBasedEndPosition + 2,
                            protein.BaseSequence.Length - spliceSite.OneBasedBeginPosition + 2,
                            $"{decoyIdentifier} {spliceSite.Description}"));
                    }
                    else
                    {
                        decoySpliceSites.Add(new SpliceSite(
                            protein.BaseSequence.Length - spliceSite.OneBasedEndPosition + 1,
                            protein.BaseSequence.Length - spliceSite.OneBasedBeginPosition + 1,
                            $"{decoyIdentifier} {spliceSite.Description}"));
                    }
                }

                // Map variants (target → decoy) with decoy-specific VCF annotations
                var decoyVariations = CreateMappedSequenceVariations(positionMapping, protein.SequenceVariations, decoyIdentifier);
                var decoyAppliedVariations = CreateMappedSequenceVariations(consensusPositionMapping, protein.AppliedSequenceVariations, decoyIdentifier);

                var decoyProtein = new Protein(
                    reversedSequence,
                    $"{decoyIdentifier}_{protein.Accession}",
                    protein.Organism,
                    protein.GeneNames.ToList(),
                    decoyModifications,
                    decoyPP,
                    protein.Name,
                    protein.FullName,
                    true,
                    protein.IsContaminant,
                    null,
                    decoyVariations,
                    decoyAppliedVariations,
                    protein.SampleNameForVariants,
                    decoyDisulfides,
                    decoySpliceSites,
                    protein.DatabaseFilePath,
                    dataset: protein.DatasetEntryTag,
                    created: protein.CreatedEntryTag,
                    modified: protein.ModifiedEntryTag,
                    version: protein.VersionEntryTag,
                    xmlns: protein.XmlnsEntryTag,
                    uniProtSequenceAttributes: protein.UniProtSequenceAttributes);

                lock (decoyProteins)
                {
                    decoyProteins.Add(decoyProtein);
                }
            });

            return decoyProteins.OrderBy(p => p.Accession).ToList();
        }

        private static int[] GeneratePositionMapping(string sequence, bool startsWithM)
        {
            int length = sequence.Length;
            int[] map = new int[length + 1]; // 1-based
            if (startsWithM)
            {
                map[1] = 1;
                for (int i = 2; i <= length; i++)
                {
                    map[i] = length - i + 2;
                }
            }
            else
            {
                for (int i = 1; i <= length; i++)
                {
                    map[i] = length - i + 1;
                }
            }
            return map;
        }

        // Shared helper to produce a decoy-specific VCF tag (ensures inequality vs target)
        private static string BuildDecoyVcfTag(string decoyIdentifier, SequenceVariation src)
        {
            string baseTag = $"{decoyIdentifier} VARIANT";
            if (src?.VariantCallFormatData == null)
            {
                // Target had no VCF metadata; still produce synthetic tag so decoy is not null
                return baseTag;
            }

            string raw = src.VariantCallFormatData.Description;
            if (string.IsNullOrWhiteSpace(raw))
            {
                raw = src.Description ?? src.SimpleString();
            }
            return string.IsNullOrWhiteSpace(raw) ? baseTag : $"{baseTag}: {raw}";
        }

        private static List<SequenceVariation> CreateMappedSequenceVariations(
            int[] positionMapping,
            List<SequenceVariation> originalVariations,
            string decoyIdentifier = "DECOY")
        {
            var result = new List<SequenceVariation>();
            if (originalVariations == null || originalVariations.Count == 0)
                return result;

            foreach (var ov in originalVariations)
            {
                if (ov == null)
                    continue;

                int newBegin = positionMapping[ov.OneBasedBeginPosition];
                int length = ov.OneBasedEndPosition - ov.OneBasedBeginPosition;
                int newEnd = newBegin + length;

                // Remap variant-specific modifications if any
                Dictionary<int, List<Modification>> newMods = new();
                if (ov.OneBasedModifications != null && ov.OneBasedModifications.Count > 0)
                {
                    foreach (var kv in ov.OneBasedModifications)
                    {
                        int mappedPos = positionMapping[kv.Key];
                        newMods[mappedPos] = kv.Value;
                    }
                }

                string decoyVcf = BuildDecoyVcfTag(decoyIdentifier, ov);

                var mapped = new SequenceVariation(
                    newBegin,
                    newEnd,
                    ov.OriginalSequence,
                    ov.VariantSequence,
                    ov.Description,
                    decoyVcf,
                    newMods);

                result.Add(mapped);
            }

            return result;
        }

        private static Dictionary<int, List<Modification>> GetReversedModifications(Protein protein, bool startsWithM)
        {
            var reversed = new Dictionary<int, List<Modification>>(protein.OneBasedPossibleLocalizedModifications.Count);
            foreach (var kv in protein.OneBasedPossibleLocalizedModifications)
            {
                if (startsWithM)
                {
                    if (kv.Key == 1)
                    {
                        reversed.Add(1, kv.Value);
                    }
                    else
                    {
                        reversed.Add(protein.BaseSequence.Length - kv.Key + 2, kv.Value);
                    }
                }
                else
                {
                    reversed.Add(protein.BaseSequence.Length - kv.Key + 1, kv.Value);
                }
            }
            return reversed;
        }

        private static List<Protein> GenerateSlideDecoys(List<Protein> proteins, int maxThreads = -1, string decoyIdentifier = "DECOY")
        {
            List<Protein> decoyProteins = new();
            Parallel.ForEach(proteins, new ParallelOptions { MaxDegreeOfParallelism = maxThreads }, protein =>
            {
                int numSlides = 20;
                char[] original = protein.BaseSequence.ToCharArray();
                char[] slided = protein.BaseSequence.ToCharArray();
                bool initMet = protein.BaseSequence.StartsWith("M", StringComparison.Ordinal);

                Dictionary<int, List<Modification>> decoyModifications = SlideProteinSequenceWithMods(slided, original, initMet, numSlides, protein);

                var slidedSequence = new string(slided);

                // Proteolysis products (length preserved)
                List<TruncationProduct> decoyPP = new();
                foreach (var pp in protein.TruncationProducts)
                {
                    decoyPP.Add(pp);
                }

                // Disulfides
                List<DisulfideBond> decoyDisulfides = new();
                foreach (var bond in protein.DisulfideBonds)
                {
                    decoyDisulfides.Add(new DisulfideBond(
                        GetNewSlidedIndex(bond.OneBasedBeginPosition - 1, numSlides, slidedSequence.Length, initMet) + 1,
                        GetNewSlidedIndex(bond.OneBasedEndPosition - 1, numSlides, slidedSequence.Length, initMet) + 1,
                        $"{decoyIdentifier} DISULFIDE BOND: {bond.Description}"));
                }

                // Splice sites
                List<SpliceSite> decoySpliceSites = new();
                foreach (var spliceSite in protein.SpliceSites)
                {
                    decoySpliceSites.Add(new SpliceSite(
                        GetNewSlidedIndex(spliceSite.OneBasedBeginPosition - 1, numSlides, slidedSequence.Length, initMet) + 1,
                        GetNewSlidedIndex(spliceSite.OneBasedEndPosition - 1, numSlides, slidedSequence.Length, initMet) + 1,
                        $"{decoyIdentifier} SPLICE SITE: {spliceSite.Description}"));
                }

                // Sequence variants (simple position sliding); keep initiator M logic where relevant
                List<SequenceVariation> decoyVariationsSlide = new();
                foreach (var sv in protein.SequenceVariations)
                {
                    int numSlidesHere = numSlides;
                    char[] variantSeqOriginal = sv.VariantSequence.ToCharArray();
                    char[] variantSeqSlided = sv.VariantSequence.ToCharArray();

                    if (sv.OneBasedBeginPosition == 1 && initMet)
                    {
                        if (numSlidesHere % variantSeqSlided.Length == 0) numSlidesHere++;
                        for (int i = 0; i < variantSeqSlided.Length; i++)
                        {
                            variantSeqSlided[i] = variantSeqOriginal[GetOldSlidedIndex(i, numSlidesHere, variantSeqOriginal.Length, true)];
                        }
                        decoyVariationsSlide.Add(new SequenceVariation(
                            oneBasedPosition: 1,
                            originalSequence: "M",
                            variantSequence: new string(variantSeqSlided),
                            description: sv.Description,
                            variantCallFormatDataString: $"{decoyIdentifier} VARIANT: Initiator Methionine Change in " + sv.VariantCallFormatData));
                    }
                    else
                    {
                        int decoyBegin = GetNewSlidedIndex(sv.OneBasedBeginPosition - 1, numSlidesHere, original.Length, initMet) + 1;
                        int decoyEnd = decoyBegin + (sv.OneBasedEndPosition - sv.OneBasedBeginPosition);

                        if (numSlidesHere % variantSeqSlided.Length == 0) numSlidesHere++;
                        for (int i = 0; i < variantSeqSlided.Length; i++)
                        {
                            variantSeqSlided[i] = variantSeqOriginal[GetOldSlidedIndex(i, numSlidesHere, variantSeqOriginal.Length, initMet)];
                        }

                        decoyVariationsSlide.Add(new SequenceVariation(
                            decoyBegin,
                            decoyEnd,
                            sv.OriginalSequence,
                            new string(variantSeqSlided),
                            sv.Description,
                            $"{decoyIdentifier} VARIANT: " + sv.VariantCallFormatData));
                    }
                }

                var decoyProteinSlide = new Protein(
                    slidedSequence,
                    $"{decoyIdentifier}_{protein.Accession}",
                    protein.Organism,
                    protein.GeneNames.ToList(),
                    decoyModifications,
                    decoyPP,
                    protein.Name,
                    protein.FullName,
                    true,
                    protein.IsContaminant,
                    null,
                    decoyVariationsSlide,
                    null,
                    protein.SampleNameForVariants,
                    decoyDisulfides,
                    decoySpliceSites,
                    protein.DatabaseFilePath,
                    dataset: protein.DatasetEntryTag,
                    created: protein.CreatedEntryTag,
                    modified: protein.ModifiedEntryTag,
                    version: protein.VersionEntryTag,
                    xmlns: protein.XmlnsEntryTag);

                lock (decoyProteins) { decoyProteins.Add(decoyProteinSlide); }
            });

            return decoyProteins.OrderBy(p => p.Accession).ToList();
        }

        private static Dictionary<int, List<Modification>> SlideProteinSequenceWithMods(char[] sequenceArraySlided, char[] sequenceArrayUnslided, bool initiatorMethionine, int numSlides, Protein protein)
        {
            int startIndex = initiatorMethionine ? 1 : 0;
            if (numSlides % (sequenceArraySlided.Length - startIndex) == 0) numSlides++;

            for (int i = startIndex; i < sequenceArraySlided.Length; i++)
            {
                sequenceArraySlided[i] = sequenceArrayUnslided[GetOldSlidedIndex(i, numSlides, protein.BaseSequence.Length, initiatorMethionine)];
            }

            Dictionary<int, List<Modification>> decoyMods = new(protein.OneBasedPossibleLocalizedModifications.Count);
            foreach (var kv in protein.OneBasedPossibleLocalizedModifications)
            {
                if (initiatorMethionine && kv.Key == 1)
                {
                    decoyMods.Add(1, kv.Value);
                }
                else
                {
                    decoyMods.Add(GetNewSlidedIndex(kv.Key - 1, numSlides, protein.BaseSequence.Length, initiatorMethionine) + 1, kv.Value);
                }
            }
            return decoyMods;
        }

        private static int GetOldSlidedIndex(int i, int numSlides, int sequenceLength, bool methioninePresent)
        {
            if (sequenceLength <= 1 || (i == 0 && methioninePresent))
                return i;

            if (methioninePresent)
            {
                i--;
                sequenceLength--;
            }

            bool forward = i % 2 == 0;
            int oldIndex = i;
            oldIndex += forward ? numSlides : -numSlides;

            while (true)
            {
                if (oldIndex < 0) forward = true;
                else if (oldIndex >= sequenceLength) forward = false;
                else return methioninePresent ? oldIndex + 1 : oldIndex;

                oldIndex = forward
                    ? (oldIndex * -1) - 1
                    : (sequenceLength * 2) - oldIndex - 1;
            }
        }

        private static int GetNewSlidedIndex(int i, int numSlides, int sequenceLength, bool methioninePresent)
        {
            if (sequenceLength <= 1 || (i == 0 && methioninePresent))
                return i;

            if (methioninePresent)
            {
                i--;
                sequenceLength--;
            }

            bool forward = i % 2 == 1;
            int newIndex = i;
            newIndex += forward ? numSlides : -numSlides;

            while (true)
            {
                if (newIndex < 0) forward = true;
                else if (newIndex >= sequenceLength) forward = false;
                else return methioninePresent ? newIndex + 1 : newIndex;

                newIndex = forward
                    ? (newIndex * -1) - 1
                    : (sequenceLength * 2) - newIndex - 1;
            }
        }
    }
}