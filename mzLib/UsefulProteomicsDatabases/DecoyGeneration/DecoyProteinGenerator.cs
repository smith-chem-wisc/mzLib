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
    /// <summary>
    /// Provides static methods for generating decoy protein sequences using various strategies (e.g., reverse, slide).
    /// Decoy proteins are used for false discovery rate estimation in proteomics workflows.
    /// </summary>
    public static class DecoyProteinGenerator
    {
        /// <summary>
        /// Generates decoy proteins from a list of target proteins using the specified decoy generation strategy.
        /// </summary>
        /// <param name="proteins">List of target proteins to generate decoys from.</param>
        /// <param name="decoyType">Type of decoy generation strategy to use.</param>
        /// <param name="maxThreads">Maximum number of threads to use for parallel processing. Default is -1 (no limit).</param>
        /// <param name="decoyIdentifier">String to prepend to decoy protein accessions and annotations. Default is "DECOY".</param>
        /// <returns>List of generated decoy proteins.</returns>
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

        /// <summary>
        /// Generates decoy proteins by reversing the sequence of each target protein, optionally preserving the initiator methionine.
        /// Also reverses associated annotations and modifications.
        /// </summary>
        /// <param name="proteins">List of target proteins to generate decoys from.</param>
        /// <param name="maxThreads">Maximum number of threads to use for parallel processing.</param>
        /// <param name="decoyIdentifier">String to prepend to decoy protein accessions and annotations.</param>
        /// <returns>List of reverse-sequence decoy proteins.</returns>
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

        /// <summary>
        /// Generates a mapping from original sequence positions to their positions in the reversed sequence.
        /// Handles special logic if the sequence starts with methionine.
        /// </summary>
        /// <param name="sequence">Protein sequence to map.</param>
        /// <param name="startsWithM">Indicates if the sequence starts with methionine.</param>
        /// <returns>Array mapping original 1-based positions to reversed positions.</returns>
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

        /// <summary>
        /// Builds a decoy-specific VCF (Variant Call Format) tag for a sequence variation, ensuring it differs from the target.
        /// </summary>
        /// <param name="decoyIdentifier">String to identify the decoy.</param>
        /// <param name="src">Source sequence variation.</param>
        /// <returns>Decoy-specific VCF tag string.</returns>
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

        /// <summary>
        /// Remaps sequence variations from the target protein to the decoy protein using a position mapping.
        /// Updates variant-specific modifications and VCF tags for the decoy.
        /// </summary>
        /// <param name="positionMapping">Mapping from original to decoy sequence positions.</param>
        /// <param name="originalVariations">List of original sequence variations.</param>
        /// <param name="decoyIdentifier">String to identify the decoy.</param>
        /// <returns>List of remapped sequence variations for the decoy.</returns>
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

        /// <summary>
        /// Reverses the positions of possible localized modifications for a protein, accounting for initiator methionine if present.
        /// </summary>
        /// <param name="protein">Protein whose modifications are to be reversed.</param>
        /// <param name="startsWithM">Indicates if the sequence starts with methionine.</param>
        /// <returns>Dictionary mapping new positions to lists of modifications.</returns>
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

        /// <summary>
        /// Generates decoy proteins by sliding the sequence of each target protein by a fixed number of positions.
        /// Modifications and annotations are adjusted accordingly.
        /// </summary>
        /// <param name="proteins">List of target proteins to generate decoys from.</param>
        /// <param name="maxThreads">Maximum number of threads to use for parallel processing.</param>
        /// <param name="decoyIdentifier">String to prepend to decoy protein accessions and annotations.</param>
        /// <returns>List of slide-sequence decoy proteins.</returns>
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

        /// <summary>
        /// Slides the sequence of a protein and its modifications by a specified number of positions.
        /// Handles initiator methionine logic and updates modification positions.
        /// </summary>
        /// <param name="sequenceArraySlided">Array to store the slided sequence.</param>
        /// <param name="sequenceArrayUnslided">Original sequence array.</param>
        /// <param name="initiatorMethionine">Indicates if the sequence starts with methionine.</param>
        /// <param name="numSlides">Number of positions to slide the sequence.</param>
        /// <param name="protein">Protein whose sequence and modifications are being slided.</param>
        /// <returns>Dictionary mapping new positions to lists of modifications after sliding.</returns>
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

        /// <summary>
        /// Calculates the original index in the unslided sequence for a given index in the slided sequence.
        /// Handles initiator methionine and sequence wrapping logic.
        /// </summary>
        /// <param name="i">Index in the slided sequence.</param>
        /// <param name="numSlides">Number of positions the sequence was slided.</param>
        /// <param name="sequenceLength">Length of the sequence.</param>
        /// <param name="methioninePresent">Indicates if the sequence starts with methionine.</param>
        /// <returns>Corresponding index in the original sequence.</returns>
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

        /// <summary>
        /// Calculates the new index in the slided sequence for a given index in the original sequence.
        /// Handles initiator methionine and sequence wrapping logic.
        /// </summary>
        /// <param name="i">Index in the original sequence.</param>
        /// <param name="numSlides">Number of positions to slide the sequence.</param>
        /// <param name="sequenceLength">Length of the sequence.</param>
        /// <param name="methioninePresent">Indicates if the sequence starts with methionine.</param>
        /// <returns>Corresponding index in the slided sequence.</returns>
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