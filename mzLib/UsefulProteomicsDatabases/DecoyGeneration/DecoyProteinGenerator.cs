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
        /// <summary>
        /// Generates decoys for a list of proteins
        /// </summary>
        /// <param name="proteins"></param>
        /// <param name="decoyType"></param>
        /// <param name="digestionParams"></param>
        /// <param name="randomSeed">Used when decoy type is shuffle for shuffling the peptides</param>
        /// <returns></returns>
        public static List<Protein> GenerateDecoys(List<Protein> proteins, DecoyType decoyType, int maxThreads = -1, string decoyIdentifier = "DECOY")
        {
            return decoyType switch
            {
                DecoyType.None => new List<Protein>(),
                DecoyType.Reverse => GenerateReverseDecoys(proteins, maxThreads, decoyIdentifier),
                DecoyType.Slide => GenerateSlideDecoys(proteins, maxThreads, decoyIdentifier),
                _ => throw new ArgumentException("Decoy type " + decoyType.ToString() + " is not implemented.")
            };
        }

        /// <summary>
        /// Generates a reverse decoy sequence
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        private static List<Protein> GenerateReverseDecoys(List<Protein> proteins, int maxThreads = -1, string decoyIdentifier = "DECOY")
        {
            List<Protein> decoyProteins = new List<Protein>();
            Parallel.ForEach(proteins, new ParallelOptions { MaxDegreeOfParallelism = maxThreads }, protein =>
            {
                // reverse sequence
                // Do not include the initiator methionine in reversal!!!
                char[] sequenceArray = protein.BaseSequence.ToCharArray();
                bool startsWithM = protein.BaseSequence.StartsWith("M", StringComparison.Ordinal);
                int[] positionMapping = GeneratePositionMapping(protein.BaseSequence, startsWithM);
                if (startsWithM)
                {
                    Array.Reverse(sequenceArray, 1, protein.BaseSequence.Length - 1);
                }
                else
                {
                    Array.Reverse(sequenceArray);
                }
                string reversedSequence = new string(sequenceArray);

                // reverse nonvariant sequence
                // Do not include the initiator methionine in reversal!!!
                char[] nonVariantSequenceArray = protein.ConsensusVariant.BaseSequence.ToCharArray();
                int[] consensusPositionMapping = GeneratePositionMapping(protein.ConsensusVariant.BaseSequence, startsWithM);
                if (protein.ConsensusVariant.BaseSequence.StartsWith("M", StringComparison.Ordinal))
                {
                    Array.Reverse(nonVariantSequenceArray, 1, protein.ConsensusVariant.BaseSequence.Length - 1);
                }
                else
                {
                    Array.Reverse(nonVariantSequenceArray);
                }
                string reversedNonVariantSequence = new string(nonVariantSequenceArray);

                // reverse modifications
                Dictionary<int, List<Modification>> decoyModifications = GetReversedModifications(protein, startsWithM);

                // reverse proteolysis products
                List<TruncationProduct> decoyPP = new List<TruncationProduct>();
                foreach (TruncationProduct pp in protein.TruncationProducts)
                {
                    // maintain lengths and approx position
                    if (startsWithM)
                    {
                        decoyPP.Add(new TruncationProduct(pp.OneBasedBeginPosition, pp.OneBasedEndPosition, $"{decoyIdentifier} {pp.Type}"));
                    }
                    else
                    {
                        decoyPP.Add(new TruncationProduct(protein.BaseSequence.Length - pp.OneBasedEndPosition + 1, protein.BaseSequence.Length - pp.OneBasedBeginPosition + 1, $"{decoyIdentifier} {pp.Type}"));
                    }
                }

                List<DisulfideBond> decoyDisulfides = new List<DisulfideBond>();
                foreach (DisulfideBond disulfideBond in protein.DisulfideBonds)
                {
                    // maintain the cysteine localizations
                    if (startsWithM)
                    {
                        decoyDisulfides.Add(new DisulfideBond(disulfideBond.OneBasedBeginPosition == 1 ? 1 : protein.BaseSequence.Length - disulfideBond.OneBasedEndPosition + 2, protein.BaseSequence.Length - disulfideBond.OneBasedBeginPosition + 2, $"{decoyIdentifier} {disulfideBond.Description}"));
                    }
                    else
                    {
                        decoyDisulfides.Add(new DisulfideBond(protein.BaseSequence.Length - disulfideBond.OneBasedEndPosition + 1, protein.BaseSequence.Length - disulfideBond.OneBasedBeginPosition + 1, $"{decoyIdentifier} {disulfideBond.Description}"));
                    }
                }

                // reverse splice sites
                List<SpliceSite> spliceSites = new List<SpliceSite>();
                foreach (SpliceSite spliceSite in protein.SpliceSites)
                {
                    // maintain the starting methionine localization
                    if (startsWithM && spliceSite.OneBasedBeginPosition == 1 && spliceSite.OneBasedEndPosition == 1)
                    {
                        spliceSites.Add(new SpliceSite(1, 1, $"{decoyIdentifier} {spliceSite.Description}"));
                    }
                    // maintain length, can't maintain localization to starting methionine in this case
                    else if (startsWithM && spliceSite.OneBasedBeginPosition == 1)
                    {
                        int end = protein.BaseSequence.Length - spliceSite.OneBasedBeginPosition + 1;
                        int begin = end - spliceSite.OneBasedEndPosition + spliceSite.OneBasedBeginPosition;
                        spliceSites.Add(new SpliceSite(begin, end, $"{decoyIdentifier} {spliceSite.Description}"));
                    }
                    else if (startsWithM)
                    {
                        spliceSites.Add(new SpliceSite(protein.BaseSequence.Length - spliceSite.OneBasedEndPosition + 2, protein.BaseSequence.Length - spliceSite.OneBasedBeginPosition + 2, $"{decoyIdentifier} {spliceSite.Description}"));
                    }
                    // maintain length and localization
                    else
                    {
                        spliceSites.Add(new SpliceSite(protein.BaseSequence.Length - spliceSite.OneBasedEndPosition + 1, protein.BaseSequence.Length - spliceSite.OneBasedBeginPosition + 1, $"{decoyIdentifier} {spliceSite.Description}"));
                    }
                }
                List<SequenceVariation> decoyVariations = CreateMappedSequenceVariations(positionMapping, protein.SequenceVariations);
                List<SequenceVariation> decoyAppliedVariations = CreateMappedSequenceVariations(consensusPositionMapping, protein.AppliedSequenceVariations);
                //List<SequenceVariation> decoyVariations = ReverseSequenceVariations(protein.SequenceVariations, protein.ConsensusVariant, reversedNonVariantSequence);
                //List<SequenceVariation> decoyAppliedVariations = ReverseSequenceVariations(protein.AppliedSequenceVariations, protein, reversedSequence);

                var decoyProtein = new Protein(
                    reversedSequence,
                    $"{decoyIdentifier}_" + protein.Accession,
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
                    spliceSites,
                    protein.DatabaseFilePath,
                    dataset: protein.DatasetEntryTag,
                    created: protein.CreatedEntryTag,
                    modified: protein.ModifiedEntryTag,
                    version: protein.VersionEntryTag,
                    xmlns: protein.XmlnsEntryTag,
                    uniProtSequenceAttributes: protein.UniProtSequenceAttributes);

                lock (decoyProteins) { decoyProteins.Add(decoyProtein); }
            });
            decoyProteins = decoyProteins.OrderBy(p => p.Accession).ToList();
            return decoyProteins;
        }
        /// <summary>
        /// Generates a mapping of amino acid positions resulting from the sequence transformation.
        /// </summary>
        /// <param name="sequence">The original sequence of the protein.</param>
        /// <param name="startsWithM">Indicates if the sequence starts with methionine.</param>
        /// <returns>An integer array mapping the original positions to the transformed positions.</returns>
        private static int[] GeneratePositionMapping(string sequence, bool startsWithM)
        {
            int length = sequence.Length;
            int[] positionMapping = new int[length + 1]; // 1-based indexing

            if (startsWithM)
            {
                // Preserve the first position (M), reverse the rest
                positionMapping[1] = 1; // M stays at position 1
                for (int i = 2; i <= length; i++)
                {
                    positionMapping[i] = length - i + 2;
                }
            }
            else
            {
                // Reverse the entire sequence
                for (int i = 1; i <= length; i++)
                {
                    positionMapping[i] = length - i + 1;
                }
            }

            return positionMapping;
        }
        /// <summary>
        /// Creates a new list of sequence variations based on the provided position mapping.
        /// </summary>
        /// <param name="positionMapping">The position mapping array (1-based indexing).</param>
        /// <param name="originalVariations">The original list of sequence variations.</param>
        /// <returns>A new list of sequence variations with updated positions and modifications.</returns>
        private static List<SequenceVariation> CreateMappedSequenceVariations(
            int[] positionMapping,
            List<SequenceVariation> originalVariations)
        {
            var newVariations = new List<SequenceVariation>();

            foreach (var originalVariation in originalVariations)
            {
                // Map the begin position using the position mapping
                int newBeginPosition = positionMapping[originalVariation.OneBasedBeginPosition];

                // Calculate the new end position
                int variationLength = originalVariation.OneBasedEndPosition - originalVariation.OneBasedBeginPosition;
                int newEndPosition = newBeginPosition + variationLength;

                // Adjust the modification dictionary
                var newModifications = new Dictionary<int, List<Modification>>();
                if (originalVariation.OneBasedModifications != null)
                {
                    foreach (var kvp in originalVariation.OneBasedModifications)
                    {
                        int newPosition = positionMapping[kvp.Key];
                        newModifications[newPosition] = kvp.Value;
                    }
                }

                // Create the new sequence variation
                var newVariation = new SequenceVariation(
                    newBeginPosition,
                    newEndPosition,
                    originalVariation.OriginalSequence,
                    originalVariation.VariantSequence,
                    originalVariation.Description,
                    originalVariation.VariantCallFormatData?.Description,
                    newModifications
                );

                // Add the new variation to the list
                newVariations.Add(newVariation);
            }

            return newVariations;
        }
        /// <summary>
        /// Extracted method to reverse modifications for a protein.
        /// </summary>
        /// <param name="protein">The protein whose modifications are being reversed.</param>
        /// <param name="startsWithM">Indicates if the protein sequence starts with methionine.</param>
        /// <returns>A dictionary of reversed modifications.</returns>
        private static Dictionary<int, List<Modification>> GetReversedModifications(Protein protein, bool startsWithM)
        {
            var reversedModifications = new Dictionary<int, List<Modification>>(protein.OneBasedPossibleLocalizedModifications.Count);

            foreach (var kvp in protein.OneBasedPossibleLocalizedModifications)
            {
                if (startsWithM)
                {
                    if (kvp.Key > 1)
                    {
                        reversedModifications.Add(protein.BaseSequence.Length - kvp.Key + 2, kvp.Value);
                    }
                    else if (kvp.Key == 1)
                    {
                        reversedModifications.Add(1, kvp.Value);
                    }
                }
                else
                {
                    reversedModifications.Add(protein.BaseSequence.Length - kvp.Key + 1, kvp.Value);
                }
            }

            return reversedModifications;
        }

        private static List<SequenceVariation> ReverseSequenceVariations(
    IEnumerable<SequenceVariation> forwardVariants,
    IBioPolymer protein,
    string reversedSequence,
    string decoyIdentifier = "DECOY")
{
    List<SequenceVariation> decoyVariations = new List<SequenceVariation>();

    static string BuildDecoyVcfTag(string decoyIdentifier, SequenceVariation src)
    {
        var baseTag = $"{decoyIdentifier} VARIANT";
        if (src?.VariantCallFormatData == null)
            return baseTag;
        var raw = src.VariantCallFormatData.Description;
        if (string.IsNullOrWhiteSpace(raw))
            raw = src.SearchableAnnotation;
        return string.IsNullOrWhiteSpace(raw) ? baseTag : $"{baseTag}: {raw}";
    }

    bool startsWithM = protein.BaseSequence.StartsWith("M", StringComparison.Ordinal);
    int seqLen = protein.BaseSequence.Length;

    foreach (SequenceVariation sv in forwardVariants)
    {
        if (sv == null)
            continue;

        string decoyVcfTag = BuildDecoyVcfTag(decoyIdentifier, sv);

        // Reverse modifications as before (not shown for brevity)
        Dictionary<int, List<Modification>> decoyVariantModifications = new Dictionary<int, List<Modification>>(sv.OneBasedModifications.Count);
        int variantSeqLength = seqLen + sv.VariantSequence.Length - sv.OriginalSequence.Length;
        bool stopGain = sv.VariantSequence.EndsWith("*", StringComparison.Ordinal);

        foreach (var kvp in sv.OneBasedModifications)
        {
            if (stopGain)
            {
                decoyVariantModifications.Add(kvp.Key, kvp.Value);
            }
            else if (startsWithM && kvp.Key > 1)
            {
                decoyVariantModifications.Add(variantSeqLength - kvp.Key + 2, kvp.Value);
            }
            else if (sv.VariantSequence.StartsWith("M", StringComparison.Ordinal) && kvp.Key == 1)
            {
                decoyVariantModifications.Add(1, kvp.Value);
            }
            else if (kvp.Key == 1)
            {
                decoyVariantModifications.Add(seqLen, kvp.Value);
            }
            else
            {
                decoyVariantModifications.Add(variantSeqLength - kvp.Key + 1, kvp.Value);
            }
        }

        char[] originalArray = sv.OriginalSequence.ToCharArray();
        char[] variationArray = sv.VariantSequence.ToCharArray();
        Array.Reverse(originalArray);
        Array.Reverse(variationArray);

        // Special handling for initiator methionine variant at position 1
        if (startsWithM && sv.OneBasedBeginPosition == 1 && sv.OriginalSequence.StartsWith("M", StringComparison.Ordinal))
        {
            decoyVariations.Add(new SequenceVariation(
                1,
                1,
                new string(originalArray),
                new string(variationArray),
                sv.Description,
                decoyVcfTag,
                decoyVariantModifications));
            continue;
        }

        // Special handling for variant at last position (C-term) in target, when M is kept at position 1 in decoy
        if (startsWithM && sv.OneBasedBeginPosition == seqLen && sv.OneBasedEndPosition == seqLen)
        {
            // Map to position 2 in decoy (since position 1 is still M)
            decoyVariations.Add(new SequenceVariation(
                2,
                2,
                new string(originalArray),
                new string(variationArray),
                sv.Description,
                decoyVcfTag,
                decoyVariantModifications));
            continue;
        }

        // All other cases: adjust as before, but account for preserved M
        int decoyBegin, decoyEnd;
        if (startsWithM)
        {
            decoyEnd = seqLen - sv.OneBasedBeginPosition + 2;
            decoyBegin = decoyEnd - originalArray.Length + 1;
        }
        else
        {
            decoyEnd = seqLen - sv.OneBasedBeginPosition + 1;
            decoyBegin = decoyEnd - originalArray.Length + 1;
        }

        decoyVariations.Add(new SequenceVariation(
            decoyBegin,
            decoyEnd,
            new string(originalArray),
            new string(variationArray),
            sv.Description,
            decoyVcfTag,
            decoyVariantModifications));
    }

    return decoyVariations;
}

        /// <summary>
        /// Generates a "slided" decoy sequence
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        private static List<Protein> GenerateSlideDecoys(List<Protein> proteins, int maxThreads = -1, string decoyIdentifier = "DECOY")
        {
            List<Protein> decoyProteins = new List<Protein>();
            Parallel.ForEach(proteins, new ParallelOptions { MaxDegreeOfParallelism = maxThreads }, protein =>
            {
                int numSlides = 20;
                char[] sequenceArrayUnslided = protein.BaseSequence.ToCharArray();
                char[] sequenceArraySlided = protein.BaseSequence.ToCharArray();

                List<DisulfideBond> decoy_disulfides_slide = new List<DisulfideBond>();
                List<SpliceSite> spliceSitesSlide = new List<SpliceSite>();
                bool initMet = protein.BaseSequence.StartsWith("M", StringComparison.Ordinal);
                Dictionary<int, List<Modification>> decoyModifications = SlideProteinSequenceWithMods(sequenceArraySlided, sequenceArrayUnslided, initMet, numSlides, protein);

                var slided_sequence = new string(sequenceArraySlided);

                List<TruncationProduct> decoyPPSlide = new List<TruncationProduct>();
                foreach (TruncationProduct pp in protein.TruncationProducts)  //can't keep all aa like you can with reverse, just keep it the same length
                {
                    decoyPPSlide.Add(pp);
                }
                foreach (DisulfideBond disulfideBond in protein.DisulfideBonds) //these actually need the same cysteines...
                {
                    decoy_disulfides_slide.Add(new DisulfideBond(GetNewSlidedIndex(disulfideBond.OneBasedBeginPosition - 1, numSlides, slided_sequence.Length, initMet) + 1, GetNewSlidedIndex(disulfideBond.OneBasedEndPosition - 1, numSlides, slided_sequence.Length, initMet) + 1, $"{decoyIdentifier} DISULFIDE BOND: " + disulfideBond.Description));
                }
                foreach (SpliceSite spliceSite in protein.SpliceSites)
                {
                    spliceSitesSlide.Add(new SpliceSite(GetNewSlidedIndex(spliceSite.OneBasedBeginPosition - 1, numSlides, slided_sequence.Length, initMet) + 1, GetNewSlidedIndex(spliceSite.OneBasedEndPosition - 1, numSlides, slided_sequence.Length, initMet) + 1, $"{decoyIdentifier} SPLICE SITE: " + spliceSite.Description));
                }

                //TODO:
                //Variants in slided and random decoys can have long reaching consequences.
                //The simplest situation (SAAV) allows for the amino acid to be substituted, but others (e.g. splicing or insertions) create new numbers or combinations of amino acids.
                //In these more complex situations, the two targets (unmodified and variant) appear largely homologous with the exception of the variant site.
                //However, the two decoys from these targets are noticeably different when the amino acids are randomized, 
                //such that the number of unique decoy peptides produced are likely to outweight the number of unique target peptides produced.
                //These issues still need to be addressed. Notably, it will be difficult to annotate the randomized variant in the decoy protein.

                //for the below code, the SAAVs will be switched in place. The downstream effects are not controlled.
                List<SequenceVariation> decoyVariationsSlide = new List<SequenceVariation>();
                foreach (SequenceVariation sv in protein.SequenceVariations)
                {
                    int numSlidesHere = numSlides;
                    char[] variationArrayUnslided = sv.VariantSequence.ToArray();
                    char[] variationArraySlided = sv.VariantSequence.ToArray();

                    //if initiator methionine, then don't move it
                    if (sv.OneBasedBeginPosition == 1 && initMet)
                    {
                        //shuffle non initiator methionine amino acids
                        if (numSlidesHere % variationArraySlided.Length == 0)
                        {
                            numSlidesHere++;
                        }
                        for (int i = 0; i < variationArraySlided.Length; i++)
                        {
                            variationArraySlided[i] = variationArrayUnslided[GetOldSlidedIndex(i, numSlidesHere, variationArrayUnslided.Length, true)];
                        }
                        decoyVariationsSlide.Add(new SequenceVariation(
                            oneBasedPosition: 1,
                            originalSequence: "M",
                            variantSequence: new string(variationArraySlided),
                            description: sv.Description,
                            variantCallFormatDataString: $"{decoyIdentifier} VARIANT: Initiator Methionine Change in " + sv.VariantCallFormatData));
                    }
                    else
                    {
                        int decoy_begin = GetNewSlidedIndex(sv.OneBasedBeginPosition - 1, numSlidesHere, sequenceArrayUnslided.Length, initMet) + 1;
                        int decoy_end = decoy_begin + sv.OneBasedEndPosition - sv.OneBasedBeginPosition;

                        //shuffle the variant sequence
                        if (numSlidesHere % variationArraySlided.Length == 0)
                        {
                            numSlidesHere++;
                        }
                        for (int i = 0; i < variationArraySlided.Length; i++)
                        {
                            variationArraySlided[i] = variationArrayUnslided[GetOldSlidedIndex(i, numSlidesHere, variationArrayUnslided.Length, initMet)];
                        }

                        decoyVariationsSlide.Add(new SequenceVariation(decoy_begin, decoy_end, sv.OriginalSequence, new string(variationArraySlided), sv.Description, $"{decoyIdentifier} VARIANT: " + sv.VariantCallFormatData));
                    }
                }
                var decoyProteinSlide = new Protein(slided_sequence, $"{decoyIdentifier}_" + protein.Accession, protein.Organism, protein.GeneNames.ToList(), decoyModifications, decoyPPSlide,
                    protein.Name, protein.FullName, true, protein.IsContaminant, null, decoyVariationsSlide, null, protein.SampleNameForVariants, decoy_disulfides_slide, spliceSitesSlide, protein.DatabaseFilePath,
                    false, protein.DatasetEntryTag, protein.CreatedEntryTag, protein.ModifiedEntryTag, protein.VersionEntryTag, protein.XmlnsEntryTag);
                lock (decoyProteins) { decoyProteins.Add(decoyProteinSlide); }
            });
            decoyProteins = decoyProteins.OrderBy(p => p.Accession).ToList();
            return decoyProteins;
        }

        private static Dictionary<int, List<Modification>> SlideProteinSequenceWithMods (char[] sequenceArraySlided, char[] sequenceArrayUnslided, bool initiatorMethionine, int numSlides, Protein protein)
        {
            // Do not include the initiator methionine in shuffle!!!
            int startIndex = initiatorMethionine ? 1 : 0;
            if (numSlides % sequenceArraySlided.Length - startIndex == 0)
            {
                numSlides++;
            }
            for (int i = startIndex; i < sequenceArraySlided.Length; i++)
            {
                sequenceArraySlided[i] = sequenceArrayUnslided[GetOldSlidedIndex(i, numSlides, protein.BaseSequence.Length, initiatorMethionine)];
            }

            Dictionary<int, List<Modification>> decoyModifications = new Dictionary<int, List<Modification>>(protein.OneBasedPossibleLocalizedModifications.Count);
            foreach (var kvp in protein.OneBasedPossibleLocalizedModifications)
            {
                if (initiatorMethionine && kvp.Key == 1)
                {
                    decoyModifications.Add(1, kvp.Value);
                }
                else
                {
                    decoyModifications.Add(GetNewSlidedIndex(kvp.Key-1, numSlides, protein.BaseSequence.Length, initiatorMethionine)+1, kvp.Value);
                }
            }

            return decoyModifications;
        }


        /// <summary>
        /// Given a new index, i, return the index of the amino acid from the unslided array
        /// </summary>
        /// <param name="i"></param>
        /// <param name="numSlides"></param>
        /// <param name="sequenceLength"></param>
        /// <param name="methioninePresent"></param>
        /// <returns></returns>
        private static int GetOldSlidedIndex(int i, int numSlides, int sequenceLength, bool methioninePresent)
        {
            if (sequenceLength > 1 && !(i == 0 && methioninePresent)) //can't shuffle a single amino acid or the initiator methionine
            {
                if (methioninePresent)
                {
                    i--;
                    sequenceLength--;
                }
                bool positiveDirection = i % 2 == 0;
                int oldIndex = i;

                if (positiveDirection)
                {
                    oldIndex += numSlides;
                }
                else
                {
                    oldIndex -= numSlides;
                }

                while (true)
                {
                    if (oldIndex < 0)
                    {
                        positiveDirection = true;
                    }
                    else if (oldIndex >= sequenceLength)
                    {
                        positiveDirection = false;
                    }
                    else
                    {
                        return methioninePresent ? oldIndex + 1 : oldIndex;
                    }

                    if (positiveDirection)
                    {
                        oldIndex = (oldIndex * -1) - 1;
                    }
                    else
                    {
                        oldIndex = (sequenceLength * 2) - oldIndex - 1;
                    }
                }
            }
            else
            {
                return i;
            }
        }


        /// <summary>
        /// Given an old index, i, return the index of the amino acid from the slided array
        /// useful for figuring out where modifications went
        /// </summary>
        /// <param name="i"></param>
        /// <param name="numSlides"></param>
        /// <param name="sequenceLength"></param>
        /// <param name="methioninePresent"></param>
        /// <returns></returns>
        private static int GetNewSlidedIndex(int i, int numSlides, int sequenceLength, bool methioninePresent)
        {
            if (sequenceLength > 1 && !(i == 0 && methioninePresent)) //can't shuffle a single amino acid or the initiator methionine
            {
                if (methioninePresent)
                {
                    i--;
                    sequenceLength--;
                }
                bool positiveDirection = i % 2 == 1;
                int newIndex = i;

                if (positiveDirection)
                {
                    newIndex += numSlides;
                }
                else
                {
                    newIndex -= numSlides;
                }

                while (true)
                {
                    if (newIndex < 0)
                    {
                        positiveDirection = true;
                    }
                    else if (newIndex >= sequenceLength)
                    {
                        positiveDirection = false;
                    }
                    else
                    {
                        return methioninePresent ? newIndex + 1 : newIndex;
                    }

                    if (positiveDirection)
                    {
                        newIndex = (newIndex * -1) - 1;
                    }
                    else
                    {
                        newIndex = (sequenceLength * 2) - newIndex - 1;
                    }
                }
            }
            else
            {
                return i;
            }
        }
    }
}