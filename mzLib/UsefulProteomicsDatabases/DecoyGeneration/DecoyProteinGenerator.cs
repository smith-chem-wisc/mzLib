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
                Dictionary<int, List<Modification>> decoyModifications = null;
                if (startsWithM)
                {
                    decoyModifications = new Dictionary<int, List<Modification>>(protein.OneBasedPossibleLocalizedModifications.Count);
                    foreach (var kvp in protein.OneBasedPossibleLocalizedModifications)
                    {
                        if (kvp.Key > 1)
                        {
                            decoyModifications.Add(protein.BaseSequence.Length - kvp.Key + 2, kvp.Value);
                        }
                        else if (kvp.Key == 1)
                        {
                            decoyModifications.Add(1, kvp.Value);
                        }
                    }
                }
                else
                {
                    decoyModifications = new Dictionary<int, List<Modification>>(protein.OneBasedPossibleLocalizedModifications.Count);
                    foreach (var kvp in protein.OneBasedPossibleLocalizedModifications)
                    {
                        decoyModifications.Add(protein.BaseSequence.Length - kvp.Key + 1, kvp.Value);
                    }
                }

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

                List<SequenceVariation> decoyVariations = ReverseSequenceVariations(protein.SequenceVariations, protein.ConsensusVariant, reversedNonVariantSequence);
                List<SequenceVariation> decoyAppliedVariations = ReverseSequenceVariations(protein.AppliedSequenceVariations, protein, reversedSequence);

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

        private static List<SequenceVariation> ReverseSequenceVariations(IEnumerable<SequenceVariation> forwardVariants, IBioPolymer protein, string reversedSequence, string decoyIdentifier = "DECOY")
        {
            List<SequenceVariation> decoyVariations = new List<SequenceVariation>();
            foreach (SequenceVariation sv in forwardVariants)
            {
                // place reversed modifications (referencing variant sequence location)
                Dictionary<int, List<Modification>> decoyVariantModifications = new Dictionary<int, List<Modification>>(sv.OneBasedModifications.Count);
                int variantSeqLength = protein.BaseSequence.Length + sv.VariantSequence.Length - sv.OriginalSequence.Length;
                bool startsWithM = protein.BaseSequence.StartsWith("M", StringComparison.Ordinal);
                bool stopGain = sv.VariantSequence.EndsWith("*");
                foreach (var kvp in sv.OneBasedModifications)
                {
                    // keeping positions for stop gain to make decoys with same length
                    if (stopGain)
                    {
                        decoyVariantModifications.Add(kvp.Key, kvp.Value);
                    }
                    // methionine retention but rest reversed
                    if (startsWithM && kvp.Key > 1)
                    {
                        decoyVariantModifications.Add(variantSeqLength - kvp.Key + 2, kvp.Value);
                    }
                    // on starting methionine
                    else if (sv.VariantSequence.StartsWith("M", StringComparison.Ordinal) && kvp.Key == 1)
                    {
                        decoyVariantModifications.Add(1, kvp.Value);
                    }
                    // on starting non-methionine
                    else if (kvp.Key == 1)
                    {
                        decoyVariantModifications.Add(protein.BaseSequence.Length, kvp.Value);
                    }
                    else
                    {
                        decoyVariantModifications.Add(variantSeqLength - kvp.Key + 1, kvp.Value);
                    }
                }

                // reverse sequence variant
                char[] originalArray = sv.OriginalSequence.ToArray();
                char[] variationArray = sv.VariantSequence.ToArray();
                int decoyEnd = protein.BaseSequence.Length - sv.OneBasedBeginPosition + 2 + Convert.ToInt32(sv.OneBasedEndPosition == reversedSequence.Length) - Convert.ToInt32(sv.OneBasedBeginPosition == 1);
                int decoyBegin = decoyEnd - originalArray.Length + 1;
                Array.Reverse(originalArray);
                Array.Reverse(variationArray);

                bool originalInitMet = sv.OneBasedBeginPosition == 1 && sv.OriginalSequence.StartsWith("M", StringComparison.Ordinal);
                bool variantInitMet = sv.OneBasedBeginPosition == 1 && sv.VariantSequence.StartsWith("M", StringComparison.Ordinal);
                bool startLoss = originalInitMet && !variantInitMet;

                // stop gains should still produce decoys with the same length
                if (stopGain)
                {
                    decoyVariations.Add(new SequenceVariation(sv.OneBasedBeginPosition,
                        reversedSequence.Substring(sv.OneBasedBeginPosition - 1, sv.OneBasedEndPosition - sv.OneBasedBeginPosition + 1),
                        new string(variationArray).Substring(1, variationArray.Length - 1) + variationArray[0],
                        $"{decoyIdentifier} VARIANT: " + sv.Description, decoyVariantModifications));
                }
                // start loss, so the variant is at the end
                else if (startLoss)
                {
                    decoyVariations.Add(new SequenceVariation(protein.BaseSequence.Length - sv.OneBasedEndPosition + 2, protein.BaseSequence.Length, new string(originalArray).Substring(0, originalArray.Length - 1), new string(variationArray), $"{decoyIdentifier} VARIANT: " + sv.Description, decoyVariantModifications));
                }
                // both start with M, but there's more
                else if (sv.VariantSequence.StartsWith("M", StringComparison.Ordinal) && sv.OneBasedBeginPosition == 1 && (sv.OriginalSequence.Length > 1 || sv.VariantSequence.Length > 1))
                {
                    string original = new string(originalArray).Substring(0, originalArray.Length - 1);
                    string variant = new string(variationArray).Substring(0, variationArray.Length - 1);
                    decoyVariations.Add(new SequenceVariation(protein.BaseSequence.Length - sv.OneBasedEndPosition + 2, protein.BaseSequence.Length, original, variant, $"{decoyIdentifier} VARIANT: " + sv.Description, decoyVariantModifications));
                }
                // gained an initiating methionine
                else if (sv.VariantSequence.StartsWith("M", StringComparison.Ordinal) && sv.OneBasedBeginPosition == 1)
                {
                    decoyVariations.Add(new SequenceVariation(1, 1, new string(originalArray), new string(variationArray), $"{decoyIdentifier} VARIANT: " + sv.Description, decoyVariantModifications));
                }
                // starting methionine, but no variations on it
                else if (startsWithM)
                {
                    decoyVariations.Add(new SequenceVariation(protein.BaseSequence.Length - sv.OneBasedEndPosition + 2, protein.BaseSequence.Length - sv.OneBasedBeginPosition + 2, new string(originalArray), new string(variationArray), $"{decoyIdentifier} VARIANT: " + sv.Description, decoyVariantModifications));
                }
                // no starting methionine
                else
                {
                    decoyVariations.Add(new SequenceVariation(protein.BaseSequence.Length - sv.OneBasedEndPosition + 1, protein.BaseSequence.Length - sv.OneBasedBeginPosition + 1, new string(originalArray), new string(variationArray), $"{decoyIdentifier} VARIANT: " + sv.Description, decoyVariantModifications));
                }
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
                        decoyVariationsSlide.Add(new SequenceVariation(1, "M", new string(variationArraySlided), $"{decoyIdentifier} VARIANT: Initiator Methionine Change in " + sv.Description));
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

                        decoyVariationsSlide.Add(new SequenceVariation(decoy_begin, decoy_end, sv.OriginalSequence, new string(variationArraySlided), $"{decoyIdentifier} VARIANT: " + sv.Description));
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