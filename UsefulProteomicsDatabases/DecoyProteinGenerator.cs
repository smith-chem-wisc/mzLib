using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

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
        public static List<Protein> GenerateDecoys(List<Protein> proteins, DecoyType decoyType, int maxThreads = -1)
        {
            if (decoyType == DecoyType.None)
            {
                return new List<Protein>();
            }
            else if (decoyType == DecoyType.Reverse)
            {
                return GenerateReverseDecoys(proteins, maxThreads);
            }
            else if (decoyType == DecoyType.Slide)
            {
                return GenerateSlideDecoys(proteins, maxThreads);
            }
            else
            {
                throw new ArgumentException("Decoy type " + decoyType.ToString() + " is not implemented.");
            }
        }

        /// <summary>
        /// Generates a reverse decoy sequence
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        private static List<Protein> GenerateReverseDecoys(List<Protein> proteins, int maxThreads = -1)
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
                char[] nonVariantSequenceArray = protein.NonVariantProtein.BaseSequence.ToCharArray();
                if (protein.NonVariantProtein.BaseSequence.StartsWith("M", StringComparison.Ordinal))
                {
                    Array.Reverse(nonVariantSequenceArray, 1, protein.BaseSequence.Length - 1);
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
                List<ProteolysisProduct> decoyPP = new List<ProteolysisProduct>();
                foreach (ProteolysisProduct pp in protein.ProteolysisProducts)
                {
                    // maintain lengths and approx position
                    if (startsWithM)
                    {
                        decoyPP.Add(new ProteolysisProduct(pp.OneBasedBeginPosition, pp.OneBasedEndPosition, $"DECOY {pp.Type}"));
                    }
                    else
                    {
                        decoyPP.Add(new ProteolysisProduct(protein.BaseSequence.Length - pp.OneBasedEndPosition + 1, protein.BaseSequence.Length - pp.OneBasedBeginPosition + 1, $"DECOY {pp.Type}"));
                    }
                }

                List<DisulfideBond> decoyDisulfides = new List<DisulfideBond>();
                foreach (DisulfideBond disulfideBond in protein.DisulfideBonds)
                {
                    // maintain the cysteine localizations
                    if (startsWithM)
                    {
                        decoyDisulfides.Add(new DisulfideBond(disulfideBond.OneBasedBeginPosition == 1 ? 1 : protein.BaseSequence.Length - disulfideBond.OneBasedEndPosition + 2, protein.BaseSequence.Length - disulfideBond.OneBasedBeginPosition + 2, $"DECOY {disulfideBond.Description}"));
                    }
                    else
                    {
                        decoyDisulfides.Add(new DisulfideBond(protein.BaseSequence.Length - disulfideBond.OneBasedEndPosition + 1, protein.BaseSequence.Length - disulfideBond.OneBasedBeginPosition + 1, $"DECOY {disulfideBond.Description}"));
                    }
                }

                // reverse splice sites
                List<SpliceSite> spliceSites = new List<SpliceSite>();
                foreach (SpliceSite spliceSite in protein.SpliceSites)
                {
                    // maintain the starting methionine localization
                    if (startsWithM && spliceSite.OneBasedBeginPosition == 1 && spliceSite.OneBasedEndPosition == 1)
                    {
                        spliceSites.Add(new SpliceSite(1, 1, $"DECOY {spliceSite.Description}"));
                    }
                    // maintain length, can't maintain localization to starting methionine in this case
                    else if (startsWithM && spliceSite.OneBasedBeginPosition == 1)
                    {
                        int end = protein.BaseSequence.Length - spliceSite.OneBasedBeginPosition + 1;
                        int begin = end - spliceSite.OneBasedEndPosition + spliceSite.OneBasedBeginPosition;
                        spliceSites.Add(new SpliceSite(begin, end, $"DECOY {spliceSite.Description}"));
                    }
                    else if (startsWithM)
                    {
                        spliceSites.Add(new SpliceSite(protein.BaseSequence.Length - spliceSite.OneBasedEndPosition + 2, protein.BaseSequence.Length - spliceSite.OneBasedBeginPosition + 2, $"DECOY {spliceSite.Description}"));
                    }
                    // maintain length and localization
                    else
                    {
                        spliceSites.Add(new SpliceSite(protein.BaseSequence.Length - spliceSite.OneBasedEndPosition + 1, protein.BaseSequence.Length - spliceSite.OneBasedBeginPosition + 1, $"DECOY {spliceSite.Description}"));
                    }
                }

                List<SequenceVariation> decoyVariations = ReverseSequenceVariations(protein.SequenceVariations, protein.NonVariantProtein, reversedNonVariantSequence);
                List<SequenceVariation> decoyAppliedVariations = ReverseSequenceVariations(protein.AppliedSequenceVariations, protein, reversedSequence);

                var decoyProtein = new Protein(
                    reversedSequence, 
                    "DECOY_" + protein.Accession, 
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
                    protein.DatabaseFilePath);

                lock (decoyProteins) { decoyProteins.Add(decoyProtein); }
            });
            return decoyProteins;
        }

        private static List<SequenceVariation> ReverseSequenceVariations(IEnumerable<SequenceVariation> forwardVariants, Protein protein, string reversedSequence)
        {
            List<SequenceVariation> decoyVariations = new List<SequenceVariation>();
            foreach (SequenceVariation sv in forwardVariants)
            {
                // place reversed modifications (referencing variant sequence location)
                Dictionary<int, List<Modification>> decoyVariantModifications = new Dictionary<int, List<Modification>>(sv.OneBasedModifications.Count);
                int variantSeqLength = protein.BaseSequence.Length + sv.VariantSequence.Length - sv.OriginalSequence.Length;
                bool startsWithM = protein.BaseSequence.StartsWith("M", StringComparison.Ordinal);
                foreach (var kvp in sv.OneBasedModifications)
                {
                    if (startsWithM && kvp.Key > 1)
                    {
                        decoyVariantModifications.Add(variantSeqLength - kvp.Key + 2, kvp.Value);
                    }
                    else if (sv.VariantSequence.StartsWith("M", StringComparison.Ordinal) && kvp.Key == 1)
                    {
                        decoyVariantModifications.Add(1, kvp.Value);
                    }
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

                // start loss, so the variant is at the end
                bool originalInitMet = sv.OriginalSequence.StartsWith("M", StringComparison.Ordinal);
                bool variantInitMet = sv.VariantSequence.StartsWith("M", StringComparison.Ordinal);
                bool startLoss = originalInitMet && !variantInitMet;
                if (startLoss)
                {
                    decoyVariations.Add(new SequenceVariation(protein.BaseSequence.Length - sv.OneBasedEndPosition + 2, protein.BaseSequence.Length, new string(originalArray).Substring(0, originalArray.Length - 1), new string(variationArray), "DECOY VARIANT: " + sv.Description, decoyVariantModifications));
                }
                // both start with M, but there's more
                else if (sv.VariantSequence.StartsWith("M", StringComparison.Ordinal) && sv.OneBasedBeginPosition == 1 && (sv.OriginalSequence.Length > 1 || sv.VariantSequence.Length > 1))
                {
                    string original = new string(originalArray).Substring(0, originalArray.Length - 1);
                    string variant = new string(variationArray).Substring(0, variationArray.Length - 1);
                    decoyVariations.Add(new SequenceVariation(protein.BaseSequence.Length - sv.OneBasedEndPosition + 2, protein.BaseSequence.Length, original, variant, "DECOY VARIANT: " + sv.Description, decoyVariantModifications));
                }
                // gained an initiating methionine
                else if (sv.VariantSequence.StartsWith("M", StringComparison.Ordinal) && sv.OneBasedBeginPosition == 1)
                {
                    decoyVariations.Add(new SequenceVariation(1, 1, new string(originalArray), new string(variationArray), "DECOY VARIANT: " + sv.Description, decoyVariantModifications));
                }
                // starting methionine, but no variations on it
                else if (startsWithM)
                {
                    decoyVariations.Add(new SequenceVariation(protein.BaseSequence.Length - sv.OneBasedEndPosition + 2, protein.BaseSequence.Length - sv.OneBasedBeginPosition + 2, new string(originalArray), new string(variationArray), "DECOY VARIANT: " + sv.Description, decoyVariantModifications));
                }
                // no starting methionine
                else
                {
                    decoyVariations.Add(new SequenceVariation(protein.BaseSequence.Length - sv.OneBasedEndPosition + 1, protein.BaseSequence.Length - sv.OneBasedBeginPosition + 1, new string(originalArray), new string(variationArray), "DECOY VARIANT: " + sv.Description, decoyVariantModifications));
                }
            }
            return decoyVariations;
        }

        /// <summary>
        /// Generates a "slided" decoy sequence
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        private static List<Protein> GenerateSlideDecoys(List<Protein> proteins, int maxThreads = -1)
        {
            List<Protein> decoyProteins = new List<Protein>();
            Parallel.ForEach(proteins, new ParallelOptions { MaxDegreeOfParallelism = maxThreads }, protein =>
            {
                Dictionary<int, List<Modification>> decoyModifications = null;
                int numSlides = 20;
                char[] sequenceArrayUnslided = protein.BaseSequence.ToCharArray();
                char[] sequenceArraySlided = protein.BaseSequence.ToCharArray();
                decoyModifications = null;
                List<DisulfideBond> decoy_disulfides_slide = new List<DisulfideBond>();
                List<SpliceSite> spliceSitesSlide = new List<SpliceSite>();
                if (protein.BaseSequence.StartsWith("M", StringComparison.Ordinal))
                {
                    // Do not include the initiator methionine in shuffle!!!
                    if (numSlides % sequenceArraySlided.Length - 1 == 0)
                    {
                        numSlides++;
                    }
                    for (int i = 1; i < sequenceArraySlided.Length; i++)
                    {
                        sequenceArraySlided[i] = sequenceArrayUnslided[GetOldShuffleIndex(i, numSlides, protein.BaseSequence.Length, true)];
                    }

                    decoyModifications = new Dictionary<int, List<Modification>>(protein.OneBasedPossibleLocalizedModifications.Count);
                    foreach (var kvp in protein.OneBasedPossibleLocalizedModifications)
                    {
                        if (kvp.Key > 1)
                        {
                            decoyModifications.Add(GetOldShuffleIndex(kvp.Key - 1, numSlides, protein.BaseSequence.Length, true) + 1, kvp.Value);
                        }
                        else if (kvp.Key == 1)
                        {
                            decoyModifications.Add(1, kvp.Value);
                        }
                    }
                }
                else
                {
                    if (numSlides % sequenceArraySlided.Length == 0)
                    {
                        numSlides++;
                    }
                    for (int i = 0; i < sequenceArraySlided.Length; i++)
                    {
                        sequenceArraySlided[i] = sequenceArrayUnslided[GetOldShuffleIndex(i, numSlides, protein.BaseSequence.Length, false)];
                    }
                    decoyModifications = new Dictionary<int, List<Modification>>(protein.OneBasedPossibleLocalizedModifications.Count);
                    foreach (var kvp in protein.OneBasedPossibleLocalizedModifications)
                    {
                        decoyModifications.Add(GetOldShuffleIndex(kvp.Key - 1, numSlides, protein.BaseSequence.Length, false) + 1, kvp.Value);
                    }
                }
                var slided_sequence = new string(sequenceArraySlided);

                List<ProteolysisProduct> decoyPPSlide = new List<ProteolysisProduct>();
                foreach (ProteolysisProduct pp in protein.ProteolysisProducts)  //can't keep all aa like you can with reverse, just keep it the same length
                {
                    decoyPPSlide.Add(pp);
                }
                foreach (DisulfideBond disulfideBond in protein.DisulfideBonds) //these actually need the same cysteines...
                {
                    decoy_disulfides_slide.Add(new DisulfideBond(GetOldShuffleIndex(disulfideBond.OneBasedBeginPosition - 1, numSlides, slided_sequence.Length, false) + 1, GetOldShuffleIndex(disulfideBond.OneBasedEndPosition - 1, numSlides, slided_sequence.Length, false) + 1, "DECOY DISULFIDE BOND: " + disulfideBond.Description));
                }
                foreach (SpliceSite spliceSite in protein.SpliceSites) //these actually need the same cysteines...
                {
                    spliceSitesSlide.Add(new SpliceSite(GetOldShuffleIndex(spliceSite.OneBasedBeginPosition - 1, numSlides, slided_sequence.Length, false) + 1, GetOldShuffleIndex(spliceSite.OneBasedEndPosition - 1, numSlides, slided_sequence.Length, false) + 1, "DECOY SPLICE SITE: " + spliceSite.Description));
                }
                List<SequenceVariation> decoyVariationsSlide = new List<SequenceVariation>();
                foreach (SequenceVariation sv in protein.SequenceVariations) //No idea what's going on here. Review is appreciated.
                {
                    char[] originalArrayUnshuffled = sv.OriginalSequence.ToArray();
                    char[] variationArrayUnslided = sv.VariantSequence.ToArray();
                    if (sv.OneBasedBeginPosition == 1)
                    {
                        bool origInitM = sv.OriginalSequence.StartsWith("M", StringComparison.Ordinal);
                        bool varInitM = sv.VariantSequence.StartsWith("M", StringComparison.Ordinal);
                        if (origInitM && !varInitM)
                        {
                            decoyVariationsSlide.Add(new SequenceVariation(1, "M", "", "DECOY VARIANT: Initiator Methionine Change in " + sv.Description));
                        }
                        originalArrayUnshuffled = sv.OriginalSequence.Substring(Convert.ToInt32(origInitM)).ToArray();
                        variationArrayUnslided = sv.VariantSequence.Substring(Convert.ToInt32(varInitM)).ToArray();
                    }
                    int decoy_end = protein.BaseSequence.Length - sv.OneBasedBeginPosition + 2 +
                        Convert.ToInt32(sv.OneBasedEndPosition == slided_sequence.Length) - Convert.ToInt32(sv.OneBasedBeginPosition == 1);
                    int decoy_begin = decoy_end - originalArrayUnshuffled.Length + 1;
                    char[] originalArraySlided = sv.OriginalSequence.ToArray();
                    char[] variationArraySlided = sv.VariantSequence.ToArray();

                    if (numSlides % originalArraySlided.Length == 0)
                    {
                        numSlides++;
                    }
                    for (int i = 0; i < originalArraySlided.Length; i++)
                    {
                        originalArraySlided[i] = originalArrayUnshuffled[GetOldShuffleIndex(i, numSlides, originalArrayUnshuffled.Length, false)];
                    }

                    if (numSlides % variationArraySlided.Length == 0)
                    {
                        numSlides++;
                    }
                    for (int i = 0; i < variationArraySlided.Length; i++)
                    {
                        variationArraySlided[i] = variationArrayUnslided[GetOldShuffleIndex(i, numSlides, variationArrayUnslided.Length, false)];
                    }

                    decoyVariationsSlide.Add(new SequenceVariation(decoy_begin, decoy_end, new string(originalArraySlided), new string(variationArraySlided), "DECOY VARIANT: " + sv.Description));
                }
                var decoyProteinSlide = new Protein(slided_sequence, "DECOY_" + protein.Accession, protein.Organism, protein.GeneNames.ToList(), decoyModifications, decoyPPSlide,
                    protein.Name, protein.FullName, true, protein.IsContaminant, null, decoyVariationsSlide, null, protein.SampleNameForVariants, decoy_disulfides_slide, spliceSitesSlide, protein.DatabaseFilePath);
                lock (decoyProteins) { decoyProteins.Add(decoyProteinSlide); }
            });
            return decoyProteins;
        }

        /// <summary>
        /// Not sure...
        /// </summary>
        /// <param name="i"></param>
        /// <param name="numSlides"></param>
        /// <param name="sequenceLength"></param>
        /// <param name="methioninePresent"></param>
        /// <returns></returns>
        private static int GetOldShuffleIndex(int i, int numSlides, int sequenceLength, bool methioninePresent)
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
    }
}