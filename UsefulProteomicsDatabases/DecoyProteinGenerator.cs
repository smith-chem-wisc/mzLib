using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace UsefulProteomicsDatabases
{
    public static class DecoyProteinGenerator
    {
        /// <summary>
        /// Generates decoys for a list of proteins
        /// </summary>
        /// <param name="proteins"></param>
        /// <param name="decoySetting"></param>
        /// <returns></returns>
        public static List<Protein> GenerateDecoys(List<Protein> proteins, DecoySetting decoySetting)
        {
            return decoySetting == null || decoySetting.DecoyType == DecoyType.None ?
                new List<Protein>() :
                proteins.Select(p => GenerateDecoy(p, decoySetting)).ToList();
        }
        
        /// <summary>
        /// Generates a decoy protein entry for a given protein
        /// </summary>
        /// <param name="proteins"></param>
        /// <param name="decoySetting"></param>
        /// <returns></returns>
        private static Protein GenerateDecoy(Protein protein, DecoySetting decoySetting)
        {
            if (decoySetting == null || decoySetting.DecoyType == DecoyType.None)
            {
                return null;
            }
            else if (decoySetting.DecoyType == DecoyType.Reverse)
            {
                return GenerateReverseDecoy(protein);
            }
            else if (decoySetting.DecoyType == DecoyType.Slide)
            {
                return GenerateSlideDecoy(protein);
            }
            else if (decoySetting.DecoyType == DecoyType.Shuffle)
            {
                return GenerateShuffledDecoy(protein);
            }
            else
            {
                throw new ArgumentException("Decoy type " + decoySetting.DecoyType.ToString() + " is not implemented.");
            }
        }

        /// <summary>
        /// Generates a reverse decoy sequence
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        private static Protein GenerateReverseDecoy(Protein protein)
        {
            Dictionary<int, List<Modification>> decoy_modifications = null;
            char[] sequence_array = protein.BaseSequence.ToCharArray();
            List<DisulfideBond> decoy_disulfides = new List<DisulfideBond>();
            if (protein.BaseSequence.StartsWith("M", StringComparison.Ordinal))
            {
                // Do not include the initiator methionine in reversal!!!
                Array.Reverse(sequence_array, 1, protein.BaseSequence.Length - 1);
                decoy_modifications = new Dictionary<int, List<Modification>>(protein.OneBasedPossibleLocalizedModifications.Count);
                foreach (var kvp in protein.OneBasedPossibleLocalizedModifications)
                {
                    if (kvp.Key > 1)
                    {
                        decoy_modifications.Add(protein.BaseSequence.Length - kvp.Key + 2, kvp.Value);
                    }
                    else if (kvp.Key == 1)
                    {
                        decoy_modifications.Add(1, kvp.Value);
                    }
                }
            }
            else
            {
                Array.Reverse(sequence_array);
                decoy_modifications = new Dictionary<int, List<Modification>>(protein.OneBasedPossibleLocalizedModifications.Count);
                foreach (var kvp in protein.OneBasedPossibleLocalizedModifications)
                {
                    decoy_modifications.Add(protein.BaseSequence.Length - kvp.Key + 1, kvp.Value);
                }
            }
            var reversed_sequence = new string(sequence_array);

            List<ProteolysisProduct> decoyPP = new List<ProteolysisProduct>();
            foreach (ProteolysisProduct pp in protein.ProteolysisProducts)
            {
                decoyPP.Add(new ProteolysisProduct(protein.BaseSequence.Length - pp.OneBasedEndPosition + 1, protein.BaseSequence.Length - pp.OneBasedBeginPosition, pp.Type));
            }
            foreach (DisulfideBond disulfideBond in protein.DisulfideBonds)
            {
                decoy_disulfides.Add(new DisulfideBond(protein.BaseSequence.Length - disulfideBond.OneBasedBeginPosition + 2, protein.BaseSequence.Length - disulfideBond.OneBasedEndPosition + 2, "DECOY DISULFIDE BOND: " + disulfideBond.Description));
            }

            List<SequenceVariation> decoy_variations = new List<SequenceVariation>();
            foreach (SequenceVariation sv in protein.SequenceVariations)
            {
                char[] original_array = sv.OriginalSequence.ToArray();
                char[] variation_array = sv.VariantSequence.ToArray();
                if (sv.OneBasedBeginPosition == 1)
                {
                    bool orig_init_m = sv.OriginalSequence.StartsWith("M", StringComparison.Ordinal);
                    bool var_init_m = sv.VariantSequence.StartsWith("M", StringComparison.Ordinal);
                    if (orig_init_m && !var_init_m)
                    {
                        decoy_variations.Add(new SequenceVariation(1, "M", "", "DECOY VARIANT: Initiator Methionine Change in " + sv.Description));
                    }
                    original_array = sv.OriginalSequence.Substring(Convert.ToInt32(orig_init_m)).ToArray();
                    variation_array = sv.VariantSequence.Substring(Convert.ToInt32(var_init_m)).ToArray();
                }
                int decoy_end = protein.BaseSequence.Length - sv.OneBasedBeginPosition + 2 + Convert.ToInt32(sv.OneBasedEndPosition == reversed_sequence.Length) - Convert.ToInt32(sv.OneBasedBeginPosition == 1);
                int decoy_begin = decoy_end - original_array.Length + 1;
                Array.Reverse(original_array);
                Array.Reverse(variation_array);
                decoy_variations.Add(new SequenceVariation(decoy_begin, decoy_end, new string(original_array), new string(variation_array), "DECOY VARIANT: " + sv.Description));
            }
            var decoy_protein = new Protein(reversed_sequence, "DECOY_" + protein.Accession, protein.Organism, protein.GeneNames.ToList(), decoy_modifications, decoyPP,
                protein.Name, protein.FullName, true, protein.IsContaminant, null, decoy_variations, decoy_disulfides, protein.DatabaseFilePath);
            return decoy_protein;
        }

        /// <summary>
        /// Generates a "slided" decoy sequence
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        private static Protein GenerateSlideDecoy(Protein protein)
        {
            Dictionary<int, List<Modification>> decoyModifications = null;
            int numSlides = 20;
            char[] sequenceArrayUnslided = protein.BaseSequence.ToCharArray();
            char[] sequenceArraySlided = protein.BaseSequence.ToCharArray();
            decoyModifications = null;
            List<DisulfideBond> decoy_disulfides_slide = new List<DisulfideBond>();
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
                protein.Name, protein.FullName, true, protein.IsContaminant, null, decoyVariationsSlide, decoy_disulfides_slide, protein.DatabaseFilePath);
            return decoyProteinSlide;
        }

        /// <summary>
        /// Generates a reverse decoy sequence
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        private static Protein GenerateShuffledDecoy(Protein protein)
        {
            return null;
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