using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace UsefulProteomicsDatabases
{
    public static class DecoyProteinGenerator
    {
        /// <summary>
        /// Generates a reverse decoy sequence
        /// </summary>
        /// <param name="block"></param>
        /// <param name="isContaminant"></param>
        /// <param name="proteinDbLocation"></param>
        /// <returns></returns>
        public static Protein GenerateReverseDecoy(ProteinXmlEntry block, bool isContaminant, string proteinDbLocation)
        {
            Dictionary<int, List<Modification>> decoy_modifications = null;
            char[] sequence_array = block.Sequence.ToCharArray();
            List<DisulfideBond> decoy_disulfides = new List<DisulfideBond>();
            if (block.Sequence.StartsWith("M", StringComparison.Ordinal))
            {
                // Do not include the initiator methionine in reversal!!!
                Array.Reverse(sequence_array, 1, block.Sequence.Length - 1);
                decoy_modifications = new Dictionary<int, List<Modification>>(block.OneBasedModifications.Count);
                foreach (var kvp in block.OneBasedModifications)
                {
                    if (kvp.Key > 1)
                    {
                        decoy_modifications.Add(block.Sequence.Length - kvp.Key + 2, kvp.Value);
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
                decoy_modifications = new Dictionary<int, List<Modification>>(block.OneBasedModifications.Count);
                foreach (var kvp in block.OneBasedModifications)
                {
                    decoy_modifications.Add(block.Sequence.Length - kvp.Key + 1, kvp.Value);
                }
            }
            var reversed_sequence = new string(sequence_array);

            List<ProteolysisProduct> decoyPP = new List<ProteolysisProduct>();
            foreach (ProteolysisProduct pp in block.ProteolysisProducts)
            {
                decoyPP.Add(new ProteolysisProduct(block.Sequence.Length - pp.OneBasedEndPosition + 1, block.Sequence.Length - pp.OneBasedBeginPosition, pp.Type));
            }
            foreach (DisulfideBond disulfideBond in block.DisulfideBonds)
            {
                decoy_disulfides.Add(new DisulfideBond(block.Sequence.Length - disulfideBond.OneBasedBeginPosition + 2, block.Sequence.Length - disulfideBond.OneBasedEndPosition + 2, "DECOY DISULFIDE BOND: " + disulfideBond.Description));
            }

            List<SequenceVariation> decoy_variations = new List<SequenceVariation>();
            foreach (SequenceVariation sv in block.SequenceVariations)
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
                int decoy_end = block.Sequence.Length - sv.OneBasedBeginPosition + 2 + Convert.ToInt32(sv.OneBasedEndPosition == reversed_sequence.Length) - Convert.ToInt32(sv.OneBasedBeginPosition == 1);
                int decoy_begin = decoy_end - original_array.Length + 1;
                Array.Reverse(original_array);
                Array.Reverse(variation_array);
                decoy_variations.Add(new SequenceVariation(decoy_begin, decoy_end, new string(original_array), new string(variation_array), "DECOY VARIANT: " + sv.Description));
            }
            var decoy_protein = new Protein(reversed_sequence, "DECOY_" + block.Accession, block.Organism, block.GeneNames, decoy_modifications, decoyPP,
                block.Name, block.FullName, true, isContaminant, null, decoy_variations, decoy_disulfides, proteinDbLocation);
            return decoy_protein;
        }

        /// <summary>
        /// Generates a "slided" decoy sequence
        /// </summary>
        /// <param name="block"></param>
        /// <param name="isContaminant"></param>
        /// <param name="proteinDbLocation"></param>
        /// <returns></returns>
        public static Protein GenerateSlideDecoy(ProteinXmlEntry block, bool isContaminant, string proteinDbLocation)
        {
            Dictionary<int, List<Modification>> decoy_modifications = null;
            int numSlides = 20;
            char[] sequence_array_unslided = block.Sequence.ToCharArray();
            char[] sequence_array_slided = block.Sequence.ToCharArray();
            decoy_modifications = null;
            List<DisulfideBond> decoy_disulfides_slide = new List<DisulfideBond>();
            if (block.Sequence.StartsWith("M", StringComparison.Ordinal))
            {
                // Do not include the initiator methionine in shuffle!!!
                if (numSlides % sequence_array_slided.Length - 1 == 0)
                {
                    numSlides++;
                }
                for (int i = 1; i < sequence_array_slided.Length; i++)
                {
                    sequence_array_slided[i] = sequence_array_unslided[GetOldShuffleIndex(i, numSlides, block.Sequence.Length, true)];
                }

                decoy_modifications = new Dictionary<int, List<Modification>>(block.OneBasedModifications.Count);
                foreach (var kvp in block.OneBasedModifications)
                {
                    if (kvp.Key > 1)
                    {
                        decoy_modifications.Add(GetOldShuffleIndex(kvp.Key - 1, numSlides, block.Sequence.Length, true) + 1, kvp.Value);
                    }
                    else if (kvp.Key == 1)
                    {
                        decoy_modifications.Add(1, kvp.Value);
                    }
                }
            }
            else
            {
                if (numSlides % sequence_array_slided.Length == 0)
                {
                    numSlides++;
                }
                for (int i = 0; i < sequence_array_slided.Length; i++)
                {
                    sequence_array_slided[i] = sequence_array_unslided[GetOldShuffleIndex(i, numSlides, block.Sequence.Length, false)];
                }
                decoy_modifications = new Dictionary<int, List<Modification>>(block.OneBasedModifications.Count);
                foreach (var kvp in block.OneBasedModifications)
                {
                    decoy_modifications.Add(GetOldShuffleIndex(kvp.Key - 1, numSlides, block.Sequence.Length, false) + 1, kvp.Value);
                }
            }
            var slided_sequence = new string(sequence_array_slided);

            List<ProteolysisProduct> decoyPP_slide = new List<ProteolysisProduct>();
            foreach (ProteolysisProduct pp in block.ProteolysisProducts)  //can't keep all aa like you can with reverse, just keep it the same length
            {
                decoyPP_slide.Add(pp);
            }
            foreach (DisulfideBond disulfideBond in block.DisulfideBonds) //these actually need the same cysteines...
            {
                decoy_disulfides_slide.Add(new DisulfideBond(GetOldShuffleIndex(disulfideBond.OneBasedBeginPosition - 1, numSlides, slided_sequence.Length, false) + 1, GetOldShuffleIndex(disulfideBond.OneBasedEndPosition - 1, numSlides, slided_sequence.Length, false) + 1, "DECOY DISULFIDE BOND: " + disulfideBond.Description));
            }
            List<SequenceVariation> decoy_variations_slide = new List<SequenceVariation>();
            foreach (SequenceVariation sv in block.SequenceVariations) //No idea what's going on here. Review is appreciated.
            {
                char[] original_array_unshuffled = sv.OriginalSequence.ToArray();
                char[] variation_array_unslided = sv.VariantSequence.ToArray();
                if (sv.OneBasedBeginPosition == 1)
                {
                    bool orig_init_m = sv.OriginalSequence.StartsWith("M", StringComparison.Ordinal);
                    bool var_init_m = sv.VariantSequence.StartsWith("M", StringComparison.Ordinal);
                    if (orig_init_m && !var_init_m)
                    {
                        decoy_variations_slide.Add(new SequenceVariation(1, "M", "", "DECOY VARIANT: Initiator Methionine Change in " + sv.Description));
                    }
                    original_array_unshuffled = sv.OriginalSequence.Substring(Convert.ToInt32(orig_init_m)).ToArray();
                    variation_array_unslided = sv.VariantSequence.Substring(Convert.ToInt32(var_init_m)).ToArray();
                }
                int decoy_end = block.Sequence.Length - sv.OneBasedBeginPosition + 2 + Convert.ToInt32(sv.OneBasedEndPosition == slided_sequence.Length) - Convert.ToInt32(sv.OneBasedBeginPosition == 1);
                int decoy_begin = decoy_end - original_array_unshuffled.Length + 1;
                char[] original_array_slided = sv.OriginalSequence.ToArray();
                char[] variation_array_slided = sv.VariantSequence.ToArray();

                if (numSlides % original_array_slided.Length == 0)
                {
                    numSlides++;
                }
                for (int i = 0; i < original_array_slided.Length; i++)
                {
                    original_array_slided[i] = original_array_unshuffled[GetOldShuffleIndex(i, numSlides, original_array_unshuffled.Length, false)];
                }

                if (numSlides % variation_array_slided.Length == 0)
                {
                    numSlides++;
                }
                for (int i = 0; i < variation_array_slided.Length; i++)
                {
                    variation_array_slided[i] = variation_array_unslided[GetOldShuffleIndex(i, numSlides, variation_array_unslided.Length, false)];
                }

                decoy_variations_slide.Add(new SequenceVariation(decoy_begin, decoy_end, new string(original_array_slided), new string(variation_array_slided), "DECOY VARIANT: " + sv.Description));
            }
            var decoy_protein_slide = new Protein(slided_sequence, "DECOY_" + block.Accession, block.Organism, block.GeneNames, decoy_modifications, decoyPP_slide,
                block.Name, block.FullName, true, isContaminant, null, decoy_variations_slide, decoy_disulfides_slide, proteinDbLocation);
            return decoy_protein_slide;
        }

        /// <summary>
        /// Not sure...
        /// </summary>
        /// <param name="i"></param>
        /// <param name="numSlides"></param>
        /// <param name="sequenceLength"></param>
        /// <param name="methioninePresent"></param>
        /// <returns></returns>
        public static int GetOldShuffleIndex(int i, int numSlides, int sequenceLength, bool methioninePresent)
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