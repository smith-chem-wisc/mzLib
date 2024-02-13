using System;
using System.Collections.Generic;
using System.Linq;
using TorchSharp;

namespace Proteomics.RetentionTimePrediction.Chronologer;
public static class ChronologerRetentionTimeEstimator
{
    public static double PredictRetentionTime(string baseSequence, string fullPeptideSequence)
    {
        torch.Tensor tensor = Tensorize(baseSequence, fullPeptideSequence);
        var model = new Chronologer();
        var prediction = model.Predict(tensor);
        return prediction[0].ToDouble();
    }

    private static torch.Tensor Tensorize(string baseSequence, string fullPeptideSequence)
    {
        var fullSequence = fullPeptideSequence.Split(new[] { '[', ']' })
            .Where(x => !x.Equals("")).ToArray();


        //section to remoce type of mod from the sequence
        //e.g. Common Fixed:Carbamidomethyl and only stay with the target aa after the mod
        for (int i = 0; i < fullSequence.Length; i++)
        {
            if (fullSequence[i].Contains(" "))
            {
                var tempSeq = fullSequence[i].Split(':');
                fullSequence[i] = tempSeq[1];
            }
        }

        if (baseSequence.Length <= 50) //Chronologer only takes sequences of length 50 or less
        {
            var tensor = torch.zeros(1, 52, torch.ScalarType.Int64);

            tensor[0][0] = 38; //C-terminus
            var tensorCounter = 1; //skips the first element which is the C-terminus in the tensor
            char modID = ' '; //takes the target aa from inside the loop to hold it for the next loop

            bool mod = fullPeptideSequence[0] == '['; // true if the first loop is a mod

            foreach (var subString in fullSequence)
            {
                //if mod, enter
                if (mod)
                {
                    var key = (modID, subString[0].ToString());

                    mod = false; //next iteration is not a mod
                    continue;
                }

                if (!subString.Contains(" "))
                {
                    //without mods
                    for (int i = 0; i < subString.Length; i++)
                    {
                        tensor[0][tensorCounter] = ChronologerDictionary[(subString[i], "")];
                        tensorCounter = tensorCounter + 1;
                    }

                    //save target aa for next loop
                    modID = subString[subString.Length - 1];
                }

                mod = true;
            }

            tensor[0][tensorCounter] = 44; //N-terminus

            return tensor;
        }

        throw new Exception("Sequence is too long for Chronologer. Expected length less or equal to 50");

    }

    private static readonly Dictionary<(char, string), int> ChronologerDictionary =
        new Dictionary<(char, string), int>()
        {
                { ('A', ""), 1 }, //'Alanine
                { ('C', ""), 2 }, //'Cysteine
                { ('D', ""), 3 }, //'Aspartate
                { ('E', ""), 4 }, //'Glutamate
                { ('F', ""), 5 }, //'Phenylalanine
                { ('G', ""), 6 }, //'Glycine
                { ('H', ""), 7 }, //'Histidine
                { ('I', ""), 8 }, //'Isoleucine
                { ('K', ""), 9 }, //'Lysine
                { ('L', ""), 10 }, //'Leucine
                { ('M', ""), 11 }, //'Methionine
                { ('N', ""), 12 }, //'Asparagine
                { ('P', ""), 13 }, //'Proline
                { ('Q', ""), 14 }, //'Glutamine
                { ('R', ""), 15 }, //'Argenine
                { ('S', ""), 16 }, //'Serine
                { ('T', ""), 17 }, //'Threonine
                { ('V', ""), 18 }, //'Valine
                { ('W', ""), 19 }, //'Tryptophane
                { ('Y', ""), 20 }, //'Tyrosine
                { ('C', "Carbamidomethyl on C"), 21 }, //'Carbamidomethyl
                { ('M', "Oxidation on M"), 22 }, //'Oxidized
                { ('E', "Glu to PyroGlu"), 24 }, //'Pyroglutamate
                { ('S', "Phosphorylation on S"), 25 }, //'Phosphoserine
                { ('T', "Phosphorylation on T"), 26 }, //'Phosphothreonine
                { ('Y', "Phosphorylation on Y"), 27 }, //'Phosphotyrosine
                { ('K', "Accetylation on K"), 28 }, //'Acetylated
                { ('K', "Succinylation on K"), 29 }, //'Succinylated
                { ('K', "Ubiquitination on K"), 30 }, //'Ubiquitinated
                { ('K', "Methylation on K"), 31 }, //'Monomethyl
                { ('K', "Dimethylation on K"), 32 }, //'Dimethyl
                { ('K', "Trimethylation on K"), 33 }, //'Trimethyl
                { ('R', "Methylation on R"), 34 }, //'Monomethyl
                { ('R', "Dimethylation on R"), 35 }, //'Dimethyl
        };
}
