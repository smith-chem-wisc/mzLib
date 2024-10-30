using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using Tensorboard;
using TorchSharp;

namespace Proteomics.RetentionTimePrediction.Chronologer
{
    public static class ChronologerEstimator
    {
        private static ChronologerModel.Chronologer ChronologerModel = new ChronologerModel.Chronologer();

        /// <summary>
        /// Uses the Chronologer model to predict C18 retention times (reported in % ACN).
        /// Only modifications present in the Chronologer dictionary are supported.
        /// Returns null if the sequence is not valid.
        /// <code>
        /// "Carbamidomethyl on C"
        /// "Oxidation on M"
        /// "Glu to PyroGlu"
        /// "Phosphorylation on S"
        /// "Phosphorylation on T"
        /// "Phosphorylation on Y"
        /// "Acetylation on K"
        /// "Succinylation on K"
        /// "Ubiquitination on K"
        /// "Methylation on K"
        /// "Dimethylation on K"
        /// "Trimethylation on K"
        /// "Methylation on R"
        /// "Dimethylation on R"
        /// </code>
        /// </summary>
        /// <param name="baseSequence"></param>
        /// <param name="fullSequence"></param>
        /// <returns></returns>
        public static float PredictRetentionTime(string baseSequence, string fullSequence)
        {
            if (baseSequence.Length > 50)
            {
                return -1;
            }
            var tensor = Tensorize(baseSequence, fullSequence);

            TorchSharp.Utils.TensorAccessor<long> tensorAccessor = tensor.data<long>();
            long[] tensorAsArray = tensorAccessor.ToArray();
            long[] incompatibleTensor = new long[52];
            incompatibleTensor.Fill(0);

            if (tensorAsArray.SequenceEqual(incompatibleTensor))
            {
                return -1;
            }

            torch.Tensor prediction = ChronologerModel.Predict(tensor);
            return prediction.data<float>()[0];
        }

        public static float[] PredictRetentionTime(string[] baseSequences, string[] fullSequences, bool gpu = true, int batch_size = 512)
        {
            torch.Device device;
            try
            {
                torch.InitializeDeviceType(DeviceType.CUDA);
                device = new torch.Device(DeviceType.CUDA);
            }
            catch
            {
                Console.WriteLine(new Exception("This device does not support CUDA, using CPU"));
                device = new torch.Device(DeviceType.CPU);
            }

            ChronologerModel.to(device);

            //vstack tensors and then batch predict to avoid running out of vRAM, I think system RAM transfer has big latency so keep em small

            float[] predictions = new float[baseSequences.Length];

            // tensorize all
            torch.Tensor[] tensorsArray = new torch.Tensor[baseSequences.Length];
            bool[] compatibleTracker = new bool[baseSequences.Length];

            Parallel.For(0, baseSequences.Length, (i, state) =>
            {
                var baseSeq = baseSequences[i];
                var fullSeq = fullSequences[i];

                var (tensor, compatible) = Tensorize(baseSeq, fullSeq, out bool chronologerCompatible);
                if (compatible)
                {
                    tensorsArray[i] = tensor;
                    compatibleTracker[i] = compatible;
                }
                else
                {
                    tensorsArray[i] = torch.zeros(1, 52, torch.ScalarType.Int64);
                    compatibleTracker[i] = compatible;
                }
            });
            // vstack and split
            var stackedTensors = torch.vstack(tensorsArray);

            float[] preds = new float[stackedTensors.size(0)];

            // make batches
            using var dataset = torch.utils.data.TensorDataset(stackedTensors);
            using var dataLoader = torch.utils.data.DataLoader(dataset, 2048);

            var predictionHolder = new List<float>();
            foreach (var batch in dataLoader)
            {
                var output = ChronologerModel.Predict(torch.vstack(batch).to(device));
                predictionHolder.AddRange(output.to(DeviceType.CPU).data<float>());
                output.Dispose();
            }

            Parallel.For(0, preds.Length, outputIndex =>
            {
                preds[outputIndex] = predictionHolder[outputIndex];
            });

            //change to -1 if same index in compatibleTracker is false, else leave as is
            //predictions = preds.SelectMany(x => x).ToArray();
            // return vstacked tensors as a matrix => float?[]

            for (var predictionsIndex = 0; predictionsIndex < predictions.Length; predictionsIndex++)
            {
                if (compatibleTracker[predictionsIndex])
                {
                    predictions[predictionsIndex] = predictionHolder[predictionsIndex];
                    continue;
                }

                predictions[predictionsIndex] = -1;
            }


            return predictions;
        }
        /// <summary>
        /// Takes the base sequence and the full peptide sequence and returns a tensor for the Chronologer model.
        /// The base sequence is the sequence without modifications and the full peptide sequence is the sequence with modifications.
        /// The model is intended to be used with sequences of length 50 or less, and only supports modifications present in the Chronologer dictionary.
        /// Sequences not supported by chronologer will return null.
        /// </summary>
        /// <param name="baseSequence"></param>
        /// <param name="fullSequence"></param>
        /// <returns></returns>
        private static (torch.Tensor, bool) Tensorize(string baseSequence, string fullSequence, out bool chronologerCompatible)
        {
            // Chronologer does not support metals
            if (fullSequence.Contains("Metal") || baseSequence.Contains("U") || baseSequence.Length > 50)
            {
                chronologerCompatible = false;
                return (torch.zeros(1, 52, torch.ScalarType.Int64), chronologerCompatible);
            }

            var fullSeq = fullSequence.Split(new[] { '[', ']' })
                .Where(x => !x.Equals("")).ToArray();


            //section to remove type of mod from the sequence
            //e.g. Common Fixed:Carbamidomethyl and only stay with the target aa after the mod
            for (int i = 0; i < fullSeq.Length; i++)
            {
                if (fullSeq[i].Contains(" "))
                {
                    var tempSeq = fullSeq[i].Split(':');
                    fullSeq[i] = tempSeq[1];
                }
            }

            if (baseSequence.Length <= 50) //Chronologer only takes sequences of length 50 or less
            {
                var tensor = torch.zeros(1, 52, torch.ScalarType.Int64);

                var tensorCounter = 1; //skips the first element which is the C-terminus in the tensor
                char modID = ' '; //takes the target aa from inside the loop to hold it for the next loop

                bool nTerminalMod = fullSequence[0] == '['; // true if the first loop is a mod

                foreach (var subString in fullSeq)
                {
                    //if mod, enter
                    if (nTerminalMod)
                    {
                        if (subString.Contains("Acetyl"))
                            tensor[0][0] = 39;
                        else
                        {
                            tensor[0][0] = 38;
                        }
                        nTerminalMod = false; //next iteration is not a mod
                        continue;
                    }

                    if (!subString.Contains(" "))
                    {
                        //without mods
                        for (int i = 0; i < subString.Length; i++)
                        {
                            tensor[0][tensorCounter] = ChronologerDictionary[(subString[i], "")];
                            tensorCounter++;
                        }

                        //save target aa for next loop
                        modID = subString[subString.Length - 1];
                    }

                    nTerminalMod = true;
                }

                tensor[0][tensorCounter] = 44; //C-terminus

                chronologerCompatible = true;
                return (tensor, chronologerCompatible);
            }

            chronologerCompatible = false;
            return (torch.zeros(1, 52, torch.ScalarType.Int64), chronologerCompatible);

        }

        private static torch.Tensor Tensorize(string baseSequence, string fullSequence)
        {
            // Chronologer does not support metals
            if (fullSequence.Contains("Metal") || baseSequence.Contains("U"))
            {
                return torch.zeros(1, 52, torch.ScalarType.Int64);
            }

            var fullSeq = fullSequence.Split(new[] { '[', ']' })
                .Where(x => !x.Equals("")).ToArray();


            //section to remove type of mod from the sequence
            //e.g. Common Fixed:Carbamidomethyl and only stay with the target aa after the mod
            for (int i = 0; i < fullSeq.Length; i++)
            {
                if (fullSeq[i].Contains(" "))
                {
                    var tempSeq = fullSeq[i].Split(':');
                    fullSeq[i] = tempSeq[1];
                }
            }

            if (baseSequence.Length <= 50) //Chronologer only takes sequences of length 50 or less
            {
                var tensor = torch.zeros(1, 52, torch.ScalarType.Int64);

                var tensorCounter = 1; //skips the first element which is the C-terminus in the tensor
                char modID = ' '; //takes the target aa from inside the loop to hold it for the next loop

                bool nTerminalMod = fullSequence[0] == '['; // true if the first loop is a mod

                foreach (var subString in fullSeq)
                {
                    //if mod, enter
                    if (nTerminalMod)
                    {
                        if (subString.Contains("Acetyl"))
                            tensor[0][0] = 39;
                        else
                        {
                            tensor[0][0] = 38;
                        }
                        nTerminalMod = false; //next iteration is not a mod
                        continue;
                    }

                    if (!subString.Contains(" "))
                    {
                        //without mods
                        for (int i = 0; i < subString.Length; i++)
                        {
                            tensor[0][tensorCounter] = ChronologerDictionary[(subString[i], "")];
                            tensorCounter++;
                        }

                        //save target aa for next loop
                        modID = subString[subString.Length - 1];
                    }

                    nTerminalMod = true;
                }

                tensor[0][tensorCounter] = 44; //C-terminus

                return tensor;
            }

            return torch.zeros(1, 52, torch.ScalarType.Int64);

        }

        private static readonly Dictionary<(char, string), int> ChronologerDictionary = new()
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
                { ('E', "Glu to PyroGlu"), 24 }, //'Pyro-glutamate
                { ('S', "Phosphorylation on S"), 25 }, //'Phosphoserine
                { ('T', "Phosphorylation on T"), 26 }, //'Phosphothreonine
                { ('Y', "Phosphorylation on Y"), 27 }, //'Phosphotyrosine
                { ('K', "Accetylation on K"), 28 }, //acetylation
                { ('K', "Succinylation on K"), 29 }, //'Succinylated
                { ('K', "Ubiquitination on K"), 30 }, //'Ubiquitinated
                { ('K', "Methylation on K"), 31 }, //'Mono-methyl
                { ('K', "Dimethylation on K"), 32 }, //'Dimethyl
                { ('K', "Trimethylation on K"), 33 }, //'Trimethyl
                { ('R', "Methylation on R"), 34 }, //'Monomethyl
                { ('R', "Dimethylation on R"), 35 }, //'Dimethyl
                {('-', "Free N-terminal"), 38},
                {('^', "N-acetyl"), 39},
                {('(', "Pyro-carbamidomethyl on C"), 40}, //S-carbamidomethylcysteine
                {(')', "pyroglutamate"), 41},
                {('&', "N-terminal TMT0"), 42},
                {('&', "N-terminal TMT10"), 43},
                {('_', "Free C-terminal"), 44},
            };
    }
}