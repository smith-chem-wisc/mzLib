using Proteomics.PSM;
using System.Diagnostics;
using TorchSharp;
using TorchSharp.Modules;

namespace MachineLearning.RetentionTimePredictionModels
{
    /// <summary>
    /// Chronologer is a deep learning model for highly accurate prediction of peptide C18 retention times (reported in % ACN).
    /// Chronologer was trained on a new large harmonized database of > 2.6 million retention time observations
    /// (2.25M unique peptides) constructed from 11 community datasets
    /// and natively supports prediction of 17 different modification types.
    /// With only a few observations of a new modification type (> 10 peptides),
    /// Chronologer can be easily re-trained to predict up to 10 user supplied modifications.
    ///
    /// Damien Beau Wilburn, Ariana E. Shannon, Vic Spicer, Alicia L. Richards, Darien Yeung, Danielle L. Swaney, Oleg V. Krokhin, Brian C. Searle
    /// bioRxiv 2023.05.30.542978; doi: https://doi.org/10.1101/2023.05.30.542978
    /// 
    /// https://github.com/searlelab/chronologer
    ///
    /// Licensed under Apache License 2.0
    /// 
    /// </summary>
    public class Chronologer : DeepTorch
    {
        public DataLoader? TrainDataLoader { get; set; }
        public DataLoader? TestDataLoader { get; set; }


        public Chronologer(Verbosity trainingVerbosity, string preTrainedModelPath = null,
            bool evalMode = true) : base(trainingVerbosity, preTrainedModelPath, evalMode) { }

        public override void Train()
        {
            throw new NotImplementedException();
        }

        public override torch.Tensor forward(torch.Tensor input)
        {
            throw new NotImplementedException();
        }

        protected override void EarlyStop()
        {
            throw new NotImplementedException();
        }

        protected override void RegisterComponents()
        {
            throw new NotImplementedException();
        }

        public static torch.Tensor Tensorize(PsmFromTsv psm, Dictionary<(char, string), int> dictionary,
            double qValueFiltering = 0.01, string ambiguityLevel = "1")
        {
            var fullSequence = psm.FullSequence.Split(new[] { '[', ']' })
                 .Where(x => !x.Equals("")).ToArray();

            //section to remoce type of mod from the sequence
            //i.e. Common Fixed:Carbamidomethyl and only stay with the target aa after the mod
            for (int i = 0; i < fullSequence.Length; i++)
            {
                if (fullSequence[i].Contains(" "))
                {
                    var tempSeq = fullSequence[i].Split(':');
                    fullSequence[i] = tempSeq[1];
                }
            }

            if (psm.BaseSeq.Length <= 50 &&
                psm.QValue <= qValueFiltering &&
                psm.AmbiguityLevel == 1.ToString()) //Chronologer only takes sequences of length 50 or less
            {
                var tensor = torch.zeros(1, 52, torch.ScalarType.Int64);

                tensor[0][0] = 38; //C-terminus
                var tensorCounter = 1; //skips the first element which is the C-terminus in the tensor
                char modID = ' '; //takes the target aa from inside the loop to hold it for the next loop
                bool mod = false; //assumes that the next loop is not a mod 

                if (fullSequence[0].Contains(" "))
                {
                    mod = true;
                    modID = fullSequence[0][fullSequence[0].Length - 1];
                }

                foreach (var subString in fullSequence)
                {
                    //if mod, enter
                    if (mod)
                    {
                        var key = (char.ToUpper(modID), subString);
                        if (dictionary.ContainsKey(key))
                        {
                            tensor[0][tensorCounter] = dictionary[key];
                            tensorCounter = tensorCounter + 1;
                        }
                        mod = false; //next iteration is not a mod
                        continue;
                    }

                    if (!subString.Contains(" "))
                    {
                        //without mods
                        for (int i = 0; i < subString.Length - 1; i++)
                        {
                            tensor[0][tensorCounter] = dictionary[(char.ToUpper(subString[i]), "")];
                            tensorCounter = tensorCounter + 1;
                        }
                        //next loop will be a mod
                        mod = true;

                        //save target aa for next loop
                        modID = subString[subString.Length - 1];
                    }
                    else
                    {
                        for (int i = 0; i < subString.Length; i++)
                        {
                            tensor[0][tensorCounter] = dictionary[(char.ToUpper(subString[i]), "")];
                            tensorCounter = tensorCounter + 1;
                        }
                        //next loop will be a mod
                        mod = true;

                    }
                }

                tensor[0][tensorCounter] = 44; //N-terminus

                Debug.WriteLine(tensor.ToString(TensorStringStyle.Julia));
                return tensor;
            }

            //if sequence is longer than 50, return a tensor of ones, quick way to remove it later from the dataset
            return torch.ones(1, 52, torch.ScalarType.Int64);
        }
    }
}
