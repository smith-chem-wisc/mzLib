using System.Collections.Generic;
using System.Linq;
using Proteomics.Fragmentation;

namespace Proteomics.ProteolyticDigestion
{
    public class ProteinDigestion
    {
        /// <summary>
        /// Initializes digestion object
        /// </summary>
        /// <param name="digestionParams"></param>
        /// <param name="allKnownFixedModifications"></param>
        /// <param name="variableModifications"></param>
        public ProteinDigestion(DigestionParams digestionParams, IEnumerable<Modification> allKnownFixedModifications, List<Modification> variableModifications)
        {
            DigestionParams = digestionParams;
            Protease = digestionParams.Protease;
            MaximumMissedCleavages = digestionParams.MaxMissedCleavages;
            InitiatorMethionineBehavior = digestionParams.InitiatorMethionineBehavior;
            MinPeptideLength = digestionParams.MinPeptideLength;
            MaxPeptideLength = digestionParams.MaxPeptideLength;
            AllKnownFixedModifications = allKnownFixedModifications ?? new List<Modification>();
            VariableModifications = variableModifications ?? new List<Modification>();
        }

        public Protease Protease { get; set; }
        public int MaximumMissedCleavages { get; set; }
        public DigestionParams DigestionParams { get; set; }
        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }
        public int MinPeptideLength { get; set; }
        public int MaxPeptideLength { get; set; }
        public IEnumerable<Modification> AllKnownFixedModifications { get; set; }
        public List<Modification> VariableModifications { get; set; }

        /// <summary>
        /// Gets peptides for speedy semispecific digestion of a protein
        /// This generates specific peptides of maximum missed cleavages
        /// These peptides need to be digested post search to their actual sequences
        /// semi-specific search enters here...
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        public IEnumerable<PeptideWithSetModifications> SpeedySemiSpecificDigestion(Protein protein) //We are only getting fully specific peptides of the maximum cleaved residues here
        {
            List<ProteolyticPeptide> peptides = new List<ProteolyticPeptide>();
            List<int> oneBasedIndicesToCleaveAfter = Protease.GetDigestionSiteIndices(protein.BaseSequence); //get peptide bonds to cleave SPECIFICALLY (termini included)

            //it's possible not to go through this loop (maxMissedCleavages+1>number of indexes), and that's okay. It will get digested in the next loops (finish C/N termini)
            for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - MaximumMissedCleavages - 1; i++)
            {
                bool retain = Protease.Retain(i, InitiatorMethionineBehavior, protein[0]);
                if (retain) //it's okay to use i instead of oneBasedIndicesToCleaveAfter[i], because the index of zero is zero and it only checks if it's the N-terminus or not
                {
                    int peptideLength = oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1] - oneBasedIndicesToCleaveAfter[i];
                    if (peptideLength >= MinPeptideLength) //if bigger than min
                    {
                        if (peptideLength <= MaxPeptideLength) //if an acceptable length (bigger than min, smaller than max), add it
                        {
                            peptides.Add(new ProteolyticPeptide(protein, oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1],
                                MaximumMissedCleavages, CleavageSpecificity.Full, "full"));
                        }
                        else if (DigestionParams.FragmentationTerminus == FragmentationTerminus.N) //make something with the maximum length and fixed N
                        {
                            int tempIndex = oneBasedIndicesToCleaveAfter[i] + 1;
                            peptides.Add(new ProteolyticPeptide(protein, tempIndex, tempIndex + MaxPeptideLength, MaximumMissedCleavages, CleavageSpecificity.Semi, "semi"));
                        }
                        else //It has to be FragmentationTerminus.C //make something with the maximum length and fixed C
                        {
                            int tempIndex = oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1];
                            peptides.Add(new ProteolyticPeptide(protein, tempIndex - MaxPeptideLength, tempIndex, MaximumMissedCleavages, CleavageSpecificity.Semi, "semi"));
                        }
                    }
                }

                if (Protease.Cleave(i, InitiatorMethionineBehavior, protein[0]) && (DigestionParams.FragmentationTerminus == FragmentationTerminus.N || !retain)) //it's okay to use i instead of oneBasedIndicesToCleaveAfter[i], because the index of zero is zero and it only checks if it's the N-terminus or not
                {
                    int peptideLength = oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1] - 1;
                    if (peptideLength >= MinPeptideLength)
                    {
                        if (peptideLength <= MaxPeptideLength)
                        {
                            peptides.Add(new ProteolyticPeptide(protein, 2, oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1], //two is hardcoded, since M=1, so the next aa is 2 (one based)
                                MaximumMissedCleavages, CleavageSpecificity.Full, "full:M cleaved"));
                        }
                        else if (DigestionParams.FragmentationTerminus == FragmentationTerminus.N)
                        {
                            peptides.Add(new ProteolyticPeptide(protein, 2, 2 + MaxPeptideLength + MaxPeptideLength, MaximumMissedCleavages, CleavageSpecificity.Semi, "semi"));
                        }
                        else //It has to be FragmentationTerminus.C //make something with the maximum length and fixed C
                        {
                            //kinda tricky, because we'll be creating a duplication if cleavage is variable
                            if (!Protease.Retain(i, InitiatorMethionineBehavior, protein[0])) //only if cleave, because then not made earlier during retain
                            {
                                int tempIndex = oneBasedIndicesToCleaveAfter[i + MaximumMissedCleavages + 1];
                                peptides.Add(new ProteolyticPeptide(protein, tempIndex - MaxPeptideLength, tempIndex, MaximumMissedCleavages, CleavageSpecificity.Semi, "semi"));
                            }
                        }
                    }
                }
            }

            //wrap up the termini that weren't hit earlier
            int lastIndex = oneBasedIndicesToCleaveAfter.Count - 1; //last cleavage index (the c-terminus)
            int maxIndexDifference = MaximumMissedCleavages < lastIndex ? MaximumMissedCleavages : lastIndex; //the number of index differences allowed. 
            //If the protein has fewer cleavage sites than allowed missed cleavages, just use the number of cleavage sites (lastIndex)
            bool nTerminusFragmentation = DigestionParams.FragmentationTerminus == FragmentationTerminus.N;
            for (int i = 1; i <= maxIndexDifference; i++) //i is the difference (in indexes) between indexes (cleavages), so it needs to start at 1, or the peptide would have length = 0
            {
                int startIndex = nTerminusFragmentation ?
                    oneBasedIndicesToCleaveAfter[lastIndex - i] :
                    oneBasedIndicesToCleaveAfter[0];
                int endIndex = nTerminusFragmentation ?
                    oneBasedIndicesToCleaveAfter[lastIndex] :
                    oneBasedIndicesToCleaveAfter[i];

                int peptideLength = endIndex - startIndex;
                if(peptideLength>=MinPeptideLength)
                {
                    if(peptideLength<=MaxPeptideLength) //if okay length, add it up to the terminus
                    {
                        peptides.Add(new ProteolyticPeptide(protein, startIndex + 1, endIndex, i - 1, CleavageSpecificity.Full, "full"));
                    }
                    else //update so that not the end of terminus
                    {
                        if(nTerminusFragmentation)
                        {
                            endIndex = startIndex + MaxPeptideLength;
                        }
                        else
                        {
                            startIndex = endIndex - MaxPeptideLength;
                        }
                        peptides.Add(new ProteolyticPeptide(protein, startIndex, endIndex, i - 1, CleavageSpecificity.Semi, "semi"));
                    }
                }        
            }

            // Also digest using the proteolysis product start/end indices
            peptides.AddRange(
                protein.ProteolysisProducts
                    .Where(proteolysisProduct => proteolysisProduct.OneBasedEndPosition.HasValue && proteolysisProduct.OneBasedBeginPosition.HasValue &&
                    (proteolysisProduct.OneBasedBeginPosition != 1 || proteolysisProduct.OneBasedEndPosition != protein.Length))
                    .Select(proteolysisProduct => new ProteolyticPeptide(protein, proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.OneBasedEndPosition.Value,
                        0, CleavageSpecificity.Full, proteolysisProduct.Type + " start")));

            return peptides.SelectMany(peptide => peptide.GetModifiedPeptides(AllKnownFixedModifications, DigestionParams, VariableModifications));
        }

        /// <summary>
        /// Gets peptides for specific protease digestion of a protein
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        public IEnumerable<PeptideWithSetModifications> Digestion(Protein protein)
        {
            var unmodifiedPeptides = Protease.GetUnmodifiedPeptides(protein, MaximumMissedCleavages, InitiatorMethionineBehavior, MinPeptideLength, MaxPeptideLength);
            return unmodifiedPeptides.SelectMany(peptide => peptide.GetModifiedPeptides(AllKnownFixedModifications, DigestionParams, VariableModifications));
        }
    }
}