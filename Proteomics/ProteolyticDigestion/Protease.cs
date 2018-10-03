using System;
using System.Collections.Generic;
using System.Linq;
using Proteomics.Fragmentation;

namespace Proteomics.ProteolyticDigestion
{
    public class Protease
    {
        public Protease(string name, IEnumerable<Tuple<string, FragmentationTerminus>> sequencesInducingCleavage, IEnumerable<Tuple<string, FragmentationTerminus>> sequencesPreventingCleavage, CleavageSpecificity cleavageSpecificity, string psiMSAccessionNumber, string psiMSName, string siteRegexp)
        {
            Name = name;
            SequencesInducingCleavage = sequencesInducingCleavage;
            SequencesPreventingCleavage = sequencesPreventingCleavage;
            CleavageSpecificity = cleavageSpecificity;
            PsiMsAccessionNumber = psiMSAccessionNumber;
            PsiMsName = psiMSName;
            SiteRegexp = siteRegexp;
        }

        public string Name { get; }
        public FragmentationTerminus CleavageTerminus { get; }
        public IEnumerable<Tuple<string, FragmentationTerminus>> SequencesInducingCleavage { get; }
        public IEnumerable<Tuple<string, FragmentationTerminus>> SequencesPreventingCleavage { get; }
        public CleavageSpecificity CleavageSpecificity { get; }
        public string PsiMsAccessionNumber { get; }
        public string PsiMsName { get; }
        public string SiteRegexp { get; }

        public override string ToString()
        {
            return Name;
        }

        public override bool Equals(object obj)
        {
            var a = obj as Protease;
            return a != null
                && a.Name.Equals(Name);
        }

        public override int GetHashCode()
        {
            return Name.GetHashCode();
        }

        public CleavageSpecificity GetCleavageSpecificity(string proteinSequence, int startIndex, int endIndex)
        {
            List<int> indicesToCleave = GetDigestionSiteIndices(proteinSequence);
            int cleavableMatches = 0;
            //if the start index is a cleavable index (-1 because one based) OR if the start index is after a cleavable methionine
            if (indicesToCleave.Contains(startIndex - 1) || (startIndex == 2 && proteinSequence[0] == 'M'))
            {
                cleavableMatches++;
            }
            if (indicesToCleave.Contains(endIndex)) //if the end index is a cleavable index
            {
                cleavableMatches++;
            }

            if (cleavableMatches == 0  //if neither were cleavable, (or it's singleN/C) then it's nonspecific
                || CleavageSpecificity == CleavageSpecificity.SingleN
                || CleavageSpecificity == CleavageSpecificity.SingleC)
            {
                return CleavageSpecificity.None;
            }
            else if (cleavableMatches == 1) //if one index was cleavable, then it's semi specific
            {
                return CleavageSpecificity.Semi;
            }
            else //2 if both, then it's fully speific
            {
                return CleavageSpecificity.Full;
            }
        }

        /// <summary>
        /// Gets intervals of a protein sequence that will result from digestion by this protease.
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="maximumMissedCleavages"></param>
        /// <param name="initiatorMethionineBehavior"></param>
        /// <param name="minPeptidesLength"></param>
        /// <param name="maxPeptidesLength"></param>
        /// <returns></returns>
        internal List<ProteolyticPeptide> GetUnmodifiedPeptides(Protein protein, int maximumMissedCleavages, InitiatorMethionineBehavior initiatorMethionineBehavior,
            int minPeptidesLength, int maxPeptidesLength)
        {
            List<ProteolyticPeptide> peptides = new List<ProteolyticPeptide>();

            // proteolytic cleavage in one spot
            if (CleavageSpecificity == CleavageSpecificity.SingleN)
            {
                bool maxTooBig = protein.Length + maxPeptidesLength < 0; //when maxPeptidesLength is too large, it becomes negative and causes issues
                //This happens when maxPeptidesLength == int.MaxValue or something close to it
                int startIndex = initiatorMethionineBehavior == InitiatorMethionineBehavior.Cleave ? 2 : 1;
                for (int proteinStart = startIndex; proteinStart <= protein.Length; proteinStart++)
                {
                    if (OkayMinLength(protein.Length - proteinStart + 1, minPeptidesLength))
                    {
                        //need Math.Max if max length is int.MaxLength, since +proteinStart will make it negative
                        //if the max length is too big to be an int (ie infinity), just do the protein length.
                        //if it's not too big to be an int, it might still be too big. Take the minimum of the protein length or the maximum length (-1, because the index is inclusive. Without -1, peptides will be one AA too long)
                        peptides.Add(new ProteolyticPeptide(protein, proteinStart, maxTooBig ? protein.Length : Math.Min(protein.Length, proteinStart + maxPeptidesLength - 1), 0, CleavageSpecificity.SingleN, "SingleN"));
                    }
                }
            }
            else if (CleavageSpecificity == CleavageSpecificity.SingleC)
            {
                int startIndex = initiatorMethionineBehavior == InitiatorMethionineBehavior.Cleave ? 2 : 1; //where does the protein start
                int lengthDifference = startIndex - 1; //take it back one for zero based index
                for (int proteinEnd = 1; proteinEnd <= protein.Length; proteinEnd++)
                {
                    //length of peptide will be at least the start index
                    if (OkayMinLength(proteinEnd - lengthDifference, minPeptidesLength)) //is the maximum possible length longer than the minimum?
                    {
                        //use the start index as the max of the N-terminus or the c-terminus minus the max (+1 because inclusive, otherwise peptides will be one AA too long)
                        peptides.Add(new ProteolyticPeptide(protein, Math.Max(startIndex, proteinEnd - maxPeptidesLength + 1), proteinEnd, 0, CleavageSpecificity.SingleC, "SingleC"));
                    }
                }
            }
            //top-down
            else if (CleavageSpecificity == CleavageSpecificity.None)
            {
                // retain methionine
                if ((initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || protein[0] != 'M')
                    && OkayLength(protein.Length, minPeptidesLength, maxPeptidesLength))
                {
                    peptides.Add(new ProteolyticPeptide(protein, 1, protein.Length, 0, CleavageSpecificity.Full, "full"));
                }

                // cleave methionine
                if ((initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && protein[0] == 'M')
                    && OkayLength(protein.Length - 1, minPeptidesLength, maxPeptidesLength))
                {
                    peptides.Add(new ProteolyticPeptide(protein, 2, protein.Length, 0, CleavageSpecificity.Full, "full:M cleaved"));
                }

                // Also digest using the proteolysis product start/end indices
                peptides.AddRange(
                    protein.ProteolysisProducts
                    .Where(proteolysisProduct => proteolysisProduct.OneBasedEndPosition.HasValue && proteolysisProduct.OneBasedBeginPosition.HasValue
                        && OkayLength(proteolysisProduct.OneBasedEndPosition.Value - proteolysisProduct.OneBasedBeginPosition.Value + 1, minPeptidesLength, maxPeptidesLength))
                    .Select(proteolysisProduct =>
                        new ProteolyticPeptide(protein, proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.OneBasedEndPosition.Value, 0, CleavageSpecificity.Full, proteolysisProduct.Type)));
            }

            // Full proteolytic cleavage
            else if (CleavageSpecificity == CleavageSpecificity.Full)
            {
                peptides.AddRange(FullDigestion(protein, initiatorMethionineBehavior, maximumMissedCleavages, minPeptidesLength, maxPeptidesLength));
            }

            // Cleavage rules for semi-specific search
            else if (CleavageSpecificity == CleavageSpecificity.Semi)
            {
                peptides.AddRange(SemiProteolyticDigestion(protein, initiatorMethionineBehavior, maximumMissedCleavages, minPeptidesLength, maxPeptidesLength));
            }
            else
            {
                throw new NotImplementedException();
            }

            return peptides;
        }

        /// <summary>
        /// Gets the indices after which this protease will cleave a given protein sequence
        /// </summary>
        /// <param name="proteinSequence"></param>
        /// <returns></returns>
        internal List<int> GetDigestionSiteIndices(string proteinSequence)
        {
            var indices = new List<int>();
            for (int i = 0; i < proteinSequence.Length - 1; i++)
            {
                foreach (var c in SequencesInducingCleavage)
                {
                    if (SequenceInducesCleavage(proteinSequence, i, c)
                        && !SequencesPreventingCleavage.Any(nc => SequencePreventsCleavage(proteinSequence, i, nc)))
                    {
                        indices.Add(i + 1);
                        break;
                    }
                }
            }
            indices.Insert(0, 0); // The start of the protein is treated as a cleavage site to retain the n-terminal peptide
            indices.Add(proteinSequence.Length); // The end of the protein is treated as a cleavage site to retain the c-terminal peptide
            return indices;
        }

        /// <summary>
        /// Retain N-terminal residue?
        /// </summary>
        /// <param name="oneBasedCleaveAfter"></param>
        /// <param name="initiatorMethionineBehavior"></param>
        /// <param name="nTerminus"></param>
        /// <returns></returns>
        internal static bool Retain(int oneBasedCleaveAfter, InitiatorMethionineBehavior initiatorMethionineBehavior, char nTerminus)
        {
            return oneBasedCleaveAfter != 0 // this only pertains to the n-terminus
                || initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave
                || nTerminus != 'M';
        }

        /// <summary>
        /// Cleave N-terminal residue?
        /// </summary>
        /// <param name="oneBasedCleaveAfter"></param>
        /// <param name="initiatorMethionineBehavior"></param>
        /// <param name="nTerminus"></param>
        /// <returns></returns>
        internal static bool Cleave(int oneBasedCleaveAfter, InitiatorMethionineBehavior initiatorMethionineBehavior, char nTerminus)
        {
            return oneBasedCleaveAfter == 0 // this only pertains to the n-terminus
                && initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain
                && nTerminus == 'M';
        }

        /// <summary>
        /// Is length of given peptide okay, given minimum and maximum?
        /// </summary>
        /// <param name="peptideLength"></param>
        /// <param name="minPeptidesLength"></param>
        /// <param name="maxPeptidesLength"></param>
        /// <returns></returns>
        internal static bool OkayLength(int? peptideLength, int? minPeptidesLength, int? maxPeptidesLength)
        {
            return OkayMinLength(peptideLength, minPeptidesLength) && OkayMaxLength(peptideLength, maxPeptidesLength);
        }

        /// <summary>
        /// Gets protein intervals for digestion by this specific protease.
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="initiatorMethionineBehavior"></param>
        /// <param name="maximumMissedCleavages"></param>
        /// <param name="minPeptidesLength"></param>
        /// <param name="maxPeptidesLength"></param>
        /// <returns></returns>
        private IEnumerable<ProteolyticPeptide> FullDigestion(Protein protein, InitiatorMethionineBehavior initiatorMethionineBehavior,
            int maximumMissedCleavages, int minPeptidesLength, int maxPeptidesLength)
        {
            List<int> oneBasedIndicesToCleaveAfter = GetDigestionSiteIndices(protein.BaseSequence);
            for (int missedCleavages = 0; missedCleavages <= maximumMissedCleavages; missedCleavages++)
            {
                for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - missedCleavages - 1; i++)
                {
                    if (Retain(i, initiatorMethionineBehavior, protein[0])
                        && OkayLength(oneBasedIndicesToCleaveAfter[i + missedCleavages + 1] - oneBasedIndicesToCleaveAfter[i], minPeptidesLength, maxPeptidesLength))
                    {
                        yield return new ProteolyticPeptide(protein, oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + missedCleavages + 1],
                            missedCleavages, CleavageSpecificity.Full, "full");
                    }
                    if (Cleave(i, initiatorMethionineBehavior, protein[0])
                        && OkayLength(oneBasedIndicesToCleaveAfter[i + missedCleavages + 1] - 1, minPeptidesLength, maxPeptidesLength))
                    {
                        yield return new ProteolyticPeptide(protein, 2, oneBasedIndicesToCleaveAfter[i + missedCleavages + 1],
                            missedCleavages, CleavageSpecificity.Full, "full:M cleaved");
                    }
                }

                // Also digest using the proteolysis product start/end indices
                foreach (var proteolysisProduct in protein.ProteolysisProducts)
                {
                    if (proteolysisProduct.OneBasedBeginPosition != 1 || proteolysisProduct.OneBasedEndPosition != protein.Length)
                    {
                        int i = 0;
                        while (oneBasedIndicesToCleaveAfter[i] < proteolysisProduct.OneBasedBeginPosition)
                        {
                            i++;
                        }

                        bool startPeptide = i + missedCleavages < oneBasedIndicesToCleaveAfter.Count
                            && oneBasedIndicesToCleaveAfter[i + missedCleavages] <= proteolysisProduct.OneBasedEndPosition
                            && proteolysisProduct.OneBasedBeginPosition.HasValue
                            && OkayLength(oneBasedIndicesToCleaveAfter[i + missedCleavages] - proteolysisProduct.OneBasedBeginPosition.Value + 1, minPeptidesLength, maxPeptidesLength);
                        if (startPeptide)
                        {
                            yield return new ProteolyticPeptide(protein, proteolysisProduct.OneBasedBeginPosition.Value, oneBasedIndicesToCleaveAfter[i + missedCleavages],
                                missedCleavages, CleavageSpecificity.Full, proteolysisProduct.Type + " start");
                        }

                        while (oneBasedIndicesToCleaveAfter[i] < proteolysisProduct.OneBasedEndPosition)
                        {
                            i++;
                        }

                        bool end = i - missedCleavages - 1 >= 0
                            && oneBasedIndicesToCleaveAfter[i - missedCleavages - 1] + 1 >= proteolysisProduct.OneBasedBeginPosition
                            && proteolysisProduct.OneBasedEndPosition.HasValue
                            && OkayLength(proteolysisProduct.OneBasedEndPosition.Value - oneBasedIndicesToCleaveAfter[i - missedCleavages - 1] + 1 - 1, minPeptidesLength, maxPeptidesLength);
                        if (end)
                        {
                            yield return new ProteolyticPeptide(protein, oneBasedIndicesToCleaveAfter[i - missedCleavages - 1] + 1, proteolysisProduct.OneBasedEndPosition.Value,
                                missedCleavages, CleavageSpecificity.Full, proteolysisProduct.Type + " end");
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Gets the protein intervals based on semiSpecific digestion rules
        /// This is the classic, slow semi-specific digestion that generates each semi-specific peptide pre-search
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="initiatorMethionineBehavior"></param>
        /// <param name="maximumMissedCleavages"></param>
        /// <param name="minPeptidesLength"></param>
        /// <param name="maxPeptidesLength"></param>
        /// <returns></returns>
        private IEnumerable<ProteolyticPeptide> SemiProteolyticDigestion(Protein protein, InitiatorMethionineBehavior initiatorMethionineBehavior,
            int maximumMissedCleavages, int minPeptidesLength, int maxPeptidesLength)
        {
            List<ProteolyticPeptide> intervals = new List<ProteolyticPeptide>();
            List<int> oneBasedIndicesToCleaveAfter = GetDigestionSiteIndices(protein.BaseSequence);
            //it's possible not to go through this loop (maxMissedCleavages+1>number of indexes), and that's okay. It will get digested in the next loops (finish C/N termini)
            for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - maximumMissedCleavages - 1; i++)
            {
                bool retain = Retain(i, initiatorMethionineBehavior, protein[0]);
                bool cleave = Cleave(i, initiatorMethionineBehavior, protein[0]);
                int cTerminusProtein = oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1];
                HashSet<int> localOneBasedIndicesToCleaveAfter = new HashSet<int>();
                for (int j = i; j < i + maximumMissedCleavages + 1; j++)
                {
                    localOneBasedIndicesToCleaveAfter.Add(oneBasedIndicesToCleaveAfter[j]);
                }
                if (retain)
                {
                    intervals.AddRange(FixedTermini(oneBasedIndicesToCleaveAfter[i], cTerminusProtein, protein, cleave, retain, minPeptidesLength, maxPeptidesLength, localOneBasedIndicesToCleaveAfter));
                }

                if (cleave)
                {
                    intervals.AddRange(FixedTermini(1, cTerminusProtein, protein, cleave, retain, minPeptidesLength, maxPeptidesLength, localOneBasedIndicesToCleaveAfter));
                }
            }

            //finish C-term of protein caused by loop being "i < oneBasedIndicesToCleaveAfter.Count - maximumMissedCleavages - 1"
            int last = oneBasedIndicesToCleaveAfter.Count - 1;
            int maxIndexSemi = maximumMissedCleavages < last ? maximumMissedCleavages : last;
            //Fringe C-term peptides
            for (int i = 1; i <= maxIndexSemi; i++)
            {
                //fixedN
                int nTerminusProtein = oneBasedIndicesToCleaveAfter[last - i];
                int cTerminusProtein = oneBasedIndicesToCleaveAfter[last];
                HashSet<int> localOneBasedIndicesToCleaveAfter = new HashSet<int>();
                for (int j = 0; j < i; j++) //include zero, the c terminus
                {
                    localOneBasedIndicesToCleaveAfter.Add(oneBasedIndicesToCleaveAfter[last - j]);
                }
                for (int j = cTerminusProtein; j > nTerminusProtein; j--)//We are hitting the c-terminus here
                {
                    if (OkayLength(j - nTerminusProtein, minPeptidesLength, maxPeptidesLength))
                    {
                        intervals.Add(localOneBasedIndicesToCleaveAfter.Contains(j) ?
                            new ProteolyticPeptide(protein, nTerminusProtein + 1, j, j - nTerminusProtein, CleavageSpecificity.Full, "full") :
                            new ProteolyticPeptide(protein, nTerminusProtein + 1, j, j - nTerminusProtein, CleavageSpecificity.Semi, "semi"));
                    }
                }
            }

            //Fringe N-term peptides
            for (int i = 1; i <= maxIndexSemi; i++)
            {
                bool retain = initiatorMethionineBehavior == InitiatorMethionineBehavior.Retain;
                //fixedC
                int nTerminusProtein = retain ? oneBasedIndicesToCleaveAfter[0] : oneBasedIndicesToCleaveAfter[0] + 1; // +1 start after M (since already covered earlier)
                int cTerminusProtein = oneBasedIndicesToCleaveAfter[i];
                HashSet<int> localOneBasedIndicesToCleaveAfter = new HashSet<int>();
                for (int j = 1; j < i; j++)//j starts at 1, because zero is n terminus
                {
                    localOneBasedIndicesToCleaveAfter.Add(oneBasedIndicesToCleaveAfter[j]);
                }
                int start = nTerminusProtein + 1;//plus one to not doublecount the n terminus (in addition to the M term skip)
                for (int j = start; j < cTerminusProtein; j++)
                {
                    if (OkayLength(cTerminusProtein - j, minPeptidesLength, maxPeptidesLength)
                    && !localOneBasedIndicesToCleaveAfter.Contains(j))
                    {
                        intervals.Add(new ProteolyticPeptide(protein, j + 1, cTerminusProtein, cTerminusProtein - j, CleavageSpecificity.Semi, "semi"));
                    }
                }
            }

            // Also digest using the proteolysis product start/end indices
            // This should only be things where the proteolysis is not K/R and the
            foreach (var proteolysisProduct in protein.ProteolysisProducts)
            {
                if (proteolysisProduct.OneBasedEndPosition.HasValue && proteolysisProduct.OneBasedBeginPosition.HasValue
                    && (proteolysisProduct.OneBasedBeginPosition != 1 || proteolysisProduct.OneBasedEndPosition != protein.Length)) //if at least one side is not a terminus
                {
                    int i = 0;
                    while (oneBasedIndicesToCleaveAfter[i] < proteolysisProduct.OneBasedBeginPosition)//"<" to prevent additions if same index as residues
                    {
                        i++; //last position in protein is an index to cleave after
                    }

                    // Start peptide
                    for (int j = proteolysisProduct.OneBasedBeginPosition.Value; j < oneBasedIndicesToCleaveAfter[i]; j++)
                    {
                        if (OkayLength(j - proteolysisProduct.OneBasedBeginPosition + 1, minPeptidesLength, maxPeptidesLength))
                        {
                            intervals.Add(new ProteolyticPeptide(protein, proteolysisProduct.OneBasedBeginPosition.Value, j,
                                j - proteolysisProduct.OneBasedBeginPosition.Value, CleavageSpecificity.Full, proteolysisProduct.Type + " start"));
                        }
                    }
                    while (oneBasedIndicesToCleaveAfter[i] < proteolysisProduct.OneBasedEndPosition) //"<" to prevent additions if same index as residues, since i-- is below
                    {
                        i++;
                    }

                    //Now that we've obtained an index to cleave after that is past the proteolysis product
                    //we need to backtrack to get the index to cleave that is immediately before the the proteolysis product
                    //to do this, we will do i--
                    //In the nitch case that the proteolysis product is already an index to cleave
                    //no new peptides will be generated using this, so we will forgo i--
                    //this makes peptides of length 0, which are not generated due to the for loop
                    //removing this if statement will result in crashes from c-terminal proteolysis product end positions
                    if (oneBasedIndicesToCleaveAfter[i] != proteolysisProduct.OneBasedEndPosition)
                    {
                        i--;
                    }

                    // End
                    for (int j = oneBasedIndicesToCleaveAfter[i] + 1; j < proteolysisProduct.OneBasedEndPosition.Value; j++)
                    {
                        if (OkayLength(proteolysisProduct.OneBasedEndPosition - j + 1, minPeptidesLength, maxPeptidesLength))
                        {
                            intervals.Add(new ProteolyticPeptide(protein, j, proteolysisProduct.OneBasedEndPosition.Value,
                                proteolysisProduct.OneBasedEndPosition.Value - j, CleavageSpecificity.Full, proteolysisProduct.Type + " end"));
                        }
                    }
                }
            }
            return intervals;
        }

        /// <summary>
        /// Get protein intervals for fixed termini. 
        /// This is used for the classic, slow semi-proteolytic cleavage that generates each semi-specific peptides pre-search.
        /// </summary>
        /// <param name="nTerminusProtein"></param>
        /// <param name="cTerminusProtein"></param>
        /// <param name="protein"></param>
        /// <param name="cleave"></param>
        /// <param name="minPeptidesLength"></param>
        /// <param name="maxPeptidesLength"></param>
        /// <returns></returns>
        private static IEnumerable<ProteolyticPeptide> FixedTermini(int nTerminusProtein, int cTerminusProtein, Protein protein, bool cleave, bool retain, int minPeptidesLength, int maxPeptidesLength, HashSet<int> localOneBasedIndicesToCleaveAfter)
        {
            bool preventMethionineFromBeingDuplicated = nTerminusProtein == 1 && cleave && retain; //prevents duplicate sequences containing N-terminal methionine
            List<ProteolyticPeptide> intervals = new List<ProteolyticPeptide>();
            if (!preventMethionineFromBeingDuplicated && OkayLength(cTerminusProtein - nTerminusProtein, minPeptidesLength, maxPeptidesLength)) //adds the full length maximum cleavages, no semi
            {
                intervals.Add(new ProteolyticPeptide(protein, nTerminusProtein + 1, cTerminusProtein,
                    cTerminusProtein - nTerminusProtein, CleavageSpecificity.Full, "full" + (cleave ? ":M cleaved" : ""))); // Maximum sequence length
            }

            // Fixed termini at each internal index
            IEnumerable<int> internalIndices = Enumerable.Range(nTerminusProtein + 1, cTerminusProtein - nTerminusProtein - 1); //every residue between them, +1 so we don't double count the original full

            List<ProteolyticPeptide> fixedCTermIntervals = new List<ProteolyticPeptide>();
            if (!preventMethionineFromBeingDuplicated)
            {
                var indexesOfAcceptableLength = internalIndices.Where(j => OkayLength(cTerminusProtein - j, minPeptidesLength, maxPeptidesLength));
                foreach (var j in indexesOfAcceptableLength)
                {
                    if (localOneBasedIndicesToCleaveAfter.Contains(j) || (j == 1 && cleave)) //if cleaved on cleavable index or after initiator methionine, record as full
                    {
                        if (j == 1 && cleave) //check we're not doubling it up
                        {
                            fixedCTermIntervals.Add(new ProteolyticPeptide(protein, j + 1, cTerminusProtein, cTerminusProtein - j, CleavageSpecificity.Full, "full:M cleaved"));
                        }
                        //else //don't allow full unless cleaved, since they're covered by Cterm
                    }
                    else //record it as a semi
                    {
                        fixedCTermIntervals.Add(new ProteolyticPeptide(protein, j + 1, cTerminusProtein, cTerminusProtein - j, CleavageSpecificity.Semi, "semi" + (cleave ? ":M cleaved" : "")));
                    }
                }
            }
            IEnumerable<ProteolyticPeptide> fixedNTermIntervals =
                internalIndices
                .Where(j => OkayLength(j - nTerminusProtein, minPeptidesLength, maxPeptidesLength))
                .Select(j => localOneBasedIndicesToCleaveAfter.Contains(j) ?
                new ProteolyticPeptide(protein, nTerminusProtein + 1, j, j - nTerminusProtein, CleavageSpecificity.Full, "full" + (cleave ? ":M cleaved" : "")) :
                new ProteolyticPeptide(protein, nTerminusProtein + 1, j, j - nTerminusProtein, CleavageSpecificity.Semi, "semi" + (cleave ? ":M cleaved" : "")));

            return intervals.Concat(fixedCTermIntervals).Concat(fixedNTermIntervals);
        }

        /// <summary>
        /// Checks if select subsequence of protein matches sequence that induces cleavage
        /// </summary>
        /// <param name="proteinSequence"></param>
        /// <param name="proteinSequenceIndex"></param>
        /// <param name="sequenceInducingCleavage"></param>
        /// <returns></returns>
        private bool SequenceInducesCleavage(string proteinSequence, int proteinSequenceIndex, Tuple<string, FragmentationTerminus> sequenceInducingCleavage)
        {
            return (sequenceInducingCleavage.Item2 != FragmentationTerminus.N
                    && proteinSequenceIndex - sequenceInducingCleavage.Item1.Length + 1 >= 0
                    && proteinSequence.Substring(proteinSequenceIndex - sequenceInducingCleavage.Item1.Length + 1, sequenceInducingCleavage.Item1.Length)
                        .Equals(sequenceInducingCleavage.Item1, StringComparison.OrdinalIgnoreCase))
                || (sequenceInducingCleavage.Item2 == FragmentationTerminus.N
                    && proteinSequenceIndex + 1 + sequenceInducingCleavage.Item1.Length <= proteinSequence.Length
                    && proteinSequence.Substring(proteinSequenceIndex + 1, sequenceInducingCleavage.Item1.Length)
                        .Equals(sequenceInducingCleavage.Item1, StringComparison.OrdinalIgnoreCase));
        }

        /// <summary>
        /// Checks if select subsequence of protein matches sequence preventing cleavage
        /// </summary>
        /// <param name="proteinSequence"></param>
        /// <param name="proteinSequenceIndex"></param>
        /// <param name="sequencePreventingCleavage"></param>
        /// <returns></returns>
        private bool SequencePreventsCleavage(string proteinSequence, int proteinSequenceIndex, Tuple<string, FragmentationTerminus> sequencePreventingCleavage)
        {
            return (sequencePreventingCleavage.Item2 != FragmentationTerminus.N
                    && proteinSequenceIndex + 1 + sequencePreventingCleavage.Item1.Length <= proteinSequence.Length
                    && proteinSequence.Substring(proteinSequenceIndex + 1, sequencePreventingCleavage.Item1.Length)
                        .Equals(sequencePreventingCleavage.Item1, StringComparison.OrdinalIgnoreCase))
                || (SequencesInducingCleavage.First().Item2 == FragmentationTerminus.N
                    && proteinSequenceIndex - sequencePreventingCleavage.Item1.Length + 1 >= 0
                    && proteinSequence.Substring(proteinSequenceIndex - sequencePreventingCleavage.Item1.Length + 1, sequencePreventingCleavage.Item1.Length)
                        .Equals(sequencePreventingCleavage.Item1, StringComparison.OrdinalIgnoreCase));
        }

        /// <summary>
        /// Is length of given peptide okay, given minimum?
        /// </summary>
        /// <param name="peptideLength"></param>
        /// <param name="minPeptidesLength"></param>
        /// <returns></returns>
        private static bool OkayMinLength(int? peptideLength, int? minPeptidesLength)
        {
            return !minPeptidesLength.HasValue || peptideLength >= minPeptidesLength;
        }

        /// <summary>
        /// Is length of given peptide okay, given maximum?
        /// </summary>
        /// <param name="peptideLength"></param>
        /// <param name="maxPeptidesLength"></param>
        /// <returns></returns>
        private static bool OkayMaxLength(int? peptideLength, int? maxPeptidesLength)
        {
            return !maxPeptidesLength.HasValue || peptideLength <= maxPeptidesLength;
        }
    }
}