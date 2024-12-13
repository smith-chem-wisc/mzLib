using System;
using System.Collections.Generic;
using System.Linq;
using Omics.Digestion;
using Omics.Modifications;

namespace Proteomics.ProteolyticDigestion
{
    public class Protease : DigestionAgent
    {
        public Protease(string name, CleavageSpecificity cleavageSpecificity, string psiMSAccessionNumber, 
            string psiMSName, List<DigestionMotif> motifList, Modification modDetails = null) 
            : base(name, cleavageSpecificity, motifList, modDetails)
        {
            PsiMsAccessionNumber = psiMSAccessionNumber;
            PsiMsName = psiMSName;
        }

        public string PsiMsAccessionNumber { get; }
        public string PsiMsName { get; }

        public override string ToString()
        {
            return Name;
        }

        public override bool Equals(object obj)
        {
            return obj is Protease a
                && (a.Name == null && Name == null || a.Name.Equals(Name));
        }

        public override int GetHashCode()
        {
            return (Name ?? "").GetHashCode();
        }

        /// <summary>
        /// This method is used to determine cleavage specificity if the cleavage specificity is unknown
        /// This occurs in the speedy nonspecific/semispecific searches when digesting post-search
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="startIndex"></param>
        /// <param name="endIndex"></param>
        /// <param name="retainMethionine"></param>
        /// <returns></returns>
        public CleavageSpecificity GetCleavageSpecificity(Protein protein, int startIndex, int endIndex, bool retainMethionine)
        {
            int cleavableMatches = 0;
            if (CleavageSpecificity != CleavageSpecificity.SingleN && CleavageSpecificity != CleavageSpecificity.SingleC) //if it's single protease, don't bother
            {
                List<int> indicesToCleave = GetDigestionSiteIndices(protein.BaseSequence);
                //if the start index is a cleavable index (-1 because one based) OR if the start index is after a cleavable methionine
                if (indicesToCleave.Contains(startIndex - 1) ||
                    (startIndex == 2 && protein.BaseSequence[0] == 'M' && !retainMethionine) ||
                    protein.ProteolysisProducts.Any(x => x.OneBasedBeginPosition == startIndex))
                {
                    cleavableMatches++;
                }
                //if the end index is a cleavable index
                if (indicesToCleave.Contains(endIndex) ||
                    protein.ProteolysisProducts.Any(x => x.OneBasedEndPosition == endIndex))
                {
                    cleavableMatches++;
                }
            }
            if (cleavableMatches == 0) //if neither were cleavable, (or it's singleN/C) then it's nonspecific
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
        /// <param name="minPeptideLength"></param>
        /// <param name="maxPeptideLength"></param>
        /// <returns></returns>
        internal List<ProteolyticPeptide> GetUnmodifiedPeptides(Protein protein, int maximumMissedCleavages, InitiatorMethionineBehavior initiatorMethionineBehavior,
            int minPeptideLength, int maxPeptideLength, Protease specificProtease, bool topDownTruncationSearch = false)
        {
            List<ProteolyticPeptide> peptides;

            // proteolytic cleavage in one spot (N)
            if (CleavageSpecificity == CleavageSpecificity.SingleN)
            {
                peptides = SingleN_Digestion(protein, initiatorMethionineBehavior, maximumMissedCleavages, minPeptideLength, maxPeptideLength, specificProtease);
            }

            // proteolytic cleavage in one spot (C)
            else if (CleavageSpecificity == CleavageSpecificity.SingleC)
            {
                peptides = SingleC_Digestion(protein, initiatorMethionineBehavior, maximumMissedCleavages, minPeptideLength, maxPeptideLength, specificProtease);
            }

            //top-down
            else if (CleavageSpecificity == CleavageSpecificity.None)
            {
                peptides = new(20);
                if (!topDownTruncationSearch)//standard top-down
                {
                    // retain methionine
                    if ((initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || protein[0] != 'M')
                        && ValidLength(protein.Length, minPeptideLength, maxPeptideLength))
                    {
                        peptides.Add(new ProteolyticPeptide(protein, 1, protein.Length, 0, CleavageSpecificity.Full, "full"));
                    }

                    // cleave methionine
                    if ((initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && protein[0] == 'M')
                        && ValidLength(protein.Length - 1, minPeptideLength, maxPeptideLength))
                    {
                        peptides.Add(new ProteolyticPeptide(protein, 2, protein.Length, 0, CleavageSpecificity.Full, "full:M cleaved"));
                    }
                }

                // Also digest using the proteolysis product start/end indices
                peptides.AddRange(
                    protein.ProteolysisProducts
                    .Where(proteolysisProduct => proteolysisProduct.OneBasedEndPosition.HasValue && proteolysisProduct.OneBasedBeginPosition.HasValue
                        && ValidLength(proteolysisProduct.OneBasedEndPosition.Value - proteolysisProduct.OneBasedBeginPosition.Value + 1, minPeptideLength, maxPeptideLength))
                    .Select(proteolysisProduct =>
                        new ProteolyticPeptide(protein, proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.OneBasedEndPosition.Value, 0, CleavageSpecificity.None, proteolysisProduct.Type)));
            }

            // Full proteolytic cleavage
            else if (CleavageSpecificity == CleavageSpecificity.Full)
            {
                peptides = FullDigestion(protein, initiatorMethionineBehavior, maximumMissedCleavages, minPeptideLength, maxPeptideLength).ToList();
            }

            // Cleavage rules for semi-specific search
            else if (CleavageSpecificity == CleavageSpecificity.Semi)
            {
                peptides = SemiProteolyticDigestion(protein, initiatorMethionineBehavior, maximumMissedCleavages, minPeptideLength, maxPeptideLength).ToList();
            }
            else
            {
                throw new NotImplementedException();
            }

            return peptides;
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
        /// Gets protein intervals for digestion by this specific protease.
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="initiatorMethionineBehavior"></param>
        /// <param name="maximumMissedCleavages"></param>
        /// <param name="minPeptideLength"></param>
        /// <param name="maxPeptideLength"></param>
        /// <returns></returns>
        private IEnumerable<ProteolyticPeptide> FullDigestion(Protein protein, InitiatorMethionineBehavior initiatorMethionineBehavior,
            int maximumMissedCleavages, int minPeptideLength, int maxPeptideLength)
        {
            List<int> oneBasedIndicesToCleaveAfter = GetDigestionSiteIndices(protein.BaseSequence);
            char firstResidueInProtein = protein[0];

            for (int missedCleavages = 0; missedCleavages <= maximumMissedCleavages; missedCleavages++)
            {
                for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - missedCleavages - 1; i++)
                {
                    if (Retain(i, initiatorMethionineBehavior, firstResidueInProtein)
                        && ValidLength(oneBasedIndicesToCleaveAfter[i + missedCleavages + 1] - oneBasedIndicesToCleaveAfter[i], minPeptideLength, maxPeptideLength))
                    {
                        yield return new ProteolyticPeptide(protein, oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + missedCleavages + 1],
                            missedCleavages, CleavageSpecificity.Full, "full");
                    }
                    if (Cleave(i, initiatorMethionineBehavior, firstResidueInProtein) && oneBasedIndicesToCleaveAfter[1] != 1 //prevent duplicates if that bond is cleaved by the protease
                        && ValidLength(oneBasedIndicesToCleaveAfter[i + missedCleavages + 1] - 1, minPeptideLength, maxPeptideLength))
                    {
                        yield return new ProteolyticPeptide(protein, 2, oneBasedIndicesToCleaveAfter[i + missedCleavages + 1],
                            missedCleavages, CleavageSpecificity.Full, "full:M cleaved");
                    }
                }

                //TODO: Generate all the proteolytic products as distinct proteins during XML reading and delete all of the code below
                // Also digest using the proteolysis product start/end indices
                foreach (var proteolysisProduct in protein.ProteolysisProducts)
                {
                    //if the proteolysis product contains something other than just the start AND end residues of the protein
                    if (proteolysisProduct.OneBasedBeginPosition != 1 || proteolysisProduct.OneBasedEndPosition != protein.Length)
                    {
                        int cleavageIndexWithinProteolysisProduct = 0;
                        //get the first cleavage index after the start of the proteolysis product
                        while (oneBasedIndicesToCleaveAfter[cleavageIndexWithinProteolysisProduct] < proteolysisProduct.OneBasedBeginPosition)
                        {
                            cleavageIndexWithinProteolysisProduct++;
                        }

                        bool startPeptide = cleavageIndexWithinProteolysisProduct + missedCleavages < oneBasedIndicesToCleaveAfter.Count //if the current missed cleavages doesn't hit the end
                            && oneBasedIndicesToCleaveAfter[cleavageIndexWithinProteolysisProduct + missedCleavages] <= proteolysisProduct.OneBasedEndPosition //and the cleavage occurs before the proteolytic end
                            && proteolysisProduct.OneBasedBeginPosition.HasValue //and the proteolytic peptide even has a beginning
                            && !oneBasedIndicesToCleaveAfter.Contains(proteolysisProduct.OneBasedBeginPosition.Value - 1) //and we haven't already cleaved here
                            && (proteolysisProduct.OneBasedBeginPosition.Value != 1 || !Cleave(0, initiatorMethionineBehavior, firstResidueInProtein)) //and it's not the initiator methionine
                            && ValidLength(oneBasedIndicesToCleaveAfter[cleavageIndexWithinProteolysisProduct + missedCleavages] - proteolysisProduct.OneBasedBeginPosition.Value + 1, minPeptideLength, maxPeptideLength); //and it's the correct size
                        if (startPeptide)
                        {
                            yield return new ProteolyticPeptide(protein, proteolysisProduct.OneBasedBeginPosition.Value, oneBasedIndicesToCleaveAfter[cleavageIndexWithinProteolysisProduct + missedCleavages],
                                missedCleavages, CleavageSpecificity.Full, proteolysisProduct.Type + " start");
                        }

                        //get the cleavage index before the end of the proteolysis product
                        while (oneBasedIndicesToCleaveAfter[cleavageIndexWithinProteolysisProduct] < proteolysisProduct.OneBasedEndPosition)
                        {
                            cleavageIndexWithinProteolysisProduct++;
                        }

                        bool endPeptide = cleavageIndexWithinProteolysisProduct - missedCleavages - 1 >= 0 //if we're not going to go out of bounds (-1 to get in front of the end)
                            && oneBasedIndicesToCleaveAfter[cleavageIndexWithinProteolysisProduct - missedCleavages - 1] + 1 >= proteolysisProduct.OneBasedBeginPosition //and it's not before the beginning
                            && proteolysisProduct.OneBasedEndPosition.HasValue //and the proteolytic peptide even has an end
                            && !oneBasedIndicesToCleaveAfter.Contains(proteolysisProduct.OneBasedEndPosition.Value) //and we haven't already cleaved here
                            && ValidLength(proteolysisProduct.OneBasedEndPosition.Value - oneBasedIndicesToCleaveAfter[cleavageIndexWithinProteolysisProduct - missedCleavages - 1] + 1 - 1, minPeptideLength, maxPeptideLength); //and it's the correct size
                        if (endPeptide)
                        {
                            yield return new ProteolyticPeptide(protein, oneBasedIndicesToCleaveAfter[cleavageIndexWithinProteolysisProduct - missedCleavages - 1] + 1, proteolysisProduct.OneBasedEndPosition.Value,
                                missedCleavages, CleavageSpecificity.Full, proteolysisProduct.Type + " end");
                        }
                    }
                }
            }

            //add intact proteolysis products (if acceptable)
            foreach (var proteolysisProduct in protein.ProteolysisProducts)
            {
                if (proteolysisProduct.OneBasedBeginPosition.HasValue //begin has value
                    && proteolysisProduct.OneBasedEndPosition.HasValue //and end has value
                    && (proteolysisProduct.OneBasedBeginPosition.Value != 1 || !Cleave(0, initiatorMethionineBehavior, firstResidueInProtein)) //and it's not the initiator methionine
                    && !oneBasedIndicesToCleaveAfter.Contains(proteolysisProduct.OneBasedBeginPosition.Value - 1) //and we haven't already cleaved here
                    && !oneBasedIndicesToCleaveAfter.Contains(proteolysisProduct.OneBasedEndPosition.Value)) //and we haven't already cleaved there
                {
                    int firstCleavage = 0;
                    //get the first cleavage index after the start of the proteolysis product
                    while (oneBasedIndicesToCleaveAfter[firstCleavage] < proteolysisProduct.OneBasedBeginPosition)
                    {
                        firstCleavage++;
                    }

                    int lastCleavage = firstCleavage;
                    //get the last cleavage index before the end of the proteolysis product
                    while (oneBasedIndicesToCleaveAfter[lastCleavage] < proteolysisProduct.OneBasedEndPosition)
                    {
                        lastCleavage++;
                    }
                    if (lastCleavage - firstCleavage < maximumMissedCleavages && //if there aren't too many missed cleavages
                        ValidLength(proteolysisProduct.OneBasedEndPosition.Value - proteolysisProduct.OneBasedBeginPosition.Value, minPeptideLength, maxPeptideLength)) //and it's the correct size
                    {
                        yield return new ProteolyticPeptide(protein, proteolysisProduct.OneBasedBeginPosition.Value, proteolysisProduct.OneBasedEndPosition.Value,
                            lastCleavage - firstCleavage, CleavageSpecificity.Full, proteolysisProduct.Type + " end");
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
        /// <param name="minPeptideLength"></param>
        /// <param name="maxPeptideLength"></param>
        /// <returns></returns>
        private IEnumerable<ProteolyticPeptide> SemiProteolyticDigestion(Protein protein, InitiatorMethionineBehavior initiatorMethionineBehavior,
            int maximumMissedCleavages, int minPeptideLength, int maxPeptideLength)
        {
            List<ProteolyticPeptide> intervals = new List<ProteolyticPeptide>();
            List<int> oneBasedIndicesToCleaveAfter = GetDigestionSiteIndices(protein.BaseSequence);

            // It's possible not to go through this loop (maxMissedCleavages+1>number of indexes), and that's okay. It will get digested in the next loops (finish C/N termini)
            for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - maximumMissedCleavages - 1; i++)
            {
                bool retain = Retain(i, initiatorMethionineBehavior, protein[0]);
                bool cleave = Cleave(i, initiatorMethionineBehavior, protein[0]) && oneBasedIndicesToCleaveAfter[1] != 1;
                int cTerminusProtein = oneBasedIndicesToCleaveAfter[i + maximumMissedCleavages + 1];
                HashSet<int> localOneBasedIndicesToCleaveAfter = new HashSet<int>();
                for (int j = i; j < i + maximumMissedCleavages + 1; j++)
                {
                    localOneBasedIndicesToCleaveAfter.Add(oneBasedIndicesToCleaveAfter[j]);
                }
                if (retain)
                {
                    intervals.AddRange(FixedTermini(oneBasedIndicesToCleaveAfter[i], cTerminusProtein, protein, cleave, retain, minPeptideLength, maxPeptideLength, localOneBasedIndicesToCleaveAfter));
                }

                if (cleave)
                {
                    intervals.AddRange(FixedTermini(1, cTerminusProtein, protein, cleave, retain, minPeptideLength, maxPeptideLength, localOneBasedIndicesToCleaveAfter));
                }
            }

            // Finish C-term of protein caused by loop being "i < oneBasedIndicesToCleaveAfter.Count - maximumMissedCleavages - 1"
            int last = oneBasedIndicesToCleaveAfter.Count - 1;
            int maxIndexSemi = maximumMissedCleavages < last ? maximumMissedCleavages : last;
            // Fringe C-term peptides
            for (int i = 1; i <= maxIndexSemi; i++)
            {
                // FixedN
                int nTerminusProtein = oneBasedIndicesToCleaveAfter[last - i];
                int cTerminusProtein = oneBasedIndicesToCleaveAfter[last];
                HashSet<int> localOneBasedIndicesToCleaveAfter = new HashSet<int>();
                for (int j = 0; j < i; j++) //include zero, the c terminus
                {
                    localOneBasedIndicesToCleaveAfter.Add(oneBasedIndicesToCleaveAfter[last - j]);
                }
                for (int j = cTerminusProtein; j > nTerminusProtein; j--)//We are hitting the c-terminus here
                {
                    if (ValidLength(j - nTerminusProtein, minPeptideLength, maxPeptideLength))
                    {
                        intervals.Add(localOneBasedIndicesToCleaveAfter.Contains(j) ?
                            new ProteolyticPeptide(protein, nTerminusProtein + 1, j, j - nTerminusProtein, CleavageSpecificity.Full, "full") :
                            new ProteolyticPeptide(protein, nTerminusProtein + 1, j, j - nTerminusProtein, CleavageSpecificity.Semi, "semi"));
                    }
                }
            }

            // Fringe N-term peptides
            for (int i = 1; i <= maxIndexSemi; i++)
            {
                bool retain = initiatorMethionineBehavior == InitiatorMethionineBehavior.Retain;
                // FixedC
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
                    if (ValidLength(cTerminusProtein - j, minPeptideLength, maxPeptideLength)
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
                    int start = proteolysisProduct.OneBasedBeginPosition.Value;
                    int end = proteolysisProduct.OneBasedEndPosition.Value;
                    int i = 0;
                    while (oneBasedIndicesToCleaveAfter[i] < start)//"<" to prevent additions if same index as residues
                    {
                        i++; // Last position in protein is an index to cleave after
                    }

                    // Start peptide
                    for (int j = start; j < oneBasedIndicesToCleaveAfter[i]; j++)
                    {
                        if (ValidLength(j - start + 1, minPeptideLength, maxPeptideLength))
                        {
                            intervals.Add(new ProteolyticPeptide(protein, start, j, j - start, CleavageSpecificity.Full, proteolysisProduct.Type + " start"));
                        }
                    }
                    while (oneBasedIndicesToCleaveAfter[i] < end) //"<" to prevent additions if same index as residues, since i-- is below
                    {
                        i++;
                    }

                    // Now that we've obtained an index to cleave after that is past the proteolysis product
                    // we need to backtrack to get the index to cleave that is immediately before the the proteolysis product
                    // to do this, we will do i--
                    // In the nitch case that the proteolysis product is already an index to cleave
                    // no new peptides will be generated using this, so we will forgo i--
                    // this makes peptides of length 0, which are not generated due to the for loop
                    // removing this if statement will result in crashes from c-terminal proteolysis product end positions
                    if (oneBasedIndicesToCleaveAfter[i] != end)
                    {
                        i--;
                    }

                    // Fin (End)
                    for (int j = oneBasedIndicesToCleaveAfter[i] + 1; j < end; j++)
                    {
                        if (ValidLength(end - j + 1, minPeptideLength, maxPeptideLength))
                        {
                            intervals.Add(new ProteolyticPeptide(protein, j, end, end - j, CleavageSpecificity.Full, proteolysisProduct.Type + " end"));
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
        /// <param name="minPeptideLength"></param>
        /// <param name="maxPeptideLength"></param>
        /// <returns></returns>
        private static IEnumerable<ProteolyticPeptide> FixedTermini(int nTerminusProtein, int cTerminusProtein, Protein protein, bool cleave, bool retain, int minPeptideLength, int maxPeptideLength, HashSet<int> localOneBasedIndicesToCleaveAfter)
        {
            bool preventMethionineFromBeingDuplicated = nTerminusProtein == 1 && cleave && retain; //prevents duplicate sequences containing N-terminal methionine
            List<ProteolyticPeptide> intervals = new List<ProteolyticPeptide>();
            if (!preventMethionineFromBeingDuplicated && ValidLength(cTerminusProtein - nTerminusProtein, minPeptideLength, maxPeptideLength)) //adds the full length maximum cleavages, no semi
            {
                intervals.Add(new ProteolyticPeptide(protein, nTerminusProtein + 1, cTerminusProtein,
                    cTerminusProtein - nTerminusProtein, CleavageSpecificity.Full, "full" + (cleave ? ":M cleaved" : ""))); // Maximum sequence length
            }

            // Fixed termini at each internal index
            IEnumerable<int> internalIndices = Enumerable.Range(nTerminusProtein + 1, cTerminusProtein - nTerminusProtein - 1); //every residue between them, +1 so we don't double count the original full

            List<ProteolyticPeptide> fixedCTermIntervals = new List<ProteolyticPeptide>();
            if (!preventMethionineFromBeingDuplicated)
            {
                var indexesOfAcceptableLength = internalIndices.Where(j => ValidLength(cTerminusProtein - j, minPeptideLength, maxPeptideLength));
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
                .Where(j => ValidLength(j - nTerminusProtein, minPeptideLength, maxPeptideLength))
                .Select(j => localOneBasedIndicesToCleaveAfter.Contains(j) ?
                new ProteolyticPeptide(protein, nTerminusProtein + 1, j, j - nTerminusProtein, CleavageSpecificity.Full, "full" + (cleave ? ":M cleaved" : "")) :
                new ProteolyticPeptide(protein, nTerminusProtein + 1, j, j - nTerminusProtein, CleavageSpecificity.Semi, "semi" + (cleave ? ":M cleaved" : "")));

            return intervals.Concat(fixedCTermIntervals).Concat(fixedNTermIntervals);
        }

        /// <summary>
        /// Gets peptides for the singleN protease
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="initiatorMethionineBehavior"></param>
        /// <param name="maximumMissedCleavages"></param>
        /// <param name="minPeptideLength"></param>
        /// <param name="maxPeptideLength"></param>
        /// <param name="specificProtease"></param>
        /// <returns></returns>
        private List<ProteolyticPeptide> SingleN_Digestion(Protein protein, InitiatorMethionineBehavior initiatorMethionineBehavior,
            int maximumMissedCleavages, int minPeptideLength, int maxPeptideLength, Protease specificProtease)
        {
            List<ProteolyticPeptide> peptides = new List<ProteolyticPeptide>();
            int proteinStart = Retain(0, initiatorMethionineBehavior, protein[0]) ? 1 : 2; //where does the protein start?

            if (Equals(specificProtease))
            {
                bool maxTooBig = protein.Length + maxPeptideLength < 0; //when maxPeptideLength is too large, it becomes negative and causes issues
                                                                        //This happens when maxPeptideLength == int.MaxValue or something close to it
                for (; proteinStart <= protein.Length; proteinStart++)
                {
                    if (ValidMinLength(protein.Length - proteinStart + 1, minPeptideLength))
                    {
                        //need Math.Max if max length is int.MaxLength, since +proteinStart will make it negative
                        //if the max length is too big to be an int (ie infinity), just do the protein length.
                        //if it's not too big to be an int, it might still be too big. Take the minimum of the protein length or the maximum length (-1, because the index is inclusive. Without -1, peptides will be one AA too long)
                        peptides.Add(new ProteolyticPeptide(protein, proteinStart, maxTooBig ? protein.Length : Math.Min(protein.Length, proteinStart + maxPeptideLength - 1), 0, CleavageSpecificity.SingleN, "SingleN"));
                    }
                }
            }
            else //if there's a specific protease, then we need to adhere to the specified missed cleavage rules
            {
                //generate only peptides with the maximum number of missed cleavages, unless the protein has fewer than the max or we're near the unselected terminus (where we run to the end of the protein)
                List<int> oneBasedIndicesToCleaveAfter = specificProtease.GetDigestionSiteIndices(protein.BaseSequence); //get peptide bonds to cleave SPECIFICALLY (termini included)
                oneBasedIndicesToCleaveAfter[0] = proteinStart - 1;//update the first cleavage to represent the initiator methionine rules
                int maximumMissedCleavagesIndexShift = maximumMissedCleavages + 1;

                for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - maximumMissedCleavagesIndexShift; i++)
                {
                    int startIndex = oneBasedIndicesToCleaveAfter[i];
                    int endProteaseIndex = oneBasedIndicesToCleaveAfter[i + maximumMissedCleavagesIndexShift];
                    int peptideLength = endProteaseIndex - startIndex;
                    if (peptideLength >= minPeptideLength) //if bigger than min
                    {
                        int endActualIndex = endProteaseIndex;
                        if (peptideLength > maxPeptideLength) //if the next cleavage is too far away, crop it to the max length
                        {
                            endActualIndex = startIndex + maxPeptideLength;
                        }
                        int nextStartIndex = oneBasedIndicesToCleaveAfter[i + 1] + 1;

                        //make SingleN peptides until we reach the next index to cleave at or until the peptides are too small
                        for (; (startIndex + 1 < nextStartIndex) && (endActualIndex - startIndex >= minPeptideLength); startIndex++)
                        {
                            peptides.Add(new ProteolyticPeptide(protein, startIndex + 1, endActualIndex, maximumMissedCleavages, CleavageSpecificity.SingleN, "SingleN"));

                            //update endIndex if needed
                            if (endActualIndex != endProteaseIndex)
                            {
                                endActualIndex++;
                            }
                        }
                    }
                }
                //wrap up the terminus
                if (oneBasedIndicesToCleaveAfter.Count < maximumMissedCleavagesIndexShift)
                {
                    maximumMissedCleavagesIndexShift = oneBasedIndicesToCleaveAfter.Count;
                }
                int lastStartIndex = oneBasedIndicesToCleaveAfter[oneBasedIndicesToCleaveAfter.Count - maximumMissedCleavagesIndexShift] + 1;
                int proteinEndIndex = oneBasedIndicesToCleaveAfter[oneBasedIndicesToCleaveAfter.Count - 1]; //end of protein
                int lastEndIndex = Math.Min(proteinEndIndex, lastStartIndex + maxPeptideLength - 1); //end of protein
                for (; lastStartIndex + minPeptideLength - 1 <= lastEndIndex; lastStartIndex++)
                {
                    peptides.Add(new ProteolyticPeptide(protein, lastStartIndex, lastEndIndex, maximumMissedCleavages, CleavageSpecificity.SingleN, "SingleN"));

                    //update the end if needed
                    if (lastEndIndex != proteinEndIndex)
                    {
                        lastEndIndex++;
                    }
                }
            }
            return peptides;
        }

        /// <summary>
        /// Gets peptides for the singleC protease
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="initiatorMethionineBehavior"></param>
        /// <param name="maximumMissedCleavages"></param>
        /// <param name="minPeptideLength"></param>
        /// <param name="maxPeptideLength"></param>
        /// <param name="specificProtease"></param>
        /// <returns></returns>
        private List<ProteolyticPeptide> SingleC_Digestion(Protein protein, InitiatorMethionineBehavior initiatorMethionineBehavior,
            int maximumMissedCleavages, int minPeptideLength, int maxPeptideLength, Protease specificProtease)
        {
            List<ProteolyticPeptide> peptides = new List<ProteolyticPeptide>();
            int proteinStart = Retain(0, initiatorMethionineBehavior, protein[0]) ? 1 : 2; //where does the protein start?
            if (Equals(specificProtease))
            {
                int lengthDifference = proteinStart - 1; //take it back one for zero based index
                for (int proteinEnd = 1; proteinEnd <= protein.Length; proteinEnd++)
                {
                    //length of peptide will be at least the start index
                    if (ValidMinLength(proteinEnd - lengthDifference, minPeptideLength)) //is the maximum possible length longer than the minimum?
                    {
                        //use the start index as the max of the N-terminus or the c-terminus minus the max (+1 because inclusive, otherwise peptides will be one AA too long)
                        peptides.Add(new ProteolyticPeptide(protein, Math.Max(proteinStart, proteinEnd - maxPeptideLength + 1), proteinEnd, 0, CleavageSpecificity.SingleC, "SingleC"));
                    }
                }
            }
            else //if there's a specific protease, then we need to adhere to the specified missed cleavage rules
            {
                //generate only peptides with the maximum number of missed cleavages, unless the protein has fewer than the max or we're near the unselected terminus (where we run to the end of the protein)
                List<int> oneBasedIndicesToCleaveAfter = specificProtease.GetDigestionSiteIndices(protein.BaseSequence); //get peptide bonds to cleave SPECIFICALLY (termini included)
                oneBasedIndicesToCleaveAfter[0] = proteinStart - 1;//update the first cleavage to represent the initiator methionine rules
                int maximumMissedCleavagesIndexShift = maximumMissedCleavages + 1;

                for (int i = oneBasedIndicesToCleaveAfter.Count - 1; i > maximumMissedCleavagesIndexShift; i--)
                {
                    int endProteaseIndex = oneBasedIndicesToCleaveAfter[i];
                    int startProteaseIndex = oneBasedIndicesToCleaveAfter[i - maximumMissedCleavagesIndexShift];
                    int peptideLength = endProteaseIndex - startProteaseIndex;
                    if (peptideLength >= minPeptideLength) //if bigger than min
                    {
                        int startActualIndex = startProteaseIndex;
                        if (peptideLength > maxPeptideLength) //if the next cleavage is too far away, crop it to the max length
                        {
                            startActualIndex = endProteaseIndex - maxPeptideLength;
                        }
                        int nextEndIndex = oneBasedIndicesToCleaveAfter[i - 1];
                        //make SingleC peptides until we reach the next index to cleave at or until the peptides are too small
                        for (; (endProteaseIndex > nextEndIndex) && (endProteaseIndex - startActualIndex >= minPeptideLength); endProteaseIndex--)
                        {
                            peptides.Add(new ProteolyticPeptide(protein, startActualIndex + 1, endProteaseIndex, maximumMissedCleavages, CleavageSpecificity.SingleC, "SingleC"));

                            //update startIndex if needed
                            if (startActualIndex != startProteaseIndex)
                            {
                                startActualIndex--;
                            }
                        }
                    }
                }
                //wrap up the terminus
                //if there are more missed cleavages allowed than there are cleavages to cleave, change the effective number of missed cleavages to the max
                if (oneBasedIndicesToCleaveAfter.Count <= maximumMissedCleavagesIndexShift)
                {
                    maximumMissedCleavagesIndexShift = oneBasedIndicesToCleaveAfter.Count - 1;
                }
                int lastEndIndex = oneBasedIndicesToCleaveAfter[maximumMissedCleavagesIndexShift];
                int startIndex = Math.Max(proteinStart, lastEndIndex - maxPeptideLength + 1);
                int minPeptideLengthOneBasedResidueShift = minPeptideLength - 1;
                for (; lastEndIndex >= startIndex + minPeptideLengthOneBasedResidueShift; lastEndIndex--)
                {
                    peptides.Add(new ProteolyticPeptide(protein, startIndex, lastEndIndex, maximumMissedCleavages, CleavageSpecificity.SingleC, "SingleC"));

                    //update the start if needed
                    if (startIndex != proteinStart)
                    {
                        startIndex--;
                    }
                }
            }
            return peptides;
        }
    }
}