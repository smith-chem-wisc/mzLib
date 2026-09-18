using System;
using System.Collections.Generic;
using System.Linq;
using Omics;
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

        /// <summary>
        /// Gets intervals of a protein sequence that will result from digestion by this protease.
        /// </summary>
        /// <param name="protein"></param>
        /// <param name="maximumMissedCleavages"></param>
        /// <param name="initiatorMethionineBehavior"></param>
        /// <param name="minPeptideLength"></param>
        /// <param name="maxPeptideLength"></param>
        /// <returns></returns>
        internal IEnumerable<ProteolyticPeptide> GetUnmodifiedPeptides(Protein protein, int maximumMissedCleavages, InitiatorMethionineBehavior initiatorMethionineBehavior,
            int minPeptideLength, int maxPeptideLength, Protease specificProtease, bool topDownTruncationSearch = false)
        {
            bool retainMethionine = initiatorMethionineBehavior != InitiatorMethionineBehavior.Cleave || protein[0] != 'M';
            bool cleaveMethionine = initiatorMethionineBehavior != InitiatorMethionineBehavior.Retain && protein[0] == 'M';
            int initialStartResidue = retainMethionine ? 1 : 2; //where does the protein start?
            int? alternateInitialStartResidue = retainMethionine && cleaveMethionine ? 2 : null;

            return CleavageSpecificity switch
            {
                // proteolytic cleavage in one spot (N)
                CleavageSpecificity.SingleN => SingleLeftSideDigestion(protein, maximumMissedCleavages, minPeptideLength, maxPeptideLength, specificProtease, initialStartResidue).Cast<ProteolyticPeptide>(),

                // proteolytic cleavage in one spot (C)
                CleavageSpecificity.SingleC => SingleRightSideDigestion(protein, maximumMissedCleavages, minPeptideLength, maxPeptideLength, specificProtease, initialStartResidue).Cast<ProteolyticPeptide>(),

                //top-down
                CleavageSpecificity.None => TopDownDigestion(protein, minPeptideLength, maxPeptideLength, topDownTruncationSearch, initialStartResidue, alternateInitialStartResidue, CleavageSpecificity.None, "full").Cast<ProteolyticPeptide>(),

                // Full proteolytic cleavage
                CleavageSpecificity.Full => FullDigestion(protein, initiatorMethionineBehavior, maximumMissedCleavages, minPeptideLength, maxPeptideLength),

                // Cleavage rules for semi-specific search
                CleavageSpecificity.Semi => SemiProteolyticDigestion(protein, initiatorMethionineBehavior, maximumMissedCleavages, minPeptideLength, maxPeptideLength),
                _ => throw new NotImplementedException()
            };
        }

        /// <summary>
        /// Gets every semi-specific peptide this protease's cleavage motifs allow, whatever this protease's own
        /// <see cref="DigestionAgent.CleavageSpecificity"/> is. Used when a fully specific protease is combined with
        /// <see cref="DigestionParams.SearchModeType"/> = <see cref="CleavageSpecificity.Semi"/>, so that asking for a
        /// semi-specific digest through the search mode and through a Semi protease gives identical peptides.
        /// </summary>
        internal IEnumerable<ProteolyticPeptide> GetSemiSpecificUnmodifiedPeptides(Protein protein, int maximumMissedCleavages,
            InitiatorMethionineBehavior initiatorMethionineBehavior, int minPeptideLength, int maxPeptideLength)
        {
            return SemiProteolyticDigestion(protein, initiatorMethionineBehavior, maximumMissedCleavages, minPeptideLength, maxPeptideLength);
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
                foreach (var proteolysisProduct in protein.TruncationProducts)
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
            foreach (var proteolysisProduct in protein.TruncationProducts)
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

            // Every peptide below reports its missed cleavages from this table: the protease's own cleavage sites that
            // fall inside the peptide. They used to be reported as a residue distance (effectively the peptide length),
            // which is neither this protease's count nor even the non-specific convention of length - 1, and which
            // disagreed with fully specific digestion for the very same peptide. See CountCleavageSitesBefore.
            int[] cleavageSitesBefore = CountCleavageSitesBefore(oneBasedIndicesToCleaveAfter, protein.Length);

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
                    intervals.AddRange(FixedTermini(oneBasedIndicesToCleaveAfter[i], cTerminusProtein, protein, cleave, retain, minPeptideLength, maxPeptideLength, localOneBasedIndicesToCleaveAfter, cleavageSitesBefore));
                }

                if (cleave)
                {
                    intervals.AddRange(FixedTermini(1, cTerminusProtein, protein, cleave, retain, minPeptideLength, maxPeptideLength, localOneBasedIndicesToCleaveAfter, cleavageSitesBefore));
                }

                // The first window, when the initiator Met must be removed (InitiatorMethionineBehavior.Cleave) and residue 1
                // is itself a cleavage site (a protease that cleaves after M, such as CNBr). Then neither branch above runs:
                // `retain` is off because the Met is not in the sample, and `cleave` is off because the peptides starting at
                // residue 2 come from the next window, which starts at that site. But the peptides that end at THIS window's
                // C-terminus with a ragged N-terminus come from no other window or end loop, so add them here. Their starts
                // run from residue 3; a start directly after one of this window's sites is fully specific and made elsewhere.
                // Without this, those semi-specific peptides were silently missing (13 of 14 for MPEPTIDEPEPTIDE with one
                // missed cleavage). See SemiDigestion_ProteaseThatCleavesAfterTheInitiatorMet_ReturnsExactlyTheReferencePeptides.
                if (i == 0 && !retain && oneBasedIndicesToCleaveAfter[1] == 1)
                {
                    for (int j = 2; j < cTerminusProtein; j++)
                    {
                        if (!localOneBasedIndicesToCleaveAfter.Contains(j) && ValidLength(cTerminusProtein - j, minPeptideLength, maxPeptideLength))
                        {
                            intervals.Add(new ProteolyticPeptide(protein, j + 1, cTerminusProtein,
                                cleavageSitesBefore[cTerminusProtein] - cleavageSitesBefore[j + 1], CleavageSpecificity.Semi, "semi"));
                        }
                    }
                }
            }

            // Finish C-term of protein caused by loop being "i < oneBasedIndicesToCleaveAfter.Count - maximumMissedCleavages - 1"
            int last = oneBasedIndicesToCleaveAfter.Count - 1;
            int maxIndexSemi = maximumMissedCleavages < last ? maximumMissedCleavages : last;

            // The initiator-methionine rules the main loop applies at the protein N-terminus. The fringe loops below must
            // apply the same rules, because when a protein has no more cleavage sites than the missed cleavages allowed,
            // the main loop never runs and the fringe loops are the only place those peptides come from.
            //   metMayBeRemoved:  residue 1 is Met and InitiatorMethionineBehavior is Variable or Cleave, so a peptide
            //                     starting at residue 2 has a specific N-terminus.
            //   metMustBeRemoved: residue 1 is Met and InitiatorMethionineBehavior is Cleave, so no peptide may start at
            //                     residue 1 at all.
            bool metMayBeRemoved = Cleave(0, initiatorMethionineBehavior, protein[0]);
            bool metMustBeRemoved = !Retain(0, initiatorMethionineBehavior, protein[0]);

            // Fringe C-term peptides: a specific N-terminus in one of the last cleavage windows, with every C-terminus up
            // to the end of the protein.
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

                // nTerminusProtein is 0 only when this window starts at the protein N-terminus, which (see above) only
                // happens when the main loop did not run. Apply the initiator-Met rules there: skip residue 1 when the Met
                // must be removed, and also start at residue 2 when it may be removed. The residue-2 start is not added
                // when residue 1 is itself a cleavage site, because window last - i == 1 already starts there.
                bool startsAtProteinNTerminus = nTerminusProtein == 0;
                if (!startsAtProteinNTerminus || !metMustBeRemoved)
                {
                    AddFixedNTerminusPeptides(nTerminusProtein, "");
                }
                if (startsAtProteinNTerminus && metMayBeRemoved && oneBasedIndicesToCleaveAfter[1] != 1)
                {
                    AddFixedNTerminusPeptides(1, ":M cleaved");
                }

                void AddFixedNTerminusPeptides(int residueBeforePeptide, string descriptionSuffix)
                {
                    for (int j = cTerminusProtein; j > residueBeforePeptide; j--)//We are hitting the c-terminus here
                    {
                        if (ValidLength(j - residueBeforePeptide, minPeptideLength, maxPeptideLength))
                        {
                            int missedCleavages = cleavageSitesBefore[j] - cleavageSitesBefore[residueBeforePeptide + 1];
                            intervals.Add(localOneBasedIndicesToCleaveAfter.Contains(j) ?
                                new ProteolyticPeptide(protein, residueBeforePeptide + 1, j, missedCleavages, CleavageSpecificity.Full, "full" + descriptionSuffix) :
                                new ProteolyticPeptide(protein, residueBeforePeptide + 1, j, missedCleavages, CleavageSpecificity.Semi, "semi" + descriptionSuffix));
                        }
                    }
                }
            }

            // Fringe N-term peptides: a specific C-terminus at one of the first cleavage sites, with a start that is not a
            // cleavage site. Starts that ARE specific come from the loops above and are skipped here: residue 1 (the
            // protein N-terminus) always, and residue 2 when an initiator Met may be removed. Residue 2 used to be skipped
            // whenever InitiatorMethionineBehavior was not Retain, even for a protein that does not start with Met, which
            // dropped every semi peptide starting at residue 2 of such proteins.
            for (int i = 1; i <= maxIndexSemi; i++)
            {
                // FixedC
                int nTerminusProtein = metMayBeRemoved ? oneBasedIndicesToCleaveAfter[0] + 1 : oneBasedIndicesToCleaveAfter[0]; // +1 start after M (since already covered earlier)
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
                        intervals.Add(new ProteolyticPeptide(protein, j + 1, cTerminusProtein, cleavageSitesBefore[cTerminusProtein] - cleavageSitesBefore[j + 1], CleavageSpecificity.Semi, "semi"));
                    }
                }
            }

            // Also digest using the proteolysis product start/end indices
            // This should only be things where the proteolysis is not K/R and the
            foreach (var proteolysisProduct in protein.TruncationProducts)
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
                            intervals.Add(new ProteolyticPeptide(protein, start, j, cleavageSitesBefore[j] - cleavageSitesBefore[start], CleavageSpecificity.Full, proteolysisProduct.Type + " start"));
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
                            intervals.Add(new ProteolyticPeptide(protein, j, end, cleavageSitesBefore[end] - cleavageSitesBefore[j], CleavageSpecificity.Full, proteolysisProduct.Type + " end"));
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
        /// <param name="cleavageSitesBefore">The table from <see cref="CountCleavageSitesBefore"/> for this protein, used to
        /// report each peptide's true missed cleavages.</param>
        /// <returns></returns>
        private static IEnumerable<ProteolyticPeptide> FixedTermini(int nTerminusProtein, int cTerminusProtein, Protein protein, bool cleave, bool retain, int minPeptideLength, int maxPeptideLength, HashSet<int> localOneBasedIndicesToCleaveAfter, int[] cleavageSitesBefore)
        {
            // Missed cleavages of the peptide from oneBasedStart to oneBasedEnd: this protease's sites inside it.
            int MissedCleavages(int oneBasedStart, int oneBasedEnd) => cleavageSitesBefore[oneBasedEnd] - cleavageSitesBefore[oneBasedStart];

            bool preventMethionineFromBeingDuplicated = nTerminusProtein == 1 && cleave && retain; //prevents duplicate sequences containing N-terminal methionine
            List<ProteolyticPeptide> intervals = new List<ProteolyticPeptide>();
            if (!preventMethionineFromBeingDuplicated && ValidLength(cTerminusProtein - nTerminusProtein, minPeptideLength, maxPeptideLength)) //adds the full length maximum cleavages, no semi
            {
                intervals.Add(new ProteolyticPeptide(protein, nTerminusProtein + 1, cTerminusProtein,
                    MissedCleavages(nTerminusProtein + 1, cTerminusProtein), CleavageSpecificity.Full, "full" + (cleave ? ":M cleaved" : ""))); // Maximum sequence length
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
                            fixedCTermIntervals.Add(new ProteolyticPeptide(protein, j + 1, cTerminusProtein, MissedCleavages(j + 1, cTerminusProtein), CleavageSpecificity.Full, "full:M cleaved"));
                        }
                        //else //don't allow full unless cleaved, since they're covered by Cterm
                    }
                    else //record it as a semi
                    {
                        fixedCTermIntervals.Add(new ProteolyticPeptide(protein, j + 1, cTerminusProtein, MissedCleavages(j + 1, cTerminusProtein), CleavageSpecificity.Semi, "semi" + (cleave ? ":M cleaved" : "")));
                    }
                }
            }
            IEnumerable<ProteolyticPeptide> fixedNTermIntervals =
                internalIndices
                .Where(j => ValidLength(j - nTerminusProtein, minPeptideLength, maxPeptideLength))
                .Select(j => localOneBasedIndicesToCleaveAfter.Contains(j) ?
                new ProteolyticPeptide(protein, nTerminusProtein + 1, j, MissedCleavages(nTerminusProtein + 1, j), CleavageSpecificity.Full, "full" + (cleave ? ":M cleaved" : "")) :
                new ProteolyticPeptide(protein, nTerminusProtein + 1, j, MissedCleavages(nTerminusProtein + 1, j), CleavageSpecificity.Semi, "semi" + (cleave ? ":M cleaved" : "")));

            return intervals.Concat(fixedCTermIntervals).Concat(fixedNTermIntervals);
        }

        protected override IEnumerable<DigestionProduct> GetConcreteProducts(IBioPolymer protein, int startResidue, int endResidue, int missedCleavages, CleavageSpecificity specificity, string description)
        {
            yield return new ProteolyticPeptide((Protein)protein, startResidue, endResidue, missedCleavages, specificity, description);
        }
    }
}
