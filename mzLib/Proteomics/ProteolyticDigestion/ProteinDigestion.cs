using System.Collections.Generic;
using System.Linq;
using Omics.BioPolymer;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;

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
            MinPeptideLength = digestionParams.MinLength;
            MaxPeptideLength = digestionParams.MaxLength;
            AllKnownFixedModifications = allKnownFixedModifications;
            VariableModifications = variableModifications;
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
        /// True when at least one configured modification -- variable or fixed -- could abolish a
        /// cleavage this protease would otherwise perform. False makes the whole cleavage-blocking
        /// correction inert: no generation slack is bought and no peptidoform can be dropped, because
        /// nothing in the search can block anything.
        ///
        /// Both modification lists are consulted, not just the variable one. Fixed modifications are
        /// appended to the variable pattern before the drop is evaluated, so a fixed blocking
        /// modification can trigger the drop and must therefore also be able to buy the slack that
        /// replaces the dropped peptidoform.
        /// </summary>
        /// <remarks>
        /// Not memoised, deliberately: <see cref="Protease"/>, <see cref="VariableModifications"/> and
        /// <see cref="AllKnownFixedModifications"/> are all settable after construction, and a cached
        /// answer would go stale behind a caller that changed one. The cost is a handful of memoised
        /// <see cref="Modification.BlocksCleavage"/> reads against the motif list, per protein, and only
        /// on the path where the flag is already on.
        /// </remarks>
        private bool AnyConfiguredModificationCanBlockCleavage =>
            (VariableModifications ?? Enumerable.Empty<Modification>())
                .Concat(AllKnownFixedModifications ?? Enumerable.Empty<Modification>())
                .Any(modification => CleavageBlockingModifications.BlocksCleavageBy(modification, Protease));

        /// <summary>
        /// Gets the fixed-terminus "seed" peptides that MetaMorpheus's non-specific search engine uses for a fast
        /// semi-specific search. These are NOT the semi-specific peptides; use <see cref="SemiSpecificDigestion"/> for those.
        /// </summary>
        /// <remarks>
        /// With <see cref="DigestionParams.FragmentationTerminus"/> = N, each seed starts at a specific N-terminus and runs
        /// as far as the missed cleavages allow (or MaxPeptideLength, whichever is shorter). Every semi-specific peptide with
        /// that N-terminus is a prefix of the seed, and they all share its N-terminal fragment ions, so the engine scores the
        /// seed once with N-terminal ions and trims it to the length the precursor mass supports. FragmentationTerminus = C is
        /// the mirror image. A complete semi-specific search therefore needs both an N pass and a C pass.
        /// <para><b>Why only a trimming engine can use seeds.</b> Trimming works because the precursor mass leaves one
        /// unknown, where the peptide ends: the engine adds residue masses along the seed until they match the precursor.
        /// An engine that has another unknown cannot do that. A glyco search, for example, takes precursor mass minus
        /// peptide mass as the glycan mass, which is undefined for a seed, and localizes glycans with fragment ions from
        /// both termini. Those engines need the peptides themselves, which <see cref="SemiSpecificDigestion"/> gives them
        /// (FragmentationTerminus Both). Asking for seeds where peptides were needed is the silent failure fixed in #1303:
        /// glyco searches lost most of their identifications with no error.</para>
        /// <para>Only call this with N or C. <c>Protein.Digest</c>
        /// routes every other terminus to <see cref="SemiSpecificDigestion"/> (see <see cref="WantsSemiSpecificSeeds"/>); below, any
        /// terminus other than N is treated as C.</para>
        /// </remarks>
        /// <param name="protein"></param>
        /// <returns></returns>
        public IEnumerable<ProteolyticPeptide> SpeedySemiSpecificDigestion(Protein protein) //We are only getting fully specific peptides of the maximum cleaved residues here
        {
            List<ProteolyticPeptide> peptides = new List<ProteolyticPeptide>();
            List<int> oneBasedIndicesToCleaveAfter = Protease.GetDigestionSiteIndices(protein.BaseSequence); //get peptide bonds to cleave SPECIFICALLY (termini included)
            int maximumMissedCleavagesIndexShift = MaximumMissedCleavages + 1;

            //it's possible not to go through this loop (maxMissedCleavages+1>number of indexes), and that's okay. It will get digested in the next loops (finish C/N termini)
            for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - maximumMissedCleavagesIndexShift; i++)
            {
                bool retain = Protease.Retain(i, InitiatorMethionineBehavior, protein[0]);
                if (retain) //it's okay to use i instead of oneBasedIndicesToCleaveAfter[i], because the index of zero is zero and it only checks if it's the N-terminus or not
                {
                    int peptideLength = oneBasedIndicesToCleaveAfter[i + maximumMissedCleavagesIndexShift] - oneBasedIndicesToCleaveAfter[i];
                    if (peptideLength >= MinPeptideLength) //if bigger than min
                    {
                        if (peptideLength <= MaxPeptideLength) //if an acceptable length (bigger than min, smaller than max), add it
                        {
                            peptides.Add(new ProteolyticPeptide(protein, oneBasedIndicesToCleaveAfter[i] + 1, oneBasedIndicesToCleaveAfter[i + maximumMissedCleavagesIndexShift],
                                MaximumMissedCleavages, CleavageSpecificity.Full, "full"));
                        }
                        else if (DigestionParams.FragmentationTerminus == FragmentationTerminus.N) //make something with the maximum length and fixed N
                        {
                            int startIndex = oneBasedIndicesToCleaveAfter[i];
                            peptides.Add(new ProteolyticPeptide(protein, startIndex + 1, startIndex + MaxPeptideLength, MaximumMissedCleavages, CleavageSpecificity.Semi, "semi"));
                        }
                        else // FragmentationTerminus.C (Protein.Digest only sends N or C here) //make something with the maximum length and fixed C
                        {
                            int endIndex = oneBasedIndicesToCleaveAfter[i + maximumMissedCleavagesIndexShift];
                            peptides.Add(new ProteolyticPeptide(protein, endIndex - MaxPeptideLength + 1, endIndex, MaximumMissedCleavages, CleavageSpecificity.Semi, "semi"));
                        }
                    }
                }

                if (Protease.Cleave(i, InitiatorMethionineBehavior, protein[0]) && (DigestionParams.FragmentationTerminus == FragmentationTerminus.N || !retain)) //it's okay to use i instead of oneBasedIndicesToCleaveAfter[i], because the index of zero is zero and it only checks if it's the N-terminus or not
                {
                    int peptideLength = oneBasedIndicesToCleaveAfter[i + maximumMissedCleavagesIndexShift] - 1;
                    if (peptideLength >= MinPeptideLength)
                    {
                        if (peptideLength <= MaxPeptideLength)
                        {
                            peptides.Add(new ProteolyticPeptide(protein, 2, oneBasedIndicesToCleaveAfter[i + maximumMissedCleavagesIndexShift], //two is hardcoded, since M=1, so the next aa is 2 (one based)
                                MaximumMissedCleavages, CleavageSpecificity.Full, "full:M cleaved"));
                        }
                        else if (DigestionParams.FragmentationTerminus == FragmentationTerminus.N)
                        {
                            peptides.Add(new ProteolyticPeptide(protein, 2, 2 + MaxPeptideLength - 1, MaximumMissedCleavages, CleavageSpecificity.Semi, "semi"));
                        }
                        else // FragmentationTerminus.C (Protein.Digest only sends N or C here) //make something with the maximum length and fixed C
                        {
                            //kinda tricky, because we'll be creating a duplication if cleavage is variable
                            if (!Protease.Retain(i, InitiatorMethionineBehavior, protein[0])) //only if cleave, because then not made earlier during retain
                            {
                                int tempIndex = oneBasedIndicesToCleaveAfter[i + maximumMissedCleavagesIndexShift];
                                peptides.Add(new ProteolyticPeptide(protein, tempIndex - MaxPeptideLength + 1, tempIndex, MaximumMissedCleavages, CleavageSpecificity.Semi, "semi"));
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

            // The initiator-methionine rules the main loop applies at i == 0 (Protease.Retain / Protease.Cleave). They
            // matter in this loop whenever a seed starts at the protein N-terminus (startIndex == 0):
            //  - EVERY C seed made here starts there. With InitiatorMethionineBehavior.Cleave the Met is not in the sample,
            //    so the seed must start after it; a seed that keeps the Met lets post-search trimming produce peptides
            //    that cannot exist. (A C seed that keeps a Met that MAY be removed is fine: trimming reaches residue 2.)
            //  - An N seed starts there only when the protein has no more cleavage sites than the missed cleavages
            //    allowed. Then, if the Met may be removed, a second N seed must start at residue 2, or no peptide starting
            //    at residue 2 is reachable at all (an N seed is only ever trimmed at its C-terminus).
            // That extra residue-2 N seed is not needed when residue 1 is itself a cleavage site: a window already starts
            // there. Being a cleavage site does NOT let a C seed keep a Met that must be removed; the C seed still starts
            // at residue 2 (a protease that cleaves after M, such as CNBr, would otherwise emit C seeds that keep it).
            bool metMayBeRemoved = Protease.Cleave(0, InitiatorMethionineBehavior, protein[0]);
            bool metMustBeRemoved = !Protease.Retain(0, InitiatorMethionineBehavior, protein[0]);
            bool residueOneIsCleavageSite = oneBasedIndicesToCleaveAfter.Count > 1 && oneBasedIndicesToCleaveAfter[1] == 1;

            for (int i = 1; i <= maxIndexDifference; i++) //i is the difference (in indexes) between indexes (cleavages), so it needs to start at 1, or the peptide would have length = 0
            {
                int startIndex = nTerminusFragmentation ?
                    oneBasedIndicesToCleaveAfter[lastIndex - i] :
                    oneBasedIndicesToCleaveAfter[0];
                int endIndex = nTerminusFragmentation ?
                    oneBasedIndicesToCleaveAfter[lastIndex] :
                    oneBasedIndicesToCleaveAfter[i];
                bool startsAtProteinNTerminus = startIndex == 0;

                if (nTerminusFragmentation)
                {
                    if (!startsAtProteinNTerminus || !metMustBeRemoved)
                    {
                        AddWrapUpSeed(startIndex, endIndex, i, "");
                    }
                    if (startsAtProteinNTerminus && metMayBeRemoved && !residueOneIsCleavageSite)
                    {
                        AddWrapUpSeed(1, endIndex, i, ":M cleaved");
                    }
                }
                else
                {
                    AddWrapUpSeed(startsAtProteinNTerminus && metMustBeRemoved ? 1 : startIndex, endIndex, i, "");
                }
            }

            // Adds the seed that starts after residue residueBeforeSeed and ends at endIndex, shortened to MaxPeptideLength
            // from its fixed terminus when it is too long.
            void AddWrapUpSeed(int residueBeforeSeed, int endIndex, int i, string descriptionSuffix)
            {
                int peptideLength = endIndex - residueBeforeSeed;
                if (peptideLength < 1 || peptideLength < MinPeptideLength)
                {
                    return;
                }
                if (peptideLength <= MaxPeptideLength) //if okay length, add it up to the terminus
                {
                    peptides.Add(new ProteolyticPeptide(protein, residueBeforeSeed + 1, endIndex, i - 1, CleavageSpecificity.Full, "full" + descriptionSuffix));
                }
                else //update so that not the end of terminus
                {
                    int startIndex = residueBeforeSeed;
                    if (nTerminusFragmentation)
                    {
                        endIndex = startIndex + MaxPeptideLength;
                    }
                    else
                    {
                        startIndex = endIndex - MaxPeptideLength;
                    }
                    peptides.Add(new ProteolyticPeptide(protein, startIndex + 1, endIndex, i - 1, CleavageSpecificity.Semi, "semi" + descriptionSuffix));
                }
            }

            // Also digest using the proteolysis product start/end indices
            foreach (TruncationProduct product in protein.TruncationProducts)
            {
                //if fixed N, we care if the start position is novel
                if (DigestionParams.FragmentationTerminus == FragmentationTerminus.N)
                {
                    //if has value and not a duplicate
                    if (product.OneBasedBeginPosition.HasValue && !oneBasedIndicesToCleaveAfter.Contains(product.OneBasedBeginPosition.Value - 1))
                    {
                        int proteaseClevageIndex = 0;

                        //get the first cleavage index after the start of the proteolysis product
                        while (oneBasedIndicesToCleaveAfter[proteaseClevageIndex] < product.OneBasedBeginPosition.Value)
                        {
                            proteaseClevageIndex++;
                        }
                        //add max missed cleavages
                        proteaseClevageIndex += MaximumMissedCleavages;

                        //set to the end if we overshot
                        if (proteaseClevageIndex >= oneBasedIndicesToCleaveAfter.Count)
                        {
                            proteaseClevageIndex = oneBasedIndicesToCleaveAfter.Count - 1;
                        }
                        int endIndex = oneBasedIndicesToCleaveAfter[proteaseClevageIndex];

                        //set to product end value if cleavages extend past
                        if (product.OneBasedEndPosition.HasValue && product.OneBasedEndPosition.Value < endIndex)
                        {
                            endIndex = product.OneBasedEndPosition.Value;
                        }

                        //limit length to the maximum allowed if necessary
                        if (endIndex - product.OneBasedBeginPosition.Value >= MaxPeptideLength)
                        {
                            endIndex = product.OneBasedBeginPosition.Value + MaxPeptideLength - 1;
                        }

                        //if it's bigger than the minimum allowed, then add it
                        if (endIndex - product.OneBasedBeginPosition.Value + 1 >= MinPeptideLength)
                        {
                            peptides.Add(new ProteolyticPeptide(protein, product.OneBasedBeginPosition.Value, endIndex, MaximumMissedCleavages, CleavageSpecificity.Full, product.Type + " start"));
                        }
                    }
                }
                else //if fixed C, we care if the end position is novel
                {
                    //if has value and not a duplicate
                    if (product.OneBasedEndPosition.HasValue && !oneBasedIndicesToCleaveAfter.Contains(product.OneBasedEndPosition.Value))
                    {
                        int proteaseClevageIndex = 0;

                        //get the first cleavage index after the start of the proteolysis product
                        while (oneBasedIndicesToCleaveAfter[proteaseClevageIndex] < product.OneBasedEndPosition.Value)
                        {
                            proteaseClevageIndex++;
                        }
                        //subtract max missed cleavages
                        proteaseClevageIndex -= (MaximumMissedCleavages + 1); //+1 because we overshot in the while loop

                        //set to the beginning if we overshot
                        if (proteaseClevageIndex < 0)
                        {
                            proteaseClevageIndex = 0;
                        }
                        int beginIndex = oneBasedIndicesToCleaveAfter[proteaseClevageIndex] + 1;

                        //set to product end value if cleavages extend past
                        if (product.OneBasedBeginPosition.HasValue && product.OneBasedBeginPosition.Value > beginIndex)
                        {
                            beginIndex = product.OneBasedBeginPosition.Value;
                        }

                        //limit length to the maximum allowed if necessary
                        if (product.OneBasedEndPosition.Value - beginIndex >= MaxPeptideLength)
                        {
                            beginIndex = product.OneBasedEndPosition.Value - MaxPeptideLength + 1;
                        }
                        //if it's bigger than the minimum allowed, then add it
                        if (product.OneBasedEndPosition.Value - beginIndex + 1 >= MinPeptideLength)
                        {
                            peptides.Add(new ProteolyticPeptide(protein, beginIndex, product.OneBasedEndPosition.Value, MaximumMissedCleavages, CleavageSpecificity.Full, product.Type + " start"));
                        }
                    }
                }
            }

            return peptides;
        }

        /// <summary>
        /// Gets every semi-specific peptide of a protein: every peptide with at least one terminus made by the protease
        /// (or a protein terminus, or a removed initiator methionine), within the missed-cleavage and length limits.
        /// Fully specific peptides are included and labelled <see cref="CleavageSpecificity.Full"/>; the rest are
        /// labelled <see cref="CleavageSpecificity.Semi"/>.
        /// </summary>
        /// <remarks>
        /// This is what <see cref="DigestionParams.SearchModeType"/> = <see cref="CleavageSpecificity.Semi"/> asks for
        /// when the caller does not want fixed-terminus seeds (see <see cref="WantsSemiSpecificSeeds"/>). It produces the
        /// same peptides as a protease whose own specificity is Semi, because both go through the protease's
        /// semi-specific enumeration with this protease's cleavage motifs.
        /// </remarks>
        public IEnumerable<ProteolyticPeptide> SemiSpecificDigestion(Protein protein)
        {
            return Protease.GetSemiSpecificUnmodifiedPeptides(protein, MaximumMissedCleavages, InitiatorMethionineBehavior, MinPeptideLength, MaxPeptideLength);
        }

        /// <summary>
        /// True when a Semi search mode is asking for <see cref="SpeedySemiSpecificDigestion"/>'s fixed-terminus seeds
        /// rather than for the semi-specific peptides themselves.
        /// </summary>
        /// <remarks>
        /// Seeds exist for MetaMorpheus's non-specific search engine, which fixes one terminus per pass, scores ions from
        /// that terminus only, and trims each seed to the length the precursor mass supports. That engine always asks
        /// with <see cref="FragmentationTerminus.N"/> or <see cref="FragmentationTerminus.C"/>, so those two values mean
        /// "seeds". Every other terminus (Both, the default) means "the semi-specific peptides", because nothing trims
        /// them afterwards. Proteases that have no cleavage motifs to be semi-specific about (top-down, singleN, singleC)
        /// keep the seed path they have always taken.
        /// </remarks>
        public static bool WantsSemiSpecificSeeds(DigestionParams digestionParams)
        {
            bool fixedTerminusRequested = digestionParams.FragmentationTerminus is FragmentationTerminus.N or FragmentationTerminus.C;
            bool proteaseCanBeSemiSpecific = digestionParams.Protease.CleavageSpecificity is CleavageSpecificity.Full or CleavageSpecificity.Semi;
            return fixedTerminusRequested || !proteaseCanBeSemiSpecific;
        }

        /// <summary>
        /// Gets peptides for specific protease digestion of a protein
        /// </summary>
        /// <param name="protein"></param>
        /// <returns></returns>
        public IEnumerable<ProteolyticPeptide> Digestion(Protein protein, bool topDownTruncationSearch = false)
        {
            // Generation uses MaximumMissedCleavages -- the instance property SpeedySemiSpecificDigestion
            // also reads, so the two digestion paths stay in step even if a caller mutates it after
            // construction. When cleavage-blocking modifications are respected in full-specificity mode we
            // add slack so the read-through form of a blocked cleavage can be generated (it costs one extra
            // missed cleavage per blocked site); the surplus is trimmed again by the open-site filter in
            // ProteolyticPeptide.GetModifiedPeptides. The slack is MaxMods, not a fixed 2: a peptidoform
            // can carry at most MaxMods variable modifications and therefore at most that many blocked
            // sites, so this reaches every variable-mod read-through no matter how many co-occur. (A fixed
            // blocking modification is unbounded and remains a documented limitation.)
            //
            // The slack is bought with enumeration, so it is only granted when it can actually be spent:
            // a search whose configured modifications cannot block a cleavage of THIS protease pays
            // nothing. That gate matters at MetaMorpheus defaults, where MaxMissedCleavages 2 plus
            // MaxMods 2 enumerates at 4 -- roughly 1.7x the unmodified peptides before modification
            // combinatorics, and it carries into the fragment index.
            int generationMaxMissedCleavages = MaximumMissedCleavages;
            if (DigestionParams.RespectCleavageBlockingModifications
                && DigestionParams.SearchModeType == CleavageSpecificity.Full
                && AnyConfiguredModificationCanBlockCleavage)
                generationMaxMissedCleavages += DigestionParams.MaxMods;

            return Protease.GetUnmodifiedPeptides(protein, generationMaxMissedCleavages, InitiatorMethionineBehavior, MinPeptideLength, MaxPeptideLength, DigestionParams.SpecificProtease, topDownTruncationSearch);
        }
    }
}