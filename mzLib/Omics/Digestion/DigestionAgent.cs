using MzLibUtil;
using Omics.Modifications;

namespace Omics.Digestion
{
    public abstract class DigestionAgent
    {
        protected static readonly HashSetPool<int> HashSetPool = new HashSetPool<int>(8);

        protected DigestionAgent(string name, CleavageSpecificity cleavageSpecificity, List<DigestionMotif> motifList, Modification cleavageMod)
        {
            Name = name;
            CleavageSpecificity = cleavageSpecificity;
            DigestionMotifs = motifList ?? new List<DigestionMotif>();
            CleavageMod = cleavageMod;
        }

        public readonly string Name;
        public CleavageSpecificity CleavageSpecificity { get; init; }
        public List<DigestionMotif> DigestionMotifs { get; init; }
        public Modification CleavageMod { get; set; }




        protected abstract IEnumerable<DigestionProduct> GetConcreteProducts(IBioPolymer parent, int oneBasedStartResidue, int oneBasedEndResidue, int missedCleavages, CleavageSpecificity specificity, string description);

        #region Digestion Methods


        /// <summary>
        /// Gets peptides for the singleN protease
        /// </summary>
        protected IEnumerable<DigestionProduct> SingleLeftSideDigestion(IBioPolymer protein, int maximumMissedCleavages, int minLength, int maxLength, DigestionAgent specificAgent, int startResidueIndex)
        {
            if (Equals(specificAgent))
            {
                bool maxTooBig = protein.Length + maxLength < 0; //when maxPeptideLength is too large, it becomes negative and causes issues
                                                                 //This happens when maxPeptideLength == int.MaxValue or something close to it
                for (; startResidueIndex <= protein.Length; startResidueIndex++)
                {
                    if (ValidMinLength(protein.Length - startResidueIndex + 1, minLength))
                    {
                        //need Math.Max if max length is int.MaxLength, since +proteinStart will make it negative
                        //if the max length is too big to be an int (ie infinity), just do the protein length.
                        //if it's not too big to be an int, it might still be too big. Take the minimum of the protein length or the maximum length (-1, because the index is inclusive. Without -1, peptides will be one AA too long)
                        int endResidue = maxTooBig ? protein.Length : Math.Min(protein.Length, startResidueIndex + maxLength - 1);
                        foreach (var product in GetConcreteProducts(protein, startResidueIndex, endResidue, 0, CleavageSpecificity.SingleN, "SingleN"))
                        {
                            yield return product;
                        }
                    }
                }
            }
            else //if there's a specific protease, then we need to adhere to the specified missed cleavage rules
            {
                //generate only peptides with the maximum number of missed cleavages, unless the protein has fewer than the max or we're near the unselected terminus (where we run to the end of the protein)
                List<int> oneBasedIndicesToCleaveAfter = specificAgent.GetDigestionSiteIndices(protein.BaseSequence); //get peptide bonds to cleave SPECIFICALLY (termini included)
                oneBasedIndicesToCleaveAfter[0] = startResidueIndex - 1;//update the first cleavage to represent the initiator methionine rules
                int maximumMissedCleavagesIndexShift = maximumMissedCleavages + 1;

                for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - maximumMissedCleavagesIndexShift; i++)
                {
                    int startIndex = oneBasedIndicesToCleaveAfter[i];
                    int endProteaseIndex = oneBasedIndicesToCleaveAfter[i + maximumMissedCleavagesIndexShift];
                    int peptideLength = endProteaseIndex - startIndex;
                    if (peptideLength >= minLength) //if bigger than min
                    {
                        int endActualIndex = endProteaseIndex;
                        if (peptideLength > maxLength) //if the next cleavage is too far away, crop it to the max length
                        {
                            endActualIndex = startIndex + maxLength;
                        }
                        int nextStartIndex = oneBasedIndicesToCleaveAfter[i + 1] + 1;

                        //make SingleN peptides until we reach the next index to cleave at or until the peptides are too small
                        for (; (startIndex + 1 < nextStartIndex) && (endActualIndex - startIndex >= minLength); startIndex++)
                        {
                            foreach (var product in GetConcreteProducts(protein, startIndex + 1, endActualIndex, maximumMissedCleavages, CleavageSpecificity.SingleN, "SingleN"))
                            {
                                yield return product;
                            }

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
                int lastEndIndex = Math.Min(proteinEndIndex, lastStartIndex + maxLength - 1); //end of protein
                for (; lastStartIndex + minLength - 1 <= lastEndIndex; lastStartIndex++)
                {
                    foreach (var product in GetConcreteProducts(protein, lastStartIndex, lastEndIndex, maximumMissedCleavages, CleavageSpecificity.SingleN, "SingleN"))
                    {
                        yield return product;
                    }

                    //update the end if needed
                    if (lastEndIndex != proteinEndIndex)
                    {
                        lastEndIndex++;
                    }
                }
            }
        }

        /// <summary>
        /// Gets peptides for the singleC protease
        /// </summary>
        protected IEnumerable<DigestionProduct> SingleRightSideDigestion(IBioPolymer parent, int maximumMissedCleavages, int minLength, int maxLength, DigestionAgent specificAgent, int startResidueIndex)
        {
            if (Equals(specificAgent))
            {
                int lengthDifference = startResidueIndex - 1; //take it back one for zero based index
                for (int parentEnd = 1; parentEnd <= parent.Length; parentEnd++)
                {
                    //length of peptide will be at least the start index
                    if (ValidMinLength(parentEnd - lengthDifference, minLength)) //is the maximum possible length longer than the minimum?
                    {
                        //use the start index as the max of the N-terminus or the c-terminus minus the max (+1 because inclusive, otherwise peptides will be one AA too long)
                        foreach (var product in GetConcreteProducts(parent, Math.Max(startResidueIndex, parentEnd - maxLength + 1), parentEnd, 0, CleavageSpecificity.SingleC, "SingleC"))
                        {
                            yield return product;
                        }
                    }
                }
            }
            else //if there's a specific protease, then we need to adhere to the specified missed cleavage rules
            {
                //generate only peptides with the maximum number of missed cleavages, unless the protein has fewer than the max or we're near the unselected terminus (where we run to the end of the protein)
                List<int> oneBasedIndicesToCleaveAfter = specificAgent.GetDigestionSiteIndices(parent.BaseSequence); //get peptide bonds to cleave SPECIFICALLY (termini included)
                oneBasedIndicesToCleaveAfter[0] = startResidueIndex - 1;//update the first cleavage to represent the initiator methionine rules
                int maximumMissedCleavagesIndexShift = maximumMissedCleavages + 1;

                for (int i = oneBasedIndicesToCleaveAfter.Count - 1; i > maximumMissedCleavagesIndexShift; i--)
                {
                    int endProteaseIndex = oneBasedIndicesToCleaveAfter[i];
                    int startProteaseIndex = oneBasedIndicesToCleaveAfter[i - maximumMissedCleavagesIndexShift];
                    int peptideLength = endProteaseIndex - startProteaseIndex;
                    if (peptideLength >= minLength) //if bigger than min
                    {
                        int startActualIndex = startProteaseIndex;
                        if (peptideLength > maxLength) //if the next cleavage is too far away, crop it to the max length
                        {
                            startActualIndex = endProteaseIndex - maxLength;
                        }
                        int nextEndIndex = oneBasedIndicesToCleaveAfter[i - 1];
                        //make SingleC peptides until we reach the next index to cleave at or until the peptides are too small
                        for (; (endProteaseIndex > nextEndIndex) && (endProteaseIndex - startActualIndex >= minLength); endProteaseIndex--)
                        {
                            foreach (var product in GetConcreteProducts(parent, startActualIndex + 1, endProteaseIndex, maximumMissedCleavages, CleavageSpecificity.SingleC, "SingleC"))
                            {
                                yield return product;
                            }

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
                int startIndex = Math.Max(startResidueIndex, lastEndIndex - maxLength + 1);
                int minPeptideLengthOneBasedResidueShift = minLength - 1;
                for (; lastEndIndex >= startIndex + minPeptideLengthOneBasedResidueShift; lastEndIndex--)
                {
                    foreach (var product in GetConcreteProducts(parent, startIndex, lastEndIndex, maximumMissedCleavages, CleavageSpecificity.SingleC, "SingleC"))
                    {
                        yield return product;
                    }

                    //update the start if needed
                    if (startIndex != startResidueIndex)
                    {
                        startIndex--;
                    }
                }
            }

        }

        protected IEnumerable<DigestionProduct> TopDownDigestion(IBioPolymer parent, int minLength, int maxLength, bool topDownTruncationSearch, int initialStartResidue, int? alternateInitialStartResidue, CleavageSpecificity truncationSpecificity, string initialDescription)
        {
            if (!topDownTruncationSearch)
            {
                if (ValidLength(parent.Length - initialStartResidue + 1, minLength, maxLength))
                {
                    foreach (var product in GetConcreteProducts(parent, initialStartResidue, parent.Length, 0,
                                 CleavageSpecificity.Full, initialDescription))
                    {
                        yield return product;
                    }
                }

                if (alternateInitialStartResidue.HasValue
                    && ValidLength(parent.Length - alternateInitialStartResidue.Value + 1, minLength, maxLength))
                {
                    foreach (var product in GetConcreteProducts(parent, alternateInitialStartResidue.Value, parent.Length, 0,
                                 CleavageSpecificity.Full, initialDescription + ":M cleaved"))
                    {
                        yield return product;
                    }
                }
            }

            foreach (var truncationProduct in parent.TruncationProducts)
            {
                if (truncationProduct.OneBasedBeginPosition.HasValue
                    && truncationProduct.OneBasedEndPosition.HasValue
                    && ValidLength(
                        truncationProduct.OneBasedEndPosition.Value - truncationProduct.OneBasedBeginPosition.Value + 1,
                        minLength,
                        maxLength))
                {
                    foreach (var product in GetConcreteProducts(
                                 parent,
                                 truncationProduct.OneBasedBeginPosition.Value,
                                 truncationProduct.OneBasedEndPosition.Value,
                                 0,
                                 truncationSpecificity,
                                 truncationProduct.Type))
                    {
                        yield return product;
                    }
                }
            }
        }

        protected IEnumerable<DigestionProduct> FullDigestion(IBioPolymer parent, int maximumMissedCleavages, int minLength, int maxLength, int initialStartResidue, int? alternateInitialStartResidue, CleavageSpecificity truncationSpecificity, string initialDescription)
        {
            List<int> cleavageIndices = GetDigestionSiteIndices(parent.BaseSequence);

            for (int missedCleavages = 0; missedCleavages <= maximumMissedCleavages; missedCleavages++)
            {
                for (int i = 0; i < cleavageIndices.Count - missedCleavages - 1; i++)
                {
                    int endResidue = cleavageIndices[i + missedCleavages + 1];
                    int startResidue = i == 0
                        ? Math.Max(initialStartResidue, cleavageIndices[i] + 1)
                        : cleavageIndices[i] + 1;

                    if (ValidLength(endResidue - startResidue + 1, minLength, maxLength))
                        foreach (var product in GetConcreteProducts(parent, startResidue, endResidue, missedCleavages, CleavageSpecificity.Full, initialDescription))
                            yield return product;

                    if (i == 0
                        && alternateInitialStartResidue.HasValue
                        && !cleavageIndices.Contains(alternateInitialStartResidue.Value - 1))
                    {
                        if (ValidLength(endResidue - alternateInitialStartResidue.Value + 1, minLength, maxLength))
                            foreach (var product in GetConcreteProducts(parent, alternateInitialStartResidue.Value, endResidue, missedCleavages, CleavageSpecificity.Full, initialDescription + ":M cleaved"))
                                yield return product;
                    }
                }

                foreach (var truncationProduct in parent.TruncationProducts)
                {
                    if (truncationProduct.OneBasedBeginPosition != 1 || truncationProduct.OneBasedEndPosition != parent.Length)
                    {
                        int cleavageIndex = 0;
                        while (cleavageIndices[cleavageIndex] < truncationProduct.OneBasedBeginPosition)
                        {
                            cleavageIndex++;
                        }

                        bool startProduct = cleavageIndex + missedCleavages < cleavageIndices.Count
                            && truncationProduct.OneBasedBeginPosition.HasValue
                            && truncationProduct.OneBasedEndPosition.HasValue
                            && cleavageIndices[cleavageIndex + missedCleavages] <= truncationProduct.OneBasedEndPosition
                            && !cleavageIndices.Contains(truncationProduct.OneBasedBeginPosition.Value - 1)
                            && truncationProduct.OneBasedBeginPosition.Value >= initialStartResidue
                            && ValidLength(cleavageIndices[cleavageIndex + missedCleavages]
                                - truncationProduct.OneBasedBeginPosition.Value + 1, minLength, maxLength);

                        if (startProduct)
                        {
                            foreach (var product in GetConcreteProducts(parent, truncationProduct.OneBasedBeginPosition.Value,
                                         cleavageIndices[cleavageIndex + missedCleavages], missedCleavages,
                                         CleavageSpecificity.Full, truncationProduct.Type + " start"))
                            {
                                yield return product;
                            }
                        }

                        while (cleavageIndices[cleavageIndex] < truncationProduct.OneBasedEndPosition)
                        {
                            cleavageIndex++;
                        }

                        bool endProduct = cleavageIndex - missedCleavages - 1 >= 0
                            && truncationProduct.OneBasedBeginPosition.HasValue
                            && truncationProduct.OneBasedEndPosition.HasValue
                            && cleavageIndices[cleavageIndex - missedCleavages - 1] + 1 >= truncationProduct.OneBasedBeginPosition
                            && !cleavageIndices.Contains(truncationProduct.OneBasedEndPosition.Value)
                            && ValidLength(truncationProduct.OneBasedEndPosition.Value
                                - cleavageIndices[cleavageIndex - missedCleavages - 1], minLength, maxLength);

                        if (endProduct)
                        {
                            foreach (var product in GetConcreteProducts(parent,
                                         cleavageIndices[cleavageIndex - missedCleavages - 1] + 1,
                                         truncationProduct.OneBasedEndPosition.Value, missedCleavages,
                                         CleavageSpecificity.Full, truncationProduct.Type + " end"))
                            {
                                yield return product;
                            }
                        }
                    }
                }
            }

            foreach (var truncationProduct in parent.TruncationProducts)
            {
                if (!truncationProduct.OneBasedBeginPosition.HasValue
                    || !truncationProduct.OneBasedEndPosition.HasValue
                    || truncationProduct.OneBasedBeginPosition.Value < initialStartResidue
                    || cleavageIndices.Contains(truncationProduct.OneBasedBeginPosition.Value - 1)
                    || cleavageIndices.Contains(truncationProduct.OneBasedEndPosition.Value)
                    || !ValidLength(truncationProduct.OneBasedEndPosition.Value - truncationProduct.OneBasedBeginPosition.Value,
                        minLength, maxLength))
                {
                    continue;
                }

                int firstCleavage = 0;
                while (cleavageIndices[firstCleavage] < truncationProduct.OneBasedBeginPosition)
                {
                    firstCleavage++;
                }

                int lastCleavage = firstCleavage;
                while (cleavageIndices[lastCleavage] < truncationProduct.OneBasedEndPosition)
                {
                    lastCleavage++;
                }

                if (lastCleavage - firstCleavage < maximumMissedCleavages)
                {
                    foreach (var product in GetConcreteProducts(parent, truncationProduct.OneBasedBeginPosition.Value,
                                 truncationProduct.OneBasedEndPosition.Value, lastCleavage - firstCleavage,
                                 truncationSpecificity, truncationProduct.Type + " end"))
                    {
                        yield return product;
                    }
                }
            }
        }

        protected IEnumerable<DigestionProduct> SpeedySemiSpecificDigestion(
            IBioPolymer parent,
            int maximumMissedCleavages,
            int minLength,
            int maxLength,
            bool fixedLeftTerminus,
            int initialStartResidue,
            int? alternateInitialStartResidue)
        {
            List<int> cleavageIndices = GetDigestionSiteIndices(parent.BaseSequence);
            int cleavageWindowSize = maximumMissedCleavages + 1;

            for (int i = 0; i < cleavageIndices.Count - cleavageWindowSize; i++)
            {
                int startBoundary = cleavageIndices[i];
                int endBoundary = cleavageIndices[i + cleavageWindowSize];

                if (i == 0)
                {
                    startBoundary = Math.Max(startBoundary, initialStartResidue - 1);
                }

                foreach (var product in GetSpeedyProducts(parent, startBoundary, endBoundary, maximumMissedCleavages,
                             fixedLeftTerminus, minLength, maxLength, ""))
                {
                    yield return product;
                }

                    if (i == 0
                        && alternateInitialStartResidue.HasValue
                        && alternateInitialStartResidue.Value != initialStartResidue
                        && !cleavageIndices.Contains(alternateInitialStartResidue.Value - 1)
                        && (fixedLeftTerminus || initialStartResidue > 1))
                {
                    foreach (var product in GetSpeedyProducts(parent, alternateInitialStartResidue.Value - 1, endBoundary,
                                 maximumMissedCleavages, fixedLeftTerminus, minLength, maxLength, ":M cleaved"))
                    {
                        yield return product;
                    }
                }
            }

            int lastIndex = cleavageIndices.Count - 1;
            int maximumIndexDifference = Math.Min(maximumMissedCleavages, lastIndex);
            bool methionineMustBeRemoved = initialStartResidue > 1;
            bool methionineMayBeRemoved = alternateInitialStartResidue.HasValue;
            bool residueOneIsCleavageSite = cleavageIndices.Count > 1 && cleavageIndices[1] == 1;

            for (int i = 1; i <= maximumIndexDifference; i++)
            {
                int startBoundary = fixedLeftTerminus ? cleavageIndices[lastIndex - i] : cleavageIndices[0];
                int endBoundary = fixedLeftTerminus ? cleavageIndices[lastIndex] : cleavageIndices[i];
                bool startsAtParentBeginning = startBoundary == 0;

                if (fixedLeftTerminus)
                {
                    if (!startsAtParentBeginning || !methionineMustBeRemoved)
                    {
                        foreach (var product in GetSpeedyProducts(parent, startBoundary, endBoundary, i - 1,
                                     true, minLength, maxLength, ""))
                        {
                            yield return product;
                        }
                    }

                    if (startsAtParentBeginning && methionineMayBeRemoved && !residueOneIsCleavageSite)
                    {
                        foreach (var product in GetSpeedyProducts(parent, 1, endBoundary, i - 1,
                                     true, minLength, maxLength, ":M cleaved"))
                        {
                            yield return product;
                        }
                    }
                }
                else
                {
                    int effectiveStartBoundary = startsAtParentBeginning && methionineMustBeRemoved ? 1 : startBoundary;
                    foreach (var product in GetSpeedyProducts(parent, effectiveStartBoundary, endBoundary, i - 1,
                                 false, minLength, maxLength, ""))
                    {
                        yield return product;
                    }
                }
            }

            foreach (var truncationProduct in parent.TruncationProducts)
            {
                if (fixedLeftTerminus)
                {
                    if (!truncationProduct.OneBasedBeginPosition.HasValue
                        || cleavageIndices.Contains(truncationProduct.OneBasedBeginPosition.Value - 1))
                    {
                        continue;
                    }

                    int cleavageIndex = 0;
                    while (cleavageIndices[cleavageIndex] < truncationProduct.OneBasedBeginPosition.Value)
                    {
                        cleavageIndex++;
                    }

                    cleavageIndex = Math.Min(cleavageIndex + maximumMissedCleavages, cleavageIndices.Count - 1);
                    int startResidue = Math.Max(truncationProduct.OneBasedBeginPosition.Value, initialStartResidue);
                    int endResidue = cleavageIndices[cleavageIndex];
                    if (truncationProduct.OneBasedEndPosition.HasValue)
                    {
                        endResidue = Math.Min(endResidue, truncationProduct.OneBasedEndPosition.Value);
                    }
                    endResidue = Math.Min(endResidue, startResidue + maxLength - 1);

                    if (endResidue - startResidue + 1 >= minLength)
                    {
                        foreach (var product in GetConcreteProducts(parent, startResidue, endResidue,
                                     maximumMissedCleavages, CleavageSpecificity.Full, truncationProduct.Type + " start"))
                        {
                            yield return product;
                        }
                    }
                }
                else
                {
                    if (!truncationProduct.OneBasedEndPosition.HasValue
                        || cleavageIndices.Contains(truncationProduct.OneBasedEndPosition.Value))
                    {
                        continue;
                    }

                    int cleavageIndex = 0;
                    while (cleavageIndices[cleavageIndex] < truncationProduct.OneBasedEndPosition.Value)
                    {
                        cleavageIndex++;
                    }

                    cleavageIndex = Math.Max(cleavageIndex - maximumMissedCleavages - 1, 0);
                    int startResidue = Math.Max(cleavageIndices[cleavageIndex] + 1, initialStartResidue);
                    if (truncationProduct.OneBasedBeginPosition.HasValue)
                    {
                        startResidue = Math.Max(startResidue, truncationProduct.OneBasedBeginPosition.Value);
                    }

                    int endResidue = truncationProduct.OneBasedEndPosition.Value;
                    startResidue = Math.Max(startResidue, endResidue - maxLength + 1);
                    if (endResidue - startResidue + 1 >= minLength)
                    {
                        foreach (var product in GetConcreteProducts(parent, startResidue, endResidue,
                                     maximumMissedCleavages, CleavageSpecificity.Full, truncationProduct.Type + " start"))
                        {
                            yield return product;
                        }
                    }
                }
            }
        }

        private IEnumerable<DigestionProduct> GetSpeedyProducts(
            IBioPolymer parent,
            int startBoundary,
            int endBoundary,
            int missedCleavages,
            bool fixedLeftTerminus,
            int minLength,
            int maxLength,
            string descriptionSuffix)
        {
            int length = endBoundary - startBoundary;
            if (length < minLength)
            {
                yield break;
            }

            int startResidue = startBoundary + 1;
            int endResidue = endBoundary;
            CleavageSpecificity specificity = CleavageSpecificity.Full;
            string description = "full" + descriptionSuffix;

            if (length > maxLength)
            {
                specificity = CleavageSpecificity.Semi;
                description = "semi" + descriptionSuffix;
                if (fixedLeftTerminus)
                {
                    endResidue = startBoundary + maxLength;
                }
                else
                {
                    startResidue = endBoundary - maxLength + 1;
                }
            }

            foreach (var product in GetConcreteProducts(parent, startResidue, endResidue, missedCleavages, specificity, description))
            {
                yield return product;
            }
        }

        #endregion

        public override string ToString()
        {
            return Name;
        }

        public override bool Equals(object? obj)
        {
            return obj is DigestionAgent agent && agent.Name == Name;
        }

        public override int GetHashCode()
        {
            return Name.GetHashCode();
        }

        #region Digestion Helpers

        /// <summary>
        /// Is length of given peptide okay, given minimum and maximum?
        /// </summary>
        /// <param name="length"></param>
        /// <param name="minLength"></param>
        /// <param name="maxLength"></param>
        /// <returns></returns>
        protected static bool ValidLength(int length, int minLength, int maxLength)
        {
            return ValidMinLength(length, minLength) && ValidMaxLength(length, maxLength);
        }

        /// <summary>
        /// Is length of given peptide okay, given minimum?
        /// </summary>
        /// <param name="length"></param>
        /// <param name="minLength"></param>
        /// <returns></returns>
        protected static bool ValidMinLength(int length, int minLength)
        {
            return length >= minLength;
        }

        /// <summary>
        /// Is length of given peptide okay, given maximum?
        /// </summary>
        /// <param name="length"></param>
        /// <param name="maxLength"></param>
        /// <returns></returns>
        protected static bool ValidMaxLength(int? length, int maxLength)
        {
            return !length.HasValue || length <= maxLength;
        }


        /// <summary>
        /// This method is used to determine cleavage specificity if the cleavage specificity is unknown
        /// This occurs in the speedy nonspecific/semispecific searches when digesting post-search
        /// </summary>
        /// <returns></returns>
        public CleavageSpecificity GetCleavageSpecificity(IBioPolymer bioPolymer, int startIndex, int endIndex, bool retainMethionine)
        {
            int cleavableMatches = 0;
            if (CleavageSpecificity != CleavageSpecificity.SingleN && CleavageSpecificity != CleavageSpecificity.SingleC) //if it's single protease, don't bother
            {
                List<int> indicesToCleave = GetDigestionSiteIndices(bioPolymer.BaseSequence);
                //if the start index is a cleavable index (-1 because one based) OR if the start index is after a cleavable methionine
                if (indicesToCleave.Contains(startIndex - 1) ||
                    (startIndex == 2 && bioPolymer.BaseSequence[0] == 'M' && !retainMethionine) ||
                    bioPolymer.TruncationProducts.Any(x => x.OneBasedBeginPosition == startIndex))
                {
                    cleavableMatches++;
                }
                //if the end index is a cleavable index
                if (indicesToCleave.Contains(endIndex) ||
                    bioPolymer.TruncationProducts.Any(x => x.OneBasedEndPosition == endIndex))
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
        /// Gets the indices after which this protease will cleave a given protein sequence
        /// </summary>
        /// <param name="sequence"></param>
        /// <returns></returns>
        public List<int> GetDigestionSiteIndices(string sequence)
        {
            var indices = HashSetPool.Get(); // use hash set to ensure no duplicates
            try // Try block is to ensure that, even if an error gets thrown, the hashset is returned to the pool
            {
                indices.Add(0); // The start of the protein is treated as a cleavage site to retain the n-terminal peptide

                for (int r = 0; r < sequence.Length; r++)
                {
                    var cutSiteIndex = -1;
                    bool cleavagePrevented = false;

                    foreach (DigestionMotif motif in DigestionMotifs)
                    {
                        var motifResults = motif.Fits(sequence, r);
                        bool motifFits = motifResults.Item1;
                        bool motifPreventsCleavage = motifResults.Item2;

                        if (motifFits && r + motif.CutIndex < sequence.Length)
                        {
                            cutSiteIndex = Math.Max(r + motif.CutIndex, cutSiteIndex);
                        }

                        if (motifPreventsCleavage) // if any motif prevents cleave
                        {
                            cleavagePrevented = true;
                        }
                    }

                    // if no motif prevents cleave
                    if (!cleavagePrevented && cutSiteIndex != -1)
                    {
                        indices.Add(cutSiteIndex);
                    }
                }

                indices.Add(sequence.Length); // The end of the protein is treated as a cleavage site to retain the c-terminal peptide
                return indices.ToList(); // convert the hashset to a list for return. 
            }
            finally
            {
                // return hashset to pool. This clears it and gets it ready for the next time it is needed from the pool.
                HashSetPool.Return(indices);
            }
        }

        /// <summary>
        /// Builds a lookup that turns "how many missed cleavages does this peptide have?" into one subtraction.
        /// </summary>
        /// <remarks>
        /// A peptide's missed cleavages are this protease's cleavage sites that fall INSIDE it: a site after residue
        /// <c>k</c> with <c>start &lt;= k &lt; end</c>. The site after the peptide's own last residue is where it was cut,
        /// not a missed cleavage, and the protein's start and end (indices 0 and length) are never inside a peptide.
        /// <para>Element <c>x</c> of the returned array is the number of sites after residues <c>1..x-1</c>, so the
        /// peptide from <c>start</c> to <c>end</c> (one-based, inclusive) has
        /// <c>table[end] - table[start]</c> missed cleavages. This is the same count fully specific digestion reports,
        /// so a peptide gets the same number whichever digestion produced it.</para>
        /// <para>For the non-specific protease every residue is a site, which gives <c>length - 1</c>; for trypsin it
        /// is the number of internal K/R. Building the table once per protein keeps semi-specific digestion, which makes
        /// many peptides per protein, from walking the site list for each one; it also does not assume that
        /// <paramref name="oneBasedIndicesToCleaveAfter"/> is sorted.</para>
        /// </remarks>
        /// <param name="oneBasedIndicesToCleaveAfter">Cleavage sites from <see cref="DigestionAgent.GetDigestionSiteIndices"/>.</param>
        /// <param name="proteinLength">Number of residues in the protein.</param>
        protected static int[] CountCleavageSitesBefore(List<int> oneBasedIndicesToCleaveAfter, int proteinLength)
        {
            var isSiteAfterResidue = new bool[proteinLength + 1];
            foreach (int site in oneBasedIndicesToCleaveAfter)
            {
                if (site >= 1 && site < proteinLength)
                {
                    isSiteAfterResidue[site] = true;
                }
            }

            var sitesBefore = new int[proteinLength + 1];
            for (int x = 1; x <= proteinLength; x++)
            {
                sitesBefore[x] = sitesBefore[x - 1] + (x - 1 >= 1 && isSiteAfterResidue[x - 1] ? 1 : 0);
            }
            return sitesBefore;
        }

        #endregion
    }
}
