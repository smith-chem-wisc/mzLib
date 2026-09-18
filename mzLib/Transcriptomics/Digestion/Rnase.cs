using Chemistry;
using Omics;
using Omics.Digestion;
using Omics.Modifications;

namespace Transcriptomics.Digestion
{
    public class Rnase : DigestionAgent, IEquatable<Rnase>
    {
        public static IHasChemicalFormula DefaultThreePrimeTerminus = ChemicalFormula.ParseFormula("H2O4P"); // Makes 3' Phosphate
        public static IHasChemicalFormula DefaultFivePrimeTerminus = ChemicalFormula.ParseFormula("O-3P-1"); // Makes 5' -OH by removing phosphate

        public IList<IHasChemicalFormula> ThreePrimeTerminusRemainder { get; set; }
        public IList<IHasChemicalFormula> FivePrimeTerminusRemainder { get; set; }

        public Rnase(string name, CleavageSpecificity cleaveSpecificity, List<DigestionMotif> motifList, Modification cleavageMod = null, IList<IHasChemicalFormula>? threePrimeTerminusRemainder = null, IList<IHasChemicalFormula>? fivePrimeTerminusRemainder = null) :
            base(name, cleaveSpecificity, motifList, cleavageMod)
        {
            CleavageSpecificity = cleaveSpecificity;
            DigestionMotifs = motifList;
            ThreePrimeTerminusRemainder = threePrimeTerminusRemainder ?? new List<IHasChemicalFormula> { DefaultThreePrimeTerminus };
            FivePrimeTerminusRemainder = fivePrimeTerminusRemainder ?? new List<IHasChemicalFormula> { DefaultFivePrimeTerminus };
        }

        public IEnumerable<NucleolyticOligo> GetUnmodifiedOligos(NucleicAcid nucleicAcid, int maxMissedCleavages, int minLength,
            int maxLength, Rnase? specificRnase = null, bool topDownTruncationSearch = false)
        {
            specificRnase ??= this;
            return CleavageSpecificity switch
            {
                // top down
                CleavageSpecificity.None => TopDownDigestion(nucleicAcid, minLength, maxLength, topDownTruncationSearch, 1, null, CleavageSpecificity.Full, "full").Cast<NucleolyticOligo>(),

                // full cleavage
                CleavageSpecificity.Full => FullDigestion(nucleicAcid, maxMissedCleavages, minLength, maxLength),

                // non-specific, anchored at one terminus
                CleavageSpecificity.SingleN => SingleLeftSideDigestion(nucleicAcid, maxMissedCleavages, minLength, maxLength, specificRnase, 1).Cast<NucleolyticOligo>(),

                CleavageSpecificity.SingleC => SingleRightSideDigestion(nucleicAcid, maxMissedCleavages, minLength, maxLength, specificRnase, 1).Cast<NucleolyticOligo>(),
                _ => throw new ArgumentException(
                    "Cleave Specificity not defined for Rna digestion, currently supports Full, None, SingleN and SingleC")
            };
        }

        protected override IEnumerable<DigestionProduct> GetConcreteProducts(IBioPolymer rna, int startResidue, int endResidue, int missedCleavages, CleavageSpecificity specificity, string description)
        {
            foreach (var (threePrimeTerminus, fivePrimeTerminus) in GetDigestedTermini(startResidue, endResidue, (NucleicAcid)rna, ThreePrimeTerminusRemainder, FivePrimeTerminusRemainder))
            {
                yield return new NucleolyticOligo((NucleicAcid)rna, startResidue, endResidue,
                    missedCleavages, specificity, fivePrimeTerminus, threePrimeTerminus, description);
            }
        }


        private IEnumerable<NucleolyticOligo> FullDigestion(NucleicAcid nucleicAcid, int maxMissedCleavages,
            int minLength, int maxLength)
        {
            List<int> oneBasedIndicesToCleaveAfter = GetDigestionSiteIndices(nucleicAcid.BaseSequence);
            for (int missedCleavages = 0; missedCleavages <= maxMissedCleavages; missedCleavages++)
            {
                for (int i = 0; i < oneBasedIndicesToCleaveAfter.Count - missedCleavages - 1; i++)
                {
                    if (ValidLength(oneBasedIndicesToCleaveAfter[i + missedCleavages + 1] - oneBasedIndicesToCleaveAfter[i],
                            minLength, maxLength))
                    {
                        int oneBasedStartResidue = oneBasedIndicesToCleaveAfter[i] + 1;
                        int oneBasedEndResidue = oneBasedIndicesToCleaveAfter[i + missedCleavages + 1];

                        foreach (var (threePrimeTerminus, fivePrimeTerminus) in GetDigestedTermini(oneBasedStartResidue, oneBasedEndResidue, nucleicAcid, ThreePrimeTerminusRemainder, FivePrimeTerminusRemainder))
                        {
                            yield return new NucleolyticOligo(nucleicAcid, oneBasedStartResidue, oneBasedEndResidue,
                                missedCleavages, CleavageSpecificity.Full, fivePrimeTerminus, threePrimeTerminus);
                        }
                    }
                }

                // Also digest using the truncation products start/end indices
                foreach (var truncation in nucleicAcid.TruncationProducts)
                {
                    if (truncation.OneBasedBeginPosition == 1 && truncation.OneBasedEndPosition == nucleicAcid.Length)
                        continue;

                        int cleavageIndexWithinTruncation = 0;
                        //get the first cleavage index after the start of the truncation
                        while (oneBasedIndicesToCleaveAfter[cleavageIndexWithinTruncation] < truncation.OneBasedBeginPosition)
                        {
                            cleavageIndexWithinTruncation++;
                        }

                        bool startPeptide = cleavageIndexWithinTruncation + missedCleavages < oneBasedIndicesToCleaveAfter.Count //if the current missed cleavages doesn't hit the end
                                && oneBasedIndicesToCleaveAfter[cleavageIndexWithinTruncation + missedCleavages] <= truncation.OneBasedEndPosition //and the cleavage occurs before the proteolytic end
                                && truncation.OneBasedBeginPosition.HasValue //and the proteolytic peptide even has a beginning
                                && !oneBasedIndicesToCleaveAfter.Contains(truncation.OneBasedBeginPosition.Value - 1) //and we haven't already cleaved here
                                && ValidLength(oneBasedIndicesToCleaveAfter[cleavageIndexWithinTruncation + missedCleavages] - truncation.OneBasedBeginPosition.Value + 1, minLength, maxLength); //and it's the correct size
                        if (startPeptide)
                        {
                            foreach (var (threePrimeTerminus, fivePrimeTerminus) in GetDigestedTermini(truncation.OneBasedBeginPosition, truncation.OneBasedEndPosition, nucleicAcid, ThreePrimeTerminusRemainder, FivePrimeTerminusRemainder))
                            {
                                yield return new NucleolyticOligo(nucleicAcid, truncation.OneBasedBeginPosition.Value, oneBasedIndicesToCleaveAfter[cleavageIndexWithinTruncation + missedCleavages],
                                missedCleavages, CleavageSpecificity.Full, fivePrimeTerminus, threePrimeTerminus, truncation.Type + " start");
                            }
                        }

                        //get the cleavage index before the end of the proteolysis product
                        while (oneBasedIndicesToCleaveAfter[cleavageIndexWithinTruncation] < truncation.OneBasedEndPosition)
                        {
                            cleavageIndexWithinTruncation++;
                        }

                        bool endPeptide = cleavageIndexWithinTruncation - missedCleavages - 1 >= 0 //if we're not going to go out of bounds (-1 to get in front of the end)
                                          && oneBasedIndicesToCleaveAfter[cleavageIndexWithinTruncation - missedCleavages - 1] + 1 >= truncation.OneBasedBeginPosition //and it's not before the beginning
                                          && truncation.OneBasedEndPosition.HasValue //and the proteolytic peptide even has an end
                                          && !oneBasedIndicesToCleaveAfter.Contains(truncation.OneBasedEndPosition.Value) //and we haven't already cleaved here
                                          && ValidLength(truncation.OneBasedEndPosition.Value - oneBasedIndicesToCleaveAfter[cleavageIndexWithinTruncation - missedCleavages - 1] + 1 - 1, minLength, maxLength); //and it's the correct size
                        if (endPeptide)
                        {
                            foreach (var (threePrimeTerminus, fivePrimeTerminus) in GetDigestedTermini(truncation.OneBasedBeginPosition, truncation.OneBasedEndPosition, nucleicAcid, ThreePrimeTerminusRemainder, FivePrimeTerminusRemainder))
                            {
                                yield return new NucleolyticOligo(nucleicAcid, oneBasedIndicesToCleaveAfter[cleavageIndexWithinTruncation - missedCleavages - 1] + 1, truncation.OneBasedEndPosition.Value,
                                    missedCleavages, CleavageSpecificity.Full, fivePrimeTerminus, threePrimeTerminus, truncation.Type + " end");
                            }
                        }
                    }
            }

            //add intact truncation (if acceptable)
            foreach (var truncation in nucleicAcid.TruncationProducts)
            {
                if (!truncation.OneBasedBeginPosition.HasValue 
                    || !truncation.OneBasedEndPosition.HasValue 
                    || !ValidLength(truncation.OneBasedEndPosition.Value - truncation.OneBasedBeginPosition.Value, minLength, maxLength) //if it's not the correct size
                    || oneBasedIndicesToCleaveAfter.Contains(truncation.OneBasedBeginPosition.Value - 1) //or we have already cleaved here
                    || oneBasedIndicesToCleaveAfter.Contains(truncation.OneBasedEndPosition.Value)) //or we have already cleaved there
                    continue; 

                int firstCleavage = 0;
                //get the first cleavage index after the start of the proteolysis product
                while (oneBasedIndicesToCleaveAfter[firstCleavage] < truncation.OneBasedBeginPosition)
                {
                    firstCleavage++;
                }

                int lastCleavage = firstCleavage;
                //get the last cleavage index before the end of the proteolysis product
                while (oneBasedIndicesToCleaveAfter[lastCleavage] < truncation.OneBasedEndPosition)
                {
                    lastCleavage++;
                }

                //if there are too many missed cleavages
                if (lastCleavage - firstCleavage >= maxMissedCleavages) 
                    continue; 

                foreach (var (threePrimeTerminus, fivePrimeTerminus) in GetDigestedTermini(truncation.OneBasedBeginPosition.Value, truncation.OneBasedEndPosition.Value, nucleicAcid, ThreePrimeTerminusRemainder, FivePrimeTerminusRemainder))
                {
                    yield return new NucleolyticOligo(nucleicAcid, truncation.OneBasedBeginPosition.Value, truncation.OneBasedEndPosition.Value,
                        lastCleavage - firstCleavage, CleavageSpecificity.Full, fivePrimeTerminus, threePrimeTerminus, truncation.Type + " end");

                }
            }
        }

        private static IEnumerable<(IHasChemicalFormula ThreePrime, IHasChemicalFormula FivePrime)> GetDigestedTermini(int? oligoStartIndex, int? oligoEndIndex, NucleicAcid nucleicAcid, IList<IHasChemicalFormula> threePrimeTerminusRemainder, IList<IHasChemicalFormula> fivePrimeTerminusRemainder)
        {
            // contains original 5' terminus ? keep it : use all rnase-specific remainders
            bool isOriginalFivePrimeTerminus = oligoStartIndex == 1;

            // contains original 3' terminus ? keep it : use all rnase-specific remainders
            bool isOriginalThreePrimeTerminus = oligoEndIndex == nucleicAcid.Length;

            if (isOriginalThreePrimeTerminus && isOriginalFivePrimeTerminus)
            {
                yield return (nucleicAcid.ThreePrimeTerminus, nucleicAcid.FivePrimeTerminus);
            }
            else if (isOriginalThreePrimeTerminus)
            {
                foreach (var fivePrime in fivePrimeTerminusRemainder)
                    yield return (nucleicAcid.ThreePrimeTerminus, fivePrime);
            }
            else if (isOriginalFivePrimeTerminus)
            {
                foreach (var threePrime in threePrimeTerminusRemainder)
                    yield return (threePrime, nucleicAcid.FivePrimeTerminus);
            }
            else
            {
                foreach (var threePrime in threePrimeTerminusRemainder)
                    foreach (var fivePrime in fivePrimeTerminusRemainder)
                        yield return (threePrime, fivePrime);
            }
        }

        public bool Equals(Rnase? other)
        {
            if (ReferenceEquals(null, other)) return false;
            if (ReferenceEquals(this, other)) return true;
            return Name == other.Name;
        }

        public override bool Equals(object? obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;
            return Equals((Rnase)obj);
        }

        public override int GetHashCode()
        {
            return Name.GetHashCode();
        }

        public override string ToString()
        {
            return Name;
        }
    }
}
