using Chemistry;
using Omics.Digestion;
using Omics.Modifications;

namespace Transcriptomics.Digestion
{
    public class Rnase : DigestionAgent, IEquatable<Rnase>
    {
        public static IHasChemicalFormula DefaultThreePrimeTerminus = ChemicalFormula.ParseFormula("H2O4P"); // Makes 3' Phosphate
        public static IHasChemicalFormula DefaultFivePrimeTerminus = ChemicalFormula.ParseFormula("O-3P-1"); // Makes 5' -OH by removing phosphate
        public Rnase(string name, CleavageSpecificity cleaveSpecificity, List<DigestionMotif> motifList, Modification cleavageMod = null) :
            base(name, cleaveSpecificity, motifList, cleavageMod)
        {
            CleavageSpecificity = cleaveSpecificity;
            DigestionMotifs = motifList;
        }

        public IEnumerable<NucleolyticOligo> GetUnmodifiedOligos(NucleicAcid nucleicAcid, int maxMissedCleavages, int minLength,
            int maxLength)
        {
            return CleavageSpecificity switch
            {
                // top down
                CleavageSpecificity.None => TopDownDigestion(nucleicAcid, minLength, maxLength),
                // full cleavage
                CleavageSpecificity.Full => FullDigestion(nucleicAcid, maxMissedCleavages, minLength, maxLength),
                _ => throw new ArgumentException(
                    "Cleave Specificity not defined for Rna digestion, currently supports Full and None")
            };
        }

        private IEnumerable<NucleolyticOligo> TopDownDigestion(NucleicAcid nucleicAcid, int minLength, int maxLength)
        {
            if (ValidLength(nucleicAcid.Length, minLength, maxLength))
                yield return new NucleolyticOligo(nucleicAcid, 1, nucleicAcid.Length,
                    0, CleavageSpecificity.Full, nucleicAcid.FivePrimeTerminus, nucleicAcid.ThreePrimeTerminus);

            // Also digest using the proteolysis product start/end indices
            foreach (var truncationProduct in nucleicAcid.TruncationProducts)
            {
                if (truncationProduct is { OneBasedEndPosition: not null, OneBasedBeginPosition: not null })
                {
                    int length = truncationProduct.OneBasedEndPosition.Value - truncationProduct.OneBasedBeginPosition.Value + 1;
                    if (!ValidLength(length, minLength, maxLength)) continue;

                    var (threePrimeTerminus, fivePrimeTerminus) = GetDigestedTermini(truncationProduct.OneBasedBeginPosition.Value,
                        truncationProduct.OneBasedEndPosition.Value, nucleicAcid);

                    yield return new NucleolyticOligo(nucleicAcid, truncationProduct.OneBasedBeginPosition.Value, truncationProduct.OneBasedEndPosition.Value,
                        0, CleavageSpecificity.Full, fivePrimeTerminus, threePrimeTerminus, truncationProduct.Type);
                }
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

                        var (threePrimeTerminus, fivePrimeTerminus) = GetDigestedTermini(oneBasedStartResidue, oneBasedEndResidue, nucleicAcid);
                        yield return new NucleolyticOligo(nucleicAcid, oneBasedStartResidue, oneBasedEndResidue,
                            missedCleavages, CleavageSpecificity.Full, fivePrimeTerminus, threePrimeTerminus);
                    }
                }

                // Also digest using the truncation products start/end indices
                foreach (var truncation in nucleicAcid.TruncationProducts)
                {
                    if (truncation.OneBasedBeginPosition == 1 && truncation.OneBasedEndPosition == nucleicAcid.Length)
                        continue;
            
                    var (threePrimeTerminus, fivePrimeTerminus) = GetDigestedTermini(truncation.OneBasedBeginPosition,
                        truncation.OneBasedEndPosition, nucleicAcid);

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
                        yield return new NucleolyticOligo(nucleicAcid, truncation.OneBasedBeginPosition.Value, oneBasedIndicesToCleaveAfter[cleavageIndexWithinTruncation + missedCleavages],
                            missedCleavages, CleavageSpecificity.Full, fivePrimeTerminus, threePrimeTerminus, truncation.Type + " start");
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
                        yield return new NucleolyticOligo(nucleicAcid, oneBasedIndicesToCleaveAfter[cleavageIndexWithinTruncation - missedCleavages - 1] + 1, truncation.OneBasedEndPosition.Value,
                            missedCleavages, CleavageSpecificity.Full, fivePrimeTerminus, threePrimeTerminus, truncation.Type + " end");
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

                var (threePrimeTerminus, fivePrimeTerminus) = GetDigestedTermini(truncation.OneBasedBeginPosition.Value, truncation.OneBasedEndPosition.Value, nucleicAcid);
                yield return new NucleolyticOligo(nucleicAcid, truncation.OneBasedBeginPosition.Value, truncation.OneBasedEndPosition.Value,
                    lastCleavage - firstCleavage, CleavageSpecificity.Full, fivePrimeTerminus, threePrimeTerminus, truncation.Type + " end");
            }
        }

        private static (IHasChemicalFormula ThreePrime, IHasChemicalFormula FivePrime) GetDigestedTermini(int? oligoStartIndex, int? oligoEndIndex, NucleicAcid nucleicAcid)
        {
            // contains original 5' terminus ? keep it : set to OH
            IHasChemicalFormula fivePrimeTerminus = oligoStartIndex == 1 ? nucleicAcid.FivePrimeTerminus : DefaultFivePrimeTerminus;

            // contains original 3' terminus ? keep it : set to phosphate
            IHasChemicalFormula threePrimeTerminus = oligoEndIndex == nucleicAcid.Length ? nucleicAcid.ThreePrimeTerminus : DefaultThreePrimeTerminus;

            return (threePrimeTerminus, fivePrimeTerminus);
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
