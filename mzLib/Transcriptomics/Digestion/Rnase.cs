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
                // non-specific, anchored at one terminus (see SingleFivePrimeDigestion)
                CleavageSpecificity.SingleN => SingleFivePrimeDigestion(nucleicAcid, minLength, maxLength),
                CleavageSpecificity.SingleC => SingleThreePrimeDigestion(nucleicAcid, minLength, maxLength),
                _ => throw new ArgumentException(
                    "Cleave Specificity not defined for Rna digestion, currently supports Full, None, SingleN and SingleC")
            };
        }

        /// <summary>
        /// Every oligo that starts at a given position and runs toward the 3' end, for every position.
        /// The counterpart of the proteomic SingleN digestion, which anchors at the N terminus; for a
        /// nucleic acid the anchored terminus is the 5' one, hence the name.
        ///
        /// Only one oligo is produced per start position -- the longest allowed one. That is the same
        /// thing Protease.SingleN_Digestion does: a non-specific search scores every prefix of this
        /// oligo as it walks the fragment index, so emitting the shorter ones here would be redundant
        /// work, not extra coverage.
        ///
        /// Unlike the proteomic version there is no specific-agent variant. Protease.SingleN_Digestion
        /// has a second branch honouring the missed-cleavage rules of a named protease, driven by
        /// DigestionParams.SpecificProtease. RnaDigestionParams has no SpecificRnase, so there is
        /// nothing to honour, and a semi-specific nucleic acid search is not expressible yet.
        /// </summary>
        private IEnumerable<NucleolyticOligo> SingleFivePrimeDigestion(NucleicAcid nucleicAcid, int minLength, int maxLength)
        {
            // maxLength is int.MaxValue by default, so start + maxLength overflows to negative
            bool maxTooBig = nucleicAcid.Length + maxLength < 0;

            for (int oneBasedStartResidue = 1; oneBasedStartResidue <= nucleicAcid.Length; oneBasedStartResidue++)
            {
                // is the longest oligo available from here still long enough?
                if (!ValidMinLength(nucleicAcid.Length - oneBasedStartResidue + 1, minLength))
                {
                    continue;
                }

                int oneBasedEndResidue = maxTooBig
                    ? nucleicAcid.Length
                    : Math.Min(nucleicAcid.Length, oneBasedStartResidue + maxLength - 1);

                var (threePrimeTerminus, fivePrimeTerminus) = GetDigestedTermini(oneBasedStartResidue, oneBasedEndResidue, nucleicAcid);

                yield return new NucleolyticOligo(nucleicAcid, oneBasedStartResidue, oneBasedEndResidue,
                    0, CleavageSpecificity.SingleN, fivePrimeTerminus, threePrimeTerminus, "SingleN");
            }
        }

        /// <summary>
        /// Every oligo that ends at a given position and runs back toward the 5' end, for every
        /// position. The counterpart of the proteomic SingleC digestion; the anchored terminus here is
        /// the 3' one. See SingleFivePrimeDigestion for why only the longest oligo per position is
        /// emitted and why there is no specific-agent variant.
        /// </summary>
        private IEnumerable<NucleolyticOligo> SingleThreePrimeDigestion(NucleicAcid nucleicAcid, int minLength, int maxLength)
        {
            bool maxTooBig = nucleicAcid.Length + maxLength < 0;

            for (int oneBasedEndResidue = 1; oneBasedEndResidue <= nucleicAcid.Length; oneBasedEndResidue++)
            {
                // is the longest oligo available back to the 5' end still long enough?
                if (!ValidMinLength(oneBasedEndResidue, minLength))
                {
                    continue;
                }

                int oneBasedStartResidue = maxTooBig
                    ? 1
                    : Math.Max(1, oneBasedEndResidue - maxLength + 1);

                var (threePrimeTerminus, fivePrimeTerminus) = GetDigestedTermini(oneBasedStartResidue, oneBasedEndResidue, nucleicAcid);

                yield return new NucleolyticOligo(nucleicAcid, oneBasedStartResidue, oneBasedEndResidue,
                    0, CleavageSpecificity.SingleC, fivePrimeTerminus, threePrimeTerminus, "SingleC");
            }
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
