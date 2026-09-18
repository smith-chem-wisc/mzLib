using Chemistry;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
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
            int maxLength, Rnase? specificRnase = null, bool topDownTruncationSearch = false,
            FragmentationTerminus fragmentationTerminus = FragmentationTerminus.Both,
            CleavageSpecificity? searchModeType = null)
        {
            specificRnase ??= this;
            CleavageSpecificity requestedSpecificity = searchModeType == CleavageSpecificity.Semi
                ? CleavageSpecificity.Semi
                : CleavageSpecificity;
            return requestedSpecificity switch
            {
                // top down
                CleavageSpecificity.None => TopDownDigestion(nucleicAcid, minLength, maxLength, topDownTruncationSearch, 1, null, CleavageSpecificity.Full, "full").Cast<NucleolyticOligo>(),

                // full cleavage
                CleavageSpecificity.Full => FullDigestion(nucleicAcid, maxMissedCleavages, minLength, maxLength,
                    1, null, CleavageSpecificity.Full, null).Cast<NucleolyticOligo>(),

                // non-specific, anchored at one terminus
                CleavageSpecificity.SingleN => SingleLeftSideDigestion(nucleicAcid, maxMissedCleavages, minLength, maxLength, specificRnase, 1).Cast<NucleolyticOligo>(),

                CleavageSpecificity.SingleC => SingleRightSideDigestion(nucleicAcid, maxMissedCleavages, minLength, maxLength, specificRnase, 1).Cast<NucleolyticOligo>(),

                CleavageSpecificity.Semi when fragmentationTerminus is FragmentationTerminus.FivePrime or FragmentationTerminus.ThreePrime
                    => SpeedySemiSpecificDigestion(nucleicAcid, maxMissedCleavages, minLength, maxLength,
                        fragmentationTerminus == FragmentationTerminus.FivePrime, 1, null).Cast<NucleolyticOligo>(),

                _ => throw new ArgumentException(
                    "Cleavage specificity or terminus is not defined for RNA digestion; currently supports Full, None, SingleN and SingleC")
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
