using Chemistry;
using Omics.Digestion;
using Omics.Modifications;

namespace Transcriptomics.Digestion
{
    public class Rnase : DigestionAgent, IEquatable<Rnase>
    {
        public Rnase(string name, CleavageSpecificity cleaveSpecificity, List<DigestionMotif> motifList, Modification cleavageMod = null) :
            base(name, cleaveSpecificity, motifList, cleavageMod)
        {
            Name = name;
            CleavageSpecificity = cleaveSpecificity;
            DigestionMotifs = motifList;
        }
        public List<NucleolyticOligo> GetUnmodifiedOligos(NucleicAcid nucleicAcid, int maxMissedCleavages, int minLength,
            int maxLength)
        {
            var oligos = new List<NucleolyticOligo>();

            // top down
            if (CleavageSpecificity == CleavageSpecificity.None)
            {
                if (ValidLength(nucleicAcid.Length, minLength, maxLength))
                    oligos.Add(new NucleolyticOligo(nucleicAcid, 1, nucleicAcid.Length,
                        0, CleavageSpecificity.Full, nucleicAcid.FivePrimeTerminus, nucleicAcid.ThreePrimeTerminus));
            }
            // full cleavage
            else if (CleavageSpecificity == CleavageSpecificity.Full)
            {
                oligos.AddRange(FullDigestion(nucleicAcid, maxMissedCleavages, minLength, maxLength));
            }
            else
            {
                throw new ArgumentException(
                    "Cleave Specificity not defined for Rna digestion, currently supports Full and None");
            }

            return oligos;
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

                        // contains original 5' terminus ? keep it : set to OH
                        IHasChemicalFormula fivePrimeTerminus = oneBasedStartResidue == 1 ? nucleicAcid.FivePrimeTerminus : ChemicalFormula.ParseFormula("O-3P-1");

                        // contains original 3' terminus ? keep it : set to phosphate
                        IHasChemicalFormula threePrimeTerminus = oneBasedEndResidue == nucleicAcid.Length ? nucleicAcid.ThreePrimeTerminus : ChemicalFormula.ParseFormula("H2O4P");

                        yield return new NucleolyticOligo(nucleicAcid, oneBasedStartResidue, oneBasedEndResidue,
                            missedCleavages, CleavageSpecificity.Full, fivePrimeTerminus, threePrimeTerminus);
                    }
                }
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
