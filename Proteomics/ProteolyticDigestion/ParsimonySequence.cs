using System;

namespace Proteomics.ProteolyticDigestion
{
    public class ParsimonySequence
    {
        public ParsimonySequence(PeptideWithSetModifications pwsm, bool TreatModPeptidesAsDifferentPeptides)
        {
            Sequence = TreatModPeptidesAsDifferentPeptides ? pwsm.FullSequence : pwsm.BaseSequence;
            Protease = pwsm.DigestionParams.Protease;
        }

        public string Sequence { get; }
        public Protease Protease { get; }


        public override bool Equals(object other)
        {
            ParsimonySequence that = (ParsimonySequence)other;

            if (that == null)
            {
                return false;
            }

            return (that.Sequence == this.Sequence && that.Protease == this.Protease);
        }

        public override int GetHashCode()
        {
            return this.Sequence.GetHashCode() ^ this.Protease.GetHashCode();
        }

        public override string ToString()
        {
            return this.Sequence.ToString();
        }
    }
}