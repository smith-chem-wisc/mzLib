using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;
using MassSpectrometry;

namespace Transcriptomics
{
    /// <summary>
    /// Class Representing a fragment from a Nucleic Acid molecule
    /// </summary>
    public class RNAFragment : IFragment<FragmentType, Terminus>
    {

        public double NeutralMass { get; }
        public double MonoisotopicMass => NeutralMass;
        public int FragmentNumber { get; }
        public int ResiduePosition { get; }
        public FragmentType ProductType { get; }
        public Terminus Terminus { get; }
        public double NeutralLoss { get; }
        public FragmentType? SecondaryProductType { get; }
        public int? SecondaryFragmentNumber { get; }
        public NucleicAcid Parent { get; private set; }
        public string Annotation => (this as IFragment<FragmentType, Terminus>).Annotation;

        public RNAFragment(NucleicAcid parent, FragmentType productType, Terminus terminus, double neutralMass,
            int fragmentNumber, int residuePosition, double neutralLoss = 0, FragmentType? secondaryProductType = null,
            int? secondaryFragmentNumber = 0)
        {
            Parent = parent;
            ProductType = productType;
            Terminus = terminus;
            NeutralMass = neutralMass;
            FragmentNumber = fragmentNumber;
            ResiduePosition = residuePosition;
            NeutralLoss = neutralLoss;
            SecondaryProductType = secondaryProductType;
            SecondaryFragmentNumber = secondaryFragmentNumber;
        }

        //TODO: Chemical Formula

        public IEnumerable<IHasMass> GetModifications()
        {
            if (Parent == null)
                yield break;

            var mods = Parent.Modifications;
            if (ProductType.GetTerminus() == Terminus.FivePrime)
            {
                for (int i = 0; i <= FragmentNumber; i++)
                {
                    if (mods[i] != null)
                        yield return mods[i];
                }
            }
            else
            {
                int length = Parent.Length + 1;
                for (int i = length - FragmentNumber; i <= length; i++)
                {
                    if (mods[i] != null)
                        yield return mods[i];
                }
            }
        }

        public string GetSequence()
        {
            if (Parent == null)
                return "";

            string parentSeq = Parent.BaseSequence;
            if (ProductType.GetTerminus() == Terminus.FivePrime)
            {
                return parentSeq.Substring(0, FragmentNumber);
            }

            return parentSeq.Substring(parentSeq.Length - FragmentNumber, FragmentNumber);
        }

        #region Interface Implementaiton and Overrides

        public override string ToString() => (this as IFragment<FragmentType, Terminus>).GetStringRepresentation();
        public bool Equals(IFragment<FragmentType, Terminus>? other)
        {
            return NeutralMass.Equals(other.NeutralMass) && FragmentNumber == other.FragmentNumber &&
                   ResiduePosition == other.ResiduePosition && ProductType == other.ProductType &&
                   Terminus == other.Terminus && NeutralLoss.Equals(other.NeutralLoss) &&
                   SecondaryProductType == other.SecondaryProductType &&
                   SecondaryFragmentNumber == other.SecondaryFragmentNumber;
        }

        public override bool Equals(object? obj)
        {
            if (ReferenceEquals(null, obj)) return false;
            if (ReferenceEquals(this, obj)) return true;
            if (obj.GetType() != this.GetType()) return false;
            return Equals((RNAFragment)obj);
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(NeutralMass, FragmentNumber, ResiduePosition, (int)ProductType, (int)Terminus,
                NeutralLoss, SecondaryProductType, SecondaryFragmentNumber);
        }

        #endregion
        
    }
}
