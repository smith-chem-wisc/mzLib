using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Chemistry;

namespace Transcriptomics
{
    public class Fragment : IHasMass, IEquatable<Fragment>
    {
        public Fragment(FragmentType type, int number, double monoisotopicMass, NucleicAcid parent)
        {
            Type = type;
            Number = number;
            Parent = parent;
            MonoisotopicMass = monoisotopicMass;
        }

        public double MonoisotopicMass { get; private set; }

        public int Number { get; private set; }

        public NucleicAcid Parent { get; private set; }

        public FragmentType Type { get; private set; }

        //TODO: Chemical Formula

        public IEnumerable<IHasMass> GetModifications()
        {
            if (Parent == null)
                yield break;

            var mods = Parent.Modifications;
            if (Type.GetTerminus() == Terminus.FivePrime)
            {
                for (int i = 0; i <= Number; i++)
                {
                    if (mods[i] != null)
                        yield return mods[i];
                }
            }
            else
            {
                int length = Parent.Length + 1;
                for (int i = length - Number; i <= length; i++)
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

            string parentSeq = Parent.Sequence;
            if (Type.GetTerminus() == Terminus.FivePrime)
            {
                return parentSeq.Substring(0, Number);
            }

            return parentSeq.Substring(parentSeq.Length - Number, Number);
        }

        #region Interface Implementaiton and Overrides

        public override string ToString()
        {
            return $"{Enum.GetName(typeof(FragmentType), Type)}{Number}";
        }

        public override int GetHashCode()
        {
            unchecked
            {
                int hCode = 23;
                hCode = hCode * 31 + Number;
                hCode = hCode * 31 + (int)Type;
                hCode = hCode * 31 + Math.Round(MonoisotopicMass).GetHashCode();
                return hCode;
            }
        }

        public bool Equals(Fragment other)
        {
            return Type.Equals(other.Type) && Number.Equals(other.Number) && MonoisotopicMass.MassEquals(other.MonoisotopicMass);
        }

        #endregion
    }
}
