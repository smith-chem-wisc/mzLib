using System.Collections.Generic;

namespace Proteomics
{
    public class ModificationComparer : IEqualityComparer<Modification>
    {
        public bool Equals(Modification m1, Modification m2)
        {
            return m1.Equals(m2);
        }

        public int GetHashCode(Modification m1)
        {
            return m1.GetCustomHashCode();
        }
    }
}
