using System.Text;

namespace Proteomics
{
    public class Modification
    {
        public readonly string id;
        public readonly string modificationType;

        public Modification(string id, string modificationType)
        {
            this.id = id;
            this.modificationType = modificationType;
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("ID   " + id);
            sb.Append("MT   " + modificationType);
            return sb.ToString();
        }

        public override bool Equals(object o)
        {
            Modification m = o as Modification;
            return o != null
                && m.id == id
                && m.modificationType == modificationType;
        }

        public override int GetHashCode()
        {
            return id.GetHashCode() ^ modificationType.GetHashCode();
        }
    }
}