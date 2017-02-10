using System.Text;

namespace Proteomics
{
    public class Modification
    {
        public readonly string id;

        public Modification(string id)
        {
            this.id = id;
        }
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("ID   " + id);
            return sb.ToString();
        }
    }
}