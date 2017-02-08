using System.Text;

namespace Proteomics
{
    public class BaseModification
    {
        public readonly string id;

        public BaseModification(string id)
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