using System.Text;

namespace Proteomics
{
    public class Modification
    {

        #region Public Fields

        public readonly string id;
        public readonly string modificationType;

        #endregion Public Fields

        #region Public Constructors

        public Modification(string id, string modificationType)
        {
            this.id = id;
            this.modificationType = modificationType;
        }

        #endregion Public Constructors

        #region Public Methods

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
            return o == null ? false : this.id == m.id;
        }

        public override int GetHashCode()
        {
            return this.id == null ? 0 : this.id.GetHashCode();
        }

        #endregion Public Methods

    }
}