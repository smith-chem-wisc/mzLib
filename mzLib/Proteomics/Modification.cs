using System.Text;
using System.Text.RegularExpressions;

namespace Proteomics
{
    public class Modification
    {

        #region Public Fields

        public readonly string id;

        #endregion Public Fields

        #region Public Constructors

        public Modification(string id)
        {
            this.id = id;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.Append("ID   " + id);
            return sb.ToString();
        }

        #endregion Public Methods
    }
}