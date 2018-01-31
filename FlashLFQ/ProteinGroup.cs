using System.Text;

namespace FlashLFQ
{
    public class ProteinGroup
    {
        #region Public Fields

        public readonly string proteinGroupName;
        public double[] intensitiesByFile;
        public string[] peptidesByFile;

        #endregion Public Fields

        #region Public Constructors

        public ProteinGroup(string name)
        {
            proteinGroupName = name;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            sb.Append("" + proteinGroupName + '\t');
            sb.Append(string.Join("\t", intensitiesByFile));
            if (peptidesByFile != null)
            {
                sb.Append("\t");
                sb.Append(string.Join("\t", peptidesByFile));
            }

            return sb.ToString();
        }

        #endregion Public Methods
    }
}