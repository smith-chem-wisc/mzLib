using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    public class ProteinGroup
    {
        public readonly string proteinGroupName;
        public double[] intensitiesByFile;
        public string[] peptidesByFile;

        public ProteinGroup(string name)
        {
            proteinGroupName = name;
        }

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
    }
}
