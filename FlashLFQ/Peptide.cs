using System.Text;

namespace FlashLFQ
{
    public class Peptide
    {
        public readonly string Sequence;
        public readonly string ProteinGroup;
        public readonly double[] intensitiesByFile;
        public readonly string[] detectionType;
        public static string[] files;

        public Peptide(string baseSeq, string proteinGroup, double[] intensitiesByFile, string[] detectionType)
        {
            Sequence = baseSeq;
            ProteinGroup = proteinGroup;
            this.intensitiesByFile = intensitiesByFile;
            this.detectionType = detectionType;
        }

        public static string TabSeparatedHeader
        {
            get
            {
                var sb = new StringBuilder();
                sb.Append("Sequence" + "\t");
                sb.Append("Protein Group" + "\t");
                for (int i = 0; i < files.Length; i++)
                    sb.Append("Intensity_" + files[i] + "\t");
                for (int i = 0; i < files.Length; i++)
                    sb.Append("Detection Type_" + files[i] + "\t");
                return sb.ToString();
            }
        }

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            
            sb.Append("" + Sequence + '\t');
            sb.Append("" + ProteinGroup + '\t');
            for (int i = 0; i < intensitiesByFile.Length; i++)
            {
                if (!(detectionType[i] == "MBR" && intensitiesByFile[i] == 0))
                    sb.Append(intensitiesByFile[i] + "\t");
                else
                    sb.Append("\t");
            }
            for (int i = 0; i < intensitiesByFile.Length; i++)
            {
                if (!(detectionType[i] == "MBR" && intensitiesByFile[i] == 0))
                    sb.Append(detectionType[i] + "\t");
                else
                    sb.Append("\t");
            }

            return sb.ToString();
        }
    }
}
