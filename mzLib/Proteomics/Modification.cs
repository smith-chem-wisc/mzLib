using System.Collections.Generic;
using System;
using System.Text;

namespace Proteomics
{
    public class Modification
    {

        public readonly Tuple<string, string> ac;
        public readonly Dictionary<string, HashSet<string>> linksToOtherDbs;
        public readonly string id;
        public readonly ModificationSites position;
        public readonly string site;


        public Modification(string id, Tuple<string, string> unimodAC, string tg, ModificationSites pos)
        {
            this.id = id;
            this.ac = unimodAC;
            this.site = tg;
            this.position = pos;
        }

        public Modification(string uniprotID, Tuple<string, string> uniprotAC, string uniprotTG, ModificationSites uniprotPP, Dictionary<string, HashSet<string>> uniprotDR)
        {
            this.id = uniprotID;
            this.ac = uniprotAC;
            this.site = uniprotTG;
            this.position = uniprotPP;
            this.linksToOtherDbs = uniprotDR;
        }
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("AC   " + ac);
            if (linksToOtherDbs != null)
                foreach (var nice in linksToOtherDbs)
                    foreach (var db in nice.Value)
                        sb.AppendLine("DR   " + nice.Key + "; " + db);
            sb.AppendLine("ID   " + id);
            sb.AppendLine("PP   " + position);
            sb.Append("TG   " + site);
            return sb.ToString();
        }
    }
}