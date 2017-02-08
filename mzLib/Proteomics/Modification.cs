using System;
using System.Collections.Generic;
using System.Text;

namespace Proteomics
{
    public class Modification : BaseModification
    {
        #region Public Fields

        public readonly Tuple<string, string> ac;
        public readonly Dictionary<string, HashSet<string>> linksToOtherDbs;
        public readonly ModificationSites position;
        public readonly string site;
        public readonly string database;

        #endregion Public Fields

        #region Public Constructors
        
        public Modification(string id, Tuple<string, string> uniprotAC, string uniprotTG, ModificationSites uniprotPP, Dictionary<string, HashSet<string>> linksToOtherDbs, string database):base(id)
        {
            this.ac = uniprotAC;
            this.site = uniprotTG;
            this.position = uniprotPP;
            this.linksToOtherDbs = linksToOtherDbs;
            this.database = database;
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("AC   " + ac);
            if (linksToOtherDbs != null)
                foreach (var nice in linksToOtherDbs)
                    foreach (var db in nice.Value)
                        sb.AppendLine("DR   " + nice.Key + "; " + db);
            sb.AppendLine("PP   " + position);
            sb.AppendLine("DB   " + database);
            sb.Append("TG   " + site);
            return sb.ToString();
        }

        #endregion Public Methods
    }
}