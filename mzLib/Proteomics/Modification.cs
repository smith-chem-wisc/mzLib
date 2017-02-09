using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Proteomics
{
    public class Modification : BaseModification
    {

        #region Public Fields

        public readonly Tuple<string, string> accession;
        public readonly Dictionary<string, HashSet<string>> linksToOtherDbs;
        public readonly ModificationSites position;
        public readonly string site;
        public readonly string database;

        #endregion Public Fields

        #region Private Fields

        private static readonly Dictionary<string, ModificationSites> modificationTypeCodes;

        #endregion Private Fields

        #region Public Constructors

        static Modification()
        {
            modificationTypeCodes = new Dictionary<string, ModificationSites>();
            modificationTypeCodes.Add("N-terminal.", ModificationSites.NProt); // Implies protein only, not peptide
            modificationTypeCodes.Add("C-terminal.", ModificationSites.ProtC);
            modificationTypeCodes.Add("Peptide N-terminal.", ModificationSites.NPep);
            modificationTypeCodes.Add("Peptide C-terminal.", ModificationSites.PepC);
            modificationTypeCodes.Add("Anywhere.", ModificationSites.Any);
            modificationTypeCodes.Add("Protein core.", ModificationSites.Any);
        }

        public Modification(string id, Tuple<string, string> accession, string site, ModificationSites position, Dictionary<string, HashSet<string>> linksToOtherDbs, string database) : base(id)
        {
            this.accession = accession;
            this.site = site;
            this.position = position;
            this.linksToOtherDbs = linksToOtherDbs;
            this.database = database;
        }

        public Modification(string id, Tuple<string, string> accession, string site, string position, Dictionary<string, HashSet<string>> linksToOtherDbs, string database)
            : this(id, accession, site, modificationTypeCodes[position], linksToOtherDbs, database)
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            sb.AppendLine("PP   " + modificationTypeCodes.First(b => b.Value.Equals(position)).Key);
            sb.AppendLine("TG   " + site);
            if (linksToOtherDbs != null)
                foreach (var nice in linksToOtherDbs)
                    foreach (var db in nice.Value)
                        sb.AppendLine("DR   " + nice.Key + "; " + db);
            sb.Append("DB   " + database);
            return sb.ToString();
        }

        #endregion Public Methods

    }
}